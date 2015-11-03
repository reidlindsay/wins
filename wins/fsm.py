#!  /usr/bin/env python

"""
Implementation of finite state machine; contains `FSM` and `Timer`.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-27 22:15:57 -0500 (Tue, 27 Sep 2011) $
* $LastChangedRevision: 5167 $

:author: Ketan Mandke <kmandke@mail.utexas.edu>

:copyright:
    Copyright 2009 The University of Texas at Austin

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

:var FSM_VERBOSE:
    Constant enumeration to control verbose thresholds of `FSM`.
  
    Setting the verbose level of a `FSM` above this threshold will cause the
    corresponding output in this file to be written (or logged).

:var TIMER_VERBOSE:
    Constant enumeration to control verbose threshold of `Timer`.
  
    Setting the verbose level of a `Timer` above this threshold will cause the
    corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import Process, SimEvent, Simulation
from SimPy.Simulation import Resource, Store, Level
from SimPy.Simulation import hold, passivate, activate, reactivate, now
from wins.trace import Traceable
from wins.queue import Queue
from wins.helper import time2usec
from wins import const

from SimPy.Simulation import FatalSimerror

FSM_VERBOSE=90
TIMER_VERBOSE=95

fsmhold = 12345     # use magic number to avoid conflict

def fsmholdfunc(a):
    _fsmhold(a[0][1], a)

def _fsmhold(proc, a):
    if len(a[0]) == 3: ## yield hold,self,delay
        delay = a[0][2]
        if delay < 0:
            raise FatalSimerror('hold: delay time negative: %s, in %s' % (
                                 delay, str(a[0][1])))
    else:              ## yield hold,self     
        delay = 0
    who = a[1]
    proc.interruptLeft = delay
    proc._inInterrupt = False
    proc.interruptCause = None
    proc.sim._post(what = who, at = proc.sim._t + delay, prior=1)

def _fixreactivate(self, *args, **kwargs):
    """Fix to make all reactivate() calls prior by default."""
    if (len(args)<4) and ('prior' not in kwargs):
        # prior is not defined -> set to True
        kwargs['prior'] = True
    return Simulation.reactivate(self, *args, **kwargs)

class FSM(Traceable, Process):
    """Finite state machine.

    State Execution Methods
    =======================

    State execution methods must accept the calling `FSM` as an argument. For
    example:

        >>> def S0(fsm):
        ...     yield hold, fsm, 1.0
        ...     yield fsm.goto(S1, 1.0)
        >>> def S1(fsm, d):
        ...     yield hold, fsm, d
        ...     yield fsm.goto(S2)
        >>> def S2(fsm):
        ...     yield hold, fsm, 0  # do nothing and kill fsm
        >>> fsm = FSM()
        >>> fsm.goto(S0)
        >>> fsm.start()

    creates a finite state machine that transitions through three states (S0,
    S1, and S2). Note that state S1 accepts any additional argument which was
    passed to it by using `goto()`.

    Use `goto()` to transition between states and to set the initial state for
    the FSM. Use `start()` to start the FSM; using `activate()` may result in
    unintended consequences.

    FSM and Ports
    =============
    This class overloads the `acquired()` and `stored()` methods from SimPy
    Process. These overloaded methods cause the `FSM` to clear any stale
    SimEvents which might be present in the `eventFired` list. This usually only
    relevant when using renege clauses with get or put requests from a SimPy
    Store (see `Queue` for more on getting/putting with a renege clause).

    :CVariables:
     * `name`: Name of FSM.
     * `tracename`: Name used in trace.
     * `sem`: Property to access generator running current state execution
       method.
     * `state`: Property to access current state of FSM, i.e. current state
       execution method.
     * `statename`: Property to get readable name of current `state`.
     * `started`: Property to indicate if `FSM` has been started.

    :IVariables:
     * `__sem`: Private variable for generator running current state execution
       method.
     * `__state`: Private variable for current state of FSM, i.e. current state
       execution method. Use `goto()` to transition to new state or set state
       before starting state machine.
     * `__started`: Private boolean flag indicates if `FSM` has been started.

    :note: Since `FSM` is a SimPy Process. All normal SimPy Process commands,
           such as `cancel()` and `interrupt()` may be used.
    """
    name = "fsm"
    tracename = "FSM"
    _fixsimhold = False
    def __init__(self, initstate=None, start=False, **kwargs):
        """Constructor.

        :param initstate: Initial state to `goto`.
        :param start: Boolean; if true, `start()` immediately.
        :param kwargs: Keywords passed to `Traceable` constructor.
        """
        Traceable.__init__(self, **kwargs)
        Process.__init__(self, name=self.name)
        self.__sem = None
        self.__state = None
        self.__started = False
        self.__initstate = self.HALT, (), {}
        # patch SimPy to fix hold
        if not FSM._fixsimhold:
            self.sim._dispatch[fsmhold] = fsmholdfunc
            self.sim._commandcodes = self.sim._dispatch.keys()
            self.sim._commandwords[fsmhold] = 'fsmhold'
            self.sim.reactivate = lambda *args, **kwargs: _fixreactivate(self.sim, *args, **kwargs)
            FSM._fixsimhold = True
        # continue initializing FSM
        if initstate: self.goto(initstate)
        if start: self.start()

    sem = property(fget=lambda self: self.__sem)
    state = property(fget=lambda self: self.__state)
    statename = property(fget=lambda self: self.get_statename() )
    started = property(fget=lambda self: self.__started)

    def start(self, sem=None, **kwargs):
        """Overloaded start to go to initial state.

        :param sem: State execution method.
        :param kwargs: Keyword arguments passed to `Process.start()`.

        State execution method (SEM) is a function that takes a `FSM` as its
        first argument. This function is like a Process execution method (PEM)
        in SimPy, but meant to be used exclusively with `FSM`.
        """
        if self.started: return
        assert ('pem' not in kwargs)
        if sem: self.__initstate = sem, (), {}
        # set initial state and start FSM
        ns, sargs, swargs = self.__initstate
        self.__started = True
        self.goto(ns, *sargs, **swargs)
        Process.start(self, pem=self.sem, **kwargs)

    def goto(self, nextstate, *args, **kwargs):
        """Transition to a new state.

        :param nextstate: Generator method for next state.
        :param args: Arguments passed to next state.
        :param kwargs: Keywords passed to next state.

        :note: `nextstate` must be a valid state execution method, see above.
        """
        if self.started:
            self.__state = nextstate
            self.__sem = nextstate(self, *args, **kwargs)
            self._nextpoint = self.sem
            if self.verbose>FSM_VERBOSE: self.log(self.statename)
        else:
            self.__initstate = nextstate, args, kwargs
        return fsmhold, self, 0

    def stop(self, *args, **kwargs):
        """Stop execution and transition to `HALT`.

        :param args: Additional arguments passed to `HALT`.
        :param kwargs: Additional keywords passed to `HALT`.

        This should only be used for "normal" transitions to `HALT`:

            >>> yield fsm.stop()

        Use `halt()` to force transitions on active `FSM` to the `HALT` state,
        e.g. when an `FSM` is waiting on an event or a queue.
        """
        #if self.active():
        #    return self.halt(*args, **kwargs)
        return self.goto(self.HALT, *args, **kwargs)

    def acquired(self, buf):
        """Overload to clear `eventsFired` when successfully acquired."""
        try:
            a = Process.acquired(self, buf)
        except AttributeError:
            # get action completed AND there was no renege clause
            if isinstance(buf, Resource): test = self in buf.activeQ
            elif isinstance(buf, Store):  test = len(self.got)
            elif isinstance(buf, Level):  test = not (self.got is None)
            else:                         test = False
            errmsg = "[FSM]: Unexpected failure in acquired()!"
            assert (test), errmsg
            a = test
        if a: self.eventsFired = []
        return a

    def stored(self, buffer):
        """Overload to clear `eventsFired` when successfully stored."""
        test = self in buffer.putQ  # test = True -> put failed
        try:
            a = Process.stored(self, buffer)
        except AttributeError:
            # put action completed AND there was no renege clause
            errmsg = "[FSM]: Unexpected failure in stored()!"
            assert (not test), errmsg
            a = not test
        if a: self.eventsFired = []
        return a

    def halt(self, *args, **kwargs):
        """Halt execution of `FSM`; send to `HALT` state.

        :param args: Additional arguments passed to `HALT` state.
        :param kwargs: Additional keywords passed to `HALT` state.

        This method forces the current running `FSM` to halt execution by use of
        the SimPy directive `cancel()`. If `FSM` is to be reactivated (or woken
        up) later, use `passivate` instead of `halt()`.

            >>> yield passivate, fsm
        """
        if self.terminated():
            return self.goto(self.HALT, *args, **kwargs)
        else:
            f = FSM()
            f.goto(self._RUN_HALT, *args, **kwargs)
            f.start()
            return hold, self, 0

    def _RUN_HALT(self, fsm, *args, **kwargs):
        """Internal SEM to halt execution of `FSM`."""
        fsm.cancel(self)
        self.goto(self.HALT, *args, **kwargs)
        reactivate(self)
        yield hold, fsm, 0

    def get_statename(self):
        """Return current state name."""
        if self.sem:
            return self.sem.gi_frame.f_code.co_name

    def HALT(self, fsm):
        """Default state to halt execution."""
        yield hold, fsm, 0      # NOP then go to sleep

    @classmethod
    def launch(cls, sem, *args, **kwargs):
        """Start a new `FSM` that starts running *immediately*.

        :param sem: State execution method (SEM).
        :param args: Additional arguments passed to SEM.
        :param kwargs: Additional keywords passed to SEM.
        :return: Newly created `FSM`.

        :note: The SEM will start running immediately.
        """
        f = cls()
        f.goto(sem, *args, **kwargs)
        f.start()
        return f

    def log(self, evt=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>FSM_VERBOSE)
        if self.verbose>FSM_VERBOSE or force:
            Traceable.log(self, evt, p, *args, **kwargs)

class Timer(FSM):
    """Stopwatch timer that fires an event after waiting for a specified time.

    After starting `Timer`, use `pause()`, `resume()`, and `stop()` to control
    the execution of the timer. The `done` and `kill` signals to monitor the
    completion of the `Timer`.

    :CVariables:
     * `CMDPAUSE`: Internal constant for pause command.
     * `CMDRESUME`: Internal constant for resume command.
     * `duration`: Property to access time that `Timer` should run.
     * `timepassed`: Property to access time elapsed by `Timer`.
     * `timeleft`: Property to access time left on `Timer`.
     * `ispaused`: Property to check if timer is currently paused.
     * `fired`: Property to check if `done` occurred.
     * `running`: Property to check if `started` and currently running.
     * `stopped`: Property to check if not `fired` and no longer active.
     * `ctrlQ`: Internal `Queue` for control messages.

    :IVariables:
     * `done`: SimEvent signalled when timer finishes successfully.
     * `kill`: SimEvent signalled when timer is prematurely stopped.
    """
    name="timer"
    tracename="TIMER"
    CMDPAUSE  = "pause"
    CMDRESUME = "resume"
    def __init__(self, duration, start=False, initstate=None, **kwargs):
        """Constructor.

        :param duration: Time for which `Timer` should run.
        :param start: Boolean; if true, `start()` immediately.
        :param initstate: Depricated (do not use).
        :param kwargs: Keywords passed to `FSM` constructor.
        """
        assert (duration>0), "[TIMER]: Cannot simulate non-positive duration!"
        FSM.__init__(self, start=False, initstate=self.RUN, **kwargs)
        self.__tpassed = 0
        self.__duration  = duration
        self.__ctrlQ = Queue()
        self.__tic   = None
        self.done = SimEvent(name=self.name+".done")
        self.kill = SimEvent(name=self.name+".kill")
        if start: self.start()

    duration   = property(fget=lambda self: self.__duration)
    timepassed = property(fget=lambda self: self.__timepassed() )
    timeleft   = property(fget=lambda self: self.duration - self.timepassed)
    ispaused   = property(fget=lambda self: (self.state==self.PAUSE) )
    fired      = property(fget=lambda self: self.done.occurred)
    running    = property(fget=lambda self: \
                          self.started and (self.state==self.RUN) )
    stopped    = property(fget=lambda self: \
                          (not self.fired) and (self.state==self.HALT) )
    ctrlQ      = property(fget=lambda self: self.__ctrlQ)

    def RUN(self, fsm):
        """RUN state; timer is active.

        Call `pause()` to pause an active timer. This method will signal `done`
        upon completion.
        """
        # set parameters
        tstart = now()
        self.__tic = tstart      # temporarily store tic as start of RUN
        queue  = self.ctrlQ
        tleft  = self.duration - self.__tpassed  # time left on timer
        # wait for timer to expire (or be paused)
        yield queue.remove(fsm, 1, renege=tleft)
        self.__tic = None
        telapsed = now() - tstart
        self.__tpassed += telapsed
        if fsm.acquired(queue):
            # PAUSE command
            assert (telapsed<tleft), \
                    "[TIMER]: Elapsed time exceeded time left during RUN!"
            assert (len(fsm.got)==1), "[TIMER]: Control queue failed!"
            cmd = fsm.got[0]
            assert (cmd==self.CMDPAUSE), \
                   "[TIMER]: Invalid control command received in RUN!"
            yield fsm.goto(self.PAUSE)
        else:
            assert (abs(self.__tpassed-self.duration)<const.EPSILON), \
                    "[TIMER]: Timer failed to complete properly in RUN!"
            self.__tpassed = self.duration
            if self.verbose>TIMER_VERBOSE: self.log("DONE")
            self.done.signal()
        yield fsm.goto(self.HALT, force=False)

    def PAUSE(self, fsm):
        """PAUSE state; timer is paused.

        Call `resume()` to restart timer.
        """
        queue = self.ctrlQ
        yield queue.remove(fsm, 1)
        assert fsm.acquired(queue) and (len(fsm.got)==1), \
               "[TIMER]: PAUSE failed to dequeue control message!"
        cmd = fsm.got[0]
        assert (cmd==self.CMDRESUME), \
                "[TIMER]: Invalid control command received in PAUSE!"
        yield fsm.goto(self.RUN)

    def resume(self, proc):
        """Blocking call to restart timer if in `PAUSE` state.

        :param proc: Process to block on resume command.
        """
        if (self.state==self.PAUSE):
            return self.ctrlQ.insert(proc, [self.CMDRESUME])
        else:
            return hold, proc, 0

    def pause(self, proc):
        """Blocking call to pause timer if in `RUN` state.

        :param proc: Process to block on pause command.
        """
        if (self.state==self.RUN):
            return self.ctrlQ.insert(proc, [self.CMDPAUSE])
        else:
            return hold, proc, 0

    def __timepassed(self):
        """Private method to determine time elapsed on timer."""
        if self.__tic is None:
            return self.__tpassed
        else:
            tdelta = now() - self.__tic
            return (self.__tpassed + tdelta)

    def HALT(self, fsm, force=True):
        """Overload `HALT` state to signal `kill` if needed."""
        if (force or (self.timeleft>0)):
            timepassed = self.timepassed
            if self.verbose>TIMER_VERBOSE:
                self.log("CANCEL", timepassed=time2usec(timepassed, fmt="%.4g"), \
                                   timeleft=self.timeleft, force=force)
            self.kill.signal(timepassed)
        yield hold, fsm, 0

    def log(self, evt=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>TIMER_VERBOSE)
        if self.verbose>TIMER_VERBOSE or force:
            FSM.log(self, evt, p, *args, **kwargs)
