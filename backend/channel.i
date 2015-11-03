/*  :file: channel.i
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-10-01 14:15:29 -0500 (Sat, 01 Oct 2011) $
 *  $LastChangedRevision: 5172 $
 *
 *  :author: Ketan Mandke <kmandke@mail.utexas.edu>
 *  
 *  :copyright:
 *    Copyright 2009-2010 The University of Texas at Austin
 *    
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    
 *      http://www.apache.org/licenses/LICENSE-2.0
 *    
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

%module(docstring="Wrapper for IT++ Array of Complex Matrices") channel

%include <std_string.i>
%include <std_complex.i>

%feature("autodoc","1");
%{
#include <wins-80211n/channel.h>
%}

%inline %{
  std::string channel_py_channel__repr(Channel *H) {
    std::string s;
    s = H->to_string();
    return s;
  }

  Channel channel_py_channel__add(Channel *A, Channel *B) {
    Channel H;
    for (int k=0; k< A->size(); k++) H(k) = (*A)(k) + (*B)(k);
    return H;
  }

  Channel channel_py_channel__sub(Channel *A, Channel *B) {
    Channel H;
    for (int k=0; k< A->size(); k++) H(k) = (*A)(k) - (*B)(k);
    return H;
  }

  Channel channel_py_channel__add(Channel *H, std::complex<double> x) {
    Channel R;
    R.set_size(H->size());
    for (int k=0; k<R.size(); k++)
      R(k) = (*H)(k) + x;
    return R;
  }

  Channel channel_py_channel__sub(Channel *H, std::complex<double> x) {
    Channel R;
    R.set_size(H->size());
    for (int k=0; k<R.size(); k++)
      R(k) = (*H)(k) - x;
    return R;
  }

  Waveform channel_py_channel__conv(Channel *H, Waveform *x) {
    Waveform y = (*H) * (*x);
    return y;
  }

  Channel channel_py_channel__mul(Channel *H, std::complex<double> t) {
    Channel R;
    R = mult(*H, t);
    return R;
  }

  int channel_py_channel__len(Channel *H) {
    return H->taps();
  }

  Waveform channel_py_channel__get(Channel *H, int k) {
    Waveform h = (*H)(k);
    return h;
  }

  void channel_py_channel__set(Channel *H, int k, Waveform *v) {
    (*H)(k) = (*v);
  }

  std::string channel_py_channel__flatten(Channel *H) {
    std::string s;
    s = H->to_string(true);
    return s;
  }
%}

%rename(MimoChannel) Channel;

class Channel {
  public:
    Channel();
    virtual ~Channel() {}

    Waveform convolve(const Waveform &);
    inline Waveform conv(const Waveform &);
    int rows();
    int cols();
    int taps();
    double power();
    inline double powerdb();
    double normpower();
    inline double normpowerdb();
    void normalize();
};

%pythoncode %{
from itertools import islice

def mimochannel_py__getitem(H, key):
    """Internal method to get taps from a Channel."""
    if isinstance(key, slice):
        args = key.indices(H.taps())
        itr = islice(range(H.taps()), *args)
        return [channel_py_channel__get(H,k) for k in itr]
    elif (-1<=(key/H.taps())<1):
        return channel_py_channel__get(H, key%H.taps())
    else:
        raise IndexError, "list index out of range"

def mimochannel_py__setitem(H, key, val):
    """Internal method to set taps of a Channel."""
    if isinstance(key, slice):
        args = key.indices(H.taps())
        itr = islice(range(H.taps()), *args)
        idx = [k for k in itr]
        if (len(idx) != len(val)):
            errmsg = "attempt to assign sequency of size %d "%(len(val)) + \
                     "to extended slice of size %d"%(len(idx))
            raise ValueError, errmsg
        if (max(idx)>=H.taps):
            errmsg = "list index out of range"
            raise IndexError, errmsg
        for k in range(len(idx)):
            channel_py_channel__set(H,idx[k], val[k])
    elif (-1<=(key/H.taps())<1):
        return channel_py_channel__set(H, key%H.taps(), val)
    else:
        raise IndexError, "list index out of range"

MimoChannel.__repr__ = channel_py_channel__repr
MimoChannel.__add__  = channel_py_channel__add
MimoChannel.__radd__ = channel_py_channel__add
MimoChannel.__sub__  = channel_py_channel__sub
MimoChannel.__rsub__ = channel_py_channel__sub
MimoChannel.__mul__  = channel_py_channel__mul
MimoChannel.__rmul__ = channel_py_channel__mul
MimoChannel.__getitem__ = mimochannel_py__getitem
MimoChannel.__setitem__ = mimochannel_py__setitem
MimoChannel.__len__  = channel_py_channel__len
MimoChannel.flatten  = channel_py_channel__flatten
%}
