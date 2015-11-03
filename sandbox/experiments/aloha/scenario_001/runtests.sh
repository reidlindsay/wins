#!  /bin/bash

# simulation parameters
NPKTS="200"
SIMTIME="$NPKTS*0.001876"   # simulation time (function of MCS)
TOPOFILE=topo.dat           # topology file
SAVETOPO="--save-topo=$TOPOFILE"    # use topology file? (save by default)
# protocol parameters
MCS=7               # rate index
PACKETLEN=1500      # in bytes
FOMAX=0.0           # in ppm (default = 13.675)
USEWAVEFORM="--use-waveform"              # enable waveform? (no -> "", yes -> "--use-waveform")
CFOCORRECTION="--disable-cfo-correction"  # enable CFO correction (yes -> "", no -> "--disable-cfo-correction")
# channel parameters
TGNMODEL="--tgn-model=A"            # TGN model type (LOS -> "", TGn CM-A -> "--tgn-model=A")
ALPHA=2.0                           # pathloss exponent
USEDOPPLER=""                       # enable doppler filter (yes -> "--use-doppler", no -> "")
USEFADING="--disable-fading"        # disable fading (yes -> "--disable-fading", no -> "")

# program definitions
TIME=`which time` -p

# definitions
ALOHAEXEC=../test.py
#runtests="0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.5 1.8 "
runtests="0.2 "

# check parameters
if [ ! -x $ALOHAEXEC ]; then
  echo "Cannot find executable $ALOHAEXEC";
  exit 0;
fi

# print simulation parameters
if [ -f README ]; then
  cat README;
  echo "";
fi

# run tests
for d in $runtests; do
  # simulation parameters
  mcs=$MCS
  load=$d;
  simtime=`echo "scale=8; $SIMTIME/($mcs + 1)"|bc`
  if [ -f $TOPOFILE ]; then
    usetopo="--use-topo=$TOPOFILE";
  else
    usetopo=$SAVETOPO;
  fi
  sargs="-s$simtime -G $load $usetopo"
  # protocol parameters
  plen=$PACKETLEN;
  fomax=$FOMAX;
  pargs="--mcs=$mcs -l$plen --fomax=$fomax $USEWAVEFORM $CFOCORRECTION";
  # channel parameters
  alpha=$ALPHA
  cargs="--alpha=$alpha $TGNMODEL $USEDOPPLER $USEFADING"
  # run simulation
  simexec="$ALOHAEXEC $sargs $pargs $cargs -o aloha-$load.raw"
  echo $simexec
  $simexec
done

exit 0;
