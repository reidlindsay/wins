#!  /bin/bash

# simulation parameters
SIMTIME=125.0     # simulation time
# protocol parameters
PACKETLEN=1500    # in bytes
FOMAX=0           # in ppm (default = 13.675)
USEWAVEFORM="--use-waveform"              # enable waveform? (no -> "", yes -> "--use-waveform")
CFOCORRECTION="--disable-cfo-correction"  # enable CFO correction (yes -> "", no -> "--disable-cfo-correction")
# channel parameters                     
TGNMODEL="--tgn-model=D"                  # TGN model type (LOS -> "", TGn CM-A -> "--tgn-model=A")
ALPHA=2.0                                 # pathloss exponent
USEDOPPLER=""                             # enable doppler filter (yes -> "--use-doppler", no -> "")
USEFADING="--disable-fading"              # disable fading (yes -> "--disable-fading", no -> "")
SNRSTEP=0.5
SNRMIN=(-2  4  6  7 13 17 17 20 )
SNRMAX=(15 12 14 24 19 24 33 27 )

# definitions
TESTEXEC=../test.py
#runtests="0 1 2 3 4 5 6 7"
runtests="3 6 "
#runtests="0 "

# check parameters
if [ ! -x $TESTEXEC ]; then
  echo "Cannot find executable $TESTEXEC";
  exit 0;
fi

# print simulation parameters
if [ -f README ]; then
  cat README;
  echo "";
fi

# run tests
for k in $runtests; do
  # simulation parameters
  sargs="-s$SIMTIME";
  # protocol parameters
  mcs=$k;
  plen=$PACKETLEN;
  fomax=$FOMAX;
  pargs="--mcs=$mcs -l$plen --fomax=$fomax $USEWAVEFORM $CFOCORRECTION";
  # channel parameters
  alpha=$ALPHA
  snrmin=${SNRMIN[k]};
  snrmax=${SNRMAX[k]};
  snrstep=$SNRSTEP;
  cargs="--snrmin=$snrmin --snrmax=$snrmax --snrstep=$snrstep --alpha=$alpha $TGNMODEL $USEDOPPLER $USEFADING"
  # run simulation
  simexec="$TESTEXEC $sargs  $pargs $cargs -o test-$mcs.raw"
  echo $simexec
  $simexec
done

exit 0;
