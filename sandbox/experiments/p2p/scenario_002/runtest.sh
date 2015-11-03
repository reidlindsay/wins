#!  /bin/bash

# simulation parameters
SIMTIME=8.000     # simulation time
# protocol parameters
PACKETLEN=1500    # in bytes
FOMAX=0           # in ppm (default = 13.675)
USEWAVEFORM="--use-waveform"              # enable waveform? (no -> "", yes -> "--use-waveform")
CFOCORRECTION="--disable-cfo-correction"  # enable CFO correction (yes -> "", no -> "--disable-cfo-correction")
# channel parameters                     
TGNMODEL="--tgn-model=F"                               # TGN model type (LOS -> "", TGn CM-A -> "--tgn-model=A")
ALPHA=2.0                                 # pathloss exponent
USEDOPPLER=""                             # enable doppler filter (yes -> "--use-doppler", no -> "")
USEFADING="--disable-fading"              # disable fading (yes -> "--disable-fading", no -> "")
SNRSTEP=1.0
SNRMIN=( 1  4  6 10 13 17 19 20 )
#SNRMAX=( 5 12 14 17 19 24 25 27 )
SNRMAX=(15 12 14 20 19 24 30 27 )

# definitions
TESTEXEC=../test.py
#runtests="0 1 2 3 4 5 6 7"
runtests="0 3 6 "

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
