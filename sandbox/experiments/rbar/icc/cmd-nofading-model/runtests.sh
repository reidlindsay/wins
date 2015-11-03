#!  /bin/bash

# simulation parameters
NFRAMES=50000
SIMTIME="$NFRAMES*0.001876" # simulation time (function of MCS)
# protocol parameters
MCS=0               # rate index
PACKETLEN=1024      # in bytes
FOMAX=0.0           # in ppm (default = 13.675)
USEWAVEFORM=""                            # enable waveform? (no -> "", yes -> "--use-waveform")
CFOCORRECTION="--disable-cfo-correction"  # enable CFO correction (yes -> "", no -> "--disable-cfo-correction")
# channel parameters
ALPHA=2.0                           # pathloss exponent
#TGNMODEL=""                         # TGN model type (LOS -> "", TGn CM-A -> "--tgn-model=A")
USEDOPPLER=""                       # enable doppler filter (yes -> "--use-doppler", no -> "")
USEFADING="--disable-fading"        # use fading (yes -> "", no -> "--disable-fading")
TGNMODEL="--tgn-model=D"            # TGN model type (LOS -> "", TGn CM-A -> "--tgn-model=A")
#USEDOPPLER="--use-doppler"          # enable doppler filter (yes -> "--use-doppler", no -> "")
#USEFADING=""                        # use fading (yes -> "", no -> "--disable-fading")

# SNR parameters
SNRSTEP=1.0
SNRMIN=( 4  1  4  6  7 13 17 17 20 )
SNRMAX=(28 29 12 14 17 19 24 27 27 )

# definitions
RBAREXEC=../test.py
runtests="1 "

# check parameters
if [ ! -x $RBAREXEC ]; then
  echo "Cannot find executable $RBAREXEC";
  exit 0;
fi

# print simulation parameters
if [ -f README ]; then
  cat README;
  echo "";
fi

# run tests
for k in $runtests; do
  # protocol parameters
  mcs=$MCS
  plen=$PACKETLEN;
  fomax=$FOMAX;
  pargs="--mcs=$mcs -l$plen --fomax=$fomax $USEWAVEFORM $CFOCORRECTION";
  # channel parameters
  alpha=$ALPHA
  fargs="--alpha=$alpha $TGNMODEL $USEDOPPLER $USEFADING"
  # channel parameters
  alpha=$ALPHA
  snrmin=${SNRMIN[k]};
  snrmax=${SNRMAX[k]};
  snrstep=$SNRSTEP;
  snrargs="--snrmin=$snrmin --snrmax=$snrmax --snrstep=$snrstep"
  cargs="$snrargs --alpha=$alpha $TGNMODEL $USEDOPPLER $USEFADING"
  # simulation parameters
  simtime=`echo "scale=8; $SIMTIME*($snrmax-$snrmin)/$snrstep"|bc`
  sargs="-s$simtime"
  # run simulation
  simexec="$RBAREXEC $sargs $pargs $cargs -o rbar-model-$k.raw"
  echo $simexec
  $simexec
done

exit 0;
