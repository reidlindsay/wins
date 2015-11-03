#!  /bin/bash

# simulation parameters
VERBOSE=54
# experiment parameters
NUMNODES=8
SIMTIME=2.5                # simulation time

# topology parameters
#BORDER=(0  100.0  0  100.0)
BORDER=(0  1000.0  0  0.0)
TOPOFILE=topo.dat           # topology file
SAVETOPO="--save-topo=$TOPOFILE"    # save topology file? (no -> "", yes -> "--save-topo=$TOPOFILE")

# agent parameters
PACKETLEN=1024      # in bytes
RATE=1.0            # in packets/second
AGTMODE="CBR"       # traffic generation mode

# DSR parameters
RREQRATE="--rreqrate=3"     # specify RREQ rate (default -> "", otherwise -> "--rreqrate=RATE")
DATARATE="--datarate=3"     # specify data rate (default -> "", otherwise -> "--datarate=RATE")
ROUTEFILE=route.dat         # routing file
#SAVEROUTE="--save-route=$ROUTEFILE" # save routing?  (no -> "", yes -> "--save-route=$ROUTEFILE")

# PHY parameters
MCS=0               # mcs index
FOMAX=0.0           # in ppm (default = 13.675)
USEWAVEFORM="--use-waveform"              # enable waveform? (no -> "", yes -> "--use-waveform")
CFOCORRECTION="--disable-cfo-correction"  # enable CFO correction (yes -> "", no -> "--disable-cfo-correction")

# channel parameters
ALPHA=2.0                               # pathloss exponent
BIDIRECTIONAL="--bidirectional-channel" # bidirectional (yes -> "--bidirectional-channel", no -> "")
#TGNMODEL=""                             # TGN model type (LOS -> "", TGn CM-A -> "--tgn-model=A")
USEDOPPLER=""                           # enable doppler filter (yes -> "--use-doppler", no -> "")
USEFADING="--disable-fading"            # use fading (yes -> "", no -> "--disable-fading")
TGNMODEL="--tgn-model=A"                # TGN model type (LOS -> "", TGn CM-A -> "--tgn-model=A")
#USEDOPPLER="--use-doppler"              # enable doppler filter (yes -> "--use-doppler", no -> "")
#USEFADING=""                            # use fading (yes -> "", no -> "--disable-fading")

# definitions
DSREXEC=../test.py
# Use number of connections as independent variable
#runtests="1 2 5 10 20"
runtests="1 "

# check parameters
if [ ! -x $DSREXEC ]; then
  echo "Cannot find executable $DSREXEC";
  exit 0;
fi

# print simulation parameters
if [ -f README ]; then
  cat README;
  echo "";
fi

######################################
# Run Tests
######################################
for d in $runtests; do
  # simulation parameters
  verbose="-v $VERBOSE"
  simtime="-s $SIMTIME"
  simargs="$simtime $verbose"

  # experiment parameters
  numnodes="-n $NUMNODES"
  nconnect="-c $d"
  expargs="$numnodes $nconnect"

  # agent parameters
  rate="-r $RATE"
  plen="-l $PACKETLEN"
  mode="--agent-mode $AGTMODE"
  agtargs="$rate $plen $mode"

  # network parameters
  if [ -f $ROUTEFILE ]; then
    useroute="--use-route=$ROUTEFILE";
  fi
  saveroute=$SAVEROUTE;
  rreqrate="$RREQRATE";
  datarate="$DATARATE";
  netargs="$rreqrate $datarate $useroute $saveroute"

  # PHY parameters
  mcs="--mcs=$MCS"
  fomax="--fomax=$FOMAX"
  usewaveform=$USEWAVEFORM
  cfocorrection=$CFOCORRECTION
  phyargs="$mcs $fomax $usewaveform $cfocorrection";

  # topology parameters
  if [ -f $TOPOFILE ]; then
    usetopo="--use-topo=$TOPOFILE";
    border=""
  else
    xborder="--xmin=${BORDER[0]} --xmax=${BORDER[1]}"
    yborder="--ymin=${BORDER[2]} --ymax=${BORDER[3]}"
    usetopo=$SAVETOPO;
    border="$xborder $yborder"
  fi
  topoargs="$usetopo $border"

  # channel parameters
  alpha="--alpha=$ALPHA"
  bidirectional="$BIDIRECTIONAL"
  tgnmodel=$TGNMODEL
  usedoppler=$USEDOPPLER
  usefading=$USEFADING
  chargs="$alpha $bidirectional $tgnmodel $usedoppler $usefading";

  # run simulation
  args="$simargs $expargs $agtargs $netargs $phyargs $topoargs $chargs"
  simexec="$DSREXEC $args -o dsr-n$NUMNODES-c$d.raw"
  echo $simexec
  $simexec
done

exit 0;
