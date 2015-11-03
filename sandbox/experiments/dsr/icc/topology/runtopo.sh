#!  /bin/bash

##############################################################################
#
#   Create multiple topologies using specified parameters.
#   
#   Revision Info
#   =============
#   * $LastChangedBy: mandke $
#   * $LastChangedDate: 2011-10-01 12:07:30 -0500 (Sat, 01 Oct 2011) $
#   * $LastChangedRevision: 5171 $
#   
#   :author: Ketan Mandke <kmandke@mail.utexas.edu>
#   
#   :copyright:
#       Copyright 2009-2011 The University of Texas at Austin
#   
#       Licensed under the Apache License, Version 2.0 (the "License");
#       you may not use this file except in compliance with the License.
#       You may obtain a copy of the License at
#   
#          http://www.apache.org/licenses/LICENSE-2.0
#   
#       Unless required by applicable law or agreed to in writing, software
#       distributed under the License is distributed on an "AS IS" BASIS,
#       WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#       See the License for the specific language governing permissions and
#       limitations under the License.
#
##############################################################################

# simulation parameters
VERBOSE=54
# experiment parameters
NUMNODES=50
SIMTIME=0.0                 # simulation time
NCONNECT=1                  # number of connections

# topology parameters
BORDER=(0 1500.0  0  300.0)
#TOPOFILE=topo.dat           # topology file
#SAVETOPO="--save-topo=$TOPOFILE"    # save topology file? (no -> "", yes -> "--save-topo=$TOPOFILE")

# agent parameters
PACKETLEN=1024      # in bytes
RATE=4.0            # in packets/second
AGTMODE="CBR"       # traffic generation mode

# DSR parameters
#ROUTEFILE=route.dat         # routing file
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
#runtests="1 "
# Create multiple topologies
runtopos="`seq 1 5`"

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
for k in $runtopos; do
  # simulation parameters
  verbose="-v $VERBOSE"
  simtime="-s $SIMTIME"
  simargs="$simtime $verbose"

  # experiment parameters
  numnodes="-n $NUMNODES"
  nconnect="-c $NCONNECT"
  expargs="$numnodes $nconnect"

  # agent parameters
  rate="-r $RATE"
  plen="-l $PACKETLEN"
  mode="--agent-mode $AGTMODE"
  agtargs="$rate $plen $mode"

  # network parameters
  netargs="";

  # PHY parameters
  mcs="--mcs=$MCS"
  fomax="--fomax=$FOMAX"
  usewaveform=$USEWAVEFORM
  cfocorrection=$CFOCORRECTION
  phyargs="$mcs $fomax $usewaveform $cfocorrection";

  # topology parameters
  xborder="--xmin=${BORDER[0]} --xmax=${BORDER[1]}"
  yborder="--ymin=${BORDER[2]} --ymax=${BORDER[3]}"
  savetopo="--save-topo=topo$k.dat";
  border="$xborder $yborder"
  topoargs="$savetopo $border"

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
