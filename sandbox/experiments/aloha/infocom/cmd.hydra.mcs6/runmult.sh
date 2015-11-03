#!  /bin/bash

# call runtests.sh multiple times and organize results into folders
NUMTESTS=10
RUNTESTS=./runtests.sh

if [ ! -x $RUNTESTS ]; then
  echo "Cannot find $RUNTESTS";
  exit 0;
fi

# macros and programs
MKDIR=`which mkdir`
MV=`which mv`

for k in $(seq 1 $NUMTESTS);
do
  # figure out string for test number
  if [ $k -lt 10 ]; then
    testnum="0$k";
  else
    testnum="$k";
  fi
  # call runtests.sh
  echo "Running test for data.$testnum ...";
  $RUNTESTS
  # create data directories
  if [ ! -d "data.$testnum" ]; then
    $MKDIR data.$testnum;
  fi
  # move output to appropriate directory
  $MV aloha-*.raw topo.dat data.$testnum/
done

exit 0;
