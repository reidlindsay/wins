#!  /bin/bash

# parse tput data from files files
NUMTESTS=10
NCOLLISION=3
PARSETPUT=../parse-tput.py

if [ ! -x $PARSETPUT ]; then
  echo "Cannot find $PARSETPUT";
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
  # check if files in directory exist
  if [ -n "`ls data.$testnum/aloha-*.raw`" ]; then
    # call parse-tput.py
    echo "Parsing TPUT data from data.$testnum for all ...";
    $PARSETPUT data.$testnum/aloha-*.raw > data.$testnum/tput.$testnum.all.dat
    echo "Parsing TPUT data from data.$testnum for c = 0 ...";
    $PARSETPUT -c 0 data.$testnum/aloha-*.raw > data.$testnum/tput.$testnum.c0.dat
  fi
done

exit 0;
