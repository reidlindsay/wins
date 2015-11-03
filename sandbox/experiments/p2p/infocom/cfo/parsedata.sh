#!  /bin/bash

# parse trace files to get FDR and PER data
PARSEFDR=../parse-fdr.py
PARSEPER=../parse-per.py

if [ ! -x $PARSEFDR ]; then
  echo "Cannot find $PARSEFDR";
  exit 0;
fi

if [ ! -x $PARSEPER ]; then
  echo "Cannot find $PARSEPER";
  exit 0;
fi

# check if trace files are in current directory
if [ -n "`ls test-*.raw`" ]; then
  # call parse-fdr.py
  echo "Parse FDR data from test-*.raw ... ";
  $PARSEFDR test-*.raw > fdr.dat
  # call parse-per.py
  echo "Parse PER data from test-*.raw ... ";
  $PARSEPER test-*.raw > per.dat
fi

exit 0;
