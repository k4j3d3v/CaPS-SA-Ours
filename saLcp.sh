#!/usr/bin/env bash

# compute SA and LCP array for all files

# exit immediately on error
set -e
# file containing the sha1sum of newly computed sa/lcp
sha1file=sha1.out

if [ $# -lt 2 ]
then
  echo "Usage:"
  echo "         $0 exe_with_options file1 [file2 ...]"
  echo
  echo "Compute SA and LCP for all file on the command line"
  echo "reporting running time and peak memory usage"
  echo "compute sha1sum's for all files and write them to $sha1file"
  echo
  echo "Sample usage:"
  echo "         $0 \"Capsa -t4\" corpus/* "      
  exit
fi

exe=$1
shift 1

echo ">>>>>>> deleting $sha1file"
rm -f $sha1file


for f in "$@"
do 
  echo ">>>>>>> sa and lcp computation "
  /usr/bin/time -o $f.stat -f"Command %C\n####### E(secs):%e Mem(kb):%M" $exe $f -w $f.sa -W $f.lcp  2>/dev/null
  cat $f.stat
  # compute sha1sum and delete (after the initial sha1 computation)
  sha1sum  $f $f.lcp $f.sa >> $sha1file 
  rm -f $f.sa $f.lcp $f.stat
done
