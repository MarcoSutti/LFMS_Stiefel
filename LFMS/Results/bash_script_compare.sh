#!/bin/bash
# NOTE : Quote it else use array to avoid problems #
# MS, 2023.06.24. Created
# MS, 2023.06.24. Created
FILES="*.mat"
for f in $FILES
do
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  # cat "$f"
  colordiff $f ../../../Dropbox/00_Scientific_Research/leapfrogstiefel/matlab/LFMS/Results/$f
done
