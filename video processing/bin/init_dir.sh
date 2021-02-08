#! /usr/bin/bash

mkdir -p box aligned leapout
chmod 0777 box aligned leapout

# gets the PWD folder name as header
export HEADER=`pwd |awk -F / '{print $NF}'`
echo "The current directory to run: /$HEADER"

# the list of .h5
ls *.h5 > h5List.txt
# list of raw .h5 videos
sed "s@^@${PWD}/@" < h5List.txt > ${HEADER}_List.txt
sed "s@${HEADER}/@${HEADER}/box/@ ;s@.h5\$@_box.h5@" < ${HEADER}_List.txt > BOX_List.txt
sed 's@_box.h5$@_box_PREDICTED.h5@' <BOX_List.txt > PRED_List.txt
sed 's@box/@aligned/@ ;s@_box.h5$@_box_aligned.h5@' <BOX_List.txt > ALIGNED_List.txt
sed 's@_aligned.h5$@_aligned_PREDICTED_2.h5@ ;s@/aligned/@/leapout/@' <ALIGNED_List.txt > PRED_2_List.txt

genSbatchFiles.sh
cp sbatchTemplates/plot* ./
