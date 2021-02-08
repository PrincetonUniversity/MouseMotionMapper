#! /usr/bin/bash

# generate sbatch files on the current directory
# make sure sbatchTemplate are in the currrent folder
export NUM_ARR=`ls *.h5 |wc -w`
echo "Number of h5 files to process: $NUM_ARR"

sed "s/SUB1/${NUM_ARR}/; \
s@SUB2@${PWD}/${HEADER}_List.txt@; \
s@SUB3@${PWD}/box@" \
< sbatchTemplates/1_runBox.sbTemp > 1_runBox.sb

sed "s/SUB1/${NUM_ARR}/; \
s@SUB2@${PWD}/BOX_List.txt@; \
s@SUB3@${PWD}/PRED_List.txt@" \
< sbatchTemplates/2_predAlignMouse.sbTemp > 2_predAlignMouse.sb

sed "s/SUB1/${NUM_ARR}/; \
s@SUB2@${PWD}/BOX_List.txt@; \
s@SUB3@${PWD}/aligned/@" \
< sbatchTemplates/3_runMouseAlign.sbTemp > 3_runMouseAlign.sb

sed "s/SUB1/${NUM_ARR}/; \
s@SUB2@${PWD}/ALIGNED_List.txt@; \
s@SUB3@${PWD}/PRED_2_List.txt@" \
< sbatchTemplates/4_predMultMouse.sbTemp > 4_predMultMouse.sb
