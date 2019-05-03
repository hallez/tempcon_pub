#!/bin/bash
set -eux
SUBJ=${SUBJ:?}
ASHS_OUT_DIR=<absolute-path>/tempcon_pub/analyzed_mri/$SUBJ/ashs
MPRAGE_DIR=<absolute-path>/tempcon_pub/analyzed_mri/$SUBJ/mprage
T2_DIR=<absolute-path>/tempcon_pub/analyzed_mri/$SUBJ/t2_19
ASHS_ATLASES=$ASHS_ROOT/atlases/final

if [ ! -d "$ASHS_OUT_DIR" ]; then
    mkdir "$ASHS_OUT_DIR"
fi

nohup $ASHS_ROOT/bin/ashs_main.sh \
    -I "$SUBJ" \
    -a $ASHS_ATLASES \
    -g $MPRAGE_DIR/s*.nii \
    -f $T2_DIR/s*.nii \
    -w $ASHS_OUT_DIR \
    &> $ASHS_OUT_DIR/nohup.out &
