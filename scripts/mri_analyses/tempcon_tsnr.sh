'''
Purpose: calculate TSNR in entorhinal subregions in the tempcon dataset (Dimsdale-Zucker et al., submitted 2019)

Usage (loop): for subj in [subj_list]; do ./tempcon_tsnr.sh $subj; done

N.B.: directory structure is set up to run on a local machine, and needs to be edited

Written by: Z. Reagh 08/30/2019
'''

# Take subject ID as input from terminal
SUBJ = $1

# Set the subject, data, and mask directories
TOPDIR = /Users/zreagh/Desktop/tempcon
SUBDIR = $TOPDIR/$SUBJ
MASKDIR = $SUBDIR/ROIs/
DATADIR = $SUBDIR/data/

# Navigate to the subject's directory
cd $SUBDIR

# Calculate TSNR map across rest run
3dTstat -tsnr -prefix rest_tsnr $DATADIR/rest.nii
# Calculate average tsnr for alEC L & output into a text file
3dmaskave -mask $MASKDIR/$SUBJ_tempcon_t2_segmented.nii.gz -mrange 21 21 $DATADIR/rest_tsnr >> $TOPDIR/alEC_L_tsnr.txt
# Calculate average tsnr for alEC R & output into a text file
3dmaskave -mask $MASKDIR/$SUBJ_tempcon_t2_segmented.nii.gz -mrange 22 22 $DATADIR/rest_tsnr >> $TOPDIR/alEC_R_tsnr.txt
# Calculate average tsnr for pmEC L & output into a text file
3dmaskave -mask $MASKDIR/$SUBJ_tempcon_t2_segmented.nii.gz -mrange 23 23 $DATADIR/rest_tsnr >> $TOPDIR/pmEC_L_tsnr.txt
# Calculate average tsnr for pmEC R & output into a text file
3dmaskave -mask $MASKDIR/$SUBJ_tempcon_t2_segmented.nii.gz -mrange 24 24 $DATADIR/rest_tsnr >> $TOPDIR/pmEC_R_tsnr.txt