#!/bin/sh

export FREESURFER_HOME=/opt/freesurfer/5.3
export SUBJECTS_DIR=$1

source $FREESURFER_HOME/SetUpFreeSurfer.sh

mksubjdirs $SUBJECTS_DIR/$2

cd $SUBJECTS_DIR
cp -f $2.mgz $SUBJECTS_DIR/$2/mri/
cp -f $2skullstrip.mgz $SUBJECTS_DIR/$2/mri/

cd $SUBJECTS_DIR/$2/mri
mri_convert -c -oc 0 0 0 $2.mgz orig.mgz
mri_convert -c -oc 0 0 0 $2skullstrip.mgz brainmask.mgz

recon-all -talairach -subjid $2
recon-all -nuintensitycor -subjid $2
recon-all -normalization -subjid $2

# part of autorecon2
recon-all -gcareg -subjid $2
recon-all -canorm -subjid $2
recon-all -careg -subjid $2
recon-all -careginv -subjid $2
recon-all -calabel -subjid $2
recon-all -normalization2 -subjid $2
recon-all -maskbfs -subjid $2
recon-all -segmentation -subjid $2
