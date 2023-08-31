#!/bin/bash

################################################################################

# Script: pp_anat_faces.sh 

# Author: Nate Cadicamo, Dan McCalley, August 2023

# Description: shell script for preprocessing faces data with AFNI.

################################################################################

# define home directory as AFNI_proc folder
home="/Users/mccalled/Desktop/Nate_data/AFNI_preproc"

# define data directory; change based on structure of stored data
# this file has data in the form of .../{group type}/{subject id}
dataDir="/Users/mccalled/Desktop/Nate_data/preprocessing"

# define group list
groups=('Abstinent' 'Relapse' 'Controls')

# define subjects
Asubs=('P039' 'P037' 'P019' 'P004' 'P002' 'B083' 'B068' 'B047' 'B031' 'B017')
Rsubs=('P014' 'P008' 'B070' 'B066' 'B045' 'B044' 'B018' 'B012' 'B011' 'B002')
Csubs=('C011' 'C010' 'C009' 'C007' 'C006' 'C005' 'C004' 'C003' 'C002' 'C001')

# define templates
t1_template="${home}/templates/MNI152_T1_1mm_brain.nii"
func_template="${home}/templates/MNI152_3mm_brain.nii.gz"

################################################################################

# Part 1: preprocess the anatomical

# define anatatomical preprocessing function
pp_anat() {

  # define passed-in variables from nested loop
  local group=$1
  local sub=$2

  # define subject input & output directories
  inDir="${dataDir}/${group}/${sub}"
  outDir="${home}/${group}/${sub}"

  # make outDir & cd to it
  mkdir $outDir
  cd $outDir

  # also make a "xfs" directory to house all xform files
  mkdir xfs

  # remove skull from t1 anatomical data
  3dSkullStrip -prefix ${sub}_t1_ns.nii.gz \
    -input $inDir/${sub}_anat.nii

  # estimate transform to put t1 in MNI152 space
  @auto_tlrc -no_ss -ok_notice -base TT_avg152T1+tlrc \
    -suffix _afni \
    -input ${sub}_t1_ns.nii.gz

  # the @auto_tlrc command produces a bunch of extra files; clean them up
  gzip ${sub}_t1_ns_afni.nii;
  mv ${sub}_t1_ns_afni.nii.gz ${sub}_t1_MNI.nii.gz;
  mv ${sub}_t1_ns_afni.Xat.1D xfs/t1_to_MNI_xform_afni;
  mv ${sub}_t1_ns_afni.nii_WarpDrive.log xfs/t1_to_MNI_xform_afni.log;
  rm ${sub}_t1_ns_afni.nii.Xaff12.1D

  # take first volume of raw functional data:
  3dTcat -output $inDir/${sub}_vol1_faces.nii.gz \
    $inDir/${sub}_faces.nii[0]

  # skull-strip functional vol
  3dSkullStrip -prefix ${sub}_vol1_faces_ns.nii.gz \
    -input $inDir/${sub}_vol1_faces.nii.gz

  # align the epi to anatomy
  align_epi_anat.py -epi2anat \
    -epi ${sub}_vol1_faces_ns.nii.gz \
    -anat ${sub}_t1_ns.nii.gz \
    -epi_base 0 \
    -tlrc_apar TT_avg152T1+tlrc \
    -epi_strip None \
    -anat_has_skull no

  # put in nifti format
  3dAFNItoNIFTI -prefix ${sub}_vol1_faces_MNI_afni.nii.gz \
    ${outDir}/${sub}_vol1_faces_ns_al+tlrc

  # remove redundant BRIK and HEAD files, since we have nifti files
  rm ${sub}_vol1_faces_ns_tlrc_al+tlrc*

  # clean up xform files
  mv ${sub}_t1_ns_al*aff12.1D xfs/t1_to_faces_xform_afni
  mv ${sub}_vol1_faces_ns_al_mat.aff12.1D xfs/faces_to_t1_xform_afni
  rm vol1_cue_ns_al_reg_mat.aff12.1D

}


################################################################################

# Part 2: preprocess the fMRI task

# threshold for determining which volumes should be censored for "bad" motion;
# unit roughly corresponds to mm; see "1d_tool.py help" for more info
censor_thresh=6

# filepaths to ROI masks
wmMaskFile="${home}/templates/MNI_3mm_WM_mask.nii.gz"
csfMaskFile="${home}/templates/MNI_3mm_CSF_mask.nii.gz"



# define function to preprocess faces task (can do other tasks with other fns)
pp_faces() {

  # define passed-in variables from nested loop
  local group=$1
  local sub=$2

  # define subject input & output directories
  inDir="${dataDir}/${group}/${sub}"
  outDir="${home}/${group}/${sub}"

  # cd into outDir
  cd $outDir

  # drop first vol, allow longitudinal magentization (t1) to reach steady state
  3dTcat -output faces1.nii.gz $inDir/${sub}_faces.nii[1..$]

  # correct for slice time differences
  3dTshift -prefix t_faces1.nii.gz -slice 0 -tpattern altplus faces1.nii.gz

  # motion correction & saves out the motion parameters in file, 'faces_vr.1D'
  3dvolreg -Fourier -twopass -zpad 4 \
    -dfile faces_vr.1D \
    -base 0 \
    -prefix m_t_faces1.nii.gz \
    t_faces1.nii.gz

  # create a “censor vector” that denotes bad movement volumes with a 0 and good
  # volumes with a 1 to be used later for glm estimation and making timecourses
  1d_tool.py -infile faces_vr.1D[1..6] \
    -show_censor_count \
    -censor_prev_TR \
    -censor_motion $censor_thresh \
    faces
  rm faces_CENSORTR.txt

  # transform functional time series to MNI space
  3dAllineate -base ${sub}_t1_MNI.nii.gz \
    -1Dmatrix_apply xfs/t1_to_MNI_xform_afni \
    -prefix w_m_t_faces_MNI_afni \
    -input m_t_faces1.nii.gz \
    -verb \
    -master BASE \
    -mast_dxyz 3.0 \
    -VERB \
    -warp aff \
    -twopass

  # convert to nifti
  3dAFNItoNIFTI -prefix w_m_t_faces_MNI_afni.nii.gz \
    ${outDir}/w_m_t_faces_MNI_afni+tlrc

  # smooth data with an 8 mm full width half max gaussian kernel
  3dmerge -1blur_fwhm 8 -doall -quiet \
    -prefix s_w_m_t_faces_MNI.nii.gz w_m_t_faces_MNI_afni.nii.gz

  # calculate the mean timeseries for each voxel
  3dTstat -mean -prefix mean_faces1.nii.gz s_w_m_t_faces_MNI.nii.gz

  # convert voxel values to be percent signal change
  3dcalc -a s_w_m_t_faces_MNI.nii.gz -b mean_faces1.nii.gz \
  -expr '((a-b)/b)*100' -prefix ps_s_w_m_t_faces_MNI.nii.gz -datum float

  # high-pass filter the data
  3dBandpass -prefix f_ps_s_w_m_t_faces_MNI.nii.gz \
    -ftop 0.011 \
    ps_s_w_m_t_faces_MNI.nii.gz

  # dump out WM and CSF  time series into separate files
  3dmaskave -mask $csfMaskFile -quiet -mrange 1 2 \
    ps_s_w_m_t_faces_MNI.nii.gz > faces_csf_ts.txt
  3dmaskave -mask $wmMaskFile -quiet -mrange 1 2 \
    ps_s_w_m_t_faces_MNI.nii.gz > faces_wm_ts.txt

  # remove redundant BRIK and HEAD files
  rm *.BRIK *.HEAD

}

################################################################################

# This section is the actual nested loop to run our preprocessing functions.
# Adapt this, and relevant "group" and "subject" variables, to other needs.

# iterate through groups and subjects in dataDir, preprocess
for group in ${groups[@]}; do

  # first group: abstinent
  if [ $group = 'Abstinent' ]; then

    # inner loop: iterate over abstinent subjects
    for sub in ${Asubs[@]}; do

      # call our functions
      pp_anat $group $sub
      pp_faces $group $sub

    done # move on to relapse group

  # second group: relapse
  elif [ $group = 'Relapse' ]; then

    # inner loop: iterate over relapse subjects
    for sub in ${Rsubs[@]}; do

      # call our functions
      pp_anat $group $sub
      pp_faces $group $sub

    done # move on to controls group

  # third group: controls
  else

    # inner loop: iterate over control subjects
    for sub in ${Csubs[@]}; do

      # call our functions
      pp_anat $group $sub
      pp_faces $group $sub

    done # all done

  fi

done
