#!/bin/bash

################################################################################

# Script: general_afni_preproc.sh

# Author: Nate Cadicamo, August 2023

# Description: shell script for preprocessing fMRI data with AFNI. Designed to
# be generalizable across different fMRI data tasks. Deployable for faces as it
# stands, ready to be adapted to cue, go-no-go, mid, and whatever else you need!

# How to use: this script is generic, and it won't run as it is. It's just a
# template! So, to use it, copy it into a new shell file and modify it as you
# need. Anywhere you see 'USER', you need to provide some specification or
# modification that depends on your given goals and computer setup. Since it is
# built to be adaptable, it should cover all or nearly all of your preprocessing
# needs! See AFNI functions online for more detail. Have fun preprocessing!

# To call this script, go to your terminal command line and cd into the
# directory where this script is stored. Call it as follows, where the {task}
# variable is to be filled in by you, USER, with something like faces.
# Command: bash general_afni_preproc.sh {task}
# example command line:
  # (base) computer3 Nate_data % bash general_afni_preproc.sh faces

# Also, thanks to Kelly (Hennigan) MacNiven, whose scripts I referenced as
# examples, and to Dan McCalley, who helped develop the scripts on which this
# template is based.

################################################################################

# USER: specify your home directory; see code below for functionality
home=""

# USER: specify your data directory (where you have raw fMRI scans)
dataDir=""

# USER: define group list, if applicable. eg: ('Abstinent' 'Relapse' 'Control')
groups=('')

# USER: define subjects. eg: ('P088' 'P089' 'P090' 'P094')
subjects=('')

# USER: define templates, such as MNI152 1mm for t1 and 3mm for functional
t1_template=""
func_template=""

# USER: from the command line, when you call general_afni_preproc.sh, specify
# which task you want this script to process. Must match options in nested loop!
# eg: ../ % bash general_afni_preproc faces
user_input_task=$1

################################################################################

# Part 1: preprocess the anatomical

# USER: change file paths, directories, and naming as necessary


# define an anatatomical preprocessing function
pp_anat() {

  # USER: define passed-in variables from nested loop
  local group=$1 # may be unneccessary for your needs
  local sub=$2

  # USER: define subject input & output directories. eg: ${home}/${group}/${sub}
  inDir="" # probably a subdirectory of your ${dataDir}
  outDir="" # probably a subdirectory of your ${home}

  # make outDir & cd to it
  mkdir $outDir
  cd $outDir

  # also make a "xfs" directory within ${outDir} to house all xform files
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

# USER: specify censor threshold for determining which volumes should be
# censored for "bad" motion; unit roughly corresponds to mm; see "1d_tool.py
# help" for more info
censor_thresh=6

# USER: define filepaths to ROI masks, such as MNI_3mm_WM_mask.nii.gz
wmMaskFile=""
csfMaskFile=""


# define function to preprocess faces
pp_faces() {

  # USER: define passed-in variables from nested loop
  local group=$1 # may be unneccessary for your needs
  local sub=$2

  # USER: define subject input & output directories
  inDir="" # probably a subdirectory of your ${dataDir}
  outDir="" # probably a subdirectory of your ${home}

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

  # USER: at this point, your s_w_m_t_faces_MNI.nii.gz file is smoothed,
  # normalized, and realligned. It is functionally ready for further analysis.
  # In the next few short commands, we get percent signal change, in case that
  # is helpful, as well as a mean timeseries for each voxel.

  # calculate the mean timeseries for each voxel
  3dTstat -mean -prefix mean_faces1.nii.gz s_w_m_t_faces_MNI.nii.gz

  # convert voxel values to be percent signal change
  3dcalc -a s_w_m_t_faces_MNI.nii.gz -b mean_faces1.nii.gz \
  -expr '((a-b)/b)*100' -prefix ps_s_w_m_t_faces_MNI.nii.gz -datum float

  # USER: high-pass filter the data. Optional, AFNI recommends 3dTproject
  3dBandpass -prefix f_ps_s_w_m_t_faces_MNI.nii.gz \
    -ftop 0.011 \
    ps_s_w_m_t_faces_MNI.nii.gz

  # dump out WM and CSF  time series into separate files
  # USER: note that you can use ps_s.* or just s_.* file as input here
  3dmaskave -mask $csfMaskFile -quiet -mrange 1 2 \
    s_w_m_t_faces_MNI.nii.gz > faces_csf_ts.txt
  3dmaskave -mask $wmMaskFile -quiet -mrange 1 2 \
    s_w_m_t_faces_MNI.nii.gz > faces_wm_ts.txt

  # remove redundant BRIK and HEAD files
  rm *.BRIK *.HEAD

}


# define function to preprocess cue
pp_cue() {

  # USER: fill this in as necessary! You can probably recycle almost everything
  # from pp_faces(), just changing some filenames

}


# define other functions for other tasks here!
pp_USER_SPECIFIED_TASK() {

}

################################################################################

# Part 3: bring it all together

# This section is the actual nested loop to run our preprocessing functions.
# Adapt this, and relevant "group" and "subject" variables, to other needs.

# The first inner loop is almost done for you. If you don't need groups, and
# you just have a list of subjects, copy and paste just the first 'for sub in'
# logic through its corresponding 'done', remove the group variable both there
# and in the corresponding functions, and adapt as you see fit.

# The second and third inner loop are almost totally empty, and can be more
# easily filled in by you, the USER, once you flesh out the first inner loop.
# You can add as many loops and nestings as you need.


# iterate through groups and subjects in groups, preprocess with our functions
for group in ${groups[@]}; do

  # USER: specify first group in '' below
  if [ $group = '' ]; then

    # inner loop: iterate over group1 subjects
    # USER: specify subject group in '' below
    for sub in ${''[@]}; do

      # to keep our sanity, note where we are
      echo -e "\ncurrently working on $sub... \n"

      # call our functions. first, we definitely want to preprocess anatomical
      pp_anat $group $sub

      # second, we want to preprocess the functional, specified by your given
      # user_input_task
      # USER: adapt this as necessary!

      if [ $user_input_task = 'faces' ]; then
        # do faces, as defined above
        pp_faces $group $sub


      elif [ $user_input_task = 'cue' ]; then
        # do cue, as defined above
        pp_cue $group $sub


      elif [ $user_input_task = '' ]; then
        # USER: specify additional functions


      else
        # USER: specify additional functions


      fi # done with this subject's functional preprocessing

    done # move on to next group


  # USER: specify next group, if necessary. copy and paste from above
  elif [ $group = '' ]; then

      # [for sub in $group logic here]

  # USER: specify last group, if necessary. copy and paste from above
  else

    # [for sub in $group logic here]


  fi # done with all the possible groups, no more comparisons to check

done # done with all the processing
