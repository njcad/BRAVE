#!/bin/bash

##############################################################################

# Script: mask_CSF_WM_rem_vol1

# Author: Nate Cadicamo, August 2023

# Description: cleaning fMRI data for BRAVE Lab / Stanford BioX project.
# Preprocessed data in SPM, now using AFNI functions to clean as follows.
# For each subject in each group in Nate_data,

  # 1) 3dTcat remove first volume and write new file into ../afni_tc_swr.nii

  # 2) 3dmaskave for CSF and WM, save as ../faces_{wm or csf}_ts.txt

  # 3) remove first line of rp_{subject id}_faces.txt, write new file into
  # rp_{subject id}_faces_first_line_removed.txt

##############################################################################

# define locations for white matter and CSF masks
maskDir="/Users/mccalled/Desktop/Nate_data/afni_proc_test/templates"
wmMaskFile="$maskDir/MNI_152_0.5_WM_bin.nii"
csfMaskFile="$maskDir/MNI_152_0.5_CSF_bin.nii.gz"

# define group list
groups=('Abstinent' 'Relapse' 'Controls')

# define subjects
Asubs=('P039' 'P037' 'P019' 'P004' 'P002' 'B083' 'B068' 'B047' 'B031' 'B017')
Rsubs=('P014' 'P008' 'B070' 'B066' 'B045' 'B044' 'B018' 'B012' 'B011' 'B002')
Csubs=('C011' 'C010' 'C009' 'C007' 'C006' 'C005' 'C004' 'C003' 'C002' 'C001')

##############################################################################

# cd into directory
cd "/Users/mccalled/Desktop/Nate_data/preprocessing/"

##############################################################################

# define function to do each of our three tasks
rearrange() {

  # get variables from our nested loops
  local group=$1
  local sub=$2

  # 1) remove first volume with 3dTcat
  3dTcat -output "${group}/${sub}/afni_tc_swr.nii" \
    "${group}/${sub}/swr${sub}_faces.nii[1..$]"

  # 2) save mask time series for WM and CSF
  3dmaskave -mask "$wmMaskFile" -quiet -mrange 1 2 \
    "${group}/${sub}/afni_tc_swr.nii" > "${group}/${sub}/faces_wm_ts.txt"
  3dmaskave -mask "$csfMaskFile" -quiet -mrange 1 2 \
    "${group}/${sub}/afni_tc_swr.nii" > "${group}/${sub}/faces_csf_ts.txt"

  # 3) remove first line of rp text file
  original="${group}/${sub}/rp_${sub}_faces.txt"
  line_removed="${group}/${sub}/rp_${sub}_faces_first_line_removed.txt"
  sed '1d' "$original" > "$line_removed"

}

##############################################################################

# outer loop: iterate over groups
for group in ${groups[@]}; do

  # first group: abstinent
  if [ $group = 'Abstinent' ]; then

    # inner loop: iterate over abstinent subjects
    for sub in ${Asubs[@]}; do
      # call our function
      rearrange $group $sub
    done # move on to relapse group

  # second group: relapse
  elif [ $group = 'Relapse' ]; then

    # inner loop: iterate over relapse subjects
    for sub in ${Rsubs[@]}; do
      # call our function
      rearrange $group $sub
    done # move on to controls group

  # third group: controls
  else

    # inner loop: iterate over control subjects
    for sub in ${Csubs[@]}; do
      # call our function
      rearrange $group $sub
    done # all done

  fi

done
