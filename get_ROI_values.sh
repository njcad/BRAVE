#!/bin/bash

################################################################################

# Script: get_ROI_values

# Author: Nate Cadicamo, August 2023

# Description: pulling activation level data from fMRI data for each subject in
# each group (abstinent, relapse, control) for each region of interest (ROI).
# Written for BRAVE Lab / BioX 2023 project.

################################################################################

# define the directory containing our regions of interest
ROI_Dir="/Users/mccalled/Desktop/Nate_data/ROI_Corrected"

# create an empty array to store file paths to these ROIs
ROIs=()

# loop through the NIfTI files in the directory and add their paths to the array
for file in "$ROI_Dir"/*.nii.gz; do
    if [ -f "$file" ]; then
        ROIs+=("$file")
    fi
done


# define group list
groups=('Abstinent' 'Relapse' 'Controls')


# define subjects
Asubs=('P039' 'P037' 'P019' 'P004' 'P002' 'B083' 'B068' 'B047' 'B031' 'B017')
Rsubs=('P014' 'P008' 'B070' 'B066' 'B045' 'B044' 'B018' 'B012' 'B011' 'B002')
Csubs=('C011' 'C010' 'C009' 'C007' 'C006' 'C005' 'C004' 'C003' 'C002' 'C001')

################################################################################

# define function to use in our loop
get_vals() {

  # define variable as passed in through loop below
  local region=$1
  local group=$2
  local sub=$3

  # define directory location for ease of access
  dir="/Users/mccalled/Desktop/Nate_data"

  # use AFNI 3dROIstats to get activation levels in given region
  3dROIstats \
    -mask "${region}" \
    -nzmean "${dir}/first_level_corrected/${group}/${sub}/con_0018.nii" \
    >> "${region%.nii.gz}_DATA/${group}.txt"

}


################################################################################

# outer loop: iterate through the regions
for region in "${ROIs[@]}"; do

  # make a directory that will store group data directories
  region_dir="${region%.nii.gz}_DATA"
  mkdir "$region_dir"

  # middle loop: iterate through the subject groups
  for group in ${groups[@]}; do

    # make a text file for this group in the region directory
    touch ${group}.txt

    # first group: abstinent
    if [ $group = 'Abstinent' ]; then

      # inner loop: iterate over abstinent subjects
      for sub in ${Asubs[@]}; do

        # call our function
        get_vals $region $group $sub

      done # move on to relapse group

    # second group: relapse
    elif [ $group = 'Relapse' ]; then

      # inner loop: iterate over relapse subjects
      for sub in ${Rsubs[@]}; do

        # call our function
        get_vals $region $group $sub

      done # move on to controls group

    # third group: controls
    else

      # inner loop: iterate over control subjects
      for sub in ${Csubs[@]}; do

        # call our function
        get_vals $region $group $sub

      done # all done

    fi

  done

done
