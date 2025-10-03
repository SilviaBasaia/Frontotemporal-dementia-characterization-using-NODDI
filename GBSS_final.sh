#!/bin/bash

# GBSS Processing Script
# Author: S. Pisano
# Description: This script performs Grey Matter-Based Spatial Statistics (GBSS) analysis using FSL tools. 
# Based on the original GBSS pipeline by Nazeri A.,  available at https://github.com/arash-n/GBSS
# Prerequisites:
#   - FSL must be installed and configured. FSLDIR should be set in your environment.
#   - Coregister T1, DTI, NODDI maps to MNI space using FNIRT and applywarp tools.
#   - For each subject, the following files should be available:
#     - FA_2_MNI_nlsubject_*.nii.gz (FA maps)
#     - ODI_2_MNI_nlsubject_*.nii.gz (ODI maps)
#     - ICVF_2_MNI_nlsubject_*.nii.gz (ICVF maps)
#     - c1t1_2_MNI_nlsubject_*.nii.gz (T1 maps)
#    where * is the subject ID.
# The script outputs:
#   - Skeletonized GM, ODI, ICVF, and FA images that are contained in the filling directory, 
#     ready for further analysis (nonparametric permutation inference using fsl randomise).

# Set thresholds
thresh=0.65   # Grey matter threshold
perc=0.7      # Percentage of subjects for mask creation

# Create temporary working directory
mkdir -p GBSS

# Merge all GM images into a single 4D file
fslmerge -t GBSS/all_GM c1t1_2_MNI_nl*.nii.gz

# Merge ODI and ICVF images for all subjects
fslmerge -t GBSS/all_ODI ODI_2_MNI_nlsubject_* 
fslmerge -t GBSS/all_ICVF ICVF_2_MNI_nlsubject_*

cd GBSS

# Create mean GM mask
fslmaths all_GM -thr ${thresh} -bin -Tmean mean_GM
fslmaths mean_GM -thr 0.2 -bin GM_mask

# Skeletonize the mean GM mask
tbss_skeleton -i mean_GM.nii.gz -o GM_skel

# Create binary skeleton mask
fslmaths GM_skel.nii.gz -bin GM_skel_mask.nii.gz

# Select skeleton voxels above threshold
fslmaths GM_skel.nii.gz -thr ${thresh} -bin GM_skel_${thresh}

# Create distance mask for skeletonization
fslmaths GM_mask -mul -1 -add -1 -add GM_skel_${thresh} GM_mean_skeleton_mask_dst

# Compute distance map
distancemap -i GM_mean_skeleton_mask_dst.nii.gz -o GM_mean_skeleton_mask_dst

# Create a zero image for masking
fslmaths ${FSLDIR}/data/standard/LowerCingulum_1mm -mul 0 zero

# Skeletonize all GM images
tbss_skeleton -i mean_GM -p ${thresh} GM_mean_skeleton_mask_dst zero all_GM all_GM_skeletonise

# Create general GM skeleton mask based on subject percentage
fslmaths all_GM_skeletonise.nii.gz -thr ${thresh} -bin -Tmean -thr $perc -bin mean_GM_skeleton_mask_general

# Skeletonize ODI, ICVF, and FA with threshold
tbss_skeleton -i mean_GM -p $thresh GM_mean_skeleton_mask_dst zero all_GM all_ODI_skeletonised -a all_ODI
tbss_skeleton -i mean_GM -p $thresh GM_mean_skeleton_mask_dst zero all_GM all_ICVF_skeletonised -a all_ICVF
tbss_skeleton -i mean_GM -p $thresh GM_mean_skeleton_mask_dst zero all_GM all_FA_skeletonised -a all_FA

# Create lesion masks for GM, ICVF, and FA
fslmaths all_GM_skeletonise -mul mean_GM_skeleton_mask_general -uthr $thresh -bin all_lesion_GM
fslmaths all_ICVF_skeletonised -mul mean_GM_skeleton_mask_general -thr 0.65 -bin all_lesion_ICVF
fslmaths all_FA_skeletonised -mul mean_GM_skeleton_mask_general -thr 0.65 -bin all_lesion_FA

# Invert skeleton mask threshold
fslmaths GM_skel_mask.nii.gz -sub GM_skel_0.65.nii.gz GM_skel_thresh_inv

# Combine lesion masks
fslmaths all_lesion_GM -add GM_skel_thresh_inv.nii.gz -bin all_lesion

# Compute mean lesion mask
fslmaths all_lesion -Tmean lesion_mean

# --- Filling missing data in lesions ---

# Prepare filling directory
mkdir filling
cp all_lesion.nii.gz filling/
cp all_ICVF_skeletonised.nii.gz all_ODI_skeletonised.nii.gz filling/
cp all_FA_skeletonised.nii.gz filling/
cd filling

# Fill ICVF lesions
fslmaths all_lesion.nii.gz -sub 1 -mul -1 -mul all_ICVF_skeletonised.nii.gz -thr 0 fIC_non_lesion
fslmaths fIC_non_lesion -s 2 fIC_non_lesion_s_2
fslmaths fIC_non_lesion -bin -s 2 fIC_non_lesion_bin_s_2
fslmaths fIC_non_lesion_s_2 -div fIC_non_lesion_bin_s_2 fIC_filler
fslmaths all_lesion.nii.gz -mul fIC_filler.nii.gz all_lesion_filled_fIC
fslmaths fIC_non_lesion_bin_s_2 -thr 0.05 -bin -mul all_lesion.nii.gz -mul fIC_filler -add fIC_non_lesion -add all_lesion_filled_fIC all_fIC_filled

# Fill ODI lesions
fslmaths all_lesion.nii.gz -sub 1 -mul -1 -mul all_ODI_skeletonised.nii.gz -thr 0 ODI_non_lesion
fslmaths ODI_non_lesion -s 2 ODI_non_lesion_s_2
fslmaths ODI_non_lesion -bin -s 2 ODI_non_lesion_bin_s_2
fslmaths ODI_non_lesion_s_2 -div ODI_non_lesion_bin_s_2 ODI_filler
fslmaths all_lesion.nii.gz -mul ODI_filler.nii.gz all_lesion_filled_ODI
fslmaths ODI_non_lesion_bin_s_2 -thr 0.05 -bin -mul all_lesion.nii.gz -mul ODI_filler -add ODI_non_lesion -add all_lesion_filled_ODI all_ODI_filled

# Fill FA lesions

fslmaths all_lesion.nii.gz -sub 1 -mul -1 -mul all_FA_skeletonised.nii.gz -thr 0 FA_non_lesion
fslmaths FA_non_lesion -s 2 FA_non_lesion_s_2
fslmaths FA_non_lesion -bin -s 2 FA_non_lesion_bin_s_2
fslmaths FA_non_lesion_s_2 -div FA_non_lesion_bin_s_2 FA_filler
fslmaths all_lesion.nii.gz -mul FA_filler.nii.gz all_lesion_filled_FA
fslmaths FA_non_lesion_bin_s_2 -thr 0.05 -bin -mul all_lesion.nii.gz -mul FA_filler -add FA_non_lesion -add all_lesion_filled_FA all_FA_filled

# End of script
# Run  nonparametric permutation inference with fsl randomise on filled FA, ICVF and ODI data