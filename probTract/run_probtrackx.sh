#!/bin/sh

##################################
### Define subject identifiers ###
##################################
bblid=$1
dateid=$2
scanid=$3 

echo $bblid
echo $dateid
echo $scanid

subj="${bblid}"/"${dateid}"x"${scanid}"

## Define path to local repo scripts
scripts_dir=/home/gbaum/baumDtiMotion/probTract

## FreeSurfer subjects directory
SUBJECTS_DIR=/data/joy/BBL/studies/pnc/processedData/structural/freesurfer53

#####################################
### Define input and output paths ###
#####################################

## Bedpostx directory containing "merged_f1samples.nii.gz" and other output
bedpostx_dir=/data/joy/BBL/studies/pnc/processedData/diffusion/pncDTI_2016_04/"${subj}"/bedpostx_64_output

## COPY BEDPOSTX OUTPUT TO TMPDIR
mkdir -p ${TMPDIR}/${subj}
cp -r ${bedpostx_dir}/ ${TMPDIR}/${subj}

tmp_bedpostx_dir=${TMPDIR}/${subj}/bedpostx_64_output

## Probtrackx output directory for
ptx_outdir=/data/joy/BBL/projects/pncBaumDti/Motion_paper/replication/subjectData/probabilistic/"${subj}"/wmEdge_p1000_pialTerm
mkdir -p "${ptx_outdir}"/input

## Subject log directory
log_dir="${ptx_outdir}"/logfiles
mkdir -p ${log_dir}

## Sym-link to dti2xcp files
ln -s /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/dti2xcp/"${bblid}"_"${dateid}"x"${scanid}"_referenceVolume.nii.gz "${ptx_outdir}"/input/"${bblid}"_"${dateid}"x"${scanid}"_referenceVolume.nii.gz

ln -s /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/coreg/"${bblid}"_"${dateid}"x"${scanid}"_struct2seq.txt "${ptx_outdir}"/input/"${bblid}"_"${dateid}"x"${scanid}"_struct2seq.txt
	
#################################
### Define WM "Waypoint Mask" ###
#################################
	
## antsCT WM segmentation in T1 space
fslmaths /data/joy/BBL/studies/pnc/processedData/structural/antsCorticalThickness/${subj}/BrainSegmentation.nii.gz -thr 3 -uthr 3 "${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_WMseg03.nii.gz

## Move into subject diffusion space
WMseg_t1="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_WMseg03.nii.gz
WMseg_dti="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_WMseg03_dti.nii.gz
	
antsApplyTransforms -d 3 -e 0 -i "${WMseg_t1}" -r /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/dti2xcp/"${bblid}"_"${dateid}"x"${scanid}"_referenceVolume.nii.gz -t /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/coreg/"${bblid}"_"${dateid}"x"${scanid}"_struct2seq.txt -o "${WMseg_dti}" -n MultiLabel
waypoint_mask="${WMseg_dti}"

echo ""
echo "Probtrackx2 White Matter Waypoints Mask"
echo ${waypoint_mask}
echo ""

###################################
### Define CSF "Avoidance Mask" ###
###################################
fslmaths /data/joy/BBL/studies/pnc/processedData/structural/antsCorticalThickness/${subj}/BrainSegmentation.nii.gz -thr 1 -uthr 1 "${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_CSFseg01.nii.gz
	
# Move into subject diffusion space
csf_t1="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_CSFseg01.nii.gz
csf_dti="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_CSFseg01_dti.nii.gz
	
antsApplyTransforms -d 3 -e 0 -i "${csf_t1}" -r /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/dti2xcp/"${bblid}"_"${dateid}"x"${scanid}"_referenceVolume.nii.gz -t /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/coreg/"${bblid}"_"${dateid}"x"${scanid}"_struct2seq.txt -o "${csf_dti}" -n MultiLabel
avoid_mask="${csf_dti}"

echo ""
echo "Probtrackx2 CSF Avoidance Mask"
echo ${avoid_mask}
echo ""

####################################
### Create Pial Termination Mask ###
####################################
mri_convert ${SUBJECTS_DIR}/${subj}/mri/orig.mgz ${SUBJECTS_DIR}/${subj}/mri/orig.nii.gz

## Remove old FreeSurfer output generated from previous runs
# rm ${SUBJECTS_DIR}/${subj}/mri/*pial*
# rm ${SUBJECTS_DIR}/${subj}/mri/*Pial*
# rm ${SUBJECTS_DIR}/${subj}/label/*PialMask*

## Convert Pial surface to GIFTI format
mris_convert ${SUBJECTS_DIR}/${subj}/surf/rh.pial ${SUBJECTS_DIR}/${subj}/surf/rh.pial.gii
mris_convert ${SUBJECTS_DIR}/${subj}/surf/lh.pial ${SUBJECTS_DIR}/${subj}/surf/lh.pial.gii

####################################################################################
### Convert Pial surface .gii files to nifti volumes (still in FS conformed space) ###
####################################################################################
surf2volume ${SUBJECTS_DIR}/${subj}/surf/rh.pial.gii ${SUBJECTS_DIR}/${subj}/mri/orig.nii.gz ${SUBJECTS_DIR}/${subj}/mri/rh_pial.nii.gz freesurfer
surf2volume ${SUBJECTS_DIR}/${subj}/surf/lh.pial.gii ${SUBJECTS_DIR}/${subj}/mri/orig.nii.gz ${SUBJECTS_DIR}/${subj}/mri/lh_pial.nii.gz freesurfer

## Merge left and right hemisphere surface volumes into one Pial surface volume
fslmaths ${SUBJECTS_DIR}/${subj}/mri/rh_pial.nii.gz -add ${SUBJECTS_DIR}/${subj}/mri/lh_pial.nii.gz ${SUBJECTS_DIR}/${subj}/mri/Pial.nii.gz

## Binarize pial volume
fslmaths ${SUBJECTS_DIR}/${subj}/mri/Pial.nii.gz -bin ${SUBJECTS_DIR}/${subj}/mri/Pial.nii.gz
	
FS_pialVol=${SUBJECTS_DIR}/${subj}/mri/Pial.nii.gz

## Dilate pial mask
ImageMath 3 ${SUBJECTS_DIR}/${subj}/mri/Pial_dil1.nii.gz GD ${FS_pialVol} 1 

dil_pialVol=${SUBJECTS_DIR}/${subj}/mri/Pial_dil1.nii.gz
	
## Mask Lausanne ROIs by dilated Pial mask to get a strip of superficial GM near Pial surface
lausanne_atlas_img=${SUBJECTS_DIR}/${subj}/label/ROIv_scale125.nii.gz 	
	
# fslmaths ${lausanne_atlas_img} -mas ${dil_pialVol} ${SUBJECTS_DIR}/${subj}/label/ROIv_scale125_dil1PialMask.nii.gz 

###########################################
### Move pial mask into diffusion space ###
###########################################
FS_termination_mask=${SUBJECTS_DIR}/${subj}/label/ROIv_scale125_dil1PialMask.nii.gz 
RPI_termination_mask=${SUBJECTS_DIR}/${subj}/label/ROIv_scale125_dil1PialMask_RPI.nii.gz
t1_termination_mask=${SUBJECTS_DIR}/${subj}/label/ROIv_scale125_dil1PialMask_T1.nii.gz
dti_termination_mask=${SUBJECTS_DIR}/${subj}/label/ROIv_scale125_dil1PialMask_dti.nii.gz
	
## Change image orientation to RPI
mri_convert ${FS_termination_mask} ${RPI_termination_mask} --out_orientation LAS
	
## Reslice to T1
antsApplyTransforms -d 3 -e 0 -i ${RPI_termination_mask} -r /data/joy/BBL/studies/pnc/processedData/structural/antsCorticalThickness/"${bblid}"/"${dateid}"x"${scanid}"/ExtractedBrain0N4.nii.gz -o "${t1_termination_mask}" -n MultiLabel

## Move into subject diffusion space
antsApplyTransforms -d 3 -e 0 -i "${t1_termination_mask}" -r /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/dti2xcp/"${bblid}"_"${dateid}"x"${scanid}"_referenceVolume.nii.gz -t /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${subj}/coreg/"${bblid}"_"${dateid}"x"${scanid}"_struct2seq.txt -o "${dti_termination_mask}" -n MultiLabel
	
#####################################
### Create Tractography Seed ROIs ###
#####################################

## Dilate t1 atlas
gm_atlas=${SUBJECTS_DIR}/"${subj}"/label/ROI_scale125.nii.gz
dil_gm_atlas=${SUBJECTS_DIR}/"${subj}"/label/ROI_scale125_dil2.nii.gz
	
ImageMath 3 ${dil_gm_atlas} GD ${gm_atlas} 2

#################################
### Change orientation to RPI ###
#################################
mri_convert ${dil_gm_atlas} ${SUBJECTS_DIR}/"${subj}"/label/ROI_scale125_dil2_RPI.nii.gz --out_orientation LAS
	
RPI_outpath=${SUBJECTS_DIR}/"${subj}"/label/ROI_scale125_dil2_RPI.nii.gz
t1_outpath=${SUBJECTS_DIR}/"${subj}"/label/ROI_scale125_dil2_T1.nii.gz
dti_atlas_outpath="${ptx_outdir}"/input/"${bblid}"_"${dateid}"x"${scanid}"_ROI_scale125_dil2_dti.nii.gz

## Reslice to T1 space (for FreeSurfer-based atlases)
antsApplyTransforms -d 3 -e 0 -i "${RPI_outpath}" -r /data/joy/BBL/studies/pnc/processedData/structural/antsCorticalThickness/"${subj}"/ExtractedBrain0N4.nii.gz -o ${t1_outpath} -n MultiLabel

#########################################
### Move t1 Atlas to diffusion  space ###
#########################################
antsApplyTransforms -d 3 -e 0 -i "${t1_outpath}" -r /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${bblid}/"${dateid}"x"${scanid}"/dti2xcp/${bblid}_"${dateid}"x"${scanid}"_referenceVolume.nii.gz -t /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${bblid}/"${dateid}"x"${scanid}"/coreg/${bblid}_"${dateid}"x"${scanid}"_struct2seq.txt -o "${dti_atlas_outpath}" -n MultiLabel


## Dilate WM segmentation by 1 voxel to reduce regional dropout
# ImageMath 3 "${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_WMseg03_dil1.nii.gz GD ${WMseg_t1} 1
# dil1_WMseg_t1="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsCT_WMseg03_dil1.nii.gz
	
## Mask Dilated GM atlas by WM segementation 
# masked_gm_atlas=${SUBJECTS_DIR}/"${subj}"/label/WM_masked_ROI_scale125_dil2_T1.nii.gz
# fslmaths "${t1_outpath}" -mas "${dil1_WMseg_t1}" "${masked_gm_atlas}" 

###################################################
### Create white matter edge in diffusion space ###
###################################################
WMseg_edge="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsWMseg_edge_dti.nii.gz
dil_WMseg_edge="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_antsWMseg_edge_dti_dil1.nii.gz

fslmaths ${WMseg_dti} -edge ${WMseg_edge}

## Dilate WM edge by 1 voxel *This is used to get better coverage for frontal poles (regions 6 and 121)
ImageMath 3 ${dil_WMseg_edge} GD ${WMseg_edge} 1
	
## Mask dilated Lausanne atlas by WMseg_edge (in diffusion space)
fslmaths ${dti_atlas_outpath} -mas ${WMseg_edge} "${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti.nii.gz
	
dti_atlas_outpath="${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti.nii.gz

dti_seedVol=${dti_atlas_outpath}

####################################################################
### Move dilated atlas to diffusion space with dti2xcp transform ###
####################################################################
dti_atlas_outpath="${ptx_outdir}"/input/"${bblid}"_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti.nii.gz

antsApplyTransforms -d 3 -e 0 -i "${masked_gm_atlas}" -r /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${bblid}/"${dateid}"x"${scanid}"/dti2xcp/${bblid}_"${dateid}"x"${scanid}"_referenceVolume.nii.gz -t /data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${bblid}/"${dateid}"x"${scanid}"/coreg/${bblid}_"${dateid}"x"${scanid}"_struct2seq.txt -o "${dti_atlas_outpath}" -n MultiLabel

##############################
### REMOVE BRAINSTEM LABEL ###
##############################
fslmaths ${dti_atlas_outpath} -uthr 233 ${dti_atlas_outpath}
	
#######################################################
### Create Individual Seed/Target ROIs for tracking ###
#######################################################
seed_dir="${ptx_outdir}"/input/seedVols/LausanneScale125
mkdir -p "${seed_dir}"

pushd "${ptx_outdir}"/input

 	matlab -nodisplay -nodesktop -r "input_dir=dir('*wmBoundary_LausanneScale125_dti.nii.gz'); nii=load_nifti(input_dir.name); vol_orig = nii.vol; for i = 1:233; vol_roi = vol_orig; vol_roi(vol_roi ~= i) = 0; nii.vol = vol_roi; save_nifti(nii, strcat('ROIseed_', int2str(i), '.nii.gz')); end; exit"

 	mv "${ptx_outdir}"/input/ROIseed_*.nii.gz "${seed_dir}"
popd
	
### Replace Frontal Pole ROIs with dilated wmEdge ROIs ###
rm "${seed_dir}"/ROIseed_6.nii.gz
rm "${seed_dir}"/ROIseed_121.nii.gz

## Mask dilated Lausanne atlas by WMseg_edge (in diffusion space)
fslmaths ${dti_atlas_outpath} -mas ${dil_WMseg_edge} "${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_dil1_wmBoundary_LausanneScale125_dti.nii.gz

fslmaths "${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_dil1_wmBoundary_LausanneScale125_dti.nii.gz -thr 5 -uthr 7 "${seed_dir}"/ROIseed_6.nii.gz

fslmaths "${ptx_outdir}"/input/${bblid}_"${dateid}"x"${scanid}"_dil1_wmBoundary_LausanneScale125_dti.nii.gz -thr 120 -uthr 122 "${seed_dir}"/ROIseed_121.nii.gz

####################################################
### Create Seed-Target text file for Probtrackx2 ###
####################################################
seedTarget_file="${ptx_outdir}"/input/"${bblid}"_"${dateid}"x"${scanid}"_LausanneScale125_wmBoundary_ptx_seedTargets.txt

for i in {1..233}; do
	  echo "${i}"
	  echo "${seed_dir}"/ROIseed_"${i}".nii.gz >> ${seedTarget_file}
done

### Re-merge seed ROIs (including new Frontal Poles) to create seed atlas ###
binVols=$(ls "${seed_dir}"/*.nii.gz)

fslmerge -t ${ptx_outdir}/input/tmp_${bblid}_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti.nii.gz ${binVols}
fslmaths ${ptx_outdir}/input/tmp_${bblid}_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti.nii.gz -Tmax ${ptx_outdir}/input/${bblid}_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti.nii.gz

rm ${ptx_outdir}/input/tmp_${bblid}_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti.nii.gz

############################
### Calculate ROI volume ###
############################
atlas_img=${dti_seedVol}
roiVol_outpath=${ptx_outdir}/input/${bblid}_"${dateid}"x"${scanid}"_wmBoundary_LausanneScale125_dti_ROI_volume.txt

## LausanneScale125
for reg in {1..233}; do 
   	echo ${reg}
   	3dBrickStat -non-zero -count ${atlas_img}"<${reg}>" 2>>/dev/null 1>> ${roiVol_outpath}
done

#############################################
### Define variables for probtrackx2 call ###
#############################################
ptx_bin="/share/apps/fsl/5.0.9/bin/probtrackx2"

echo ""
echo "Probtrackx Output Directory"
echo ${ptx_outdir}
echo ""

## Number of streamlines propogated for each seed voxel
ncount=1000
echo ""
echo "Number of streamlines initiated in each seed voxel: ${ncount}"
echo ""

termination_mask=${SUBJECTS_DIR}/${subj}/label/ROIv_scale125_dil1PialMask_dti.nii.gz 
	
## Add sym-link for Termination mask -> input directory
ln -s ${termination_mask} ${ptx_outdir}/input/"${bblid}"_"${dateid}"x"${scanid}"_ROIv_scale125_dil1PialMask_dti.nii.gz
	
seed_paths=$(cat "${seedTarget_file}")
	
## Run probtrackx2 and setup output structure for each seed region
for seed in ${seed_paths};do
	seedName=$(basename ${seed} .nii.gz)
	echo ${seedName}

	curr_outdir=${ptx_outdir}/output/${seedName}_output
		
	mkdir -p ${curr_outdir}
	mkdir -p ${TMPDIR}/${subj}/${seedName}_output

	${ptx_bin} -s ${tmp_bedpostx_dir}/merged -m ${tmp_bedpostx_dir}/nodif_brain_mask.nii.gz -x ${seed} --seedref=${dti_seedVol} --avoid=${avoid_mask} --waypoints=${waypoint_mask} --stop=${termination_mask} --ompl --os2t --s2tastext --opd -l -c 0.2 -S 2000 --steplength=0.5 -P ${ncount} -V 1 --forcedir --dir=${TMPDIR}/${subj}/${seedName}_output --targetmasks=${seedTarget_file}

	# Matrix output seed-to-targets (text files)
	cp ${TMPDIR}/${subj}/${seedName}_output/matrix_seeds_to_all_targets ${curr_outdir}/
	cp ${TMPDIR}/${subj}/${seedName}_output/matrix_seeds_to_all_targets_lengths ${curr_outdir}/

	# Waytotal
	cp ${TMPDIR}/${subj}/${seedName}_output/waytotal ${curr_outdir}/

	# Path distributions
	cp ${TMPDIR}/${subj}/${seedName}_output/fdt_paths* ${curr_outdir}/

	# Probtrackx log file
	cp ${TMPDIR}/${subj}/${seedName}_output/probtrackx.log ${curr_outdir}/

done

####################################################
### Generate Probabilistic Connectivity Matrices ###
####################################################
conmat_outpath=${ptx_dir}/output/${bblid}_"${dateid}"x"${scanid}"_ptx_p1000_wmBoundary_LausanneScale125_conmats.mat

roiVol_path=${roiVol_outpath}

pushd ${scripts_dir}

matlab -nosplash -nodesktop -logfile ${log_dir}/gen_ptx_conmat_LausanneScale125.log -r "generate_ptx_conmat ${ptx_dir} ${conmat_outpath} ${roiVol_path} ${bblid} ${dateid} ${scanid}; exit()"
popd
