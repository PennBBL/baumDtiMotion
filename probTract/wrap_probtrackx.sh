#!/bin/sh

#################################################################
### Define subject list in 'bblid_scanid_dateidofscan' format ###
#################################################################
subject_list=$(cat /data/joy/BBL/projects/pncBaumDti/Motion_paper/subjects/subject_lists/n949_dtiQA_b0_FSexclude_LTNv2_subjList.txt)

## Define path to local repo scripts
scripts_dir=/home/gbaum/baumDtiMotion/probTract

for name in ${subject_list}; do

	bblid=$(basename ${name} | cut -d_ -f1)
	scanid=$(basename ${name} | cut -d_ -f2)
	dateid=$(basename ${name} | cut -d_ -f3)

	echo $bblid
	echo $scanid
	echo $dateid

	subj="${bblid}"/"${dateid}"x"${scanid}"

	#####################################
	### Define input and output paths ###
	#####################################

	## Bedpostx directory containing "merged_f1samples.nii.gz" and other output
	bedpostx_dir=/data/joy/BBL/studies/pnc/processedData/diffusion/pncDTI_2016_04/"${subj}"/bedpostx_64_output

	## Probtrackx output directory for
	ptx_outdir=/data/joy/BBL/projects/pncBaumDti/probtrackx_2017/Motion_Paper/"${subj}"/wmEdge_p1000_pialTerm
	mkdir -p "${ptx_outdir}"/input

	## Subject log directory
	log_dir="${ptx_outdir}"/logfiles
	mkdir -p ${log_dir}

 	var0="pushd ${scripts_dir}; ./run_probtrackx.sh ${bblid} ${dateid} ${scanid}; popd"

	echo -e "${var0}"
	echo -e "${var0}" >> ${log_dir}/runPTx_wmBoundary_LausanneScale125_p1000_pialTerm_"${bblid}"_"${dateid}"x"${scanid}".sh

	subject_script=${log_dir}/runPTx_wmBoundary_LausanneScale125_p1000_pialTerm_"${bblid}"_"${dateid}"x"${scanid}".sh
 	
	## Execute qsub job for probtrackx2 runs for each subject 
	qsub -q all.q,basic.q -m e -M grahamlbaum@gmail.com -wd ${log_dir} -l h_vmem=4G,s_vmem=3G ${subject_script}

done

