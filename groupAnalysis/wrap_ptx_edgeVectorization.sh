#!/bin/sh

subject_list=$(cat /data/joy/BBL/projects/pncBaumDti/Motion_paper/subjects/subject_lists/n949_dtiQA_b0_FSexclude_LTNv2_subjList.txt)

# /data/joy/BBL/projects/pncBaumDti/Motion_paper/subjects/subject_lists/n148_19to23_dtiQA_b0_FSexclude_LTNv2_subjList.txt)

## Define path to local repo scripts
scripts_dir=/home/gbaum/baumDtiMotion/groupAnalysis

########################################################
### Create output files for vectorized network edges ###
########################################################
outdir=/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic
mkdir -p ${outdir}

## Probabilistic Streamline Count
SC_outpath=${outdir}/ptx_edgeVec_streamlineCount_LausanneScale125_n949.txt
rm "${SC_outpath}"
touch "${SC_outpath}"

## Volume-normalized Streamline count
volNormSC_outpath=${outdir}/ptx_edgeVec_volNormSC_LausanneScale125_n949.txt
touch "${volNormSC_outpath}"

## Connectivity Probability
connProb_outpath=${outdir}/ptx_edgeVec_connectivityProbability_LausanneScale125_n949.txt
touch "${connProb_outpath}"

## Mean Streamline Length
meanLength_outpath=${outdir}/ptx_edgeVec_meanStreamlineLength_LausanneScale125_n949.txt
rm "${meanLength_outpath}"
touch "${meanLength_outpath}"

##############################################################################
### Loop through subjects to vectorize probabilistic connectivity matrices ###
##############################################################################
for name in ${subject_list}; do

	bblid=$(basename ${name} | cut -d_ -f1)
	scanid=$(basename ${name} | cut -d_ -f2)
	dateid=$(basename ${name} | cut -d_ -f3)

	echo $bblid
	echo $scanid
	echo $dateid

	subj="${bblid}"/"${dateid}"x"${scanid}"

	ptx_adjmatpath=/data/joy/BBL/projects/pncBaumDti/probtrackx_2017/Motion_Paper/${subj}/wmEdge_p1000_pialTerm/output/"${bblid}"_"${dateid}"x"${scanid}"_ptx_p1000_wmBoundary_LausanneScale125_conmats.mat

	# eucDist_path=/data/joy/BBL/projects/pncBaumDti/probtrackx_2017/Motion_Paper/${subj}/wmEdge_p1000_pialTerm/roi/LausanneScale125/distance/"${bblid}"_"${dateid}"x"${scanid}"_ptx_wmBoundary_LausanneScale125_dti.txt 

	##############################################
	### Run Edge Vectorization (MATLAB) script ###
	##############################################
	pushd ${scripts_dir}

	matlab -nosplash -nodesktop -r "ptx_edgeVectorization ${ptx_adjmatpath} ${SC_outpath} ${volNormSC_outpath} ${connProb_outpath} ${meanLength_outpath}; exit()"

	popd

done
