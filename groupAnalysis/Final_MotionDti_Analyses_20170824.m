clear all

cd /data/joy/BBL/projects/pncBaumDti/Motion_paper
% load matlab_probabilistic_analysis_workspace_20170604.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in subject demographics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demog=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/subjects/demogs/ML_Formatted_Motion_DTI_sample_n949_demogs.csv',',',1,1);

scanid=demog(:,1)
age=demog(:,2)
ageSq=demog(:,3)
ageCub=demog(:,4)
sex=demog(:,5)
motion=demog(:,6)
envSES=demog(:,7)
tsnr=demog(:,9)
clipCount=demog(:,13)
outmax=demog(:,16)
outmean=demog(:,17)

restMotion = dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/n949_restRelMeanRMSMotion.txt');

eddyAvgTrans = dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/subjects/demogs/n949_avg_vol2vol_translation_eddyOutput.txt')

eddyAvgRot = dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/subjects/demogs/n949_avg_vol2vol_rotation_eddyOutput.txt')


%% Yeo 7-system partition (Lausanne 234)
Yeo_part=dlmread('/data/joy/BBL/projects/pncBaumDti/Yeo_in_Lausanne_partitions/Yeo_7system_in_Lausanne234.txt');
S=Yeo_part;
S(234)=[];
Yeo_part(234)=[] % remove brainstem
numComms=numel(unique(Yeo_part));
%% Within/Between Module index
withinBetween_mat=~bsxfun(@eq,Yeo_part,Yeo_part');
withinBetween_edge_index=squareform(withinBetween_mat)';
within_idx=find(withinBetween_edge_index==0);
between_idx=find(withinBetween_edge_index==1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in vectorized PROBABILISTIC connectivitity matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Araw_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/edgevec/ptx_edgeVec_streamlineCount_LausanneScale125_n949.txt');
%  Araw_edgeVec([highMotion_idx],:)=[];

Araw_volNorm_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/edgevec/ptx_edgeVec_volNormSC_LausanneScale125_n949.txt');
% Araw_volNorm_edgeVec([highMotion_idx],:)=[];

Aprop_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/edgevec/ptx_edgeVec_connectivityProbability_LausanneScale125_n949.txt');
% Aprop_edgeVec([highMotion_idx],:)=[];

Alength_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/edgevec/ptx_edgeVec_meanStreamlineLength_LausanneScale125_n949.txt');
% Alength_edgeVec([highMotion_idx],:)=[];

% Euclidean Distance matrix (derived from MNI template -> FreeSurfer recon -> Lausanne)
eucDist_mat = dlmread('/data/joy/BBL/studies/pnc/processedData/structural/freesurfer53/MNI_1mm_template/label/distance/MNI1mm_LausanneScale125_FS.txt');
	eucDist_mat(:,[234])=[];
	eucDist_mat([234],:)=[];

eucDist_edgeVec = squareform(eucDist_mat)';

% eucDist_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/edgevec/ptx_edgeVec_euclideanDistance_LausanneScale125_n949.txt');
% eucDist_edgeVec([highMotion_idx],:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define number of Subjects, Regions, Edges %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsub=size(Araw_edgeVec,1)
nreg=length(squareform(Araw_edgeVec(1,:)))
nedge=length(Araw_edgeVec(1,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mean length measures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_length=mean(Alength_edgeVec)';
median_length=median(Alength_edgeVec)';
mean_eucDist=mean(eucDist_edgeVec)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in vectorized DETERMINISTIC connectivitity matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create inverse ADC matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADC_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/deterministic/edgeVec/full_EdgeVec_det_ADC_LausanneScale125_n949.txt');

inv_ADC_edgeVec=zeros(nsub,nedge);

for i = 1:nsub
	A=ADC_edgeVec(i,:);
	inv_A= 1 ./ A;
	inv_A(isinf(inv_A)) = 0;
	inv_A=squareform(inv_A);
	inv_A(:,[234])=[];
	inv_A([234],:)=[];
	inv_A=squareform(inv_A);
	inv_ADC_edgeVec(i,:)=inv_A;
end

subj_mean_ADC = mean(inv_ADC_edgeVec)';
mean_ADC_mat = squareform(subj_mean_ADC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mean Streamline Length %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_Length_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/deterministic/edgeVec/full_EdgeVec_det_mean_fiberLength_LausanneScale125_n949.txt');
Length_edgeVec=zeros(nsub, nedge);

for i=1:nsub
	A=squareform(tmp_Length_edgeVec(i,:));
	A(:,[234])=[];
	A([234],:)=[];
	Length_edgeVec(i,:)=squareform(A);
end
subj_mean_length=mean(Length_edgeVec)';
subj_median_length=median(Length_edgeVec)';


%%%%%%%%%%%%%%%
%%% Mean FA %%%
%%%%%%%%%%%%%%%
tmp_FA_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/deterministic/edgeVec/full_EdgeVec_det_FA_LausanneScale125_n949.txt');
FA_edgeVec=zeros(nsub, nedge);

%%% REMOVE BRAINSTEM %%%
for i=1:nsub
	A=squareform(tmp_FA_edgeVec(i,:));
	A(:,[234])=[];
	A([234],:)=[];
	FA_edgeVec(i,:)=squareform(A);
end

subj_mean_FA=mean(FA_edgeVec)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Volume-normalized Streamline Count %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_volNormSC_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/deterministic/edgeVec/full_EdgeVec_det_volNormSC_LausanneScale125_n949.txt');
volNormSC_edgeVec=zeros(nsub, nedge);

for i=1:nsub
	A=squareform(tmp_volNormSC_edgeVec(i,:));
	A(:,[234])=[];
	A([234],:)=[];
	volNormSC_edgeVec(i,:)=squareform(A);
end

subj_mean_volNormSC=mean(volNormSC_edgeVec)';


%%%%%%%%%%%%%%%%%%%%%%%%
%%% Streamline Count %%%
%%%%%%%%%%%%%%%%%%%%%%%%
tmp_SC_edgeVec=dlmread('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/deterministic/edgeVec/full_EdgeVec_det_streamlineCount_LausanneScale125_n949.txt');

SC_edgeVec=zeros(nsub, nedge);

for i=1:nsub
	A=squareform(tmp_SC_edgeVec(i,:));
	A(:,[234])=[];
	A([234],:)=[];
	SC_edgeVec(i,:)=squareform(A);
end

subj_mean_SC=mean(SC_edgeVec)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Network Density %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
det_netDensity=zeros(nsub,1);
for i = 1:nsub
	A=squareform(inv_ADC_edgeVec(i,:));
	kden = density_und(A);
	det_netDensity(i) = kden;
end

mean(det_netDensity)
std(det_netDensity)


prob_netDensity=zeros(nsub,1);
for i = 1:nsub
	A=squareform(Araw_edgeVec(i,:));
	kden = density_und(A);
	prob_netDensity(i) = kden;
end

mean(prob_netDensity)
std(prob_netDensity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BrainNet Viewer: Edgewise motion effect Renderings %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Connectivity Probability
[Edgewise_Aprop_motion_r,Edgewise_Aprop_motion_p]=partialcorr(Aprop_edgeVec, motion, [age ageSq sex]);
[Edgewise_Aprop_pthr,Edgewise_Aprop_pcor,Edgewise_Aprop_padj]=fdr(Edgewise_Aprop_motion_p,0.05);
NOTsig_edge_idx=find(Edgewise_Aprop_motion_p > Edgewise_Aprop_pthr);
Edgewise_Aprop_motion_p(NOTsig_edge_idx)=0;
Edgewise_Aprop_motion_r(NOTsig_edge_idx)=0;
Aprop_r_mat=squareform(Edgewise_Aprop_motion_r);
density_und(Aprop_r_mat)
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/BrainNet/probabilistic/Edgewise_Aprop_motion_partialR_matrix_n949_FDRthresh.edge', Aprop_r_mat, 'delimiter','\t')

%% Streamline Count
[Edgewise_Araw_motion_r,Edgewise_Araw_motion_p]=partialcorr(Araw_edgeVec,motion,[age ageSq sex]);
[Edgewise_Araw_pthr,Edgewise_Araw_pcor,Edgewise_Araw_padj]=fdr(Edgewise_Araw_motion_p,0.05);
NOTsig_edge_idx=find(Edgewise_Araw_motion_p >= Edgewise_Araw_pthr);
Edgewise_Araw_motion_p(NOTsig_edge_idx)=0;
Edgewise_Araw_motion_r(NOTsig_edge_idx)=0;
Araw_r_mat=squareform(Edgewise_Araw_motion_r);
density_und(Araw_r_mat)
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/BrainNet/probabilistic/Edgewise_Araw_motion_partialR_matrix_n949_FDRthresh.edge', Araw_r_mat, 'delimiter','\t')

%% Streamline Length
[Edgewise_Alength_motion_r,Edgewise_Alength_motion_p]=partialcorr(Alength_edgeVec,motion,[age ageSq sex]);
[Edgewise_Alength_pthr,Edgewise_Alength_pcor,Edgewise_Alength_padj]=fdr(Edgewise_Alength_motion_p,0.05);
NOTsig_edge_idx=find(Edgewise_Alength_motion_p >= Edgewise_Alength_pthr);
Edgewise_Alength_motion_p(NOTsig_edge_idx)=0;
Edgewise_Alength_motion_r(NOTsig_edge_idx)=0;

rm_idx=find(Edgewise_Alength_motion_r > -0.3 & Edgewise_Alength_motion_r < 0);
Edgewise_Alength_motion_r(rm_idx)=0;
% rm_idx=find(Edgewise_Alength_motion_r < 0.4 & Edgewise_Alength_motion_r > 0) 
% Edgewise_Alength_motion_r(rm_idx)=0;
% rm_idx=find(Edgewise_Alength_motion_r < 0.10 & Edgewise_Alength_motion_r > 0) 
% Edgewise_Alength_motion_r(rm_idx)=0;
Alength_r_mat=squareform(Edgewise_Alength_motion_r);
density_und(Alength_r_mat)
figure; hist(nonzeros(Edgewise_Alength_motion_r))

dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/BrainNet/probabilistic/THRESH_Edgewise_Alength_motion_partialR_matrix_n949_FDRthresh.edge', Alength_r_mat, 'delimiter','\t')


%%%%%%%%%
%% ADC %%
%%%%%%%%%
[Edgewise_ADC_motion_r,Edgewise_ADC_motion_p]=partialcorr(inv_ADC_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_ADC_motion_r));
Edgewise_ADC_motion_r(isnan_idx)=[];
Edgewise_ADC_motion_p(isnan_idx)=[];
% FDR correction
[Edgewise_ADC_motion_r,Edgewise_ADC_motion_p]=partialcorr(inv_ADC_edgeVec,motion,[age ageSq sex]);
% Re-calc pvals
[Edgewise_ADC_pthr,Edgewise_ADC_pcor,Edgewise_ADC_padj]=fdr(Edgewise_ADC_motion_p,0.05);
NOTsig_edge_idx=find(Edgewise_ADC_motion_p > Edgewise_ADC_pthr);
Edgewise_ADC_motion_p(NOTsig_edge_idx)=0;
Edgewise_ADC_motion_r(NOTsig_edge_idx)=0;
isnan_idx=find(isnan(Edgewise_ADC_motion_r));
Edgewise_ADC_motion_r(isnan_idx)=0;
Edgewise_ADC_motion_p(isnan_idx)=0;
ADC_r_mat=squareform(Edgewise_ADC_motion_r);
density_und(ADC_r_mat)
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/BrainNet/deterministic/Edgewise_invADC_motion_partialR_matrix_n949_FDRthresh.edge', ADC_r_mat, 'delimiter','\t')

%%%%%%%%%%%%%%%%
%% Det Length %%
%%%%%%%%%%%%%%%%
[Edgewise_length_motion_r,Edgewise_length_motion_p]=partialcorr(Length_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_length_motion_r));
Edgewise_length_motion_r(isnan_idx)=[];
Edgewise_length_motion_p(isnan_idx)=[];
% FDR correction
[Edgewise_length_pthr,Edgewise_length_pcor,Edgewise_length_padj]=fdr(Edgewise_length_motion_p,0.05);
% Re-calc pvals
[Edgewise_length_motion_r,Edgewise_length_motion_p]=partialcorr(Length_edgeVec,motion,[age ageSq sex]);
NOTsig_edge_idx=find(Edgewise_length_motion_p > Edgewise_length_pthr);
Edgewise_length_motion_p(NOTsig_edge_idx)=0;
Edgewise_length_motion_r(NOTsig_edge_idx)=0;
isnan_idx=find(isnan(Edgewise_length_motion_r));
Edgewise_length_motion_r(isnan_idx)=0;
Edgewise_length_motion_p(isnan_idx)=0;
length_r_mat=squareform(Edgewise_length_motion_r);
density_und(length_r_mat)
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/BrainNet/deterministic/Edgewise_streamlineLength_motion_partialR_matrix_n949_FDRthresh.edge', length_r_mat, 'delimiter','\t')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROBABILISTIC: Streamline Count %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Coefficient of Variation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_mats=zeros(nreg,nreg,nsub);
	for i = 1:nsub
		A=squareform(Araw_edgeVec(i,:)); 
		group_mats(:, :,i) = A;
	end
[W_thr, Wcv] = GLB_threshold_consistency(group_mats, 1);
%% Set diagonal to zero
Wcv =(Wcv - diag(diag(Wcv)));
%% Set NaN to zero
isnan_idx=find(isnan(Wcv));
Wcv(isnan_idx)=0;
Araw_Wcv=squareform(Wcv)';

%% Node strength FDR thresh
orig_Araw_nodeStrength=zeros(nsub,nreg);
	for i = 1:nsub
		A=squareform(Araw_edgeVec(i,:));
		nodeStrength=sum(A,1);
		orig_Araw_nodeStrength(i,:)=nodeStrength;
	end

	[r,p]=partialcorr(orig_Araw_nodeStrength, motion,[age ageSq sex]) ;
	[orig_Nodewise_Araw_motion_pthr, Nodewise_Araw_motion_pcor, Nodewise_Araw_motion_padj] = fdr(p,0.05);
	orig_Nodewise_Araw_motion_pthr

%% Motion partial correlation
[Edgewise_Araw_motion_r,Edgewise_Araw_motion_p]=partialcorr(Araw_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_Araw_motion_r));
Edgewise_Araw_motion_p(isnan_idx)=[];
Edgewise_Araw_motion_r(isnan_idx)=[];
tmp_mean_length=mean_length; 
tmp_mean_length(isnan_idx)=[];
tmp_Araw_Wcv=Araw_Wcv; tmp_Araw_Wcv(isnan_idx)=[];

corr(Edgewise_Araw_motion_r,tmp_mean_length)
corr(Edgewise_Araw_motion_r,-log(tmp_Araw_Wcv))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export dataframe for 3rd-level scatterplots in R -- FIGURE 2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Araw_output_df=zeros(nedge,3);
Araw_output_df(:,1)=Edgewise_Araw_motion_r;
Araw_output_df(:,2)=-log(Araw_Wcv);
Araw_output_df(:,3)=mean_length;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/Araw_edgewise_motion_effects_n949.csv',Araw_output_df)

%% 3rd-level correlation plot
corr(-log(Araw_Wcv),Edgewise_Araw_motion_r) % Consistency
corr(Edgewise_Araw_motion_r,mean_length) % mean Streamline length

%% Scatterplot
c = linspace(1,max(mean_length),length(mean_length)); figure;scatter(-log(Araw_Wcv),Edgewise_Araw_motion_r,[],c); set(gcf,'color','white')

%% Calculate FDR threshold for Motion and Age Effects
[Edgewise_Araw_motion_r,Edgewise_Araw_motion_p]=partialcorr(Araw_edgeVec,motion,[age ageSq sex]);
[Edgewise_Araw_pthr,Edgewise_Araw_pcor,Edgewise_Araw_padj]=fdr(Edgewise_Araw_motion_p,0.05);
sig_edge_idx=find(Edgewise_Araw_motion_p < Edgewise_Araw_pthr);
sig_edges=Edgewise_Araw_motion_r(sig_edge_idx);
sigPos_edges=find(sig_edges > 0);size(sigPos_edges)
sigNeg_edges=find(sig_edges < 0); size(sigNeg_edges)
perc_Araw_sig_Motion = (length(sig_edge_idx) / length(Edgewise_Araw_motion_p)) * 100;

%% Age effect with motion
[Edgewise_Araw_age_r,Edgewise_Araw_age_p]=partialcorr(Araw_edgeVec, age,[ageSq motion sex]);
isnan_idx=find(isnan(Edgewise_Araw_age_p));
Edgewise_Araw_age_p(isnan_idx)=[];
Edgewise_Araw_age_r(isnan_idx)=[];
[Araw_age_pthr,Araw_age_pcor, Araw_age_padj] = fdr(Edgewise_Araw_age_p,0.05);
Age_sig_edge_idx=find(Edgewise_Araw_age_p < Araw_age_pthr);

[Nodewise_Araw_age_r,Nodewise_Araw_age_p]=partialcorr(orig_Araw_nodeStrength, age,[ageSq motion sex]);
isnan_idx=find(isnan(Nodewise_Araw_age_p));
Nodewise_Araw_age_p(isnan_idx)=[];
Nodewise_Araw_age_r(isnan_idx)=[];
[Araw_Node_age_pthr, Araw_Node_age_pcor, Araw_Node_age_padj] = fdr(Nodewise_Araw_age_p,0.05);
Age_sig_node_idx=find(Nodewise_Araw_age_p < Araw_Node_age_pthr);

%% Age Effect without motion
[Edgewise_Araw_age_noMotion_r, Edgewise_Araw_age_noMotion_p]=partialcorr(Araw_edgeVec,age,[ageSq sex]);
isnan_idx=find(isnan(Edgewise_Araw_age_noMotion_p));
Edgewise_Araw_age_noMotion_p(isnan_idx)=[];
Edgewise_Araw_age_noMotion_r(isnan_idx)=[];
[Araw_age_noMotion_pthr, Araw_age_noMotion_pcor,age_noMotionpadj] = fdr(Edgewise_Araw_age_noMotion_p,0.05);
Age_noMotion_sig_edge_idx=find(Edgewise_Araw_age_noMotion_p < Araw_age_noMotion_pthr);

[Nodewise_Araw_age_noMotion_r, Nodewise_Araw_age_noMotion_p]=partialcorr(orig_Araw_nodeStrength, age,[ageSq sex]);
% isnan_idx=find(isnan(Nodewise_Araw_age_p));
% Nodewise_Araw_age_noMotion_p(isnan_idx)=[];
% Nodewise_Araw_age_noMotion_r(isnan_idx)=[];
[Araw_Node_age_noMotion_pthr, Araw_Node_age_noMotion_pcor, Araw_Node_age_noMotion_padj] = fdr(Nodewise_Araw_age_noMotion_p, 0.05);
Age_noMotion_sig_node_idx=find(Nodewise_Araw_age_noMotion_p < Araw_Node_age_noMotion_pthr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOTION and AGE EFFECTS ACROSS CONSISTENCY THRESHOLDS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold_range=[1.0; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3; 0.2; 0.1];

%% Edge measures
W_Araw_acrossThresh=zeros(nedge,length(threshold_range)); % Group mean connectivity matrices after thresholding
Araw_edgeVec_acrossThresh=zeros(nsub,nedge,length(threshold_range));
Edgewise_Araw_motion_partR_acrThresh=zeros(nedge,length(threshold_range));
Edgewise_Araw_motion_pval_acrThresh=zeros(nedge,length(threshold_range));

%% Node Strength
prob_Araw_nodeStrength=zeros(nsub, nreg, length(threshold_range));
prob_Araw_nodeStrength_partR_acrThresh=zeros(nreg, length(threshold_range),1);
prob_Araw_nodeStrength_pval_acrThresh=zeros(nreg, length(threshold_range),1);
perc_Araw_nodeStrength_sig_Motion=zeros(length(threshold_range),1);

%% Network Strength - Motion Effects
prob_Araw_netStrength=zeros(nsub,length(threshold_range));
prob_Araw_netStrength_partR_acrThresh=zeros(length(threshold_range),1);
prob_Araw_netStrength_pval_acrThresh=zeros(length(threshold_range),1);
perc_Araw_MotionEffect_acrossThresh=zeros(length(threshold_range),1);

%% Network Strength - Age Effects
prob_Araw_netStrength_Age_partR_acrThresh=zeros(length(threshold_range),1);
prob_Araw_netStrength_Age_noMotion_partR_acrThresh=zeros(length(threshold_range),1);

% 3rd-level correlations
thirdLevelConsistency_Araw_r=zeros(length(threshold_range),1);
thirdLevelmeanLength_Araw_r=zeros(length(threshold_range),1);

% Edgewise Age Effects
perc_Araw_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Araw_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);

Edgewise_Araw_age_partialR_acrossThresh=zeros(nedge,length(threshold_range));
Edgewise_Araw_age_noMotion_partialR_acrossThresh=zeros(nedge,length(threshold_range));

% Nodewise Age Effects
perc_Araw_nodewise_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Araw_nodewise_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);
Araw_nodewise_AgeEffect_partR_acrossThresh=zeros(nreg, length(threshold_range),1);
Araw_nodewise_AgeEffect_noMotion_partR_acrossThresh=zeros(nreg, length(threshold_range),1);

%% Within-module connectivity
withinConn_Araw_acrossThresh=zeros(nsub,length(threshold_range));
withinConn_Araw_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Araw_age_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Araw_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Global Efficiency
globEff_Araw_acrossThresh=zeros(nsub,length(threshold_range));
globEff_Araw_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Araw_age_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Araw_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Distance-dependence of Age effects (3rd-level)
Araw_age_length_thirdLevel_r=zeros(length(threshold_range),1);
Araw_age_noMotion_length_thirdLevel_r=zeros(length(threshold_range),1);

Araw_age_consistency_thirdLevel_r=zeros(length(threshold_range),1);
Araw_age_noMotion_consistency_thirdLevel_r=zeros(length(threshold_range),1);

for T=1:length(threshold_range)
	%% Define group of Probabilistic connectivity matrices for thresholding
	thresh=threshold_range(T);
	% perc_out=thresh * 100
	thresh_Araw_edgeVec=Araw_edgeVec;
	[W_thr, Wcv] = GLB_threshold_consistency(group_mats, thresh);
	sq_W=squareform(W_thr)';
	W_Araw_acrossThresh(:,T)=sq_W;
	edgeThresh_index=find(W_Araw_acrossThresh(:,T) == 0)';
	thresh_Araw_edgeVec(:,[edgeThresh_index])=0;
	Araw_edgeVec_acrossThresh(:,:,T)=thresh_Araw_edgeVec;	
	% filename=strcat('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/consistencyThresh/ptx_streamlineCount_BreakspearThresh_',num2str(perc_out),'perc_ArawSC_edgeVec_n949.txt')
	% dlmwrite(filename,thresh_Araw_edgeVec)


	%% Motion effects on Network Strength for each Consistency Threshold
	netStrength=sum(thresh_Araw_edgeVec,2);
	prob_Araw_netStrength(:,T)=netStrength;
	[r,p]=partialcorr(netStrength,motion,[age ageSq sex]);
	prob_Araw_netStrength_partR_acrThresh(T)=r;
	prob_Araw_netStrength_pval_acrThresh(T)=p;

	%% Age effects on Network Strength for each Consistency Threshold (with and without motion)
	
	% Controlling for Motion
	[r,p]=partialcorr(netStrength,age,[motion ageSq sex]);
	prob_Araw_netStrength_Age_partR_acrThresh(T)=r;
	% NOT controlling for motion
	[r,p]=partialcorr(netStrength,age,[ageSq sex]);
	prob_Araw_netStrength_Age_noMotion_partR_acrThresh(T)=r;

	%% Edgewise Motion Effects for each Consistency Threshold
	[Edgewise_Araw_motion_r,Edgewise_Araw_motion_p]=partialcorr(thresh_Araw_edgeVec,motion,[age ageSq sex]);
	Edgewise_Araw_motion_partR_acrThresh(:,T)=Edgewise_Araw_motion_r;
	Edgewise_Araw_motion_pval_acrThresh(:,T)=Edgewise_Araw_motion_p;


	%% Node strength effects
	for i = 1:nsub
		A=squareform(thresh_Araw_edgeVec(i,:));
		nodeStrength=sum(A,1);
		prob_Araw_nodeStrength(i,:,T)=nodeStrength;
	end

	[r,p]=partialcorr(prob_Araw_nodeStrength(:,:,T),motion, [age ageSq sex]);
	sig_node_idx=find(p < orig_Nodewise_Araw_motion_pthr);
	NOT_sig_node_idx=find(p > orig_Nodewise_Araw_motion_pthr);
	r(NOT_sig_node_idx)=0;
	prob_Araw_nodeStrength_partR_acrThresh(:,T)=r;
	perc_Araw_nodeStrength_sig_Motion(T)=(length(sig_node_idx) / length(p)) * 100;

	%% Percentage of Edgewise motion effects across thresholds
	[Edgewise_Araw_motion_r,Edgewise_Araw_motion_p]=partialcorr(thresh_Araw_edgeVec,motion,[age ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Araw_motion_p));
	Edgewise_Araw_motion_p(isnan_idx)=[];
	Edgewise_Araw_motion_r(isnan_idx)=[];
	tmp_consistency=Araw_Wcv;
	tmp_consistency(isnan_idx)=[];
	tmp_length=mean_length;
	tmp_length(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Araw_motion_p < Edgewise_Araw_pthr);
	perc_Araw_sig_Motion=(length(sig_edge_idx) / length(Edgewise_Araw_motion_p)) * 100;
	perc_Araw_MotionEffect_acrossThresh(T)=perc_Araw_sig_Motion;

	%% 3rd-level Consistency Effect
	[r,p]=corr(Edgewise_Araw_motion_r, -log(tmp_consistency));
	thirdLevelConsistency_r(T)=r;

	[r,p]=corr(Edgewise_Araw_motion_r, tmp_length);
	thirdLevelmeanLength_r(T)=r;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Edgewise Age effects with/without controlling for motion %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% WITH Motion: Save age effects across thresholds before removing NaNs
	[Edgewise_Araw_age_r, Edgewise_Araw_age_p] = partialcorr(thresh_Araw_edgeVec, age,[motion ageSq sex]);
	NOTsig_edge_idx=find(Edgewise_Araw_age_p >= Araw_age_pthr);
	isnan_idx=find(isnan(Edgewise_Araw_age_p));
	Edgewise_Araw_age_r(isnan_idx)=0;
	Edgewise_Araw_age_r(NOTsig_edge_idx)=0;
	Edgewise_Araw_age_partialR_acrossThresh(:,T)=Edgewise_Araw_age_r;
	%% Now remove NaNs to calculate percentage of sig edges
	[Edgewise_Araw_age_r, Edgewise_Araw_age_p]=partialcorr(thresh_Araw_edgeVec,age,[motion ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Araw_age_p));
	Edgewise_Araw_age_p(isnan_idx)=[];
	Edgewise_Araw_age_r(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Araw_age_p < Araw_age_pthr);
	perc_Araw_sig_Age=(length(sig_edge_idx) / length(Edgewise_Araw_age_p)) * 100;
	perc_Araw_AgeEffect_acrossThresh(T)=perc_Araw_sig_Age;	

	%% Distance-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Araw_age_r, tmp_length);
	Araw_age_length_thirdLevel_r(T)=r;

	%% Consistency-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Araw_age_r, -log(tmp_consistency));
	Araw_age_consistency_thirdLevel_r(T)=r;

	%% WITHOUT Motion: Save age effects across thresholds before removing NaNs
	[Edgewise_Araw_age_noMotion_r,Edgewise_Araw_age_noMotion_p]=partialcorr(thresh_Araw_edgeVec,age,[ageSq sex]);
	NOTsig_edge_idx=find(Edgewise_Araw_age_noMotion_p >= Araw_age_noMotion_pthr);
	isnan_idx=find(isnan(Edgewise_Araw_age_noMotion_p));
	Edgewise_Araw_age_noMotion_r(isnan_idx)=0;
	Edgewise_Araw_age_noMotion_partialR_acrossThresh(:,T)=Edgewise_Araw_age_noMotion_r;
	%% Now remove NaNs to calculate percentage of sig edges
	[Edgewise_Araw_age_noMotion_r,Edgewise_Araw_age_noMotion_p]=partialcorr(thresh_Araw_edgeVec,age,[ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Araw_age_noMotion_p));
	Edgewise_Araw_age_noMotion_p(isnan_idx)=[];
	Edgewise_Araw_age_noMotion_r(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Araw_age_noMotion_p < Araw_age_noMotion_pthr);
	perc_Araw_sig_Age_noMotion=(length(sig_edge_idx) / length(Edgewise_Araw_age_noMotion_p)) * 100;
	perc_Araw_AgeEffect_noMotion_acrossThresh(T)=perc_Araw_sig_Age_noMotion;	

	%% Distance-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Araw_age_noMotion_r, tmp_length);
	Araw_age_noMotion_length_thirdLevel_r(T)=r;

	%% Consistency-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Araw_age_noMotion_r, -log(tmp_consistency));
	Araw_age_noMotion_consistency_thirdLevel_r(T)=r;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Nodewise Age effects with/without controlling for motion %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[Nodewise_Araw_age_r, Nodewise_Araw_age_p]=partialcorr(prob_Araw_nodeStrength(:,:,T), age,[motion ageSq sex]);
	% isnan_idx=find(isnan(Nodewise_Araw_age_p));
	% Nodewise_Araw_age_p(isnan_idx)=[];
	% Nodewise_Araw_age_r(isnan_idx)=[];
	sig_node_idx=find(Nodewise_Araw_age_p < Araw_Node_age_pthr); 
	NOTsig_node_idx=find(Nodewise_Araw_age_p >= Araw_Node_age_pthr);
	perc_node_Araw_sig_Age=(length(sig_node_idx) / length(Nodewise_Araw_age_p)) * 100;
	perc_Araw_nodewise_AgeEffect_acrossThresh(T)=perc_node_Araw_sig_Age;
	Nodewise_Araw_age_r(NOTsig_node_idx)=0;
	Araw_nodewise_AgeEffect_partR_acrossThresh(:,T)=Nodewise_Araw_age_r;


	[Nodewise_Araw_age_noMotion_r, Nodewise_Araw_age_noMotion_p]=partialcorr(prob_Araw_nodeStrength(:,:,T), age,[ageSq sex]);
	% isnan_idx=find(isnan(Nodewise_Araw_age_noMotion_p));
	% Nodewise_Araw_age_noMotion_p(isnan_idx)=[];
	% Nodewise_Araw_age_noMotion_r(isnan_idx)=[];
	sig_node_idx=find(Nodewise_Araw_age_noMotion_p < Araw_Node_age_noMotion_pthr);
	NOTsig_node_idx=find(Nodewise_Araw_age_noMotion_p >= Araw_Node_age_noMotion_pthr);
	perc_node_Araw_sig_Age_noMotion=(length(sig_node_idx) / length(Nodewise_Araw_age_noMotion_p)) * 100;
	perc_Araw_nodewise_AgeEffect_noMotion_acrossThresh(T)= perc_node_Araw_sig_Age_noMotion;
	Nodewise_Araw_age_noMotion_r(NOTsig_node_idx)=0;
	Araw_nodewise_AgeEffect_noMotion_partR_acrossThresh(:,T)=Nodewise_Araw_age_noMotion_r;


	%% WITHIN-MODULE CONNECTIVITY and GLOBAL EFFICIENCY
	
	for i = 1:nsub
		A=thresh_Araw_edgeVec(i,:)';
		curr_within=A(within_idx);
		% curr_within=nonzeros(curr_within);
		withinConn=mean(curr_within);
		withinConn_Araw_acrossThresh(i,T)=withinConn;

		% Eglob=efficiency_wei(squareform(A));
		% globEff_Araw_acrossThresh(i,T)=Eglob;
	end

	% Motion effects on withinMod
	[r,p]=partialcorr(withinConn_Araw_acrossThresh(:,T), motion, [prob_Araw_netStrength(:,T) age ageSq sex]);
	withinConn_Araw_motion_partialR_acrossThresh(T)=r;
	
	% Age effects on withinMod (controlling for motion)
	[r,p]=partialcorr(withinConn_Araw_acrossThresh(:,T), age, [prob_Araw_netStrength(:,T) motion ageSq sex]);
	withinConn_Araw_age_partialR_acrossThresh(T)=r;

	% Age effects on withinMod (controlling for motion)
	[r,p]=partialcorr(withinConn_Araw_acrossThresh(:,T), age, [prob_Araw_netStrength(:,T) ageSq sex]);
	withinConn_Araw_age_noMotion_partialR_acrossThresh(T)=r;

	% Motion effects on Global Efficiency
	[r,p]=partialcorr(globEff_Araw_acrossThresh(:,T), motion, [age ageSq sex]); % prob_Araw_netStrength(:,T)
	globEff_Araw_motion_partialR_acrossThresh(T)=r;
	
	% % Age effects on GlobEff (controlling for motion)
	% [r,p]=partialcorr(globEff_Araw_acrossThresh(:,T), age, [motion ageSq sex]); % prob_Araw_netStrength(:,T)
	% globEff_Araw_age_partialR_acrossThresh(T)=r;

	% % Age effects on GlobEff (controlling for motion)
	% [r,p]=partialcorr(globEff_Araw_acrossThresh(:,T), age, [ageSq sex]); % prob_Araw_netStrength(:,T)
	% globEff_Araw_age_noMotion_partialR_acrossThresh(T)=r;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROBABILISTIC: Streamline Count %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Coefficient of Variation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_mats=zeros(nreg,nreg,nsub);
	for i = 1:nsub
		A=squareform(Alength_edgeVec(i,:)); 
		group_mats(:, :,i) = A;
	end
[W_thr, Wcv] = GLB_threshold_consistency(group_mats, 1);
%% Set diagonal to zero
Wcv =(Wcv - diag(diag(Wcv)));
%% Set NaN to zero
isnan_idx=find(isnan(Wcv));
Wcv(isnan_idx)=0;
Alength_Wcv=squareform(Wcv)';

%% Node strength FDR thresh
orig_Alength_nodeStrength=zeros(nsub,nreg);
	for i = 1:nsub
		A=squareform(Alength_edgeVec(i,:));
		nodeStrength=sum(A,1);
		orig_Alength_nodeStrength(i,:)=nodeStrength;
	end

	[r,p]=partialcorr(orig_Alength_nodeStrength, motion,[age ageSq sex]) ;
	[orig_Nodewise_Alength_motion_pthr, Nodewise_Alength_motion_pcor, Nodewise_Alength_motion_padj] = fdr(p,0.05);
	orig_Nodewise_Alength_motion_pthr

%% Motion partial correlation
[Edgewise_Alength_motion_r,Edgewise_Alength_motion_p]=partialcorr(Alength_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_Alength_motion_r));
Edgewise_Alength_motion_p(isnan_idx)=[];
Edgewise_Alength_motion_r(isnan_idx)=[];
tmp_mean_length=mean_length; 
tmp_mean_length(isnan_idx)=[];
tmp_Alength_Wcv=Alength_Wcv; tmp_Alength_Wcv(isnan_idx)=[];

corr(Edgewise_Alength_motion_r,tmp_mean_length)
corr(Edgewise_Alength_motion_r,-log(tmp_Alength_Wcv))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export dataframe for 3rd-level scatterplots in R -- FIGURE 2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Alength_output_df=zeros(nedge,3);
Alength_output_df(:,1)=Edgewise_Alength_motion_r;
Alength_output_df(:,2)=-log(Alength_Wcv);
Alength_output_df(:,3)=mean_length;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/Alength_edgewise_motion_effects_n949.csv',Alength_output_df)

%% 3rd-level correlation plot
corr(-log(Alength_Wcv),Edgewise_Alength_motion_r) % Consistency
corr(Edgewise_Alength_motion_r,mean_length) % mean Streamline length

%% Scatterplot
c = linspace(1,max(mean_length),length(mean_length)); figure;scatter(-log(Alength_Wcv),Edgewise_Alength_motion_r,[],c); set(gcf,'color','white')

%% Calculate FDR threshold for Motion and Age Effects
[Edgewise_Alength_motion_r,Edgewise_Alength_motion_p]=partialcorr(Alength_edgeVec,motion,[age ageSq sex]);
[Edgewise_Alength_pthr,Edgewise_Alength_pcor,Edgewise_Alength_padj]=fdr(Edgewise_Alength_motion_p,0.05);
sig_edge_idx=find(Edgewise_Alength_motion_p < Edgewise_Alength_pthr);
sig_edges=Edgewise_Alength_motion_r(sig_edge_idx);
sigPos_edges=find(sig_edges > 0);size(sigPos_edges)
sigNeg_edges=find(sig_edges < 0); size(sigNeg_edges)
perc_Alength_sig_Motion = (length(sig_edge_idx) / length(Edgewise_Alength_motion_p)) * 100;

%% Age effect with motion
[Edgewise_Alength_age_r,Edgewise_Alength_age_p]=partialcorr(Alength_edgeVec, age,[ageSq motion sex]);
isnan_idx=find(isnan(Edgewise_Alength_age_p));
Edgewise_Alength_age_p(isnan_idx)=[];
Edgewise_Alength_age_r(isnan_idx)=[];
[Alength_age_pthr,Alength_age_pcor, Alength_age_padj] = fdr(Edgewise_Alength_age_p,0.05);
Age_sig_edge_idx=find(Edgewise_Alength_age_p < Alength_age_pthr);

[Nodewise_Alength_age_r,Nodewise_Alength_age_p]=partialcorr(orig_Alength_nodeStrength, age,[ageSq motion sex]);
isnan_idx=find(isnan(Nodewise_Alength_age_p));
Nodewise_Alength_age_p(isnan_idx)=[];
Nodewise_Alength_age_r(isnan_idx)=[];
[Alength_Node_age_pthr, Alength_Node_age_pcor, Alength_Node_age_padj] = fdr(Nodewise_Alength_age_p,0.05);
Age_sig_node_idx=find(Nodewise_Alength_age_p < Alength_Node_age_pthr);

%% Age Effect without motion
[Edgewise_Alength_age_noMotion_r, Edgewise_Alength_age_noMotion_p]=partialcorr(Alength_edgeVec,age,[ageSq sex]);
isnan_idx=find(isnan(Edgewise_Alength_age_noMotion_p));
Edgewise_Alength_age_noMotion_p(isnan_idx)=[];
Edgewise_Alength_age_noMotion_r(isnan_idx)=[];
[Alength_age_noMotion_pthr, Alength_age_noMotion_pcor,age_noMotionpadj] = fdr(Edgewise_Alength_age_noMotion_p,0.05);
Age_noMotion_sig_edge_idx=find(Edgewise_Alength_age_noMotion_p < Alength_age_noMotion_pthr);

[Nodewise_Alength_age_noMotion_r, Nodewise_Alength_age_noMotion_p]=partialcorr(orig_Alength_nodeStrength, age,[ageSq sex]);
% isnan_idx=find(isnan(Nodewise_Alength_age_p));
% Nodewise_Alength_age_noMotion_p(isnan_idx)=[];
% Nodewise_Alength_age_noMotion_r(isnan_idx)=[];
[Alength_Node_age_noMotion_pthr, Alength_Node_age_noMotion_pcor, Alength_Node_age_noMotion_padj] = fdr(Nodewise_Alength_age_noMotion_p, 0.05);
Age_noMotion_sig_node_idx=find(Nodewise_Alength_age_noMotion_p < Alength_Node_age_noMotion_pthr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOTION and AGE EFFECTS ACROSS CONSISTENCY THRESHOLDS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold_range=[1.0; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3; 0.2; 0.1];

%% Edge measures
W_Alength_acrossThresh=zeros(nedge,length(threshold_range)); % Group mean connectivity matrices after thresholding
Alength_edgeVec_acrossThresh=zeros(nsub,nedge,length(threshold_range));
Edgewise_Alength_motion_partR_acrThresh=zeros(nedge,length(threshold_range));
Edgewise_Alength_motion_pval_acrThresh=zeros(nedge,length(threshold_range));

%% Node Strength
prob_Alength_nodeStrength=zeros(nsub, nreg, length(threshold_range));
prob_Alength_nodeStrength_partR_acrThresh=zeros(nreg, length(threshold_range),1);
prob_Alength_nodeStrength_pval_acrThresh=zeros(nreg, length(threshold_range),1);
perc_Alength_nodeStrength_sig_Motion=zeros(length(threshold_range),1);

%% Network Strength - Motion Effects
prob_Alength_netStrength=zeros(nsub,length(threshold_range));
prob_Alength_netStrength_partR_acrThresh=zeros(length(threshold_range),1);
prob_Alength_netStrength_pval_acrThresh=zeros(length(threshold_range),1);
perc_Alength_MotionEffect_acrossThresh=zeros(length(threshold_range),1);

%% Network Strength - Age Effects
prob_Alength_netStrength_Age_partR_acrThresh=zeros(length(threshold_range),1);
prob_Alength_netStrength_Age_noMotion_partR_acrThresh=zeros(length(threshold_range),1);

% 3rd-level correlations
thirdLevelConsistency_Alength_r=zeros(length(threshold_range),1);
thirdLevelmeanLength_Alength_r=zeros(length(threshold_range),1);

% Edgewise Age Effects
perc_Alength_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Alength_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);

Edgewise_Alength_age_partialR_acrossThresh=zeros(nedge,length(threshold_range));
Edgewise_Alength_age_noMotion_partialR_acrossThresh=zeros(nedge,length(threshold_range));

% Nodewise Age Effects
perc_Alength_nodewise_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Alength_nodewise_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);
Alength_nodewise_AgeEffect_partR_acrossThresh=zeros(nreg, length(threshold_range),1);
Alength_nodewise_AgeEffect_noMotion_partR_acrossThresh=zeros(nreg, length(threshold_range),1);

%% Within-module connectivity
withinConn_Alength_acrossThresh=zeros(nsub,length(threshold_range));
withinConn_Alength_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Alength_age_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Alength_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Global Efficiency
globEff_Alength_acrossThresh=zeros(nsub,length(threshold_range));
globEff_Alength_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Alength_age_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Alength_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Distance-dependence of Age effects (3rd-level)
Alength_age_length_thirdLevel_r=zeros(length(threshold_range),1);
Alength_age_noMotion_length_thirdLevel_r=zeros(length(threshold_range),1);

Alength_age_consistency_thirdLevel_r=zeros(length(threshold_range),1);
Alength_age_noMotion_consistency_thirdLevel_r=zeros(length(threshold_range),1);

for T=1:length(threshold_range)
	%% Define group of Probabilistic connectivity matrices for thresholding
	thresh=threshold_range(T);
	% perc_out=thresh * 100
	thresh_Alength_edgeVec=Alength_edgeVec;
	[W_thr, Wcv] = GLB_threshold_consistency(group_mats, thresh);
	sq_W=squareform(W_thr)';
	W_Alength_acrossThresh(:,T)=sq_W;
	edgeThresh_index=find(W_Alength_acrossThresh(:,T) == 0)';
	thresh_Alength_edgeVec(:,[edgeThresh_index])=0;
	Alength_edgeVec_acrossThresh(:,:,T)=thresh_Alength_edgeVec;	
	% filename=strcat('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/consistencyThresh/ptx_streamlineCount_BreakspearThresh_',num2str(perc_out),'perc_AlengthSC_edgeVec_n949.txt')
	% dlmwrite(filename,thresh_Alength_edgeVec)


	%% Motion effects on Network Strength for each Consistency Threshold
	netStrength=sum(thresh_Alength_edgeVec,2);
	prob_Alength_netStrength(:,T)=netStrength;
	[r,p]=partialcorr(netStrength,motion,[age ageSq sex]);
	prob_Alength_netStrength_partR_acrThresh(T)=r;
	prob_Alength_netStrength_pval_acrThresh(T)=p;

	%% Age effects on Network Strength for each Consistency Threshold (with and without motion)
	
	% Controlling for Motion
	[r,p]=partialcorr(netStrength,age,[motion ageSq sex]);
	prob_Alength_netStrength_Age_partR_acrThresh(T)=r;
	% NOT controlling for motion
	[r,p]=partialcorr(netStrength,age,[ageSq sex]);
	prob_Alength_netStrength_Age_noMotion_partR_acrThresh(T)=r;

	%% Edgewise Motion Effects for each Consistency Threshold
	[Edgewise_Alength_motion_r,Edgewise_Alength_motion_p]=partialcorr(thresh_Alength_edgeVec,motion,[age ageSq sex]);
	Edgewise_Alength_motion_partR_acrThresh(:,T)=Edgewise_Alength_motion_r;
	Edgewise_Alength_motion_pval_acrThresh(:,T)=Edgewise_Alength_motion_p;


	%% Node strength effects
	for i = 1:nsub
		A=squareform(thresh_Alength_edgeVec(i,:));
		nodeStrength=sum(A,1);
		prob_Alength_nodeStrength(i,:,T)=nodeStrength;
	end

	[r,p]=partialcorr(prob_Alength_nodeStrength(:,:,T),motion, [age ageSq sex]);
	sig_node_idx=find(p < orig_Nodewise_Alength_motion_pthr);
	NOT_sig_node_idx=find(p > orig_Nodewise_Alength_motion_pthr);
	r(NOT_sig_node_idx)=0;
	prob_Alength_nodeStrength_partR_acrThresh(:,T)=r;
	perc_Alength_nodeStrength_sig_Motion(T)=(length(sig_node_idx) / length(p)) * 100;

	%% Percentage of Edgewise motion effects across thresholds
	[Edgewise_Alength_motion_r,Edgewise_Alength_motion_p]=partialcorr(thresh_Alength_edgeVec,motion,[age ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Alength_motion_p));
	Edgewise_Alength_motion_p(isnan_idx)=[];
	Edgewise_Alength_motion_r(isnan_idx)=[];
	tmp_consistency=Alength_Wcv;
	tmp_consistency(isnan_idx)=[];
	tmp_length=mean_length;
	tmp_length(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Alength_motion_p < Edgewise_Alength_pthr);
	perc_Alength_sig_Motion=(length(sig_edge_idx) / length(Edgewise_Alength_motion_p)) * 100;
	perc_Alength_MotionEffect_acrossThresh(T)=perc_Alength_sig_Motion;

	%% 3rd-level Consistency Effect
	[r,p]=corr(Edgewise_Alength_motion_r, -log(tmp_consistency));
	thirdLevelConsistency_r(T)=r;

	[r,p]=corr(Edgewise_Alength_motion_r, tmp_length);
	thirdLevelmeanLength_r(T)=r;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Edgewise Age effects with/without controlling for motion %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% WITH Motion: Save age effects across thresholds before removing NaNs
	[Edgewise_Alength_age_r, Edgewise_Alength_age_p] = partialcorr(thresh_Alength_edgeVec, age,[motion ageSq sex]);
	NOTsig_edge_idx=find(Edgewise_Alength_age_p >= Alength_age_pthr);
	isnan_idx=find(isnan(Edgewise_Alength_age_p));
	Edgewise_Alength_age_r(isnan_idx)=0;
	Edgewise_Alength_age_r(NOTsig_edge_idx)=0;
	Edgewise_Alength_age_partialR_acrossThresh(:,T)=Edgewise_Alength_age_r;
	%% Now remove NaNs to calculate percentage of sig edges
	[Edgewise_Alength_age_r, Edgewise_Alength_age_p]=partialcorr(thresh_Alength_edgeVec,age,[motion ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Alength_age_p));
	Edgewise_Alength_age_p(isnan_idx)=[];
	Edgewise_Alength_age_r(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Alength_age_p < Alength_age_pthr);
	perc_Alength_sig_Age=(length(sig_edge_idx) / length(Edgewise_Alength_age_p)) * 100;
	perc_Alength_AgeEffect_acrossThresh(T)=perc_Alength_sig_Age;	

	%% Distance-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Alength_age_r, tmp_length);
	Alength_age_length_thirdLevel_r(T)=r;

	%% Consistency-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Alength_age_r, -log(tmp_consistency));
	Alength_age_consistency_thirdLevel_r(T)=r;

	%% WITHOUT Motion: Save age effects across thresholds before removing NaNs
	[Edgewise_Alength_age_noMotion_r,Edgewise_Alength_age_noMotion_p]=partialcorr(thresh_Alength_edgeVec,age,[ageSq sex]);
	NOTsig_edge_idx=find(Edgewise_Alength_age_noMotion_p >= Alength_age_noMotion_pthr);
	isnan_idx=find(isnan(Edgewise_Alength_age_noMotion_p));
	Edgewise_Alength_age_noMotion_r(isnan_idx)=0;
	Edgewise_Alength_age_noMotion_partialR_acrossThresh(:,T)=Edgewise_Alength_age_noMotion_r;
	%% Now remove NaNs to calculate percentage of sig edges
	[Edgewise_Alength_age_noMotion_r,Edgewise_Alength_age_noMotion_p]=partialcorr(thresh_Alength_edgeVec,age,[ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Alength_age_noMotion_p));
	Edgewise_Alength_age_noMotion_p(isnan_idx)=[];
	Edgewise_Alength_age_noMotion_r(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Alength_age_noMotion_p < Alength_age_noMotion_pthr);
	perc_Alength_sig_Age_noMotion=(length(sig_edge_idx) / length(Edgewise_Alength_age_noMotion_p)) * 100;
	perc_Alength_AgeEffect_noMotion_acrossThresh(T)=perc_Alength_sig_Age_noMotion;	

	%% Distance-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Alength_age_noMotion_r, tmp_length);
	Alength_age_noMotion_length_thirdLevel_r(T)=r;

	%% Consistency-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Alength_age_noMotion_r, -log(tmp_consistency));
	Alength_age_noMotion_consistency_thirdLevel_r(T)=r;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Nodewise Age effects with/without controlling for motion %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[Nodewise_Alength_age_r, Nodewise_Alength_age_p]=partialcorr(prob_Alength_nodeStrength(:,:,T), age,[motion ageSq sex]);
	% isnan_idx=find(isnan(Nodewise_Alength_age_p));
	% Nodewise_Alength_age_p(isnan_idx)=[];
	% Nodewise_Alength_age_r(isnan_idx)=[];
	sig_node_idx=find(Nodewise_Alength_age_p < Alength_Node_age_pthr); 
	NOTsig_node_idx=find(Nodewise_Alength_age_p >= Alength_Node_age_pthr);
	perc_node_Alength_sig_Age=(length(sig_node_idx) / length(Nodewise_Alength_age_p)) * 100;
	perc_Alength_nodewise_AgeEffect_acrossThresh(T)=perc_node_Alength_sig_Age;
	Nodewise_Alength_age_r(NOTsig_node_idx)=0;
	Alength_nodewise_AgeEffect_partR_acrossThresh(:,T)=Nodewise_Alength_age_r;


	[Nodewise_Alength_age_noMotion_r, Nodewise_Alength_age_noMotion_p]=partialcorr(prob_Alength_nodeStrength(:,:,T), age,[ageSq sex]);
	% isnan_idx=find(isnan(Nodewise_Alength_age_noMotion_p));
	% Nodewise_Alength_age_noMotion_p(isnan_idx)=[];
	% Nodewise_Alength_age_noMotion_r(isnan_idx)=[];
	sig_node_idx=find(Nodewise_Alength_age_noMotion_p < Alength_Node_age_noMotion_pthr);
	NOTsig_node_idx=find(Nodewise_Alength_age_noMotion_p >= Alength_Node_age_noMotion_pthr);
	perc_node_Alength_sig_Age_noMotion=(length(sig_node_idx) / length(Nodewise_Alength_age_noMotion_p)) * 100;
	perc_Alength_nodewise_AgeEffect_noMotion_acrossThresh(T)= perc_node_Alength_sig_Age_noMotion;
	Nodewise_Alength_age_noMotion_r(NOTsig_node_idx)=0;
	Alength_nodewise_AgeEffect_noMotion_partR_acrossThresh(:,T)=Nodewise_Alength_age_noMotion_r;


	%% WITHIN-MODULE CONNECTIVITY and GLOBAL EFFICIENCY
	
	for i = 1:nsub
		A=thresh_Alength_edgeVec(i,:)';
		curr_within=A(within_idx);
		% curr_within=nonzeros(curr_within);
		withinConn=mean(curr_within);
		withinConn_Alength_acrossThresh(i,T)=withinConn;

		% Eglob=efficiency_wei(squareform(A));
		% globEff_Alength_acrossThresh(i,T)=Eglob;
	end

	% Motion effects on withinMod
	[r,p]=partialcorr(withinConn_Alength_acrossThresh(:,T), motion, [prob_Alength_netStrength(:,T) age ageSq sex]);
	withinConn_Alength_motion_partialR_acrossThresh(T)=r;
	
	% Age effects on withinMod (controlling for motion)
	[r,p]=partialcorr(withinConn_Alength_acrossThresh(:,T), age, [prob_Alength_netStrength(:,T) motion ageSq sex]);
	withinConn_Alength_age_partialR_acrossThresh(T)=r;

	% Age effects on withinMod (controlling for motion)
	[r,p]=partialcorr(withinConn_Alength_acrossThresh(:,T), age, [prob_Alength_netStrength(:,T) ageSq sex]);
	withinConn_Alength_age_noMotion_partialR_acrossThresh(T)=r;

	% Motion effects on Global Efficiency
	[r,p]=partialcorr(globEff_Alength_acrossThresh(:,T), motion, [age ageSq sex]); % prob_Alength_netStrength(:,T)
	globEff_Alength_motion_partialR_acrossThresh(T)=r;
	
	% % Age effects on GlobEff (controlling for motion)
	% [r,p]=partialcorr(globEff_Alength_acrossThresh(:,T), age, [motion ageSq sex]); % prob_Alength_netStrength(:,T)
	% globEff_Alength_age_partialR_acrossThresh(T)=r;

	% % Age effects on GlobEff (controlling for motion)
	% [r,p]=partialcorr(globEff_Alength_acrossThresh(:,T), age, [ageSq sex]); % prob_Alength_netStrength(:,T)
	% globEff_Alength_age_noMotion_partialR_acrossThresh(T)=r;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROBABILISTIC: Streamline Count %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Coefficient of Variation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_mats=zeros(nreg,nreg,nsub);
	for i = 1:nsub
		A=squareform(Aprop_edgeVec(i,:)); 
		group_mats(:, :,i) = A;
	end
[W_thr, Wcv] = GLB_threshold_consistency(group_mats, 1);
%% Set diagonal to zero
Wcv =(Wcv - diag(diag(Wcv)));
%% Set NaN to zero
isnan_idx=find(isnan(Wcv));
Wcv(isnan_idx)=0;
Aprop_Wcv=squareform(Wcv)';

%% Node strength FDR thresh
orig_Aprop_nodeStrength=zeros(nsub,nreg);
	for i = 1:nsub
		A=squareform(Aprop_edgeVec(i,:));
		nodeStrength=sum(A,1);
		orig_Aprop_nodeStrength(i,:)=nodeStrength;
	end

	[r,p]=partialcorr(orig_Aprop_nodeStrength, motion,[age ageSq sex]) ;
	[orig_Nodewise_Aprop_motion_pthr, Nodewise_Aprop_motion_pcor, Nodewise_Aprop_motion_padj] = fdr(p,0.05);
	orig_Nodewise_Aprop_motion_pthr

%% Motion partial correlation
[Edgewise_Aprop_motion_r,Edgewise_Aprop_motion_p]=partialcorr(Aprop_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_Aprop_motion_r));
Edgewise_Aprop_motion_p(isnan_idx)=[];
Edgewise_Aprop_motion_r(isnan_idx)=[];
tmp_mean_length=mean_length; 
tmp_mean_length(isnan_idx)=[];
tmp_Aprop_Wcv=Aprop_Wcv; tmp_Aprop_Wcv(isnan_idx)=[];

corr(Edgewise_Aprop_motion_r,tmp_mean_length)
corr(Edgewise_Aprop_motion_r,-log(tmp_Aprop_Wcv))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export dataframe for 3rd-level scatterplots in R -- FIGURE 2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aprop_output_df=zeros(nedge,3);
Aprop_output_df(:,1)=Edgewise_Aprop_motion_r;
Aprop_output_df(:,2)=-log(Aprop_Wcv);
Aprop_output_df(:,3)=mean_length;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/Aprop_edgewise_motion_effects_n949.csv',Aprop_output_df)

%% 3rd-level correlation plot
corr(-log(Aprop_Wcv),Edgewise_Aprop_motion_r) % Consistency
corr(Edgewise_Aprop_motion_r,mean_length) % mean Streamline length

%% Scatterplot
c = linspace(1,max(mean_length),length(mean_length)); figure;scatter(-log(Aprop_Wcv),Edgewise_Aprop_motion_r,[],c); set(gcf,'color','white')

%% Calculate FDR threshold for Motion and Age Effects
[Edgewise_Aprop_motion_r,Edgewise_Aprop_motion_p]=partialcorr(Aprop_edgeVec,motion,[age ageSq sex]);
[Edgewise_Aprop_pthr,Edgewise_Aprop_pcor,Edgewise_Aprop_padj]=fdr(Edgewise_Aprop_motion_p,0.05);
sig_edge_idx=find(Edgewise_Aprop_motion_p < Edgewise_Aprop_pthr);
sig_edges=Edgewise_Aprop_motion_r(sig_edge_idx);
sigPos_edges=find(sig_edges > 0);size(sigPos_edges)
sigNeg_edges=find(sig_edges < 0); size(sigNeg_edges)
perc_Aprop_sig_Motion = (length(sig_edge_idx) / length(Edgewise_Aprop_motion_p)) * 100;

%% Age effect with motion
[Edgewise_Aprop_age_r,Edgewise_Aprop_age_p]=partialcorr(Aprop_edgeVec, age,[ageSq motion sex]);
isnan_idx=find(isnan(Edgewise_Aprop_age_p));
Edgewise_Aprop_age_p(isnan_idx)=[];
Edgewise_Aprop_age_r(isnan_idx)=[];
[Aprop_age_pthr,Aprop_age_pcor, Aprop_age_padj] = fdr(Edgewise_Aprop_age_p,0.05);
Age_sig_edge_idx=find(Edgewise_Aprop_age_p < Aprop_age_pthr);

[Nodewise_Aprop_age_r,Nodewise_Aprop_age_p]=partialcorr(orig_Aprop_nodeStrength, age,[ageSq motion sex]);
isnan_idx=find(isnan(Nodewise_Aprop_age_p));
Nodewise_Aprop_age_p(isnan_idx)=[];
Nodewise_Aprop_age_r(isnan_idx)=[];
[Aprop_Node_age_pthr, Aprop_Node_age_pcor, Aprop_Node_age_padj] = fdr(Nodewise_Aprop_age_p,0.05);
Age_sig_node_idx=find(Nodewise_Aprop_age_p < Aprop_Node_age_pthr);

%% Age Effect without motion
[Edgewise_Aprop_age_noMotion_r, Edgewise_Aprop_age_noMotion_p]=partialcorr(Aprop_edgeVec,age,[ageSq sex]);
isnan_idx=find(isnan(Edgewise_Aprop_age_noMotion_p));
Edgewise_Aprop_age_noMotion_p(isnan_idx)=[];
Edgewise_Aprop_age_noMotion_r(isnan_idx)=[];
[Aprop_age_noMotion_pthr, Aprop_age_noMotion_pcor,age_noMotionpadj] = fdr(Edgewise_Aprop_age_noMotion_p,0.05);
Age_noMotion_sig_edge_idx=find(Edgewise_Aprop_age_noMotion_p < Aprop_age_noMotion_pthr);

[Nodewise_Aprop_age_noMotion_r, Nodewise_Aprop_age_noMotion_p]=partialcorr(orig_Aprop_nodeStrength, age,[ageSq sex]);
% isnan_idx=find(isnan(Nodewise_Aprop_age_p));
% Nodewise_Aprop_age_noMotion_p(isnan_idx)=[];
% Nodewise_Aprop_age_noMotion_r(isnan_idx)=[];
[Aprop_Node_age_noMotion_pthr, Aprop_Node_age_noMotion_pcor, Aprop_Node_age_noMotion_padj] = fdr(Nodewise_Aprop_age_noMotion_p, 0.05);
Age_noMotion_sig_node_idx=find(Nodewise_Aprop_age_noMotion_p < Aprop_Node_age_noMotion_pthr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOTION and AGE EFFECTS ACROSS CONSISTENCY THRESHOLDS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold_range=[1.0; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3; 0.2; 0.1];

%% Edge measures
W_Aprop_acrossThresh=zeros(nedge,length(threshold_range)); % Group mean connectivity matrices after thresholding
Aprop_edgeVec_acrossThresh=zeros(nsub,nedge,length(threshold_range));
Edgewise_Aprop_motion_partR_acrThresh=zeros(nedge,length(threshold_range));
Edgewise_Aprop_motion_pval_acrThresh=zeros(nedge,length(threshold_range));

%% Node Strength
prob_Aprop_nodeStrength=zeros(nsub, nreg, length(threshold_range));
prob_Aprop_nodeStrength_partR_acrThresh=zeros(nreg, length(threshold_range),1);
prob_Aprop_nodeStrength_pval_acrThresh=zeros(nreg, length(threshold_range),1);
perc_Aprop_nodeStrength_sig_Motion=zeros(length(threshold_range),1);

%% Network Strength - Motion Effects
prob_Aprop_netStrength=zeros(nsub,length(threshold_range));
prob_Aprop_netStrength_partR_acrThresh=zeros(length(threshold_range),1);
prob_Aprop_netStrength_pval_acrThresh=zeros(length(threshold_range),1);
perc_Aprop_MotionEffect_acrossThresh=zeros(length(threshold_range),1);

%% Network Strength - Age Effects
prob_Aprop_netStrength_Age_partR_acrThresh=zeros(length(threshold_range),1);
prob_Aprop_netStrength_Age_noMotion_partR_acrThresh=zeros(length(threshold_range),1);

% 3rd-level correlations
thirdLevelConsistency_Aprop_r=zeros(length(threshold_range),1);
thirdLevelmeanLength_Aprop_r=zeros(length(threshold_range),1);

% Edgewise Age Effects
perc_Aprop_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Aprop_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);

Edgewise_Aprop_age_partialR_acrossThresh=zeros(nedge,length(threshold_range));
Edgewise_Aprop_age_noMotion_partialR_acrossThresh=zeros(nedge,length(threshold_range));

% Nodewise Age Effects
perc_Aprop_nodewise_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Aprop_nodewise_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);
Aprop_nodewise_AgeEffect_partR_acrossThresh=zeros(nreg, length(threshold_range),1);
Aprop_nodewise_AgeEffect_noMotion_partR_acrossThresh=zeros(nreg, length(threshold_range),1);

%% Within-module connectivity
withinConn_Aprop_acrossThresh=zeros(nsub,length(threshold_range));
withinConn_Aprop_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Aprop_age_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Aprop_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Global Efficiency
globEff_Aprop_acrossThresh=zeros(nsub,length(threshold_range));
globEff_Aprop_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Aprop_age_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Aprop_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Distance-dependence of Age effects (3rd-level)
Aprop_age_length_thirdLevel_r=zeros(length(threshold_range),1);
Aprop_age_noMotion_length_thirdLevel_r=zeros(length(threshold_range),1);

Aprop_age_consistency_thirdLevel_r=zeros(length(threshold_range),1);
Aprop_age_noMotion_consistency_thirdLevel_r=zeros(length(threshold_range),1);

for T=1:length(threshold_range)
	%% Define group of Probabilistic connectivity matrices for thresholding
	thresh=threshold_range(T);
	% perc_out=thresh * 100
	thresh_Aprop_edgeVec=Aprop_edgeVec;
	[W_thr, Wcv] = GLB_threshold_consistency(group_mats, thresh);
	sq_W=squareform(W_thr)';
	W_Aprop_acrossThresh(:,T)=sq_W;
	edgeThresh_index=find(W_Aprop_acrossThresh(:,T) == 0)';
	thresh_Aprop_edgeVec(:,[edgeThresh_index])=0;
	Aprop_edgeVec_acrossThresh(:,:,T)=thresh_Aprop_edgeVec;	
	% filename=strcat('/data/joy/BBL/projects/pncBaumDti/Motion_paper/network_measures/n949/probabilistic/consistencyThresh/ptx_streamlineCount_BreakspearThresh_',num2str(perc_out),'perc_ApropSC_edgeVec_n949.txt')
	% dlmwrite(filename,thresh_Aprop_edgeVec)


	%% Motion effects on Network Strength for each Consistency Threshold
	netStrength=sum(thresh_Aprop_edgeVec,2);
	prob_Aprop_netStrength(:,T)=netStrength;
	[r,p]=partialcorr(netStrength,motion,[age ageSq sex]);
	prob_Aprop_netStrength_partR_acrThresh(T)=r;
	prob_Aprop_netStrength_pval_acrThresh(T)=p;

	%% Age effects on Network Strength for each Consistency Threshold (with and without motion)
	
	% Controlling for Motion
	[r,p]=partialcorr(netStrength,age,[motion ageSq sex]);
	prob_Aprop_netStrength_Age_partR_acrThresh(T)=r;
	% NOT controlling for motion
	[r,p]=partialcorr(netStrength,age,[ageSq sex]);
	prob_Aprop_netStrength_Age_noMotion_partR_acrThresh(T)=r;

	%% Edgewise Motion Effects for each Consistency Threshold
	[Edgewise_Aprop_motion_r,Edgewise_Aprop_motion_p]=partialcorr(thresh_Aprop_edgeVec,motion,[age ageSq sex]);
	Edgewise_Aprop_motion_partR_acrThresh(:,T)=Edgewise_Aprop_motion_r;
	Edgewise_Aprop_motion_pval_acrThresh(:,T)=Edgewise_Aprop_motion_p;


	%% Node strength effects
	for i = 1:nsub
		A=squareform(thresh_Aprop_edgeVec(i,:));
		nodeStrength=sum(A,1);
		prob_Aprop_nodeStrength(i,:,T)=nodeStrength;
	end

	[r,p]=partialcorr(prob_Aprop_nodeStrength(:,:,T),motion, [age ageSq sex]);
	sig_node_idx=find(p < orig_Nodewise_Aprop_motion_pthr);
	NOT_sig_node_idx=find(p > orig_Nodewise_Aprop_motion_pthr);
	r(NOT_sig_node_idx)=0;
	prob_Aprop_nodeStrength_partR_acrThresh(:,T)=r;
	perc_Aprop_nodeStrength_sig_Motion(T)=(length(sig_node_idx) / length(p)) * 100;

	%% Percentage of Edgewise motion effects across thresholds
	[Edgewise_Aprop_motion_r,Edgewise_Aprop_motion_p]=partialcorr(thresh_Aprop_edgeVec,motion,[age ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Aprop_motion_p));
	Edgewise_Aprop_motion_p(isnan_idx)=[];
	Edgewise_Aprop_motion_r(isnan_idx)=[];
	tmp_consistency=Aprop_Wcv;
	tmp_consistency(isnan_idx)=[];
	tmp_length=mean_length;
	tmp_length(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Aprop_motion_p < Edgewise_Aprop_pthr);
	perc_Aprop_sig_Motion=(length(sig_edge_idx) / length(Edgewise_Aprop_motion_p)) * 100;
	perc_Aprop_MotionEffect_acrossThresh(T)=perc_Aprop_sig_Motion;

	%% 3rd-level Consistency Effect
	[r,p]=corr(Edgewise_Aprop_motion_r, -log(tmp_consistency));
	thirdLevelConsistency_r(T)=r;

	[r,p]=corr(Edgewise_Aprop_motion_r, tmp_length);
	thirdLevelmeanLength_r(T)=r;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Edgewise Age effects with/without controlling for motion %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% WITH Motion: Save age effects across thresholds before removing NaNs
	[Edgewise_Aprop_age_r, Edgewise_Aprop_age_p] = partialcorr(thresh_Aprop_edgeVec, age,[motion ageSq sex]);
	NOTsig_edge_idx=find(Edgewise_Aprop_age_p >= Aprop_age_pthr);
	isnan_idx=find(isnan(Edgewise_Aprop_age_p));
	Edgewise_Aprop_age_r(isnan_idx)=0;
	Edgewise_Aprop_age_r(NOTsig_edge_idx)=0;
	Edgewise_Aprop_age_partialR_acrossThresh(:,T)=Edgewise_Aprop_age_r;
	%% Now remove NaNs to calculate percentage of sig edges
	[Edgewise_Aprop_age_r, Edgewise_Aprop_age_p]=partialcorr(thresh_Aprop_edgeVec,age,[motion ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Aprop_age_p));
	Edgewise_Aprop_age_p(isnan_idx)=[];
	Edgewise_Aprop_age_r(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Aprop_age_p < Aprop_age_pthr);
	perc_Aprop_sig_Age=(length(sig_edge_idx) / length(Edgewise_Aprop_age_p)) * 100;
	perc_Aprop_AgeEffect_acrossThresh(T)=perc_Aprop_sig_Age;	

	%% Distance-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Aprop_age_r, tmp_length);
	Aprop_age_length_thirdLevel_r(T)=r;

	%% Consistency-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Aprop_age_r, -log(tmp_consistency));
	Aprop_age_consistency_thirdLevel_r(T)=r;

	%% WITHOUT Motion: Save age effects across thresholds before removing NaNs
	[Edgewise_Aprop_age_noMotion_r,Edgewise_Aprop_age_noMotion_p]=partialcorr(thresh_Aprop_edgeVec,age,[ageSq sex]);
	NOTsig_edge_idx=find(Edgewise_Aprop_age_noMotion_p >= Aprop_age_noMotion_pthr);
	isnan_idx=find(isnan(Edgewise_Aprop_age_noMotion_p));
	Edgewise_Aprop_age_noMotion_r(isnan_idx)=0;
	Edgewise_Aprop_age_noMotion_partialR_acrossThresh(:,T)=Edgewise_Aprop_age_noMotion_r;
	%% Now remove NaNs to calculate percentage of sig edges
	[Edgewise_Aprop_age_noMotion_r,Edgewise_Aprop_age_noMotion_p]=partialcorr(thresh_Aprop_edgeVec,age,[ageSq sex]);
	isnan_idx=find(isnan(Edgewise_Aprop_age_noMotion_p));
	Edgewise_Aprop_age_noMotion_p(isnan_idx)=[];
	Edgewise_Aprop_age_noMotion_r(isnan_idx)=[];
	sig_edge_idx=find(Edgewise_Aprop_age_noMotion_p < Aprop_age_noMotion_pthr);
	perc_Aprop_sig_Age_noMotion=(length(sig_edge_idx) / length(Edgewise_Aprop_age_noMotion_p)) * 100;
	perc_Aprop_AgeEffect_noMotion_acrossThresh(T)=perc_Aprop_sig_Age_noMotion;	

	%% Distance-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Aprop_age_noMotion_r, tmp_length);
	Aprop_age_noMotion_length_thirdLevel_r(T)=r;

	%% Consistency-dependence of Edgewise Age effects (3rd-level)
	[r,p]=corr(Edgewise_Aprop_age_noMotion_r, -log(tmp_consistency));
	Aprop_age_noMotion_consistency_thirdLevel_r(T)=r;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Nodewise Age effects with/without controlling for motion %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[Nodewise_Aprop_age_r, Nodewise_Aprop_age_p]=partialcorr(prob_Aprop_nodeStrength(:,:,T), age,[motion ageSq sex]);
	% isnan_idx=find(isnan(Nodewise_Aprop_age_p));
	% Nodewise_Aprop_age_p(isnan_idx)=[];
	% Nodewise_Aprop_age_r(isnan_idx)=[];
	sig_node_idx=find(Nodewise_Aprop_age_p < Aprop_Node_age_pthr); 
	NOTsig_node_idx=find(Nodewise_Aprop_age_p >= Aprop_Node_age_pthr);
	perc_node_Aprop_sig_Age=(length(sig_node_idx) / length(Nodewise_Aprop_age_p)) * 100;
	perc_Aprop_nodewise_AgeEffect_acrossThresh(T)=perc_node_Aprop_sig_Age;
	Nodewise_Aprop_age_r(NOTsig_node_idx)=0;
	Aprop_nodewise_AgeEffect_partR_acrossThresh(:,T)=Nodewise_Aprop_age_r;


	[Nodewise_Aprop_age_noMotion_r, Nodewise_Aprop_age_noMotion_p]=partialcorr(prob_Aprop_nodeStrength(:,:,T), age,[ageSq sex]);
	% isnan_idx=find(isnan(Nodewise_Aprop_age_noMotion_p));
	% Nodewise_Aprop_age_noMotion_p(isnan_idx)=[];
	% Nodewise_Aprop_age_noMotion_r(isnan_idx)=[];
	sig_node_idx=find(Nodewise_Aprop_age_noMotion_p < Aprop_Node_age_noMotion_pthr);
	NOTsig_node_idx=find(Nodewise_Aprop_age_noMotion_p >= Aprop_Node_age_noMotion_pthr);
	perc_node_Aprop_sig_Age_noMotion=(length(sig_node_idx) / length(Nodewise_Aprop_age_noMotion_p)) * 100;
	perc_Aprop_nodewise_AgeEffect_noMotion_acrossThresh(T)= perc_node_Aprop_sig_Age_noMotion;
	Nodewise_Aprop_age_noMotion_r(NOTsig_node_idx)=0;
	Aprop_nodewise_AgeEffect_noMotion_partR_acrossThresh(:,T)=Nodewise_Aprop_age_noMotion_r;


	%% WITHIN-MODULE CONNECTIVITY and GLOBAL EFFICIENCY
	
	for i = 1:nsub
		A=thresh_Aprop_edgeVec(i,:)';
		curr_within=A(within_idx);
		% curr_within=nonzeros(curr_within);
		withinConn=mean(curr_within);
		withinConn_Aprop_acrossThresh(i,T)=withinConn;

		% Eglob=efficiency_wei(squareform(A));
		% globEff_Aprop_acrossThresh(i,T)=Eglob;
	end

	% Motion effects on withinMod
	[r,p]=partialcorr(withinConn_Aprop_acrossThresh(:,T), motion, [prob_Aprop_netStrength(:,T) age ageSq sex]);
	withinConn_Aprop_motion_partialR_acrossThresh(T)=r;
	
	% Age effects on withinMod (controlling for motion)
	[r,p]=partialcorr(withinConn_Aprop_acrossThresh(:,T), age, [prob_Aprop_netStrength(:,T) motion ageSq sex]);
	withinConn_Aprop_age_partialR_acrossThresh(T)=r;

	% Age effects on withinMod (controlling for motion)
	[r,p]=partialcorr(withinConn_Aprop_acrossThresh(:,T), age, [prob_Aprop_netStrength(:,T) ageSq sex]);
	withinConn_Aprop_age_noMotion_partialR_acrossThresh(T)=r;

	% Motion effects on Global Efficiency
	[r,p]=partialcorr(globEff_Aprop_acrossThresh(:,T), motion, [age ageSq sex]); % prob_Aprop_netStrength(:,T)
	globEff_Aprop_motion_partialR_acrossThresh(T)=r;
	
	% % Age effects on GlobEff (controlling for motion)
	% [r,p]=partialcorr(globEff_Aprop_acrossThresh(:,T), age, [motion ageSq sex]); % prob_Aprop_netStrength(:,T)
	% globEff_Aprop_age_partialR_acrossThresh(T)=r;

	% % Age effects on GlobEff (controlling for motion)
	% [r,p]=partialcorr(globEff_Aprop_acrossThresh(:,T), age, [ageSq sex]); % prob_Aprop_netStrength(:,T)
	% globEff_Aprop_age_noMotion_partialR_acrossThresh(T)=r;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPORT PROBABILISTIC DATA FOR FIGURE 5 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

perc=threshold_range * 100;

% Define X-axis for plots
x= 100 - perc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 4 - Consistency, Length, Age, and Motion Effects %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Araw
fig4_Araw_df = zeros(nedge, 5);
[Edgewise_Araw_age_r, Edgewise_Araw_age_p] = partialcorr(Araw_edgeVec, age, [motion ageSq sex]);
[Edgewise_Araw_motion_r, Edgewise_Araw_motion_p] = partialcorr(Araw_edgeVec, motion, [age ageSq sex]);
fig4_Araw_df(:,1)= Edgewise_Araw_age_r;
fig4_Araw_df(:,2)= Edgewise_Araw_motion_r;
fig4_Araw_df(:,3)= -log(Araw_Wcv);
fig4_Araw_df(:,4)= mean_length;
fig4_Araw_df(:,5)= median_length;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/Araw_edgewise_Age_Motion_Consistency_Length_effects_n949.csv', fig4_Araw_df)

% % Aprop
% fig4_Aprop_df = zeros(nedge, 4);
% [Edgewise_Aprop_age_r, Edgewise_Aprop_age_p] = partialcorr(Aprop_edgeVec, age, [motion ageSq sex]);
% [Edgewise_Aprop_motion_r, Edgewise_Aprop_motion_p] = partialcorr(Aprop_edgeVec, motion, [age ageSq sex]);
% fig4_Aprop_df(:,1)=Edgewise_Aprop_age_r;
% fig4_Aprop_df(:,2)=Edgewise_Aprop_motion_r;
% fig4_Aprop_df(:,3)=-log(Aprop_Wcv);
% fig4_Aprop_df(:,4)=mean_length;
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/Aprop_edgewise_Age_Motion_Consistency_Length_effects_n949.csv', fig4_Aprop_df)

% Alength
% fig4_Alength_df = zeros(nedge, 4);
% [Edgewise_Alength_age_r, Edgewise_Alength_age_p] = partialcorr(Alength_edgeVec, age, [motion ageSq sex]);
% [Edgewise_Alength_motion_r, Edgewise_Alength_motion_p] = partialcorr(Alength_edgeVec, motion, [age ageSq sex]);
% fig4_Alength_df(:,1)=Edgewise_Alength_age_r;
% fig4_Alength_df(:,2)=Edgewise_Alength_motion_r;
% fig4_Alength_df(:,3)=-log(Alength_Wcv);
% fig4_Alength_df(:,4)=mean_length;
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/Alength_edgewise_Age_Motion_Consistency_Length_effects_n949.csv', fig4_Alength_df)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Partial r threshold for Figure 4 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimate motion Effects
[Edgewise_Araw_motion_r, Edgewise_Araw_motion_p]=partialcorr(Araw_edgeVec, motion,[ageSq age sex]);
isnan_idx=find(isnan(Edgewise_Araw_motion_p));
Edgewise_Araw_motion_p(isnan_idx)=[];
Edgewise_Araw_motion_r(isnan_idx)=[];

%% FDR correction
[Araw_motion_pthr, Araw_motion_pcor, Araw_motion_padj] = fdr(Edgewise_Araw_motion_p, 0.05);
motion_sig_edge_idx=find(Edgewise_Araw_motion_p < Araw_motion_pthr);

%% Subset significant motion effects
sig_idx = find(Edgewise_Araw_motion_p < Araw_motion_pthr);
sig_Araw_motion_effects = Edgewise_Araw_motion_r(sig_idx);

%% Define partial r threshold
min(abs(sig_Araw_motion_effects))

%%%%%%%%%%%%%%%%%%
%%%%%% ADC %%%%%%%
%%%%%%%%%%%%%%%%%%

%% Estimate motion Effects
[Edgewise_inv_ADC_motion_r, Edgewise_inv_ADC_motion_p]=partialcorr(inv_ADC_edgeVec, motion,[ageSq age sex]);
isnan_idx=find(isnan(Edgewise_inv_ADC_motion_p));
Edgewise_inv_ADC_motion_p(isnan_idx)=[];
Edgewise_inv_ADC_motion_r(isnan_idx)=[];

%% FDR correction
[inv_ADC_motion_pthr, inv_ADC_motion_pcor, inv_ADC_motion_padj] = fdr(Edgewise_inv_ADC_motion_p, 0.05);
motion_sig_edge_idx=find(Edgewise_inv_ADC_motion_p < inv_ADC_motion_pthr);

%% Subset significant motion effects
sig_idx = find(Edgewise_inv_ADC_motion_p < inv_ADC_motion_pthr);
sig_inv_ADC_motion_effects = Edgewise_inv_ADC_motion_r(sig_idx);

%% Define partial r threshold
min(abs(sig_inv_ADC_motion_effects))

%%%%%%%%%%%%%%%%
%%%%%% FA %%%%%%
%%%%%%%%%%%%%%%%
%% Estimate motion Effects
[Edgewise_FA_motion_r, Edgewise_FA_motion_p]=partialcorr(FA_edgeVec, motion,[ageSq age sex]);
isnan_idx=find(isnan(Edgewise_FA_motion_p));
Edgewise_FA_motion_p(isnan_idx)=[];
Edgewise_FA_motion_r(isnan_idx)=[];

%% FDR correction
[FA_motion_pthr, FA_motion_pcor, FA_motion_padj] = fdr(Edgewise_FA_motion_p, 0.05);
motion_sig_edge_idx=find(Edgewise_FA_motion_p < FA_motion_pthr);

%% Subset significant motion effects
sig_idx = find(Edgewise_FA_motion_p < FA_motion_pthr);
sig_FA_motion_effects = Edgewise_FA_motion_r(sig_idx);

%% Define partial r threshold
min(abs(sig_FA_motion_effects))


%%%%%%%%%%%%%%%%
%%% FIGURE 5 %%%
%%%%%%%%%%%%%%%%

%% Percent Edges with Sig MOTION Effect 
fig5_Araw_df = zeros(length(threshold_range),2);
fig5_Araw_df(:,1) = x;
fig5_Araw_df(:,2) = perc_Araw_MotionEffect_acrossThresh;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Araw_percEdge_sigMotion_df', fig5_Araw_df)

fig5_Alength_df = zeros(length(threshold_range),2);
fig5_Alength_df(:,1) = x;
fig5_Alength_df(:,2) = perc_Alength_MotionEffect_acrossThresh;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Alength_percEdge_sigMotion_df', fig5_Alength_df)

% fig5_Aprop_df = zeros(length(threshold_range),2);
% fig5_Aprop_df(:,1) = x;
% fig5_Aprop_df(:,2) = perc_Aprop_MotionEffect_acrossThresh;
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Aprop_percEdge_sigMotion_df', fig5_Aprop_df)

%% Fig 5 -- Node strength at 50th percentile
Araw_fiftyPerc_effect = prob_Araw_nodeStrength_partR_acrThresh(:,6);
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Araw_sig_Motion_nodeStrength_50perc.txt', Araw_fiftyPerc_effect )

%% Fig 5 -- Motion effect on NetStrength across thresh
fig5_Araw_netStrength_df = zeros(length(threshold_range),2);
fig5_Araw_netStrength_df(:,1) = x;
fig5_Araw_netStrength_df(:,2) = prob_Araw_netStrength_partR_acrThresh;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Araw_sig_Motion_NetStrength_acrThresh.txt', fig5_Araw_netStrength_df)


% %% Motion effect on Node Strength 
% Aprop_noThresh_effect = prob_Aprop_nodeStrength_partR_acrThresh(:,1);
% Aprop_fiftyPerc_effect = prob_Aprop_nodeStrength_partR_acrThresh(:,6);
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/probabilistic/nodeStrength/fig5_ptx_Aprop_sig_Motion_nodeStrength_noThresh.txt', Aprop_noThresh_effect)
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/probabilistic/nodeStrength/fig5_ptx_Aprop_sig_Motion_nodeStrength_50perc.txt', Aprop_fiftyPerc_effect )

% Araw_noThresh_effect = prob_Araw_nodeStrength_partR_acrThresh(:,1);
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Araw_sig_Motion_nodeStrength_noThresh.txt', Araw_noThresh_effect)

% Alength_noThresh_effect = prob_Alength_nodeStrength_partR_acrThresh(:,1);
% Alength_fiftyPerc_effect = prob_Alength_nodeStrength_partR_acrThresh(:,6);
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/probabilistic/nodeStrength/fig5_ptx_Alength_sig_Motion_nodeStrength_noThresh.txt', Alength_noThresh_effect)
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/probabilistic/nodeStrength/fig5_ptx_Alength_sig_Motion_nodeStrength_50perc.txt', Alength_fiftyPerc_effect )


% fig5_Aprop_netStrength_df = zeros(length(threshold_range),2);
% fig5_Aprop_netStrength_df(:,1) = x;
% fig5_Aprop_netStrength_df(:,2) = prob_Aprop_netStrength_partR_acrThresh;
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Aprop_sig_Motion_NetStrength_acrThresh.txt', fig5_Aprop_netStrength_df)

% fig5_Alength_netStrength_df = zeros(length(threshold_range),2);
% fig5_Alength_netStrength_df(:,1) = x;
% fig5_Alength_netStrength_df(:,2) = prob_Alength_netStrength_partR_acrThresh;
% dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Alength_sig_Motion_NetStrength_acrThresh.txt', fig5_Alength_netStrength_df)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% DETERMINISTIC - FA %
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Breakspear Consistency for Deterministic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nreg=233
group_mats=zeros(nreg,nreg,nsub);
	for i = 1:nsub
		A=squareform(FA_edgeVec(i,:)); 
		group_mats(:, :,i) = A;
	end
[W_thr, Wcv] = GLB_threshold_consistency(group_mats, 1);
%% Set diagonal to zero
Wcv =(Wcv - diag(diag(Wcv)));
%% Set NaN to zero
isnan_idx=find(isnan(Wcv));
Wcv(isnan_idx)=0;
FA_Wcv=squareform(Wcv)';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Binary Consistency %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Binary_edgeVec=zeros(nsub,nedge);

for i = 1:nsub
	A=FA_edgeVec(i,:);
	edge_idx=find(A>0);
	A(edge_idx)=1; 
	Binary_edgeVec(i,:)=A;
end

bin_sum=sum(Binary_edgeVec)';
bin_consistency=bin_sum ./ 949;
bin_consistency=bin_consistency .* 100;

[Edgewise_FA_motion_r,Edgewise_FA_motion_p]=partialcorr(FA_edgeVec, motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_FA_motion_r))
Edgewise_FA_motion_r(isnan_idx)=[];
tmp_bin_consistency=bin_consistency;
tmp_bin_consistency(isnan_idx)=[];
tmp_mean_length=subj_mean_length
tmp_mean_length(isnan_idx)=[];
tmp_FA_Wcv=FA_Wcv;
tmp_FA_Wcv(isnan_idx)=[];

tmp_eucDist = eucDist_edgeVec;
tmp_eucDist(isnan_idx)=[];

corr(tmp_bin_consistency, Edgewise_FA_motion_r) % Consistency-dependence
corr(tmp_mean_length, Edgewise_FA_motion_r) % Length-dependence


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export Data for 3rd-level R plot (Figure 3) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FA_output_df=zeros(length(Edgewise_FA_motion_r),5);
FA_output_df(:,1)= Edgewise_FA_motion_r; % Edgewise motion effect
FA_output_df(:,2)= tmp_bin_consistency; % Binary edge consistency
FA_output_df(:,3)= -log(tmp_FA_Wcv); % Weighted-Breakspear consistency
FA_output_df(:,4)= tmp_mean_length; % mean streamline length (across subj)
FA_output_df(:,5)= tmp_eucDist; % mean Euclidean distance (MNI-Lausanne template)

dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_3/det_FA_edgewise_motion_effects_n949.csv', FA_output_df)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE FDR THRESHOLDS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate FDR threshold for Motion effects
[Edgewise_FA_motion_r,Edgewise_FA_motion_p]=partialcorr(FA_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_FA_motion_r))
Edgewise_FA_motion_r(isnan_idx)=[]; size(Edgewise_FA_motion_r)
Edgewise_FA_motion_p(isnan_idx)=[];
[Edgewise_FA_motion_pthr, Edgewise_FA_pcor,Edgewise_FA_padj]=fdr(Edgewise_FA_motion_p,0.05);
sig_edge_idx=find(Edgewise_FA_motion_p < Edgewise_FA_motion_pthr);
sig_edges=Edgewise_FA_motion_r(sig_edge_idx); size(sig_edges)
sigPos_edges=find(sig_edges > 0);size(sigPos_edges)
sigNeg_edges=find(sig_edges < 0); size(sigNeg_edges)
perc_FA_sig_Motion = (length(sig_edge_idx) / length(Edgewise_FA_motion_p)) * 100;

%% Calculate FDR threshold for Age Effects (with motion)
[Edgewise_FA_age_r,Edgewise_FA_age_p]=partialcorr(FA_edgeVec,age,[motion ageSq sex]);
Edgewise_FA_age_r(isnan_idx)=[];
Edgewise_FA_age_p(isnan_idx)=[];
[Edgewise_FA_age_pthr,Edgewise_FA_age_pcor,Edgewise_FA_age_padj] = fdr(Edgewise_FA_age_p,0.05);

%% Calculate FDR threshold for Age Effects (without motion)
[Edgewise_FA_age_noMotion_r, Edgewise_FA_age_noMotion_p]=partialcorr(FA_edgeVec,age,[ageSq sex]);
Edgewise_FA_age_noMotion_r(isnan_idx)=[];
Edgewise_FA_age_noMotion_p(isnan_idx)=[];
[Edgewise_FA_age_noMotion_pthr, Edgewise_FA_age_noMotion_pcor, Edgewise_FA_age_noMotion_padj] = fdr(Edgewise_FA_age_noMotion_p,0.05);

%% Calculate Percentage of Edges with Significant Motion Effect
sig_edge_idx=find(Edgewise_FA_motion_p < Edgewise_FA_motion_pthr);
perc_FA_sig_Motion= (length(sig_edge_idx) / length(Edgewise_FA_motion_p)) * 100

%%%%%%%%%%%%
%%% NODE %%%
%%%%%%%%%%%%
%% Node strength FDR thresh
orig_FA_nodeStrength=zeros(nsub,nreg);
	for i = 1:nsub
		A=squareform(FA_edgeVec(i,:));
		nodeStrength=sum(A,1);
		orig_FA_nodeStrength(i,:)=nodeStrength;
	end

% Motion
[r,p]=partialcorr(orig_FA_nodeStrength, motion,[age ageSq sex]) ;
[Nodewise_FA_motion_pthr, Nodewise_FA_motion_pcor, Nodewise_FA_motion_padj] = fdr(p,0.05);
Nodewise_FA_motion_pthr;
Motion_sig_node_idx=find(r < Nodewise_FA_motion_pthr);

% Age (with motion)
[Nodewise_FA_age_r, Nodewise_FA_age_p]=partialcorr(orig_FA_nodeStrength, age,[ageSq motion sex]);
[Node_FA_age_pthr, Node_FA_age_pcor, Node_FA_age_padj] = fdr(Nodewise_FA_age_p,0.05);
Age_sig_node_idx=find(Nodewise_FA_age_p < Nodewise_FA_age_p);

% Age (without motion)
[Nodewise_FA_age_noMotion_r, Nodewise_FA_age_noMotion_p]=partialcorr(orig_FA_nodeStrength, age,[ageSq sex]);
[Node_FA_age_noMotion_pthr, Node_FA_age_noMotion_pcor, Node_age_noMotion_padj] = fdr(Nodewise_FA_age_noMotion_p, 0.05);
Age_noMotion_sig_node_idx=find(Nodewise_FA_age_noMotion_p < Node_FA_age_noMotion_pthr);

tmp_length=subj_mean_length;
tmp_length(isnan_idx)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edgewise Motion Effects across Consistency Thresholds %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

det_consistency_percentiles = prctile(bin_consistency,[0 10 20 30 40 50 60 70 80 90])'

threshold_range =[0; 10; 20; 30; 40; 50; 60; 70; 80; 90; 100]
% threshold_range = det_consistency_percentiles;
nsub=949;
nedge=size(FA_edgeVec,2);

Binary_edgeVec=zeros(nsub,nedge);

for i = 1:nsub
	A=FA_edgeVec(i,:);
	edge_idx=find(A>0);
	A(edge_idx)=1; 
	Binary_edgeVec(i,:)=A;
end

bin_sum=sum(Binary_edgeVec)';
bin_consistency=bin_sum ./ 949;
bin_consistency=bin_consistency .* 100;

%%%%%%%%%%
% FA %
%%%%%%%%%%
[Edgewise_FA_motion_r,Edgewise_FA_motion_p]=partialcorr(FA_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_FA_motion_p));
Edgewise_FA_motion_p(isnan_idx)=[];
tmp_bin_consistency=bin_consistency;
tmp_bin_consistency(isnan_idx)=[];
tmp_length=subj_mean_length;
tmp_length(isnan_idx)=[];
% subj_mean_eucDist(isnan_idx)=[]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-allocate output measures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Edgewise_FA_motion_partR_acrThresh=zeros(nedge,length(threshold_range));
Edgewise_FA_motion_pval_acrThresh=zeros(nedge,length(threshold_range));
perc_FA_MotionEffect_acrossThresh=zeros(length(threshold_range),1);
perc_FA_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_FA_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);
thirdLevelbinConsistency_FA_r=zeros(length(threshold_range),1);
thirdLevelBreakspearConsistency_FA_r=zeros(length(threshold_range),1);
% FA_Modularity_across_thresh=zeros(nsub,length(threshold_range));
% FA_GlobalEfficiency_across_thresh=zeros(nsub,length(threshold_range));
FA_netStrength_acrossThresh=zeros(nsub,length(threshold_range));
FA_netStrength_motionEffect_acrossThresh_partialR=zeros(length(threshold_range),1);
FA_netStrength_motionEffect_acrossThresh_pval=zeros(length(threshold_range),1);

%% Node strength
FA_nodeStrength=zeros(nsub, nreg, length(threshold_range));
FA_nodeStrength_partR_acrThresh=zeros(nreg, length(threshold_range),1);
FA_nodeStrength_pval_acrThresh=zeros(nreg, length(threshold_range),1);
perc_FA_nodeStrength_sig_Motion=zeros(length(threshold_range),1);

%% Within-module connectivity
withinConn_FA_acrossThresh=zeros(nsub,length(threshold_range));
withinConn_FA_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_FA_age_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_FA_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Global Efficiency
globEff_FA_acrossThresh=zeros(nsub,length(threshold_range));
globEff_FA_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_FA_age_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_FA_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

for T=1:length(threshold_range)
	% Create subset of edges after thresholding for each subject
		thresh_group_FA_edges=FA_edgeVec;
		thresh_value=threshold_range(T);
		retain_edge_idx=find(bin_consistency >= thresh_value);
		thresh_idx=setdiff(1:numel(bin_consistency),retain_edge_idx);
		thresh_group_FA_edges(:,[thresh_idx])=0;
		% outThresh=thresh_value * 100;
		% tmp_mean_Length=subj_mean_length;
		% tmp_mean_Length(thresh_idx)=0;
		% mean_fibFA_across_thresh(:,T)=tmp_mean_Length;

		%% Percentage of Edgewise motion effects across thresholds
		[Edgewise_FA_motion_r,Edgewise_FA_motion_p]=partialcorr(thresh_group_FA_edges,motion,[age ageSq sex]);
		isnan_idx=find(isnan(Edgewise_FA_motion_p));
		Edgewise_FA_motion_p(isnan_idx)=[];
		Edgewise_FA_motion_r(isnan_idx)=[];
		tmp_bin_consistency=bin_consistency;
		tmp_bin_consistency(isnan_idx)=[];
		sig_edge_idx=find(Edgewise_FA_motion_p < Edgewise_FA_motion_pthr);
		perc_FA_sig_Motion=(length(sig_edge_idx) / length(Edgewise_FA_motion_p)) * 100;
		perc_FA_MotionEffect_acrossThresh(T)=perc_FA_sig_Motion;

		%% 3rd-level Consistency Effect
		[r,p]=corr(Edgewise_FA_motion_r,tmp_bin_consistency);
		thirdLevelbinConsistency_FA_r(T)=r;

		tmp_FA_Wcv=FA_Wcv;
		tmp_FA_Wcv(isnan_idx)=[];

		[r,p]=corr(Edgewise_FA_motion_r, -log(tmp_FA_Wcv));
		thirdLevelBreakspearConsistency_FA_r(T)=r;

		%% Network strength across thresholds
		netStrength=sum(thresh_group_FA_edges,2);
		FA_netStrength_acrossThresh(:,T)=netStrength;
		[r,p]=partialcorr(netStrength, motion,[age ageSq sex]);
		FA_netStrength_motionEffect_acrossThresh_partialR(T)=r;
		FA_netStrength_motionEffect_acrossThresh_pval(T)=p;

		%% Age effects with/without controlling for motion
		[Edgewise_FA_age_r, Edgewise_FA_age_p]=partialcorr(thresh_group_FA_edges, age,[motion ageSq sex]);
		isnan_idx=find(isnan(Edgewise_FA_age_p));
		Edgewise_FA_age_p(isnan_idx)=[];
		Edgewise_FA_age_r(isnan_idx)=[];
		% [age_pthr,age_pcor,age_padj] = fdr(Edgewise_FA_age_p,0.05);
		sig_edge_idx=find(Edgewise_FA_age_p < Edgewise_FA_age_pthr);
		perc_FA_sig_Age=(length(sig_edge_idx) / length(Edgewise_FA_age_p)) * 100;
		perc_FA_AgeEffect_acrossThresh(T)=perc_FA_sig_Age;

		[Edgewise_FA_age_noMotion_r, Edgewise_FA_age_noMotion_p]=partialcorr(thresh_group_FA_edges, age,[ageSq sex]);
		isnan_idx=find(isnan(Edgewise_FA_age_noMotion_p));
		Edgewise_FA_age_noMotion_p(isnan_idx)=[];
		Edgewise_FA_age_noMotion_r(isnan_idx)=[];
		% [age_pthr,age_pcor,age_padj] = fdr(Edgewise_FA_age_noMotion_p,0.05);
		sig_edge_idx=find(Edgewise_FA_age_noMotion_p < Edgewise_FA_age_noMotion_pthr);
		perc_FA_sig_Age_noMotion=(length(sig_edge_idx) / length(Edgewise_FA_age_noMotion_p)) * 100;
		perc_FA_AgeEffect_noMotion_acrossThresh(T)=perc_FA_sig_Age_noMotion;

	%%% Node Strength %%%
	for i = 1:nsub
		A=squareform(thresh_group_FA_edges(i,:));
		nodeStrength=sum(A,1);
		FA_nodeStrength(i,:,T)=nodeStrength;
	end

	[r,p]=partialcorr(FA_nodeStrength(:,:,T),motion, [age ageSq sex]);
	sig_node_idx=find(p < Nodewise_FA_motion_pthr);
	NOT_sig_node_idx=find(p >= Nodewise_FA_motion_pthr);
	r(NOT_sig_node_idx)=0;
	FA_nodeStrength_partR_acrThresh(:,T)=r;
	perc_FA_nodeStrength_sig_Motion(T)=(length(sig_node_idx) / length(p)) * 100;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXPORT DATA FOR FIGURES IN R %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%% FIGURE 4 - FA %%
%%%%%%%%%%%%%%%%%%%
fig4_FA_df = zeros(nedge, 6);
[Edgewise_FA_age_r, Edgewise_FA_age_p] = partialcorr(FA_edgeVec, age, [motion ageSq sex]);
[Edgewise_FA_motion_r, Edgewise_FA_motion_p] = partialcorr(FA_edgeVec, motion, [age ageSq sex]);
nan_idx= find(isnan(Edgewise_FA_motion_r));
tmp_eucDist=eucDist_edgeVec;
fig4_FA_df(:,1)= Edgewise_FA_age_r;
fig4_FA_df(:,2)= Edgewise_FA_motion_r;
fig4_FA_df(:,3)= bin_consistency
fig4_FA_df(:,4)= subj_mean_length;
fig4_FA_df(:,5)= subj_median_length;
fig4_FA_df(:,6)= tmp_eucDist

fig4_FA_df([nan_idx],:) = []; %% Remove edges where all subjects have 0
size(fig4_FA_df)
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/det_FA_edgewise_Age_Motion_Consistency_Length_effects_n949.csv', fig4_FA_df)

%%%%%%%%%%%%%%%%%%%
%% FIGURE 5 - FA %%
%%%%%%%%%%%%%%%%%%%
%% Percent edges with significant MOTION effects (across Thresh)
fig5_FA_df = zeros(length(threshold_range),2);
fig5_FA_df(:,1) = threshold_range;
fig5_FA_df(:,2) = perc_FA_MotionEffect_acrossThresh;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_FA_percEdge_sigMotion_df.csv', fig5_FA_df)

%% Motion effect on Node Strength (50th percentile consistency)
FA_fiftyPerc_effect = FA_nodeStrength_partR_acrThresh(:,6);
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_FA_sig_Motion_nodeStrength_50perc_consistency.txt', FA_fiftyPerc_effect )

%% Motion effect on Total Network Strength (across Thresh)
fig5_FA_netStrength_df = zeros(length(threshold_range),2);
fig5_FA_netStrength_df(:,1) = threshold_range;
fig5_FA_netStrength_df(:,2) = FA_netStrength_motionEffect_acrossThresh_partialR;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_FA_sig_Motion_NetStrength_acrThresh.txt', fig5_FA_netStrength_df)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINISTIC - inv_ADC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Breakspear Consistency for Deterministic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nreg=233
group_mats=zeros(nreg,nreg,nsub);
	for i = 1:nsub
		A=squareform(inv_ADC_edgeVec(i,:)); 
		group_mats(:, :,i) = A;
	end
[W_thr, Wcv] = GLB_threshold_consistency(group_mats, 1);
%% Set diagonal to zero
Wcv =(Wcv - diag(diag(Wcv)));
%% Set NaN to zero
isnan_idx=find(isnan(Wcv));
Wcv(isnan_idx)=0;
inv_ADC_Wcv=squareform(Wcv)';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Binary Consistency %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Binary_edgeVec=zeros(nsub,nedge);

for i = 1:nsub
	A=inv_ADC_edgeVec(i,:);
	edge_idx=find(A>0);
	A(edge_idx)=1; 
	Binary_edgeVec(i,:)=A;
end

bin_sum=sum(Binary_edgeVec)';
bin_consistency=bin_sum ./ 949;
bin_consistency=bin_consistency .* 100;

[Edgewise_inv_ADC_motion_r,Edgewise_inv_ADC_motion_p]=partialcorr(inv_ADC_edgeVec, motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_inv_ADC_motion_r))
Edgewise_inv_ADC_motion_r(isnan_idx)=[];
tmp_bin_consistency=bin_consistency;
tmp_bin_consistency(isnan_idx)=[];
tmp_mean_length=subj_mean_length
tmp_mean_length(isnan_idx)=[];
tmp_inv_ADC_Wcv=inv_ADC_Wcv;
tmp_inv_ADC_Wcv(isnan_idx)=[];

tmp_eucDist = eucDist_edgeVec;
tmp_eucDist(isnan_idx)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export Data for 3rd-level R plot (Figure 3) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inv_ADC_output_df=zeros(length(Edgewise_inv_ADC_motion_r),5);
inv_ADC_output_df(:,1)= Edgewise_inv_ADC_motion_r; % Edgewise motion effect
inv_ADC_output_df(:,2)= tmp_bin_consistency; % Binary edge consistency
inv_ADC_output_df(:,3)= -log(tmp_inv_ADC_Wcv); % Weighted-Breakspear consistency
inv_ADC_output_df(:,4)= tmp_mean_length; % mean streamline length (across subj)
inv_ADC_output_df(:,5)= tmp_eucDist; % mean Euclidean distance (MNI-Lausanne template)

dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_3/det_inv_ADC_edgewise_motion_effects_n949.csv', inv_ADC_output_df)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE FDR THRESHOLDS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate FDR threshold for Motion effects
[Edgewise_inv_ADC_motion_r,Edgewise_inv_ADC_motion_p]=partialcorr(inv_ADC_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_inv_ADC_motion_r))
Edgewise_inv_ADC_motion_r(isnan_idx)=[]; size(Edgewise_inv_ADC_motion_r)
Edgewise_inv_ADC_motion_p(isnan_idx)=[];
[Edgewise_inv_ADC_motion_pthr, Edgewise_inv_ADC_pcor,Edgewise_inv_ADC_padj]=fdr(Edgewise_inv_ADC_motion_p,0.05);
sig_edge_idx=find(Edgewise_inv_ADC_motion_p < Edgewise_inv_ADC_motion_pthr);
sig_edges=Edgewise_inv_ADC_motion_r(sig_edge_idx); size(sig_edges)
sigPos_edges=find(sig_edges > 0);size(sigPos_edges)
sigNeg_edges=find(sig_edges < 0); size(sigNeg_edges)
perc_inv_ADC_sig_Motion = (length(sig_edge_idx) / length(Edgewise_inv_ADC_motion_p)) * 100;

%% Calculate FDR threshold for Age Effects (with motion)
[Edgewise_inv_ADC_age_r,Edgewise_inv_ADC_age_p]=partialcorr(inv_ADC_edgeVec,age,[motion ageSq sex]);
Edgewise_inv_ADC_age_r(isnan_idx)=[];
Edgewise_inv_ADC_age_p(isnan_idx)=[];
[Edgewise_inv_ADC_age_pthr,Edgewise_inv_ADC_age_pcor,Edgewise_inv_ADC_age_padj] = fdr(Edgewise_inv_ADC_age_p,0.05);

%% Calculate FDR threshold for Age Effects (without motion)
[Edgewise_inv_ADC_age_noMotion_r, Edgewise_inv_ADC_age_noMotion_p]=partialcorr(inv_ADC_edgeVec,age,[ageSq sex]);
Edgewise_inv_ADC_age_noMotion_r(isnan_idx)=[];
Edgewise_inv_ADC_age_noMotion_p(isnan_idx)=[];
[Edgewise_inv_ADC_age_noMotion_pthr, Edgewise_inv_ADC_age_noMotion_pcor, Edgewise_inv_ADC_age_noMotion_padj] = fdr(Edgewise_inv_ADC_age_noMotion_p,0.05);

%% Calculate Percentage of Edges with Significant Motion Effect
sig_edge_idx=find(Edgewise_inv_ADC_motion_p < Edgewise_inv_ADC_motion_pthr);
perc_inv_ADC_sig_Motion= (length(sig_edge_idx) / length(Edgewise_inv_ADC_motion_p)) * 100

%%%%%%%%%%%%
%%% NODE %%%
%%%%%%%%%%%%
%% Node strength FDR thresh
orig_inv_ADC_nodeStrength=zeros(nsub,nreg);
	for i = 1:nsub
		A=squareform(inv_ADC_edgeVec(i,:));
		nodeStrength=sum(A,1);
		orig_inv_ADC_nodeStrength(i,:)=nodeStrength;
	end

% Motion
[r,p]=partialcorr(orig_inv_ADC_nodeStrength, motion,[age ageSq sex]) ;
[Nodewise_inv_ADC_motion_pthr, Nodewise_inv_ADC_motion_pcor, Nodewise_inv_ADC_motion_padj] = fdr(p,0.05);
Nodewise_inv_ADC_motion_pthr;
Motion_sig_node_idx=find(r < Nodewise_inv_ADC_motion_pthr);

% Age (with motion)
[Nodewise_inv_ADC_age_r, Nodewise_inv_ADC_age_p]=partialcorr(orig_inv_ADC_nodeStrength, age,[ageSq motion sex]);
[Node_inv_ADC_age_pthr, Node_inv_ADC_age_pcor, Node_inv_ADC_age_padj] = fdr(Nodewise_inv_ADC_age_p,0.05);
Age_sig_node_idx=find(Nodewise_inv_ADC_age_p < Nodewise_inv_ADC_age_p);

% Age (without motion)
[Nodewise_inv_ADC_age_noMotion_r, Nodewise_inv_ADC_age_noMotion_p]=partialcorr(orig_inv_ADC_nodeStrength, age,[ageSq sex]);
[Node_inv_ADC_age_noMotion_pthr, Node_inv_ADC_age_noMotion_pcor, Node_age_noMotion_padj] = fdr(Nodewise_inv_ADC_age_noMotion_p, 0.05);
Age_noMotion_sig_node_idx=find(Nodewise_inv_ADC_age_noMotion_p < Node_inv_ADC_age_noMotion_pthr);

tmp_length=subj_mean_length;
tmp_length(isnan_idx)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edgewise Motion Effects across Consistency Thresholds %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

det_consistency_percentiles = prctile(bin_consistency,[0 10 20 30 40 50 60 70 80 90])'

threshold_range =[0; 10; 20; 30; 40; 50; 60; 70; 80; 90; 100]
% threshold_range = det_consistency_percentiles;
nsub=949;
nedge=size(inv_ADC_edgeVec,2);

Binary_edgeVec=zeros(nsub,nedge);

for i = 1:nsub
	A=inv_ADC_edgeVec(i,:);
	edge_idx=find(A>0);
	A(edge_idx)=1; 
	Binary_edgeVec(i,:)=A;
end

bin_sum=sum(Binary_edgeVec)';
bin_consistency=bin_sum ./ 949;
bin_consistency=bin_consistency .* 100;

%%%%%%%%%%
% inv_ADC %
%%%%%%%%%%
[Edgewise_inv_ADC_motion_r,Edgewise_inv_ADC_motion_p]=partialcorr(inv_ADC_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_inv_ADC_motion_p));
Edgewise_inv_ADC_motion_p(isnan_idx)=[];
tmp_bin_consistency=bin_consistency;
tmp_bin_consistency(isnan_idx)=[];
tmp_length=subj_mean_length;
tmp_length(isnan_idx)=[];
% subj_mean_eucDist(isnan_idx)=[]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-allocate output measures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Edgewise_inv_ADC_motion_partR_acrThresh=zeros(nedge,length(threshold_range));
Edgewise_inv_ADC_motion_pval_acrThresh=zeros(nedge,length(threshold_range));
perc_inv_ADC_MotionEffect_acrossThresh=zeros(length(threshold_range),1);
perc_inv_ADC_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_inv_ADC_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);
thirdLevelbinConsistency_inv_ADC_r=zeros(length(threshold_range),1);
thirdLevelBreakspearConsistency_inv_ADC_r=zeros(length(threshold_range),1);
% inv_ADC_Modularity_across_thresh=zeros(nsub,length(threshold_range));
% inv_ADC_GlobalEfficiency_across_thresh=zeros(nsub,length(threshold_range));
inv_ADC_netStrength_acrossThresh=zeros(nsub,length(threshold_range));
inv_ADC_netStrength_motionEffect_acrossThresh_partialR=zeros(length(threshold_range),1);
inv_ADC_netStrength_motionEffect_acrossThresh_pval=zeros(length(threshold_range),1);

%% Node strength
inv_ADC_nodeStrength=zeros(nsub, nreg, length(threshold_range));
inv_ADC_nodeStrength_partR_acrThresh=zeros(nreg, length(threshold_range),1);
inv_ADC_nodeStrength_pval_acrThresh=zeros(nreg, length(threshold_range),1);
perc_inv_ADC_nodeStrength_sig_Motion=zeros(length(threshold_range),1);

%% Within-module connectivity
withinConn_inv_ADC_acrossThresh=zeros(nsub,length(threshold_range));
withinConn_inv_ADC_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_inv_ADC_age_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_inv_ADC_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Global Efficiency
globEff_inv_ADC_acrossThresh=zeros(nsub,length(threshold_range));
globEff_inv_ADC_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_inv_ADC_age_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_inv_ADC_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

for T=1:length(threshold_range)
	% Create subset of edges after thresholding for each subject
		thresh_group_inv_ADC_edges=inv_ADC_edgeVec;
		thresh_value=threshold_range(T);
		retain_edge_idx=find(bin_consistency >= thresh_value);
		thresh_idx=setdiff(1:numel(bin_consistency),retain_edge_idx);
		thresh_group_inv_ADC_edges(:,[thresh_idx])=0;
		% outThresh=thresh_value * 100;
		% tmp_mean_Length=subj_mean_length;
		% tmp_mean_Length(thresh_idx)=0;
		% mean_fibinv_ADC_across_thresh(:,T)=tmp_mean_Length;

		%% Percentage of Edgewise motion effects across thresholds
		[Edgewise_inv_ADC_motion_r,Edgewise_inv_ADC_motion_p]=partialcorr(thresh_group_inv_ADC_edges,motion,[age ageSq sex]);
		isnan_idx=find(isnan(Edgewise_inv_ADC_motion_p));
		Edgewise_inv_ADC_motion_p(isnan_idx)=[];
		Edgewise_inv_ADC_motion_r(isnan_idx)=[];
		tmp_bin_consistency=bin_consistency;
		tmp_bin_consistency(isnan_idx)=[];
		sig_edge_idx=find(Edgewise_inv_ADC_motion_p < Edgewise_inv_ADC_motion_pthr);
		perc_inv_ADC_sig_Motion=(length(sig_edge_idx) / length(Edgewise_inv_ADC_motion_p)) * 100;
		perc_inv_ADC_MotionEffect_acrossThresh(T)=perc_inv_ADC_sig_Motion;

		%% 3rd-level Consistency Effect
		[r,p]=corr(Edgewise_inv_ADC_motion_r,tmp_bin_consistency);
		thirdLevelbinConsistency_inv_ADC_r(T)=r;

		tmp_inv_ADC_Wcv=inv_ADC_Wcv;
		tmp_inv_ADC_Wcv(isnan_idx)=[];

		[r,p]=corr(Edgewise_inv_ADC_motion_r, -log(tmp_inv_ADC_Wcv));
		thirdLevelBreakspearConsistency_inv_ADC_r(T)=r;

		%% Network strength across thresholds
		netStrength=sum(thresh_group_inv_ADC_edges,2);
		inv_ADC_netStrength_acrossThresh(:,T)=netStrength;
		[r,p]=partialcorr(netStrength, motion,[age ageSq sex]);
		inv_ADC_netStrength_motionEffect_acrossThresh_partialR(T)=r;
		inv_ADC_netStrength_motionEffect_acrossThresh_pval(T)=p;

		%% Age effects with/without controlling for motion
		[Edgewise_inv_ADC_age_r, Edgewise_inv_ADC_age_p]=partialcorr(thresh_group_inv_ADC_edges, age,[motion ageSq sex]);
		isnan_idx=find(isnan(Edgewise_inv_ADC_age_p));
		Edgewise_inv_ADC_age_p(isnan_idx)=[];
		Edgewise_inv_ADC_age_r(isnan_idx)=[];
		% [age_pthr,age_pcor,age_padj] = fdr(Edgewise_inv_ADC_age_p,0.05);
		sig_edge_idx=find(Edgewise_inv_ADC_age_p < Edgewise_inv_ADC_age_pthr);
		perc_inv_ADC_sig_Age=(length(sig_edge_idx) / length(Edgewise_inv_ADC_age_p)) * 100;
		perc_inv_ADC_AgeEffect_acrossThresh(T)=perc_inv_ADC_sig_Age;

		[Edgewise_inv_ADC_age_noMotion_r, Edgewise_inv_ADC_age_noMotion_p]=partialcorr(thresh_group_inv_ADC_edges, age,[ageSq sex]);
		isnan_idx=find(isnan(Edgewise_inv_ADC_age_noMotion_p));
		Edgewise_inv_ADC_age_noMotion_p(isnan_idx)=[];
		Edgewise_inv_ADC_age_noMotion_r(isnan_idx)=[];
		% [age_pthr,age_pcor,age_padj] = fdr(Edgewise_inv_ADC_age_noMotion_p,0.05);
		sig_edge_idx=find(Edgewise_inv_ADC_age_noMotion_p < Edgewise_inv_ADC_age_noMotion_pthr);
		perc_inv_ADC_sig_Age_noMotion=(length(sig_edge_idx) / length(Edgewise_inv_ADC_age_noMotion_p)) * 100;
		perc_inv_ADC_AgeEffect_noMotion_acrossThresh(T)=perc_inv_ADC_sig_Age_noMotion;

	%%% Node Strength %%%
	for i = 1:nsub
		A=squareform(thresh_group_inv_ADC_edges(i,:));
		nodeStrength=sum(A,1);
		inv_ADC_nodeStrength(i,:,T)=nodeStrength;
	end

	[r,p]=partialcorr(inv_ADC_nodeStrength(:,:,T),motion, [age ageSq sex]);
	sig_node_idx=find(p < Nodewise_inv_ADC_motion_pthr);
	NOT_sig_node_idx=find(p >= Nodewise_inv_ADC_motion_pthr);
	r(NOT_sig_node_idx)=0;
	inv_ADC_nodeStrength_partR_acrThresh(:,T)=r;
	perc_inv_ADC_nodeStrength_sig_Motion(T)=(length(sig_node_idx) / length(p)) * 100;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXPORT DATA FOR FIGURES IN R %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%% FIGURE 4 - inv_ADC %%
%%%%%%%%%%%%%%%%%%%
fig4_inv_ADC_df = zeros(nedge, 6);
[Edgewise_inv_ADC_age_r, Edgewise_inv_ADC_age_p] = partialcorr(inv_ADC_edgeVec, age, [motion ageSq sex]);
[Edgewise_inv_ADC_motion_r, Edgewise_inv_ADC_motion_p] = partialcorr(inv_ADC_edgeVec, motion, [age ageSq sex]);
nan_idx= find(isnan(Edgewise_inv_ADC_motion_r));
tmp_eucDist=eucDist_edgeVec;
fig4_inv_ADC_df(:,1)= Edgewise_inv_ADC_age_r;
fig4_inv_ADC_df(:,2)= Edgewise_inv_ADC_motion_r;
fig4_inv_ADC_df(:,3)= bin_consistency
fig4_inv_ADC_df(:,4)= subj_mean_length;
fig4_inv_ADC_df(:,5)= subj_median_length;
fig4_inv_ADC_df(:,6)= tmp_eucDist

fig4_inv_ADC_df([nan_idx],:) = []; %% Remove edges where all subjects have 0
size(fig4_inv_ADC_df)
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/det_inv_ADC_edgewise_Age_Motion_Consistency_Length_effects_n949.csv', fig4_inv_ADC_df)

%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 5 - inv_ADC %%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Percent edges with significant MOTION effects (across Thresh)
fig5_inv_ADC_df = zeros(length(threshold_range),2);
fig5_inv_ADC_df(:,1) = threshold_range;
fig5_inv_ADC_df(:,2) = perc_inv_ADC_MotionEffect_acrossThresh;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_inv_ADC_percEdge_sigMotion_df.csv', fig5_inv_ADC_df)

%% Motion effect on Node Strength (50th percentile consistency)
inv_ADC_fiftyPerc_effect = inv_ADC_nodeStrength_partR_acrThresh(:,6);
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_inv_ADC_sig_Motion_nodeStrength_50perc_consistency.txt', inv_ADC_fiftyPerc_effect )

%% Motion effect on Total Network Strength (across Thresh)
fig5_inv_ADC_netStrength_df = zeros(length(threshold_range),2);
fig5_inv_ADC_netStrength_df(:,1) = threshold_range;
fig5_inv_ADC_netStrength_df(:,2) = inv_ADC_netStrength_motionEffect_acrossThresh_partialR;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_inv_ADC_sig_Motion_NetStrength_acrThresh.txt', fig5_inv_ADC_netStrength_df)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINISTIC - Length %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Breakspear Consistency for Deterministic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nreg=233
group_mats=zeros(nreg,nreg,nsub);
	for i = 1:nsub
		A=squareform(Length_edgeVec(i,:)); 
		group_mats(:, :,i) = A;
	end
[W_thr, Wcv] = GLB_threshold_consistency(group_mats, 1);
%% Set diagonal to zero
Wcv =(Wcv - diag(diag(Wcv)));
%% Set NaN to zero
isnan_idx=find(isnan(Wcv));
Wcv(isnan_idx)=0;
Length_Wcv=squareform(Wcv)';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Binary Consistency %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Binary_edgeVec=zeros(nsub,nedge);

for i = 1:nsub
	A=Length_edgeVec(i,:);
	edge_idx=find(A>0);
	A(edge_idx)=1; 
	Binary_edgeVec(i,:)=A;
end

bin_sum=sum(Binary_edgeVec)';
bin_consistency=bin_sum ./ 949;
bin_consistency=bin_consistency .* 100;

[Edgewise_Length_motion_r,Edgewise_Length_motion_p]=partialcorr(Length_edgeVec, motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_Length_motion_r))
Edgewise_Length_motion_r(isnan_idx)=[];
tmp_bin_consistency=bin_consistency;
tmp_bin_consistency(isnan_idx)=[];
tmp_mean_length=subj_mean_length
tmp_mean_length(isnan_idx)=[];
tmp_Length_Wcv=Length_Wcv;
tmp_Length_Wcv(isnan_idx)=[];

tmp_mean_FA = subj_mean_FA;
tmp_mean_FA(isnan_idx)=[];


tmp_eucDist = eucDist_edgeVec;
tmp_eucDist(isnan_idx)=[];


corr(tmp_bin_consistency, Edgewise_Length_motion_r)% Consistency-dependence
corr(tmp_mean_length, Edgewise_Length_motion_r) % Length-dependence
corr(tmp_mean_FA, Edgewise_Length_motion_r) % FA-dependence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export Data for 3rd-level R plot (Figure 3) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Length_output_df=zeros(length(Edgewise_Length_motion_r),6);
Length_output_df(:,1)= Edgewise_Length_motion_r; % Edgewise motion effect
Length_output_df(:,2)= tmp_bin_consistency; % Binary edge consistency
Length_output_df(:,3)= -log(tmp_Length_Wcv); % Weighted-Breakspear consistency
Length_output_df(:,4)= tmp_mean_length; % mean streamline length (across subj)
Length_output_df(:,5)= tmp_eucDist; % mean Euclidean distance (MNI-Lausanne template)
Length_output_df(:,6)= tmp_mean_FA; % mean FA (across subj)

dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_3/det_Length_edgewise_motion_effects_n949.csv', Length_output_df)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE FDR THRESHOLDS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate FDR threshold for Motion effects
[Edgewise_Length_motion_r,Edgewise_Length_motion_p]=partialcorr(Length_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_Length_motion_r))
Edgewise_Length_motion_r(isnan_idx)=[]; size(Edgewise_Length_motion_r)
Edgewise_Length_motion_p(isnan_idx)=[];
[Edgewise_Length_motion_pthr, Edgewise_Length_pcor,Edgewise_Length_padj]=fdr(Edgewise_Length_motion_p,0.05);
sig_edge_idx=find(Edgewise_Length_motion_p < Edgewise_Length_motion_pthr);
sig_edges=Edgewise_Length_motion_r(sig_edge_idx); size(sig_edges)
sigPos_edges=find(sig_edges > 0);size(sigPos_edges)
sigNeg_edges=find(sig_edges < 0); size(sigNeg_edges)
perc_Length_sig_Motion = (length(sig_edge_idx) / length(Edgewise_Length_motion_p)) * 100;

%% Calculate FDR threshold for Age Effects (with motion)
[Edgewise_Length_age_r,Edgewise_Length_age_p]=partialcorr(Length_edgeVec,age,[motion ageSq sex]);
Edgewise_Length_age_r(isnan_idx)=[];
Edgewise_Length_age_p(isnan_idx)=[];
[Edgewise_Length_age_pthr,Edgewise_Length_age_pcor,Edgewise_Length_age_padj] = fdr(Edgewise_Length_age_p,0.05);

%% Calculate FDR threshold for Age Effects (without motion)
[Edgewise_Length_age_noMotion_r, Edgewise_Length_age_noMotion_p]=partialcorr(Length_edgeVec,age,[ageSq sex]);
Edgewise_Length_age_noMotion_r(isnan_idx)=[];
Edgewise_Length_age_noMotion_p(isnan_idx)=[];
[Edgewise_Length_age_noMotion_pthr, Edgewise_Length_age_noMotion_pcor, Edgewise_Length_age_noMotion_padj] = fdr(Edgewise_Length_age_noMotion_p,0.05);

%% Calculate Percentage of Edges with Significant Motion Effect
sig_edge_idx=find(Edgewise_Length_motion_p < Edgewise_Length_motion_pthr);
perc_Length_sig_Motion= (length(sig_edge_idx) / length(Edgewise_Length_motion_p)) * 100

%%%%%%%%%%%%
%%% NODE %%%
%%%%%%%%%%%%
%% Node strength FDR thresh
orig_Length_nodeStrength=zeros(nsub,nreg);
	for i = 1:nsub
		A=squareform(Length_edgeVec(i,:));
		nodeStrength=sum(A,1);
		orig_Length_nodeStrength(i,:)=nodeStrength;
	end

% Motion
[r,p]=partialcorr(orig_Length_nodeStrength, motion,[age ageSq sex]) ;
[Nodewise_Length_motion_pthr, Nodewise_Length_motion_pcor, Nodewise_Length_motion_padj] = fdr(p,0.05);
Nodewise_Length_motion_pthr;
Motion_sig_node_idx=find(r < Nodewise_Length_motion_pthr);

% Age (with motion)
[Nodewise_Length_age_r, Nodewise_Length_age_p]=partialcorr(orig_Length_nodeStrength, age,[ageSq motion sex]);
[Node_Length_age_pthr, Node_Length_age_pcor, Node_Length_age_padj] = fdr(Nodewise_Length_age_p,0.05);
Age_sig_node_idx=find(Nodewise_Length_age_p < Nodewise_Length_age_p);

% Age (without motion)
[Nodewise_Length_age_noMotion_r, Nodewise_Length_age_noMotion_p]=partialcorr(orig_Length_nodeStrength, age,[ageSq sex]);
[Node_Length_age_noMotion_pthr, Node_Length_age_noMotion_pcor, Node_age_noMotion_padj] = fdr(Nodewise_Length_age_noMotion_p, 0.05);
Age_noMotion_sig_node_idx=find(Nodewise_Length_age_noMotion_p < Node_Length_age_noMotion_pthr);

tmp_length=subj_mean_length;
tmp_length(isnan_idx)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edgewise Motion Effects across Consistency Thresholds %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

det_consistency_percentiles = prctile(bin_consistency,[0 10 20 30 40 50 60 70 80 90])'

threshold_range =[0; 10; 20; 30; 40; 50; 60; 70; 80; 90; 100]
% threshold_range = det_consistency_percentiles;
nsub=949;
nedge=size(Length_edgeVec,2);

Binary_edgeVec=zeros(nsub,nedge);

for i = 1:nsub
	A=Length_edgeVec(i,:);
	edge_idx=find(A>0);
	A(edge_idx)=1; 
	Binary_edgeVec(i,:)=A;
end

bin_sum=sum(Binary_edgeVec)';
bin_consistency=bin_sum ./ 949;
bin_consistency=bin_consistency .* 100;

%%%%%%%%%%
% Length %
%%%%%%%%%%
[Edgewise_Length_motion_r,Edgewise_Length_motion_p]=partialcorr(Length_edgeVec,motion,[age ageSq sex]);
isnan_idx=find(isnan(Edgewise_Length_motion_p));
Edgewise_Length_motion_p(isnan_idx)=[];
tmp_bin_consistency=bin_consistency;
tmp_bin_consistency(isnan_idx)=[];
tmp_length=subj_mean_length;
tmp_length(isnan_idx)=[];
% subj_mean_eucDist(isnan_idx)=[]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-allocate output measures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Edgewise_Length_motion_partR_acrThresh=zeros(nedge,length(threshold_range));
Edgewise_Length_motion_pval_acrThresh=zeros(nedge,length(threshold_range));
perc_Length_MotionEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Length_AgeEffect_acrossThresh=zeros(length(threshold_range),1);
perc_Length_AgeEffect_noMotion_acrossThresh=zeros(length(threshold_range),1);
thirdLevelbinConsistency_Length_r=zeros(length(threshold_range),1);
thirdLevelBreakspearConsistency_Length_r=zeros(length(threshold_range),1);
% Length_Modularity_across_thresh=zeros(nsub,length(threshold_range));
% Length_GlobalEfficiency_across_thresh=zeros(nsub,length(threshold_range));
Length_netStrength_acrossThresh=zeros(nsub,length(threshold_range));
Length_netStrength_motionEffect_acrossThresh_partialR=zeros(length(threshold_range),1);
Length_netStrength_motionEffect_acrossThresh_pval=zeros(length(threshold_range),1);

%% Node strength
Length_nodeStrength=zeros(nsub, nreg, length(threshold_range));
Length_nodeStrength_partR_acrThresh=zeros(nreg, length(threshold_range),1);
Length_nodeStrength_pval_acrThresh=zeros(nreg, length(threshold_range),1);
perc_Length_nodeStrength_sig_Motion=zeros(length(threshold_range),1);

%% Within-module connectivity
withinConn_Length_acrossThresh=zeros(nsub,length(threshold_range));
withinConn_Length_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Length_age_partialR_acrossThresh=zeros(length(threshold_range),1);
withinConn_Length_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

%% Global Efficiency
globEff_Length_acrossThresh=zeros(nsub,length(threshold_range));
globEff_Length_motion_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Length_age_partialR_acrossThresh=zeros(length(threshold_range),1);
globEff_Length_age_noMotion_partialR_acrossThresh=zeros(length(threshold_range),1);

for T=1:length(threshold_range)
	% Create subset of edges after thresholding for each subject
		thresh_group_Length_edges=Length_edgeVec;
		thresh_value=threshold_range(T);
		retain_edge_idx=find(bin_consistency >= thresh_value);
		thresh_idx=setdiff(1:numel(bin_consistency),retain_edge_idx);
		thresh_group_Length_edges(:,[thresh_idx])=0;
		% outThresh=thresh_value * 100;
		% tmp_mean_Length=subj_mean_length;
		% tmp_mean_Length(thresh_idx)=0;
		% mean_fibLength_across_thresh(:,T)=tmp_mean_Length;

		%% Percentage of Edgewise motion effects across thresholds
		[Edgewise_Length_motion_r,Edgewise_Length_motion_p]=partialcorr(thresh_group_Length_edges,motion,[age ageSq sex]);
		isnan_idx=find(isnan(Edgewise_Length_motion_p));
		Edgewise_Length_motion_p(isnan_idx)=[];
		Edgewise_Length_motion_r(isnan_idx)=[];
		tmp_bin_consistency=bin_consistency;
		tmp_bin_consistency(isnan_idx)=[];
		sig_edge_idx=find(Edgewise_Length_motion_p < Edgewise_Length_motion_pthr);
		perc_Length_sig_Motion=(length(sig_edge_idx) / length(Edgewise_Length_motion_p)) * 100;
		perc_Length_MotionEffect_acrossThresh(T)=perc_Length_sig_Motion;

		%% 3rd-level Consistency Effect
		[r,p]=corr(Edgewise_Length_motion_r,tmp_bin_consistency);
		thirdLevelbinConsistency_Length_r(T)=r;

		tmp_Length_Wcv=Length_Wcv;
		tmp_Length_Wcv(isnan_idx)=[];

		[r,p]=corr(Edgewise_Length_motion_r, -log(tmp_Length_Wcv));
		thirdLevelBreakspearConsistency_Length_r(T)=r;

		%% Network strength across thresholds
		netStrength=sum(thresh_group_Length_edges,2);
		Length_netStrength_acrossThresh(:,T)=netStrength;
		[r,p]=partialcorr(netStrength, motion,[age ageSq sex]);
		Length_netStrength_motionEffect_acrossThresh_partialR(T)=r;
		Length_netStrength_motionEffect_acrossThresh_pval(T)=p;

		%% Age effects with/without controlling for motion
		[Edgewise_Length_age_r, Edgewise_Length_age_p]=partialcorr(thresh_group_Length_edges, age,[motion ageSq sex]);
		isnan_idx=find(isnan(Edgewise_Length_age_p));
		Edgewise_Length_age_p(isnan_idx)=[];
		Edgewise_Length_age_r(isnan_idx)=[];
		% [age_pthr,age_pcor,age_padj] = fdr(Edgewise_Length_age_p,0.05);
		sig_edge_idx=find(Edgewise_Length_age_p < Edgewise_Length_age_pthr);
		perc_Length_sig_Age=(length(sig_edge_idx) / length(Edgewise_Length_age_p)) * 100;
		perc_Length_AgeEffect_acrossThresh(T)=perc_Length_sig_Age;

		[Edgewise_Length_age_noMotion_r, Edgewise_Length_age_noMotion_p]=partialcorr(thresh_group_Length_edges, age,[ageSq sex]);
		isnan_idx=find(isnan(Edgewise_Length_age_noMotion_p));
		Edgewise_Length_age_noMotion_p(isnan_idx)=[];
		Edgewise_Length_age_noMotion_r(isnan_idx)=[];
		% [age_pthr,age_pcor,age_padj] = fdr(Edgewise_Length_age_noMotion_p,0.05);
		sig_edge_idx=find(Edgewise_Length_age_noMotion_p < Edgewise_Length_age_noMotion_pthr);
		perc_Length_sig_Age_noMotion=(length(sig_edge_idx) / length(Edgewise_Length_age_noMotion_p)) * 100;
		perc_Length_AgeEffect_noMotion_acrossThresh(T)=perc_Length_sig_Age_noMotion;

	%%% Node Strength %%%
	for i = 1:nsub
		A=squareform(thresh_group_Length_edges(i,:));
		nodeStrength=sum(A,1);
		Length_nodeStrength(i,:,T)=nodeStrength;
	end

	[r,p]=partialcorr(Length_nodeStrength(:,:,T),motion, [age ageSq sex]);
	sig_node_idx=find(p < Nodewise_Length_motion_pthr);
	NOT_sig_node_idx=find(p >= Nodewise_Length_motion_pthr);
	r(NOT_sig_node_idx)=0;
	Length_nodeStrength_partR_acrThresh(:,T)=r;
	perc_Length_nodeStrength_sig_Motion(T)=(length(sig_node_idx) / length(p)) * 100;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXPORT DATA FOR FIGURES IN R %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 4 - Length %%
%%%%%%%%%%%%%%%%%%%%%%%
fig4_Length_df = zeros(nedge, 6);
[Edgewise_Length_age_r, Edgewise_Length_age_p] = partialcorr(Length_edgeVec, age, [motion ageSq sex]);
[Edgewise_Length_motion_r, Edgewise_Length_motion_p] = partialcorr(Length_edgeVec, motion, [age ageSq sex]);
nan_idx= find(isnan(Edgewise_Length_motion_r));
tmp_eucDist=eucDist_edgeVec;
fig4_Length_df(:,1)= Edgewise_Length_age_r;
fig4_Length_df(:,2)= Edgewise_Length_motion_r;
fig4_Length_df(:,3)= bin_consistency
fig4_Length_df(:,4)= subj_mean_length;
fig4_Length_df(:,5)= subj_median_length;
fig4_Length_df(:,6)= tmp_eucDist

fig4_Length_df([nan_idx],:) = []; %% Remove edges where all subjects have 0
size(fig4_Length_df)
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/det_Length_edgewise_Age_Motion_Consistency_Length_effects_n949.csv', fig4_Length_df)

%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 5 - Length %%
%%%%%%%%%%%%%%%%%%%%%%%
%% Percent edges with significant MOTION effects (across Thresh)
fig5_Length_df = zeros(length(threshold_range),2);
fig5_Length_df(:,1) = threshold_range;
fig5_Length_df(:,2) = perc_Length_MotionEffect_acrossThresh;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_Length_percEdge_sigMotion_df.csv', fig5_Length_df)

%% Motion effect on Node Strength (50th percentile consistency)
Length_fiftyPerc_effect = Length_nodeStrength_partR_acrThresh(:,6);
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_Length_sig_Motion_nodeStrength_50perc_consistency.txt', Length_fiftyPerc_effect )

%% Motion effect on Total Network Strength (across Thresh)
fig5_Length_netStrength_df = zeros(length(threshold_range),2);
fig5_Length_netStrength_df(:,1) = threshold_range;
fig5_Length_netStrength_df(:,2) = Length_netStrength_motionEffect_acrossThresh_partialR;
dlmwrite('/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_Length_sig_Motion_NetStrength_acrThresh.txt', fig5_Length_netStrength_df)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('/data/joy/BBL/projects/pncBaumDti/Motion_paper/matlab_probabilistic_analysis_workspace_20170904.mat')

