###############################
### Load Relevant Libraries ###
###############################
require(ggplot2)
require(gtools)
require(Hmisc)
require(MASS)
require(mgcv)
require(ppcor)
require(stringr)
require(visreg)
require(parallel)
require(multilevel)
require(stats)
require(Formula)
require(mediation)
require(RColorBrewer)
require(colorspace)

#######################################################################################
### This script will construct a subject cohort for DTI-Motion analyses, retaining  ###
###   subjects who passed quality assurance protocols for T1, DTI, and ltnExcludev2 ###
#######################################################################################
# T1
t1_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv")
# DTI 
dti_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_dti_qa_20170301.csv")
# Resting-state QA
rest_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_RestQAData_20170714.csv")
# B0 Acquisition
protVal <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/n1601_pnc_protocol_validation_params_status_20161220.csv")
# LTN Status
health <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")
# Demographics
demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
# SES
envs <- read.csv("/data/joy/BBL/studies/pnc/n9498_dataFreeze/environment/n9498_go1_environment_factor_scores_tymoore_20150909.csv")
# Cognitive Scores
cognitive <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factorscores.csv")
# Clinical Score
clinical <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_clinical_factor_scores_20161212.csv")
# Subject identifier
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]

##############################
### Merge demographic data ###
##############################
demo <- merge(demo, health, by=c("bblid","scanid"))
demo <- merge(demo, t1_qa, by=c("bblid","scanid"))
demo <- merge(demo, dti_qa, by=c("bblid","scanid"))
demo <- merge(demo, rest_qa, by=c("bblid","scanid"))
demo <- merge(demo, protVal, by=c("bblid","scanid"))
demo <- merge(demo, clinical, by=c("bblid","scanid"))
demo <- merge(demo, cognitive, by=c("bblid","scanid"))
demo <- merge(demo, envs, by="bblid")

#########################
### DTI Motion Sample ###
#########################
Motion_df <- demo
Motion_df <- Motion_df[which(Motion_df$fsFinalExclude == 0), ]
Motion_df <- Motion_df[which(Motion_df$dti64ProtocolValidationStatus == 1), ]
Motion_df <- Motion_df[which(Motion_df$b0ProtocolValidationStatus==1), ]
Motion_df <-  Motion_df[which(Motion_df$ltnExcludev2==0), ]
Motion_df <- Motion_df[which(Motion_df$dti64Exclude == 0), ]
dim(Motion_df)

Motion_n949_df <- as.data.frame(Motion_df)
Motion_df <- NULL
Motion_n949_df$age <- Motion_n949_df$ageAtScan1 / 12

Motion_n949_df$sex <- as.factor(Motion_n949_df$sex)
Motion_n949_df$race <- as.factor(Motion_n949_df$race)
Motion_n949_df$race2 <- as.factor(Motion_n949_df$race2)

## Define quadratic and cubic age terms (demeaned)
Motion_n949_df$ageSq_demeaned <-(Motion_n949_df$age-mean(Motion_n949_df$age))^2
Motion_n949_df$ageCub_demeaned <-(Motion_n949_df$age-mean(Motion_n949_df$age))^3

hist(Motion_n949_df$age,col="navyblue",breaks=15)
axis(1, at=seq(8,23,by=1), labels=seq(8,23,by=1))

# hist(Motion_n949_df$dti64MeanRelRMS,col="navyblue")

# ageMotion_lm <- lm(age ~ dti64MeanRelRMS + sex, data=Motion_n949_df)
# visreg(ageMotion_lm, xlab="Age (years)", ylab="Mean Framewise Displacement (mm)")

###########################################
## Probabilistic Edgewise Motion Effects ##
###########################################

Araw_fig2_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/Araw_edgewise_motion_effects_n949.csv", header=FALSE)
colnames(Araw_fig2_df) <- c("motion_partial_r", "edge_consistency", "streamline_length")

# Aprop_fig2_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/edgewise_effects/Aprop_edgewise_motion_effects_n949.csv", header=FALSE)
# colnames(Aprop_fig2_df) <- c("motion_partial_r", "edge_consistency", "streamline_length")

Alength_fig2_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/Alength_edgewise_motion_effects_n949.csv", header=FALSE)
colnames(Alength_fig2_df) <- c("motion_partial_r", "edge_consistency", "streamline_length")


#################
### FIGURE 2B ###
#################
# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_2B_Araw.pdf", width=7, height=6)

# ggplot call with Edge Consistency on X-axis
Araw_plot_fig2B <- ggplot(data=Araw_fig2_df, aes(edge_consistency)) + geom_point(aes(y=motion_partial_r, col=streamline_length)) + scale_colour_gradientn(colours=c("lightpink", "orangered","red4","black")) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

Araw_plot_fig2B + ylim(-0.3, 0.4) + scale_x_continuous(name="Edge Consistency", breaks=c(-3, -2, -1, 0, 1, 2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/supplement/Araw_thirdLevel_length_suppFig.pdf", width=6, height=7)

# ggplot call with Streamline Length on X-axis
Araw_plot_fig2B <- ggplot(data=Araw_fig2_df, aes(streamline_length)) + geom_point(aes(y=motion_partial_r, alpha=0.8)) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

Araw_plot_fig2B + ylim(-0.3, 0.3) + scale_x_continuous(name="Streamline Length") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()

# 3rd-level Correlations with Consistency and Streamline Length
rcorr(Araw_fig2_df$edge_consistency, Araw_fig2_df$motion_partial_r)
rcorr(Araw_fig2_df$streamline_length, Araw_fig2_df$motion_partial_r)

#################
### FIGURE 2D ###
#################
# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_2D_Alength.pdf", width=7, height=6)

# ggplot call
Alength_plot_fig2D <- ggplot(data=Alength_fig2_df, aes(edge_consistency)) + geom_point(aes(y=motion_partial_r,col=streamline_length)) + scale_colour_gradientn(colours=c("lightpink", "orangered","red4","black")) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

Alength_plot_fig2D + ylim(-0.5, 0.2) + scale_x_continuous(name="Edge Consistency", breaks=c(-3, -2, -1, 0, 1, 2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/supplement/prob_length_motionEffects_consistency.pdf", width=6, height=7)

# ggplot call
Alength_plot_fig2D <- ggplot(data=Alength_fig2_df, aes(edge_consistency)) + geom_point(aes(y=motion_partial_r), alpha=0.8) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

Alength_plot_fig2D + ylim(-0.5, 0.2) + scale_x_continuous(name="Edge Consistency", breaks=c(-3, -2, -1, 0, 1, 2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()

# 3rd-level correlations
rcorr(Alength_fig2_df$edge_consistency, Alength_fig2_df$motion_partial_r)
rcorr(Alength_fig2_df$streamline_length, Alength_fig2_df$motion_partial_r)
#################

##############
## FIGURE 3 ##
##############

############
## invADC ##
############
# ADC_fig2_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_2/edgewise_effects/deterministic/det_inv_ADC_edgewise_motion_effects_n949.csv", header=FALSE)
# colnames(ADC_fig2_df) <- c("motion_partial_r", "binary_edge_consistency", "breakspear_edge_consistency", "streamline_length")

#################
### FIGURE 3B ###
#################
# Export pdf
# pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/supplement/SuppFig_det_invADC_edgewiseMotionEffects.pdf", width=7, height=6)

# ggplot call
# invADC_plot_fig3B <- ggplot(data= ADC_fig2_df, aes(binary_edge_consistency)) + geom_point(aes(y=motion_partial_r, col=streamline_length)) + scale_colour_gradientn(colours=c("lightpink", "orangered","red4","black")) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

# invADC_plot_fig3B + ylim(-0.3, 0.3) + scale_x_continuous(name="Edge Consistency", breaks=c(0, 25, 50, 75, 100)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

# dev.off()

# 3rd-level correlations
# rcorr(ADC_fig2_df$binary_edge_consistency, ADC_fig2_df$motion_partial_r)
# rcorr(ADC_fig2_df$streamline_length, ADC_fig2_df$motion_partial_r)


#############################################
## Figure 3B -- FA edgewise motion effects ##
#############################################
FA_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_3/det_FA_edgewise_motion_effects_n949.csv", header=FALSE)
colnames(FA_df) <- c("motion_partial_r", "binary_edge_consistency", "breakspear_edge_consistency", "streamline_length", "eucDist")

FA_consistencyCorr <- rcorr(FA_df$binary_edge_consistency, FA_df$motion_partial_r)
FA_consistencyCorr$r

FA_lengthCorr <- rcorr(FA_df$streamline_length, FA_df$motion_partial_r)
FA_lengthCorr$r

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_3B_det_FA.pdf", width=7, height=6)

# Consistency on X-axis
FA_plot1 <- ggplot(data= FA_df, aes(binary_edge_consistency)) + geom_point(aes(y=motion_partial_r, col=streamline_length)) + scale_colour_gradientn(colours=c("lightpink", "orangered","red4","black")) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

FA_plot1 + scale_x_continuous(name="Edge Consistency", breaks=c(0, 25, 50, 75, 100)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/supplement/det_FA_thirdLevel_length_suppFig.pdf", width=6, height=7)

# ggplot call with Streamline Length on X-axis
FA_plot2 <- ggplot(data= FA_df, aes(streamline_length)) + scale_colour_gradientn(colours=c("lightpink", "orangered","red4","black")) + geom_point(aes(y=motion_partial_r, col=streamline_length)) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

FA_plot2 + scale_x_continuous(name="Streamline Length") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()


############
## Length ##
############
Length_fig3_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_3/det_Length_edgewise_motion_effects_n949.csv", header=FALSE)

colnames(Length_fig3_df) <- c("motion_partial_r", "binary_edge_consistency", "breakspear_edge_consistency", "streamline_length", "eucDist", "mean_FA")

#################
### FIGURE 3D ###
#################
# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_3D_det_meanStreamlineLength.pdf", width=7, height=6)

# ggplot call
detLength_plot_fig3D <- ggplot(data= Length_fig3_df, aes(binary_edge_consistency)) + geom_point(aes(y=motion_partial_r, col= mean_FA)) + scale_colour_gradientn(colours=c("lightpink", "orangered","red4","black")) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

detLength_plot_fig3D  + scale_x_continuous(name="Edge Consistency", breaks=c(0, 25, 50, 75, 100)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) 

# + guides(fill = guide_colorbar(ticks = FALSE))

# + ylim(-0.3, 0.3)
dev.off()

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/supplement/det_length_motionEffects_consistency.pdf", width=6, height=7)

# ggplot call
detLength_plot_fig3D <- ggplot(data= Length_fig3_df, aes(binary_edge_consistency)) + geom_point(aes(y=motion_partial_r, col= mean_FA)) + scale_colour_gradientn(colours=c("lightpink", "orangered","red4","black")) + geom_smooth(method="lm",size=2, aes(y = motion_partial_r))

detLength_plot_fig3D  + scale_x_continuous(name="Edge Consistency", breaks=c(0, 25, 50, 75, 100)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) 

# + guides(fill = guide_colorbar(ticks = FALSE))

# + ylim(-0.3, 0.3)
dev.off()

# 3rd-level correlations
rcorr(Length_fig3_df$binary_edge_consistency,Length_fig3_df$motion_partial_r)
rcorr(Length_fig3_df$streamline_length, Length_fig3_df$motion_partial_r)
rcorr(Length_fig3_df$mean_FA, Length_fig3_df$motion_partial_r)

#################

#################
### Figure 4A ###
#################
fig4_Araw_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/Araw_edgewise_Age_Motion_Consistency_Length_effects_n949.csv", header=FALSE)
colnames(fig4_Araw_df) <- c("age_r", "motion_r", "edge_consistency", "streamline_length")

fdr_thresh <- 0.0892

tmp_pal <- brewer.pal(n=10, name="RdBu")
tmp_pal <- rev(tmp_pal)

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_4A_prob_Araw.pdf", width=7, height=6)

Araw_consistency_length_plot <- ggplot(data= subset(fig4_Araw_df, abs(fig4_Araw_df$motion_r) > fdr_thresh), aes(edge_consistency)) + geom_point(aes(y= streamline_length, col=motion_r)) + scale_colour_gradientn(colours=c(tmp_pal), guide="colorbar", limits=c(-.4, .4)) + geom_smooth(aes(y = streamline_length), color="black", method = "gam", formula = y ~ s(x, k=6), size=2, se=FALSE)
Araw_consistency_length_plot + scale_x_continuous(name="Edge Consistency", breaks=c(-3, -2, -1, 0, 1, 2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))
# + ylim(0,160)

dev.off()

#################
### Figure 4B ###
#################
fig4_FA_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_4/det_FA_edgewise_Age_Motion_Consistency_Length_effects_n949.csv", header=FALSE)

#fig4_invADC_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/det_invADC_edgewise_Age_Motion_Consistency_Length_effects_n949.csv", header=FALSE)

colnames(fig4_FA_df) <- c("age_r", "motion_r", "edge_consistency", "mean_streamline_length", "median_streamline_length", "eucDist")

# fdr_thresh <- 0.0883

# Partial R threshold (based on FDR correction)
fdr_thresh <- 0.0879
tmp_pal <- brewer.pal(n=10, name="RdBu")
tmp_pal <- rev(tmp_pal)

pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_4B_det_FA.pdf", width=7, height=6)

fig4_FA_df$log_consistency <- log(fig4_FA_df$edge_consistency)
fig4_FA_df$arctanconsistency <- atan(fig4_FA_df$edge_consistency)

FA_consistency_meanLength_plot <- ggplot(data= subset(fig4_FA_df, abs(fig4_FA_df$motion_r) > fdr_thresh) , aes(edge_consistency)) + geom_point(aes(y= mean_streamline_length, col=motion_r)) + scale_colour_gradientn(colours=c(tmp_pal), guide="colorbar", limits=c(-.4, .4)) + geom_smooth(aes(y = mean_streamline_length), color="black", method = "gam", formula = y ~ s(x, k=3), size=2, se=FALSE)

FA_consistency_meanLength_plot + scale_x_continuous(name="Edge Consistency") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()

#####################################################
### FIGURE 5A - Prob Motion effects across thresh ###
#####################################################

## Probabilistic Streamline Count
fig5_Araw_edges_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Araw_percEdge_sigMotion_df", header=FALSE)
colnames(fig5_Araw_edges_df) <- c("thresh_percentile","perc_sig_edges")

## Read in bootstrapped results
Araw_boot_percEdge_SD <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/bootstrapping/Araw_percEdge_boot_SD.txt", header=FALSE)

colnames(Araw_boot_percEdge_SD) <- "Araw_boot_percEdge_SD"
Araw_boot_percEdge_SD$SE <- Araw_boot_percEdge_SD$Araw_boot_percEdge_SD / 10

fig5_Araw_edges_df <- cbind(fig5_Araw_edges_df, Araw_boot_percEdge_SD)

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_5A_Araw_percSigEdges.pdf", width=7, height=6)

## Percent EDGES plots
Araw_plot <- ggplot(data= fig5_Araw_edges_df, aes(thresh_percentile)) +  geom_point(aes(y= perc_sig_edges), color="firebrick", size=6)

# geom_smooth(aes(y = perc_sig_edges), color="red4", method = "gam", formula = y ~ s(x, k=3), size=2, se=FALSE) # 
Araw_plot  + scale_y_continuous(breaks=c(10, 20, 30, 40)) + scale_x_continuous(breaks=c(0, 25, 50, 75)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + geom_errorbar(aes(ymin=perc_sig_edges - Araw_boot_percEdge_SD, ymax=perc_sig_edges + Araw_boot_percEdge_SD), colour="black", width=1, size=0.8) 

dev.off()

####################################################
### FIGURE 5D - Det Motion effects across thresh ###
####################################################

## Deterministic
fig5_FA_edges_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_FA_percEdge_sigMotion_df.csv", header=FALSE)
colnames(fig5_FA_edges_df) <- c("thresh_percentile","perc_sig_edges")

fig5_FA_edges_df <- fig5_FA_edges_df[-11,]

## Read in bootstrapped results
FA_boot_percEdge_SD <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/bootstrapping/FA_percEdge_boot_SD.txt", header=FALSE)

FA_boot_percEdge_SD <- FA_boot_percEdge_SD[-11,]

fig5_FA_edges_df <- cbind(fig5_FA_edges_df, FA_boot_percEdge_SD)

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_5D_FA_percSigEdges.pdf", width=7, height=6)

# FA_plot <- ggplot(data=fig5_FA_edges_df, aes(thresh_percentile)) + geom_smooth(aes(y = perc_sig_edges), color="red4", method = "gam", formula = y ~ s(x,k=6), size=2, se=FALSE) # geom_point(aes(y=perc_sig_edges), size=1)

FA_plot <- ggplot(data=fig5_FA_edges_df, aes(thresh_percentile)) + geom_point(aes(y = perc_sig_edges), color="firebrick", size=6) # geom_point(aes(y=perc_sig_edges), size=1)

FA_plot  + scale_x_continuous(breaks=c(0, 25, 50, 75)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +  theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_y_continuous(limit=c(10,70),breaks=c(20, 40, 60)) + geom_errorbar(aes(ymin=perc_sig_edges - FA_boot_percEdge_SD, ymax=perc_sig_edges + FA_boot_percEdge_SD), colour="black", width=1, size=0.8) 

dev.off()

#################
### FIGURE 5C ###
#################
fig5_Araw_netStrength_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/probabilistic/fig5_ptx_Araw_sig_Motion_NetStrength_acrThresh.txt", header=FALSE)
colnames(fig5_Araw_netStrength_df) <- c("thresh_percentile","motion_r")

## Read in bootstrapped results
Araw_boot_netStrength_SD <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/bootstrapping/Araw_netStrength_boot_SD.txt", header=FALSE)

fig5_Araw_netStrength_df <- cbind(fig5_Araw_netStrength_df, Araw_boot_netStrength_SD)

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_5C_Araw_netStrength.pdf", width=7, height=6)

Araw_netStrength_plot <- ggplot(data=fig5_Araw_netStrength_df, aes(thresh_percentile)) + geom_point(aes(y = motion_r), color="firebrick", size=6) # geom_point(aes(y=perc_sig_edges), size=1)

Araw_netStrength_plot  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_x_continuous(breaks=c(0, 25, 50, 75)) + scale_y_continuous(breaks=c(-0.1, -0.2, -0.3, -0.4)) +  theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + geom_errorbar(aes(ymin= motion_r - Araw_boot_netStrength_SD, ymax= motion_r + Araw_boot_netStrength_SD), colour="black", width=1, size=0.8)  + ylim(-0.4, -0.2)

dev.off()

#################
### FIGURE 5E ###
#################
fig5_FA_netStrength_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_5/deterministic/fig5_det_FA_sig_Motion_NetStrength_acrThresh.txt", header=FALSE)
colnames(fig5_FA_netStrength_df) <- c("thresh_percentile","motion_r")


dim(fig5_FA_netStrength_df)


## Read in bootstrapped results
FA_boot_netStrength_SD <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/bootstrapping/FA_netStrength_boot_SD.txt", header=FALSE)
colnames(FA_boot_netStrength_SD) <- "FA_boot_netStrength_SD"

# FA_boot_netStrength_SD <- FA_boot_netStrength_SD[-11,]

# Merge data for making error bars
fig5_FA_netStrength_df <- cbind(fig5_FA_netStrength_df, FA_boot_netStrength_SD)

# Remove 11th threshold (100% consistency)
fig5_FA_netStrength_df <- fig5_FA_netStrength_df[-11,]

# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_5C_det_FA_netStrength.pdf", width=7, height=6)

FA_netStrength_plot <- ggplot(data= fig5_FA_netStrength_df, aes(thresh_percentile)) + geom_point(aes(y= motion_r), color="firebrick", size=6)

# + geom_smooth(aes(y = motion_r), color="red4", method = "gam", formula = y ~ s(x, k=4), size=2, se=FALSE)  

FA_netStrength_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_x_continuous(breaks=c(0, 25, 50, 75)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + geom_errorbar(aes(ymin= motion_r - FA_boot_netStrength_SD, ymax= motion_r + FA_boot_netStrength_SD), colour="black", width=1, size=0.8) 

# + ylim(-0.4, -0.2)

dev.off()
################################################


#################
### FIGURE 6A ###
#################

## Motion vs. Age Plot 
# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_6A_MotionOverAge.pdf", width=7, height=6)

motionAge_plot <- ggplot(data=Motion_n949_df, aes(age)) + geom_point(aes(y=dti64MeanRelRMS)) + geom_smooth(aes(y =dti64MeanRelRMS), color="red3", method="lm", size=1.5)

motionAge_plot + ylim(0,2) + scale_x_continuous(name="Age (years)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

dev.off()

#############################################################################
## PROBABILISTIC: Edges with Significant Age Effect (Prob Streamline Count ##
#############################################################################
Araw_edgeStrength <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_Araw_sigAgeEffect_sexCov_edgeWeights_n949.csv", header=FALSE)
colnames(Araw_edgeStrength)[6928] <- "scanid"

demogs_and_Araw_edgeStrength <- merge(Araw_edgeStrength, Motion_n949_df, by="scanid")

#######################################################
### Loop to Sobel test across all edges for Alength ###
#	   Write out pvals and Z for indirect effect      #
#######################################################
motion_sobel_df <- as.data.frame(matrix(NA, nrow=6928, ncol=3))
colnames(motion_sobel_df) <- c("edge_index","sobel_zscore", "indirect_coeff")
# Node index
motion_sobel_df[,1] <- 0:6927
# n <- 1

# Read in sig_edges_df, made in MATLAB (Motion_mediation_setup.m)
sig_Araw_edges_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_Araw_sigAgeEffect_edgeConsistency_streamlineLength_motionAgeEffects_sexCov_df.csv", header=FALSE)
# read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/Araw_sigAgeEffect_edgeConsistency_streamlineLength_df.csv", header=FALSE)

colnames(sig_Araw_edges_df) <- c("orig_edge_idx", "edge_consistency", "mean_length", "motion_r", "age_r")

for (i in 2:6928){
   
  # Define edge weight
  currEdge <- as.data.frame(demogs_and_Araw_edgeStrength[i])
  # Run sobel test 
  Age <- as.data.frame(demogs_and_Araw_edgeStrength$age)
  Motion <- as.data.frame(demogs_and_Araw_edgeStrength$dti64MeanRelRMS)
  tmp_df <- cbind(Age, Motion,  currEdge)
  colnames(tmp_df) <- c("Age", "Motion", "currEdge")
  sob_node_test <- sobel(pred= tmp_df$Age, med= tmp_df$Motion, out= tmp_df$currEdge)
  sob_Zscore <-sob_node_test$z.value
  sob_indirectCoeff <- sob_node_test$Indirect.Effect
  # Sobel Z-score (indirect effect)
  motion_sobel_df[i,2] <- sob_Zscore
  # Sobel pval (indirect effect)
  motion_sobel_df[i,3] <- sob_indirectCoeff
}

# Remove first row of NAs
motion_sobel_df <- motion_sobel_df[-1,]

# Merge with length and consistency data
motion_sobel_df <- cbind(motion_sobel_df, sig_Araw_edges_df)

# Edge consistency -- Negative log transformation for figs
# motion_sobel_df$edge_consistency <- -log(motion_sobel_df$edge_consistency)

# Look at range of results
hist(motion_sobel_df$sobel_zscore)
min(motion_sobel_df$sobel_zscore)
max(motion_sobel_df$sobel_zscore)

## Define index for edges showing significant positive Age effects
motion_sobel_df$posNeg_Index <- 0
motion_sobel_df$posNeg_Index[which(motion_sobel_df$sobel_zscore > 0)] <- 1
motion_sobel_df$posNeg_Index[which(motion_sobel_df$sobel_zscore < 0)] <- 2
motion_sobel_df$posNeg_Index <- as.factor(motion_sobel_df$posNeg_Index)

## Define significant mediation effects using FDR correction on converted sobel Z-score
motion_sobel_df$sobel_pvals <- 2 * pnorm(-abs(motion_sobel_df$sobel_zscore))
motion_sobel_df$FDRcorr_sobel_pvals <- p.adjust(motion_sobel_df$sobel_pvals, method = "fdr")

sig_edgeStrength_df <- subset(motion_sobel_df, motion_sobel_df$FDRcorr_sobel_pvals < 0.05)
dim(sig_edgeStrength_df)

## Significant Positive Mediation Effects
sig_Pos_edgeStrength_df <- sig_edgeStrength_df[which(sig_edgeStrength_df$sobel_zscore > 0),]
dim(sig_Pos_edgeStrength_df)
# hist(sig_Pos_edgeStrength_df$mean_length)
mean(sig_Pos_edgeStrength_df$mean_length)

sum(sig_Pos_edgeStrength_df$age_r > 0)
sum(sig_Pos_edgeStrength_df$age_r < 0)

sum(sig_Pos_edgeStrength_df$motion_r > 0)
sum(sig_Pos_edgeStrength_df$motion_r < 0)

## Significant Negative Mediation Effects
sig_Neg_edgeStrength_df <- sig_edgeStrength_df[which(sig_edgeStrength_df$sobel_zscore < 0),]
dim(sig_Neg_edgeStrength_df)
# hist(sig_Neg_edgeStrength_df$mean_length)
mean(sig_Neg_edgeStrength_df$mean_length)

sum(sig_Neg_edgeStrength_df$age_r > 0)
sum(sig_Neg_edgeStrength_df$age_r < 0)

sum(sig_Neg_edgeStrength_df$motion_r > 0)
sum(sig_Neg_edgeStrength_df$motion_r < 0)

#################
### FIGURE 6E ###
#################
# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_6E_prob_Araw_mediation_consistency_vioplot.pdf", width=7, height=6)

sigMediation_consistency_plot <- ggplot(data=sig_edgeStrength_df, aes(posNeg_Index, edge_consistency))

sigMediation_consistency_plot + geom_violin(aes(col = factor(posNeg_Index), fill=factor(posNeg_Index)), trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + scale_x_discrete(name="Direction of Mediation Effects (Positve, Negative)") + scale_y_continuous(name="Edge Consistency (-log CV)") + theme(panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_colour_manual(values = c("red3", "blue3")) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_fill_manual(values = alpha(c("red3", "blue3"),0.2))
# (values = c("red3", "blue3")) 

dev.off()

## Export original edge index and posNeg index for BrainNet Viewer renderings
write.table(sig_edgeStrength_df$orig_edge_idx, "/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_sig_ptx_Araw_mediation_orig_edge_idx.txt", col.names=FALSE, row.names=FALSE)
write.table(sig_edgeStrength_df$posNeg_Index, "/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_sig_ptx_Araw_mediation_posNeg_idx.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

#########################
## PERMUTATION TESTING ##
#########################

#############################################################################################
### Permute posNeg index and build a distribution of differences in mean edge consistency ###
#############################################################################################

# Observed difference in means
obs_mean_diff_consistency <- (mean(sig_Pos_edgeStrength_df$edge_consistency) - mean(sig_Neg_edgeStrength_df$edge_consistency))

# Define empty variable for mean difference
nPerm=10000
permuted_meanDiff_edgeConsistency <- as.data.frame(matrix(NA, nrow= nPerm, ncol=1))

set.seed(66)

for (i in 1:nPerm){

tmp_df <- sig_edgeStrength_df

# Permute posNeg index
tmp_df$perm_posNeg_idx <- permute(tmp_df$posNeg_Index)

# Define groups of edges using permuted label
tmp_pos_edges <- subset(tmp_df, tmp_df$perm_posNeg_idx == 1)
tmp_neg_edges <- subset(tmp_df, tmp_df$perm_posNeg_idx == 2)

# Calculate mean consistency for edges with (permuted) positive/negative mediation effects
mean_pos_consistency <- mean(tmp_pos_edges$edge_consistency)
mean_neg_consistency <- mean(tmp_neg_edges$edge_consistency)

permuted_meanDiff_edgeConsistency[i,] <- (mean_pos_consistency - mean_neg_consistency)

}

colnames(permuted_meanDiff_edgeConsistency) <- "permuted_meanDiff_edgeConsistency"
hist(permuted_meanDiff_edgeConsistency$permuted_meanDiff_edgeConsistency, col="navyblue")

## Compare empirical to null distribution
num_sig_perms <- sum(permuted_meanDiff_edgeConsistency$permuted_meanDiff_edgeConsistency > obs_mean_diff_consistency)
num_sig_perms

##################################################################################
## DETERMINISTIC: Edges with Significant Age Effect (Deterministic inverse ADC) ##
##################################################################################
FA_edgeStrength <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_det_FA_sigAgeEffect_sexCov_edgeWeights_n949.csv", header=FALSE)
colnames(FA_edgeStrength)[1579] <- "scanid"

demogs_and_FA_edgeStrength <- merge(FA_edgeStrength, Motion_n949_df, by="scanid")  # sort=FALSE)

#######################################################
### Loop to Sobel test across all edges for Alength ###
#	   Write out pvals and Z for indirect effect      #
#######################################################
motion_sobel_df <- as.data.frame(matrix(NA, nrow=1579, ncol=3))
colnames(motion_sobel_df) <- c("edge_index","sobel_zscore", "indirect_coeff")
# Node index
motion_sobel_df[,1] <- 0:1578

for (i in 2:1579){
   
  # Define edge weight
  currEdge <- as.data.frame(demogs_and_FA_edgeStrength[i])
  # Run sobel test 
  Age <- as.data.frame(demogs_and_FA_edgeStrength$age)
  Motion <- as.data.frame(demogs_and_FA_edgeStrength$dti64MeanRelRMS)
  tmp_df <- cbind(Age, Motion,  currEdge)
  colnames(tmp_df) <- c("Age", "Motion", "currEdge")
  sob_node_test <- sobel(pred= tmp_df$Age, med= tmp_df$Motion, out= tmp_df$currEdge)
  sob_Zscore <-sob_node_test$z.value
  sob_indirectCoeff <- sob_node_test$Indirect.Effect

  # Sobel Z-score (indirect effect)
  motion_sobel_df[i,2] <- sob_Zscore
  # Sobel pval (indirect effect)
  motion_sobel_df[i,3] <- sob_indirectCoeff
}

# Remove first row of NAs
motion_sobel_df <- motion_sobel_df[-1,]

# Read in sig_edges_df, made in MATLAB (Motion_mediation_setup.m)
sig_FA_edges_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_det_FA_sigAgeEffect_edgeConsistency_streamlineLength_sexCov_df.csv", header=FALSE)
colnames(sig_FA_edges_df) <- c("orig_edge_idx", "edge_consistency", "mean_length")

# Merge consistnecy and length data
motion_sobel_df <- cbind(motion_sobel_df,sig_FA_edges_df)
# motion_sobel_df <- cbind(motion_sobel_df, edge_consistency)

# Look at range of results
hist(motion_sobel_df$sobel_zscore)
min(motion_sobel_df$sobel_zscore)
max(motion_sobel_df$sobel_zscore)

## Define index for edges showing significant positive Age effects
motion_sobel_df$posNeg_Index <- 0
motion_sobel_df$posNeg_Index[which(motion_sobel_df$sobel_zscore > 0)] <- 1
motion_sobel_df$posNeg_Index[which(motion_sobel_df$sobel_zscore < 0)] <- 2
motion_sobel_df$posNeg_Index <- as.factor(motion_sobel_df$posNeg_Index)

## Define significant mediation effects using FDR correction on converted sobel Z-score
motion_sobel_df$sobel_pvals <- 2 * pnorm(-abs(motion_sobel_df$sobel_zscore))
motion_sobel_df$FDRcorr_sobel_pvals <- p.adjust(motion_sobel_df$sobel_pvals, method = "fdr")
sig_edgeStrength_df <- subset(motion_sobel_df, motion_sobel_df$FDRcorr_sobel_pvals < 0.05)
dim(sig_edgeStrength_df)

## Significant Positive Mediation Effects
sig_Pos_edgeStrength_df <- sig_edgeStrength_df[which(sig_edgeStrength_df$sobel_zscore > 0),]
dim(sig_Pos_edgeStrength_df)
# hist(sig_Pos_edgeStrength_df$mean_length)
mean(sig_Pos_edgeStrength_df$mean_length)

## Significant Negative Mediation Effects
sig_Neg_edgeStrength_df <- sig_edgeStrength_df[which(sig_edgeStrength_df$sobel_zscore < 0),]
dim(sig_Neg_edgeStrength_df)
# hist(sig_Neg_edgeStrength_df$mean_length)
mean(sig_Neg_edgeStrength_df$mean_length)

#################
### FIGURE 6F ###
#################
# Export pdf
pdf(file="/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Routput/figure_6F_det_FA_mediation_consistency_vioplot.pdf", width=7, height=6)

sigMediation_FA_consistency_plot <- ggplot(data=sig_edgeStrength_df, aes(posNeg_Index, edge_consistency))

sigMediation_FA_consistency_plot + geom_violin(aes(col = factor(posNeg_Index), fill = factor(posNeg_Index)),trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + scale_x_discrete(name="Direction of Mediation Effects (Positve, Negative)") + theme(panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_colour_manual(values = c("red3", "blue3")) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_fill_manual(values = alpha(c("red3", "blue3"),0.2)) + scale_y_continuous(breaks=c(0, 25, 50, 75, 100))

dev.off()

#########################
## PERMUTATION TESTING ##
#########################

#############################################################################################
### Permute posNeg index and build a distribution of differences in mean edge consistency ###
#############################################################################################

# Observed difference in means
obs_mean_diff_consistency <- (mean(sig_Pos_edgeStrength_df$edge_consistency) - mean(sig_Neg_edgeStrength_df$edge_consistency))

# Define empty variable for mean difference
nPerm=10000
permuted_meanDiff_edgeConsistency <- as.data.frame(matrix(NA, nrow= nPerm, ncol=1))

set.seed(66)

for (i in 1:nPerm){
tmp_df <- sig_edgeStrength_df

# Permute posNeg index
tmp_df$perm_posNeg_idx <- permute(tmp_df$posNeg_Index)

# Define groups of edges using permuted label
tmp_pos_edges <- subset(tmp_df, tmp_df$perm_posNeg_idx == 1)
tmp_neg_edges <- subset(tmp_df, tmp_df$perm_posNeg_idx == 2)

# Calculate mean consistency for edges with (permuted) positive/negative mediation effects
mean_pos_consistency <- mean(tmp_pos_edges$edge_consistency)
mean_neg_consistency <- mean(tmp_neg_edges$edge_consistency)

permuted_meanDiff_edgeConsistency[i,] <- (mean_pos_consistency - mean_neg_consistency)
}

colnames(permuted_meanDiff_edgeConsistency) <- "permuted_meanDiff_edgeConsistency"
hist(permuted_meanDiff_edgeConsistency$permuted_meanDiff_edgeConsistency, col="navyblue")

## Compare empirical to null distribution
num_sig_perms <- sum(permuted_meanDiff_edgeConsistency$permuted_meanDiff_edgeConsistency > obs_mean_diff_consistency)
num_sig_perms

## Export original edge index for BrainNet Viewer renderings
write.table(sig_edgeStrength_df$orig_edge_idx, "/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_sig_det_FA_mediation_orig_edge_idx.txt", col.names=FALSE, row.names=FALSE)
write.table(as.numeric(sig_edgeStrength_df$posNeg_Index), "/data/joy/BBL/projects/pncBaumDti/Motion_paper/figures/Figure_6/mediation/NEW_sig_det_FA_mediation_posNeg_idx.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

save.image("/data/joy/BBL/projects/pncBaumDti/Motion_paper/scripts/Group_analysis/DTI_Motion_Workspace.Rdata")
