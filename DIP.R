#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This is a script for performing DIP analyses.
#The input data for DIP are mulitple fasta files, each of which represnts a single-locus alignment
#The related script, Make_locus_alignments.R, can be used to generate these single-locus alignments.

#Usage:
#Rscript --vanilla DIP.R <A1:Jobname> <A2:type_of_alignment> <A3:directory_containing_alignments_and_listfile> <A4:Simulation replicate> 
#<A5:Total_number_of_taxa_in_alignments> <A6:P1_taxon_string> <A7:P2_taxon_string> <A8:P3_taxon_string> <A9:Outgroup_taxon_string>

#Set working directory to the directory in which this script lives
setwd("/DIP/")

#designated the package to install/library
package_list<-c("phyclust", "parallel", "plyr", "seqinr", "reshape2", "ggplot2", "gplots", "RColorBrewer", "ape")

#Loop to install/library needed packages
for(p in 1:length(package_list)){

  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

#Set Bootstrap cutoff. This is a little bit a relic.
#By setting to 0, we essentially bypass the bootstap step.
BS_cutoff<-0.0

#Read jobname from command line arguments
#NOTE: this is used as a search term for your alignment files. 
#Should be a string that is present in the names of all alignments
jobname<-args[1]

#Read alignment type (fasta/phylip) from command line arguments
Aln_type<-args[2]

#Designate the directory in which single-locus fasta files live (read from command line args)
sl_aln_dir<-args[3]

#This is specific to simulation analyses. 
#If analyzing empirical data, make this the same as "jobname" above
sim_rep<-args[4]

#create list of single-locus alignment file names
sl_listfile<-list.files(path = sl_aln_dir, pattern = jobname)

#Extract only the alignment names containing the sim_rep
#If analyzing empirical data, this step is redundant
sl_listfile<-sl_listfile[grep(sim_rep, sl_listfile)]

#Trim down to only inlude names with .fa or .phy suffix
#This is just a step to help make sure only desired alignment files are analyzed
if(Aln_type=="phylip"){
  sl_listfile<-sl_listfile[grep("phy", sl_listfile)]
}else if(Aln_type=="fasta"){
  sl_listfile<-sl_listfile[grep("fa", sl_listfile)]
}

#The total number of taxa in each alignment
#This is just another sanity check to make sure the correct alignments are analyzed
num_taxa<-as.numeric(args[5])

#Define P1, P2, P3, and outgroup respectively
taxa_prefix<-c(args[6], args[7], args[8], args[9])

#initialize an empty matrix to store the results of the loop below
DandT_mat<-matrix(, ncol = 10, nrow =length(sl_listfile))

#Loop through all the single-locus alignments specified by sl_listfile
for(x in 1:length(sl_listfile)){

#Read the alignments (either fasta or phylip format)
  if(Aln_type=="phylip"){
    sl_aln<-read.dna(file = paste(sl_aln_dir, sl_listfile[x], sep = ""), format = "sequential")
  }else if(Aln_type=="fasta"){
  sl_aln<-read.dna(file = paste(sl_aln_dir, sl_listfile[x], sep = ""), format = "fasta")
}

#if statement to make sure the alignment has (at least) the minimum expected number of taxa
if(length(labels(sl_aln))>=as.numeric(paste(num_taxa))){
  
#Assign the taxon names
P1<-grep(taxa_prefix[1],labels(sl_aln))
P2<-grep(taxa_prefix[2],labels(sl_aln))
P3<-grep(taxa_prefix[3],labels(sl_aln))
P_out<-grep(taxa_prefix[4],labels(sl_aln))

#Prune the alignment to contain the focal species
sl_aln_pruned<-sl_aln[c(P1, P2, P3, P_out), ]

#convert to the proper class
dist_mat<-dist.dna(sl_aln_pruned, as.matrix = TRUE)

#Skip loci that contain NA in dist matrix
if(!any(is.na(dist_mat))){

#Extract divergence values
K23<-dist_mat[2,3] 
K12<-dist_mat[1,2] 
K13<-dist_mat[1,3] 

#Extract ALT divergence values (P1 and P2 are switched)
K23_alt<-dist_mat[1,3] #Formerly D1
K12_alt<-dist_mat[1,2] #Formerly D2
K13_alt<-dist_mat[2,3] #Formerly D3

#Get the neighbor joining topology of tree
tree<-nj(dist_mat)

#Root the tree
root_tree<-root.phylo(tree, outgroup = grep(taxa_prefix[4],tree$tip.label), resolve.root = TRUE)

#Get the topology of the full nj tree
if(is.monophyletic(tree, tips = c(1, 2))){
  Topology<-"12top"
}else if(is.monophyletic(tree, tips = c(2, 3))){
  Topology<-"23top"
}else if(is.monophyletic(tree, tips = c(1, 3))){
  Topology<-"13top"
}

#Get the bs support for the gene tree topology
#Note: This is a relic; the bootstrap component is no longer included 
bs_for_top<-1.0

#Get the newick tree. Store as an object
tree_newick<-write.tree(root_tree, file = "")

#Make sure there were no errors in the above data extraction
if(exists("Topology") && exists("K23") && exists("K12") && exists("K13") && 
   !is.na(Topology) && !is.na(K23) && !is.na(K12) && !is.na(K13)){
  DandT_mat[x,]<-c(paste(sl_listfile[x]), Topology, bs_for_top, K23, K12, K13, K23_alt, K12_alt, K13_alt, tree_newick)

} #End of is.na if statement
} #End loop to see if dist.dna contains NAs
} #end of number of taxa if 
} #end alignment for loop

#Convert output matrix into a dataframe
DandT_df<-as.data.frame(DandT_mat)

#Add names to the columns
names(DandT_df)<-c("Alignment_name", "Topology", "bs_for_top", "K23", "K12", "K13", "K23_alt", "K12_alt", "K13_alt", "tree_newick")

#Make sure each variable is  stored in the correct class
DandT_df_clean<-data.frame(Alignment_name=DandT_df$Alignment_name, 
                           Topology=DandT_df$Topology,
                           bs_for_top=as.numeric(paste(DandT_df$bs_for_top)),
                           K23=as.numeric(paste(DandT_df$K23)),
                           K12=as.numeric(paste(DandT_df$K12)),
                           K13=as.numeric(paste(DandT_df$K13)),
                           K23_alt=as.numeric(paste(DandT_df$K23_alt)),
                           K12_alt=as.numeric(paste(DandT_df$K12_alt)),
                           K13_alt=as.numeric(paste(DandT_df$K13_alt)),
                           tree_newick=paste(DandT_df$tree_newick)
                           )

#Write a copy of the divergence/tree data for each locus
write.csv(DandT_df_clean, file = paste(jobname, sim_rep, "_DIP_data_out", Sys.Date(), ".csv", sep = ""))

### INFER MODE ###
#Subset dataframe by Topology
SP_loci<-subset(DandT_df_clean, Topology=="12top" & bs_for_top >= BS_cutoff)
IG_loci<-subset(DandT_df_clean, Topology=="23top" & bs_for_top >= BS_cutoff)
ILS_loci<-subset(DandT_df_clean, Topology=="13top" & bs_for_top >= BS_cutoff)

#Get the deltas for the full genome
### INFER MODE ###
deltaK23_full<-(mean(SP_loci$K23)-mean(IG_loci$K23))
deltaK12_full<-(mean(IG_loci$K12)-mean(SP_loci$K12))
deltaK13_full<-(mean(SP_loci$K13)-mean(IG_loci$K13))
deltaK23_alt_full<-(mean(SP_loci$K23_alt)-mean(ILS_loci$K23_alt))
deltaK12_alt_full<-(mean(ILS_loci$K12_alt)-mean(SP_loci$K12_alt))
deltaK13_alt_full<-(mean(SP_loci$K13_alt)-mean(ILS_loci$K13_alt))
delta_delta_full<-deltaK12_full-deltaK13_full
delta_delta_alt_full<-deltaK12_alt_full-deltaK13_alt_full
deltaX3_full<-delta_delta_full-delta_delta_alt_full
deltaX3_wt_full<-delta_delta_full-((nrow(ILS_loci)/nrow(IG_loci))*delta_delta_alt_full)
### END INFER MODE ###

# ### OMNISCIENT MODE ###
# #Subset dataframe by Topology
# SP_loci_omn<-DandT_df_clean[1:(nrow(DandT_df_clean)/2),]
# IG_loci_omn<-DandT_df_clean[(nrow(DandT_df_clean)/2):nrow(DandT_df_clean),]
# 
# #Get the deltas for the full genome
# deltaK23_full_omn<-(mean(SP_loci_omn$K23)-mean(IG_loci_omn$K23))
# deltaK12_full_omn<-(mean(IG_loci_omn$K12)-mean(SP_loci_omn$K12))
# deltaK13_full_omn<-(mean(SP_loci_omn$K13)-mean(IG_loci_omn$K13))
# delta_delta_full_omn<-deltaK12_full_omn-deltaK13_full_omn
# ### END OMNISCIENT MODE ###
# 
# ### OMNISCIENT PLUS INFER MODE ###
# #Subset dataframe by Topology
# SP_loci_oi<-subset(DandT_df_clean[1:(nrow(DandT_df_clean)/2),], Topology=="12top" & bs_for_top >= BS_cutoff)
# IG_loci_oi<-subset(DandT_df_clean[(nrow(DandT_df_clean)/2):nrow(DandT_df_clean),], Topology=="23top" & bs_for_top >= BS_cutoff)
# 
# #Get the deltas for the full genome (no resampled)
# deltaK23_full_oi<-(mean(SP_loci_oi$K23)-mean(IG_loci_oi$K23))
# deltaK12_full_oi<-(mean(IG_loci_oi$K12)-mean(SP_loci_oi$K12))
# deltaK13_full_oi<-(mean(SP_loci_oi$K13)-mean(IG_loci_oi$K13))
# delta_delta_full_oi<-deltaK12_full_oi-deltaK13_full_oi
# ### END OMNISCIENT PLUS INFER MODE ###


#Set number of resampling reps
#1000 is used by default
delta_reps<-1000

#Creat DIP function (to be replicated <delta_reps> times)
delta_delta_func<-function(){
  #Bootstrap resample the test data
  temp_df<-DandT_df_clean[sample(1:nrow(DandT_df_clean), nrow(DandT_df_clean), replace = TRUE), ]

  #Subset SP and IG loci
    #subset based on gene tree topology
    IG_loci_temp<-subset(temp_df,
                         Topology=="23top" &
                           bs_for_top >= BS_cutoff)
    SP_loci_temp<-subset(temp_df,
                         Topology=="12top" &
                           bs_for_top >= BS_cutoff)
    ILS_loci_temp<-subset(temp_df,
                         Topology=="13top" &
                           bs_for_top >= BS_cutoff)

  #Get the deltas
  deltaK23<-(mean(SP_loci_temp$K23)-mean(IG_loci_temp$K23))
  deltaK12<-(mean(IG_loci_temp$K12)-mean(SP_loci_temp$K12))
  deltaK13<-(mean(SP_loci_temp$K13)-mean(IG_loci_temp$K13))
  deltaK23_alt<-(mean(SP_loci_temp$K23_alt)-mean(ILS_loci_temp$K23_alt))
  deltaK12_alt<-(mean(ILS_loci_temp$K12_alt)-mean(SP_loci_temp$K12_alt))
  deltaK13_alt<-(mean(SP_loci_temp$K13_alt)-mean(ILS_loci_temp$K13_alt))
  delta_delta<-deltaK12-deltaK13
  delta_delta_alt<-deltaK12_alt-deltaK13_alt
  deltaX3<-delta_delta-delta_delta_alt
  deltaX3_wt<-delta_delta-((nrow(ILS_loci_temp)/nrow(IG_loci_temp))*delta_delta_alt)

  #Return the results when function is called
  return(c(deltaK23, deltaK12, deltaK13, deltaK23_alt, deltaK12_alt, deltaK13_alt, delta_delta, delta_delta_alt, deltaX3, deltaX3_wt))
}

#Replicate the above function
delta_delta_out<-replicate(delta_reps, delta_delta_func())

#Convert function output to dataframe
delta_delta_out_df<-data.frame(
  deltaK23_reps=delta_delta_out[1, ],
  deltaK12_reps=delta_delta_out[2, ],
  deltaK13_reps=delta_delta_out[3, ],
  deltaK23_alt_reps=delta_delta_out[4, ],
  deltaK12_alt_reps=delta_delta_out[5, ],
  deltaK13_alt_reps=delta_delta_out[6, ],
  delta_delta_reps=delta_delta_out[7, ],
  delta_delta_alt_reps=delta_delta_out[8, ],
  deltaX3_reps=delta_delta_out[9, ],
  deltaX3_wt_reps=delta_delta_out[10, ]
)

### MAKE PLOTS AND WRITE OUTPUT ###

###1xDIP FRAMEWORK

#Caculate the three p-values for 1xDIP
prof_new_pvalK23<-length(which(delta_delta_out_df$deltaK23_reps<=0))/length(delta_delta_out_df$deltaK23_reps)
prof_new_pvalK12<-length(which(delta_delta_out_df$deltaK12_reps<=0))/length(delta_delta_out_df$deltaK12_reps)
prof_new_pvalK13<-length(which(delta_delta_out_df$deltaK13_reps<=0))/length(delta_delta_out_df$deltaK13_reps)

# #Get the 1xDIP results (results are coded as 1 (P3->P2), -1 (P2->P3), or 0 (bidirectional))
# #This is only used for graphing large parameter scans
# if(prof_new_pvalK23<=0.05 && prof_new_pvalK12<=0.05 && prof_new_pvalK13>=0.05){
#   DIPX1_result<-1
# }else if(prof_new_pvalK23<=0.05 && prof_new_pvalK12>=0.05 && prof_new_pvalK13<=0.05){
#   DIPX1_result<-(-1)
# }else{
#   DIPX1_result<-0
# }
   
#Create Violin Plot for 1xDIP
#Rearrange data to make it suitable to plot with ggplot
violin_data<-rbind(
data.frame(Divergence_type="K23", Divergence_value=delta_delta_out_df$deltaK23_reps),
data.frame(Divergence_type="K12", Divergence_value=delta_delta_out_df$deltaK12_reps),
data.frame(Divergence_type="K13", Divergence_value=delta_delta_out_df$deltaK13_reps)
)

#get rid of NAs. There shouldn't be any of these, just a precaution
violin_data<-subset(violin_data, violin_data$Divergence_value!="NA")

#Print pdf of 1xDIP test
pdf(paste(jobname, "_1xDIP_", Sys.Date(), ".pdf", sep = ""),width=6,height=6)

#Plot violit plot
ggplot(violin_data, aes(x = Divergence_type, y = Divergence_value)) +
  geom_violin()+
  geom_hline(yintercept = 0, linetype="dashed", size=1)+
theme(axis.text.x = element_text(face="italic", size=10, angle=270, hjust = 0, vjust = 0.5),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) +
      labs(title=paste("K23 p-val: ",prof_new_pvalK23, "\n",
                       "K12 p-val: ",prof_new_pvalK12, "\n",
                       "K13 p-val: ",prof_new_pvalK13, sep = "" ),
       x ="Delta comparison", y = "Delta")

dev.off()

### 2xDIP
#P-value from 2xDIP
if(length(which(delta_delta_out_df$delta_delta_reps>0))>0.5*delta_reps){
  dd_direction<-"P3P2"
  dd_pval<-2*(length(which(delta_delta_out_df$delta_delta_reps<0)))/length(delta_delta_out_df$delta_delta_reps)
}else if(length(which(delta_delta_out_df$delta_delta_reps<0))>0.5*delta_reps){
  dd_direction<-"P2P3"
  dd_pval<-2*(length(which(delta_delta_out_df$delta_delta_reps>0)))/length(delta_delta_out_df$delta_delta_reps)
}

# #Get 2xDIP results (coded as numeric as for 1xDIP)
# #This is only used for plotting parameter scans
# if(dd_direction=="P3P2" && dd_pval<=0.05){
#   DIPX2_result<-1
# }else if(dd_direction=="P2P3" && dd_pval<=0.05){
#   DIPX2_result<-(-1)
# }else{
#   DIPX2_result<-0
# }

#### WRITE THE DIP RESULT FILE
#write.csv(c(DIPX1, DIPX2), file = paste(jobname, sim_rep, "DIP_results.csv", sep = ""))

#Print pdf of 2xDIP
pdf(paste(jobname, "_2xDIP_", Sys.Date(), ".pdf", sep = ""),width=6,height=6)

#Plot the distribution of delta delta replicates
dens<-density(delta_delta_out_df$delta_delta_reps, na.rm = TRUE, adjust = 5)
plot(dens, main = paste("P-value:", dd_pval)
    #, xlim = c(-0.02, 0.03) #This can be set manually to make the plot look better
     )
#Add a line at 0
abline(v = 0)

#This sets up coloring the area under the curve
x1 <- min(which(dens$x >= 0))
x2 <- max(which(dens$x <  10000000)) #arbitrary large number
x3 <- min(which(dens$x >= -10000000)) #arbitrary largely negative number
x4 <- max(which(dens$x <  0))

#Color the area under the curve
if(is.finite(x1) & is.finite(x2) & is.finite(x3) & is.finite(x4)){
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="red"))
with(dens, polygon(x=c(x[c(x3,x3:x4,x4)]), y= c(0, y[x3:x4], 0), col="gray"))
}

dev.off()

## 3xDIP PLOT
#Get p-values
if(length(which(delta_delta_out_df$deltaX3_wt_rep>0))>0.5*delta_reps){
  ddx3_direction<-"P3P2"
  ddx3_pval<-2*(length(which(delta_delta_out_df$deltaX3_wt_rep<0)))/length(delta_delta_out_df$deltaX3_wt_rep)
}else if(length(which(delta_delta_out_df$deltaX3_wt_rep<0))>0.5*delta_reps){
  ddx3_direction<-"P2P3"
  ddx3_pval<-2*(length(which(delta_delta_out_df$deltaX3_wt_rep>0)))/length(delta_delta_out_df$deltaX3_wt_rep)
}

#Print pdf of 3xDIP
pdf(paste(jobname, "_3xDIP_", Sys.Date(), ".pdf", sep = ""),width=6,height=6)

#Plot the distribution of 3xDIP replicates
dens<-density(delta_delta_out_df$deltaX3_wt_reps, na.rm = TRUE, adjust = 5)
plot(dens, main = paste("P-value:", ddx3_pval)
     , xlim = c(-0.03, 0.03)
)
#Add a line at 0
abline(v = 0)

#This sets up coloring the area under the curve
x1 <- min(which(dens$x >= 0))
x2 <- max(which(dens$x <  10000000))
x3 <- min(which(dens$x >= -10000000))
x4 <- max(which(dens$x <  0))

if(is.finite(x1) & is.finite(x2) & is.finite(x3) & is.finite(x4)){
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="red"))
  with(dens, polygon(x=c(x[c(x3,x3:x4,x4)]), y= c(0, y[x3:x4], 0), col="gray"))
}

dev.off()
