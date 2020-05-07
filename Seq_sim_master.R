#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Designate the directory where this script lives (be sure to end the path with a "/")
work_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/IG_direction/DIP/DIP_revisions/NEW_EVERYTHING/NEW_VERSION_FOR_GITHUB/"
  
#Set the working directory
setwd(work_dir)

#Desgnate which packages to load and library
package_list<-c("phyclust", "parallel", "plyr", "seqinr", "reshape2", "ggplot2", "gplots", "RColorBrewer", "ape")

#Loop to check if package is installed and libraried
for(p in 1:length(package_list)){
  
  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

#Check if output directory exists (create one if not)
if(!dir.exists(paste(work_dir, "Simulated_alignments/", sep = ""))){
  system(paste("mkdir ", paste(work_dir, "Simulated_alignments/", sep = ""), sep = ""))
}

#Designate a directory to store simulated alignments
aln_dir<-paste(work_dir, "Simulated_alignments/", sep = "")

#Directory to store DIP output
out_csv_dir<-work_dir

#Check if output directory exists for single-locus alignments (create one if not)
if(!dir.exists(paste(work_dir, "Simulated_alignments/Windows/", sep = ""))){
  system(paste("mkdir ", paste(work_dir, "Simulated_alignments/Windows/", sep = ""), sep = ""))
}

#Directory where single-locus (window) alignments will output
window_dir<-paste(work_dir, "Simulated_alignments/Windows/", sep = "")

#How many loci in simulated genome?
reps<-as.numeric(paste(args[1]))

#How long are loci?
loci_ln<-as.numeric(paste(args[2]))

#Set pIG
pIG<-as.numeric(paste(args[3]))


#Set p32 range
p32<-seq(0, 1, by = 0.05)

#For iterating through different scaling factors, give range.
#scale_fact<-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)

#For standard seqsim analysis use scaling factor of 1
scale_fact<-1

#If performing replicates, number the reps
#genome_reps<-1:10

#If only doing one rep
genome_reps<-1

#Start for loop to loop through scaling factors
for(s in 1:length(scale_fact)){
  #Start loop for different values of p32
  for(p in 1:length(p32)){
  
  #Also loop through the genome simulation reps
  for(r in 1: length(genome_reps)){
  
    #Check if file exists
    if(
      !length(
        grep(paste("pIG_",pIG,"p32_", p32[p], "_5000_",scale_fact[s],"_g_rep", genome_reps[r], sep = ""), 
             list.files(path = aln_dir, pattern = paste("_", scale_fact[s], "_g", sep = "")))
      )>0
    ){  
      
#create an empty alignment matrix
aln_mat<-matrix(, nrow = 4, ncol = reps*loci_ln)

#Start a ticker for filling the matrix
win_ticker<-0

#Simulate noIG gene tree
for(n in 1:floor(reps*(1-pIG))){
ret.ms <- ms(nsam = 4, nreps = 1,
             opts = paste("-T -t 50 -I 4 1 1 1 1 -ej ",4*scale_fact[s]," 2 1 -ej ",8*scale_fact[s]," 3 1 -ej ",12*scale_fact[s]," 4 1 -r 5 ", loci_ln, sep = "")
)

#Use seqgen to simulate sequences
seqs<-seqgen(opts = "-mHKY -l5000 -s 0.01", rooted.tree = read.tree(text = ret.ms[3]))

#Extract each sequence
seq1<-strsplit(c2s(substring(text = (seqs[grep("s1", seqs)]), 11, 5010)), split = "")
seq2<-strsplit(c2s(substring(text = (seqs[grep("s2", seqs)]), 11, 5010)), split = "")
seq3<-strsplit(c2s(substring(text = (seqs[grep("s3", seqs)]), 11, 5010)), split = "")
seq4<-strsplit(c2s(substring(text = (seqs[grep("s4", seqs)]), 11, 5010)), split = "")

#Populate a portion of the matrix
aln_mat[1,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq1)
aln_mat[2,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq2)
aln_mat[3,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq3)
aln_mat[4,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq4)
#Advance the window ticker
win_ticker<-win_ticker+1

}#End loop

if(floor(reps*(pIG*p32[p]))>0){  #Are there any P3->P2 loci?
#Simulate P3->P2 IG tree
for(n in 1:floor(reps*(pIG*p32[p]))){
  ret.ms <- ms(nsam = 4, nreps = 1,
               opts = paste("-T -t 50 -I 4 1 1 1 1 -ej ",4*scale_fact[s]," 2 1 -ej ",8*scale_fact[s]," 3 1 -ej ",12*scale_fact[s]," 4 1 -es ",1*scale_fact[s]," 2 0 -ej ",1*scale_fact[s]," 5 3 -r 5 ", loci_ln, sep = "")
  )

  #Use seqgen to simulate sequences
  seqs<-seqgen(opts = "-mHKY -l5000 -s 0.01", rooted.tree = read.tree(text = ret.ms[3]))

  #Extract each sequence
  seq1<-strsplit(c2s(substring(text = (seqs[grep("s1", seqs)]), 11, 5010)), split = "")
  seq2<-strsplit(c2s(substring(text = (seqs[grep("s2", seqs)]), 11, 5010)), split = "")
  seq3<-strsplit(c2s(substring(text = (seqs[grep("s3", seqs)]), 11, 5010)), split = "")
  seq4<-strsplit(c2s(substring(text = (seqs[grep("s4", seqs)]), 11, 5010)), split = "")

  #Populate a portion of the matrix
  aln_mat[1,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq1)
  aln_mat[2,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq2)
  aln_mat[3,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq3)
  aln_mat[4,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq4)
  #Advance the windo ticker
  win_ticker<-win_ticker+1

}#End loop
}#End if

if(floor(reps*(pIG*(1-p32[p])))>0){ #Are there any P2->P3 loci?
#Simulate P2->P3 IG
for(n in 1:floor(reps*(pIG*(1-p32[p])))){
  ret.ms <- ms(nsam = 4, nreps = 1,
               opts = paste("-T -t 50 -I 4 1 1 1 1 -ej ",4*scale_fact[s]," 2 1 -ej ",8*scale_fact[s]," 3 1 -ej 6 4 1 -es ",1*scale_fact[s]," 3 0 -ej ",1*scale_fact[s]," 5 2 -r 5 ", loci_ln, sep = "")
  )

  #Use seqgen to simulate sequences
  seqs<-seqgen(opts = "-mHKY -l5000 -s 0.01", rooted.tree = read.tree(text = ret.ms[3]))

  #Extract each sequence
  seq1<-strsplit(c2s(substring(text = (seqs[grep("s1", seqs)]), 11, 5010)), split = "")
  seq2<-strsplit(c2s(substring(text = (seqs[grep("s2", seqs)]), 11, 5010)), split = "")
  seq3<-strsplit(c2s(substring(text = (seqs[grep("s3", seqs)]), 11, 5010)), split = "")
  seq4<-strsplit(c2s(substring(text = (seqs[grep("s4", seqs)]), 11, 5010)), split = "")

  #Populate a portion of the matrix
  aln_mat[1,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq1)
  aln_mat[2,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq2)
  aln_mat[3,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq3)
  aln_mat[4,((win_ticker*loci_ln)+1):((win_ticker*loci_ln)+loci_ln)]<-unlist(seq4)
  #Advance the windo ticker
  win_ticker<-win_ticker+1

}#End loop
}#End if

#Write the full alignment to fasta file
aln_mat_bin<-as.DNAbin(aln_mat)
write.FASTA(aln_mat_bin, file =
              paste(aln_dir,"pIG_",pIG,"p32_", p32[p],"_", reps,"_", scale_fact[s],"_g_rep",genome_reps[r],".fa", sep = "")
            )

    }#End of if statement to see if files exist
}#End of genome rep loop
  
  ### Run Make_locus_alignments.R script
  #See if it's been done already
  if(!length(grep(paste("pIG_",pIG,"p32_", p32[p], "_", sep = ""),
    list.files(path = window_dir,
                pattern = paste("_", scale_fact[s], "_g", sep = "")))) > 0
  ){
  system(paste(
    paste("Rscript --vanilla ", work_dir, "Make_locus_alignments.R", sep = ""),
    aln_dir,
    window_dir,
    paste("pIG_",pIG,"p32_",p32[p],"_5000_",scale_fact[s], "_g", sep = ""),
    "4",
    5000,
    sep = " "
  ))
  }#end if statement to see if it's been done already
    
  #Run DIP.R
  for(r in 1: length(genome_reps)){
    system(paste(
      paste("Rscript --vanilla ", work_dir, "DIP.R", sep = ""),
      paste("pIG_",pIG,"p32_",p32[p],"_5000_", scale_fact[s], sep = ""), #arg1
      "fasta", #arg2
      window_dir, #arg3
      paste("rep", genome_reps[r], ".fa", sep = ""), #arg4
      "4", "1", "2", "3", "4", #arg5-9
      sep = " "
    ))
  }#End of local g_reps loop
  }#End p32 loop
}#End of scaling factor loop

# ###### Everything below here is used to make parameter scan plots ######
# scale_fact<-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
# genome_reps<-1:10
# 
# #Initialize a matrix and a ticker
# SF_out_mat<-matrix(, nrow = length(scale_fact)*length(genome_reps), ncol = 10)
# row_ticker<-0
# out_csv_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/IG_direction/DIP/Output_0.6/"
# 
# #Start for loop to loop through scaling factors
# for(s in 1:length(scale_fact)){
#   #Also loop through the genome reps (different simulated genomes)
#   for(r in 1: length(genome_reps)){
# 
# row_ticker<-row_ticker+1
# SF_out_mat[row_ticker,]<-c(scale_fact[s],
#                            genome_reps[r],
#                            read.csv(file = paste(out_csv_dir, "p32_0.6_5000_", scale_fact[s], "rep", genome_reps[r], ".fa_full_genome_stats.csv", sep = ""))$x
#                            )
# }}#End the nested loops
# 
# #Convert to df and rename cols
# SF_out_df<-as.data.frame(SF_out_mat)
# names(SF_out_df)<-c("SF", "Rep", "SP_trees", "IG_trees", "ILS_trees", "DIPX2", "DIPX3", "DIPX3_wt", "DIPX2_omni", "DIPX2_omni_inf")
# 
# #Make some plots
# plot(jitter(SF_out_df$SF), SF_out_df$SP_trees, ylim = c(0, 2500), col = "orange")
# points(jitter(SF_out_df$SF), SF_out_df$IG_trees, col = "dark green")
# points(jitter(SF_out_df$SF), SF_out_df$ILS_trees, col = "purple")
# 
# plot(jitter(SF_out_df$SF), SF_out_df$DIPX2, pch = 1, cex=1, ylim = c(-0.035, 0.035))
# abline(h=0)
# 
# plot(jitter(SF_out_df$SF), SF_out_df$DIPX3_wt, pch = 1, cex=1, ylim = c(-0.035, 0.035))
# abline(h=0)
# 
# plot(jitter(SF_out_df$SF), SF_out_df$DIPX2_omni, pch = 1, cex=1, ylim = c(-0.035, 0.035))
# abline(h=0)
# 
# #Initialize a matrix and a ticker
# p32<-c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
# scale_fact<-c(0.1, 0.2, 0.3)
# genome_reps<-1:10
# p32scan_out_mat<-matrix(, nrow = length(p32)*length(scale_fact)*length(genome_reps), ncol = 11)
# row_ticker_p32scan<-0
# out_csv_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/IG_direction/DIP/Output_P32scan/"
# 
# #Loop through p32
# for(p in 1:length(p32)){
#   #Loop through scale factor
#   for(s in 1:length(scale_fact)){
#   #Loop through the genome reps (different simulated genomes)
#   for(r in 1: length(genome_reps)){
# 
# row_ticker_p32scan<-row_ticker_p32scan+1
# p32scan_out_mat[row_ticker_p32scan,]<-c(p32[p], scale_fact[s],
#                            genome_reps[r],
#                            read.csv(file = paste(out_csv_dir, "p32_", p32[p],"_5000_", scale_fact[s], "rep", genome_reps[r], ".fa_full_genome_stats.csv", sep = ""))$x
#                            )
# }}}#End the nested loops
# 
# #Convert to df and rename cols
# p32scan_out_df<-as.data.frame(p32scan_out_mat)
# names(p32scan_out_df)<-c("p32", "SF", "Rep", "SP_trees", "IG_trees", "ILS_trees", "DIPX2", "DIPX3", "DIPX3_wt", "DIPX2_omni", "DIPX2_omni_inf")
# 
# #DIPX2
# plot(jitter(p32scan_out_df$p32), p32scan_out_df$DIPX2, col= as.factor(p32scan_out_df$SF), pch = 1, cex=1, ylim = c(-0.035, 0.035))
# abline(h=0)
# abline(coef(glm(p32scan_out_df$DIPX2[which(p32scan_out_df$SF==0.1)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.1)])), col="black")
# abline(coef(glm(p32scan_out_df$DIPX2[which(p32scan_out_df$SF==0.2)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.2)])), col="red")
# abline(coef(glm(p32scan_out_df$DIPX2[which(p32scan_out_df$SF==0.3)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.3)])), col="green")
# 
# #DIPX3
# plot(jitter(p32scan_out_df$p32), p32scan_out_df$DIPX3_wt, col= as.factor(p32scan_out_df$SF), pch = 1, cex=1, ylim = c(-0.035, 0.035))
# abline(h=0)
# abline(coef(glm(p32scan_out_df$DIPX3_wt[which(p32scan_out_df$SF==0.1)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.1)])), col="black")
# abline(coef(glm(p32scan_out_df$DIPX3_wt[which(p32scan_out_df$SF==0.2)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.2)])), col="red")
# abline(coef(glm(p32scan_out_df$DIPX3_wt[which(p32scan_out_df$SF==0.3)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.3)])), col="green")
# 
# #DIPX2 with omniscence
# plot(jitter(p32scan_out_df$p32), p32scan_out_df$DIPX2_omni, col= as.factor(p32scan_out_df$SF), pch = 1, cex=1, ylim = c(-0.035, 0.035))
# abline(h=0)
# abline(coef(glm(p32scan_out_df$DIPX2_omni[which(p32scan_out_df$SF==0.1)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.1)])), col="black")
# abline(coef(glm(p32scan_out_df$DIPX2_omni[which(p32scan_out_df$SF==0.2)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.2)])), col="red")
# abline(coef(glm(p32scan_out_df$DIPX2_omni[which(p32scan_out_df$SF==0.3)]~p32scan_out_df$p32[which(p32scan_out_df$SF==0.2)])), col="green")
# 
# 
# #### LARGE PARAMETER SCAN ####
# outSCAN_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/IG_direction/DIP/Output_SCAN/"
# scan_rep<-"Rep05/"
# SCAN_files<-list.files(path = paste(outSCAN_dir, scan_rep, sep = ""), pattern = ".csv")
# 
# #Set parameter values to look at
# p32<-seq(0, 1, by = 0.05)
# pIG<-seq(0.05, 0.95, by = 0.05)
# #p<-1
# #i<-1
# 
# #Creat empty matrix
# SCAN_mat<-matrix(, nrow = length(p32)*length(pIG), ncol = 4)
# #Ticker
# SCAN_ticker<-0
# 
# for(p in 1:length(p32)){
#   for(i in 1:length(pIG)){
#     SCAN_ticker<-SCAN_ticker+1
#     SCAN_mat[SCAN_ticker, ]<-c(pIG[i], p32[p],
#       read.csv(file = paste(outSCAN_dir, scan_rep, "pIG_",pIG[i],"p32_",p32[p],"_5000_1repSCAN.faDIP_results.csv", sep = ""))$x)
#   }
# }
# 
# #Convert to df and clean
# SCAN_df<-as.data.frame(SCAN_mat)
# names(SCAN_df)<-c("pIG", "p32", "DIPX1_result", "DIPX2_result")
# 
# #Create custom colors so that all graphs have consistent colors:
# myColors <- c("red", "grey", "white")
# names(myColors) <- c(1, -1, 0)
# colScale <- scale_fill_manual(name = "DIP_SCAN",values = myColors)
# 
# ggplot(data = SCAN_df, aes(x=p32, y = pIG, fill=as.factor(SCAN_df$DIPX1_result)), colour=colScale) +
#   geom_tile()+ colScale +
# labs(title="DIPX1",
#      x="p(P3=>P2)", y="p(IG)")
# 
# ggplot(data = SCAN_df, aes(x=p32, y = pIG, fill=as.factor(SCAN_df$DIPX2_result)), colour=colScale) +
#   geom_tile()+ colScale +
#   labs(title="DIPX2",
#        x="p(P3=>P2)", y="p(IG)")
# 
# #Combine the SCAN reps
# #SCAN_df_rep05<-SCAN_df
# 
# SCAN_df_ALL<-SCAN_df_rep01+SCAN_df_rep02+SCAN_df_rep03+SCAN_df_rep04+SCAN_df_rep05
# 
# #Create custom colors so that all graphs have consistent colors:
# #myColors <- c("#333333", #666666", "#999999", "#CCCCCC", "#CCCCCC", "white", "#FF9999", "#FF6666", "#FF3333", "#CC0000")
# myColors <- c(rev(brewer.pal(n = 6, name = "Greys")[1:5]), "white", brewer.pal(n = 6, name = "Reds")[1:5])
# names(myColors) <- seq(-5, 5, by = 1)
# colScale <- scale_fill_manual(name = "DIP_SCAN",values = myColors)
# 
# 
# ggplot(data = SCAN_df_ALL, aes(x=p32, y = pIG, fill=as.factor(SCAN_df_ALL$DIPX1_result)), colour=colScale) +
#   geom_tile()+ colScale +
#   labs(title="DIPX1",
#        x="p(P3=>P2)", y="p(IG)")
# 
# ggplot(data = SCAN_df_ALL, aes(x=p32, y = pIG, fill=as.factor(SCAN_df_ALL$DIPX2_result)), colour=colScale) +
#   geom_tile()+ colScale +
#   labs(title="DIPX2",
#        x="p(P3=>P2)", y="p(IG)")
# 
