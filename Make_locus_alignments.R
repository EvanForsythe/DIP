#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This is a script for processing whole-chromosome alignments into many single-locus alignments. 

#Set working directory on HPC
setwd("/DIP/")

#Load packages
package_list<-c("ape")

#Loop to check if package is installed and libraried
for(p in 1:length(package_list)){
  
  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

#Read in all the arguments passed by the user
aln_dir<-args[1]

#output directory
out_dir<-args[2]
  
#Name of the listfile
listfile<-list.files(path = aln_dir, pattern = args[3])

#Number of taxa
n_tax<-as.numeric(args[4])

#set the size of windows
win_length<-as.numeric(args[5])

###Create a loop to look through alignments###
##############################################
#loop through each full genome/chromosome alignment file
for(x in 1:length(listfile)){

#Read the alignment
aln<-read.dna(file = paste(aln_dir, listfile[x], sep = ""), format = "fasta")

#define the starting sites of each window
starts <- seq(from =0, to =length(aln[1,])-win_length, by = win_length)

#Make the first window start at 1 instead of 0
starts[1]<-1

#Define the number of windows
n_wins<-length(starts)

###Create a loop to look through windows in alignment###
########################################################
#Run loop
for(y in 1:n_wins){
  #Get the window
  window<-aln[1:n_tax, starts[y]:(starts[y]+(win_length-1))]

  #Is the window filled with "N"s?
  if(length(which(dist.dna(window)=="NaN"))==0){
    #Write a single-locus alignment for each window
    write.dna(window, 
              file = paste(out_dir, listfile[x],"__", sprintf("%08d", starts[y]), ".fa", sep = ""),
              format = "fasta"
              )
    } #if statement
  } #window loop
rm(aln) #Remove the aln object from the R environment
}#Chromosome file loop

