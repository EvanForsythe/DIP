import os
import sys
from subprocess import call

#save the chromosome that is to be downloaded
chrom= str(sys.argv[1])

#make sure a valid chromosome was chosen
if 1 <= int(chrom) <= 23:
    # Verify that the chromosome was interpreted and is valid
    print("Preparing to download data for chrom: " + str(chrom))
else:
    print("Usage: python Download_vcfs.py <interger indicating which chrom to study>\nNot a valid chromosome. Exiting...")
    exit()  ## exit python script

#Check
#Define a function for running commands
def CheckAndRun():
    #check if input file (specified as file_name_temp) exists and is not empty
    if os.path.isfile(file_name_temp) and os.stat(file_name_temp).st_size > 0:
        print(file_name_temp+" exists and is not empty\nrunning the following command...\n"+cmd_temp)
        os.system(cmd_temp)
    else:
        print("ERROR:"+file_name_temp+" does not exist or is empty. The following command could NOT be run...\n"+cmd_temp+"\nmoving on...")
        #sys.exit()

#store command for downloading reference genome files

ref_cmd=("curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr"+str(chrom)+".fa.gz> hg19_chrom"+str(chrom)+".fa.gz")

#Print message
print("downloading reference with the following command...")
print(ref_cmd)

#call the command to download the file
os.system(ref_cmd)

###Perfrom unzip
file_name_temp=("hg19_chrom"+str(chrom)+".fa.gz")
cmd_temp=("gunzip "+file_name_temp)
CheckAndRun()

### Perform find and replace of reference genome

#replace ">chromX" with ">X" in fasta file
file_name_temp=str("hg19_chrom"+str(chrom)+".fa")
cmd_temp=("sed 's/chr"+str(chrom)+"/"+str(chrom)+"/g' "+file_name_temp+ ">hg19_chrom"+str(chrom)+"clean.fa")
CheckAndRun()

#Delete the original fasta file
file_name_temp=str("hg19_chrom"+str(chrom)+"clean.fa")
cmd_temp=("rm hg19_chrom"+str(chrom)+".fa")
CheckAndRun()

####################################
###Download VCF files for Vindija###
####################################

#create commands
vindi_cmd=("curl http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr"+str(chrom)+"_mq25_mapab100.vcf.gz> Vindija_chr"+str(chrom)+".vcf.gz")
vindi_cmd_tbi=("curl http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr"+str(chrom)+"_mq25_mapab100.vcf.gz.tbi> Vindija_chr"+str(chrom)+".vcf.gz.tbi")

print("Preparing to download Vindija files with the following commands...")
print(vindi_cmd)
print(vindi_cmd_tbi)

#call the command to download the file
os.system(vindi_cmd)
os.system(vindi_cmd_tbi)

#Remove indels from Vindija file
file_name_temp=str("Vindija_chr"+str(chrom)+".vcf.gz")
cmd_temp=("vcftools --gzvcf "+file_name_temp+" --remove-indels --recode --recode-INFO-all --out Vindija_chr"+str(chrom)+"_SNPs")
CheckAndRun()

#Rezip and index cleaned/"recoded" (i.e. no indels) vcf file
file_name_temp="Vindija_chr"+str(chrom)+"_SNPs.recode.vcf"
cmd_temp=("bgzip "+file_name_temp+"; tabix -p vcf "+file_name_temp+".gz")
CheckAndRun()

#Extract fasta file
file_name_temp="Vindija_chr"+str(chrom)+"_SNPs.recode.vcf.gz"
cmd_temp="cat hg19_chrom"+str(chrom)+"clean.fa | bcftools consensus "+file_name_temp+" >Vindija_chr"+str(chrom)+".fa"
CheckAndRun()

#Print message if fasta file was successfully generated and remove all the temporary files
if os.path.isfile("Vindija_chr"+str(chrom)+".fa") and os.stat("Vindija_chr"+str(chrom)+".fa").st_size > 0:
    print("Vindija_chr"+str(chrom)+".fa"+ " was successfully created!!\n")
    os.system(str("rm Vindija_chr"+str(chrom)+".vcf.gz"))
    os.system(str("rm Vindija_chr" + str(chrom) + ".vcf.gz.tbi"))
    os.system(str("rm Vindija_chr"+str(chrom)+"_SNPs.recode.vcf"))
    os.system(str("rm Vindija_chr"+str(chrom)+"_SNPs.recode.vcf.gz"))
    os.system(str("rm Vindija_chr"+str(chrom)+"_SNPs.recode.vcf.gz.tbi"))
else:
    print("ERROR: Vindija_chr"+str(chrom)+".fa"+" was NOT created\nSorry about that :(\nMoving on...")


####################################
###Download VCF files for Altai ####
####################################
#http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz

#create commands
altai_cmd=("curl http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr"+str(chrom)+"_mq25_mapab100.vcf.gz> Altai_chr"+str(chrom)+".vcf.gz")
altai_cmd_tbi=("curl http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr"+str(chrom)+"_mq25_mapab100.vcf.gz> Altai_chr"+str(chrom)+".vcf.gz.tbi")

print("Preparing to download Altai files with the following commands...")
print(altai_cmd)
print(altai_cmd_tbi)

#call the command to download the file
os.system(altai_cmd)
os.system(altai_cmd_tbi)

#Remove indels from Altai file
file_name_temp=str("Altai_chr"+str(chrom)+".vcf.gz")
cmd_temp=("vcftools --gzvcf "+file_name_temp+" --remove-indels --recode --recode-INFO-all --out Altai_chr"+str(chrom)+"_SNPs")
CheckAndRun()

#Rezip and index cleaned/"recoded" (i.e. no indels) vcf file
file_name_temp="Altai_chr"+str(chrom)+"_SNPs.recode.vcf"
cmd_temp=("bgzip "+file_name_temp+"; tabix -p vcf "+file_name_temp+".gz")
CheckAndRun()

#Extract fasta file
file_name_temp="Altai_chr"+str(chrom)+"_SNPs.recode.vcf.gz"
cmd_temp="cat hg19_chrom"+str(chrom)+"clean.fa | bcftools consensus "+file_name_temp+" >Altai_chr"+str(chrom)+".fa"
CheckAndRun()

#Print message if fasta file was successfully generated
if os.path.isfile("Altai_chr"+str(chrom)+".fa") and os.stat("Altai_chr"+str(chrom)+".fa").st_size > 0:
    print("Altai_chr"+str(chrom)+".fa"+ " was successfully created!!\n")
    os.system(str("rm Altai_chr"+str(chrom)+".vcf.gz"))
    os.system(str("rm Altai_chr" + str(chrom) + ".vcf.gz.tbi"))
    os.system(str("rm Altai_chr"+str(chrom)+"_SNPs.recode.vcf"))
    os.system(str("rm Altai_chr"+str(chrom)+"_SNPs.recode.vcf.gz"))
    os.system(str("rm Altai_chr"+str(chrom)+"_SNPs.recode.vcf.gz.tbi"))
else:
    print("ERROR: Altai_chr"+str(chrom)+".fa"+" was NOT created\nSorry about that :(\nMoving on...")

#########################################
### Download VCF files for Denisovan ####
#########################################
#http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz
#create commands
deni_cmd=("curl http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr"+str(chrom)+"_mq25_mapab100.vcf.gz> Deni_chr"+str(chrom)+".vcf.gz")
deni_cmd_tbi=("curl http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr"+str(chrom)+"_mq25_mapab100.vcf.gz> Deni_chr"+str(chrom)+".vcf.gz.tbi")

print("Preparing to download Denisovan files with the following commands...")
print(deni_cmd)
print(deni_cmd_tbi)

#call the command to download the file
os.system(deni_cmd)
os.system(deni_cmd_tbi)

#Remove indels from Deni file
file_name_temp=str("Deni_chr"+str(chrom)+".vcf.gz")
cmd_temp=("vcftools --gzvcf "+file_name_temp+" --remove-indels --recode --recode-INFO-all --out Deni_chr"+str(chrom)+"_SNPs")
CheckAndRun()

#Rezip and index cleaned/"recoded" (i.e. no indels) vcf file
file_name_temp="Deni_chr"+str(chrom)+"_SNPs.recode.vcf"
cmd_temp=("bgzip "+file_name_temp+"; tabix -p vcf "+file_name_temp+".gz")
CheckAndRun()

#Extract fasta file
file_name_temp="Deni_chr"+str(chrom)+"_SNPs.recode.vcf.gz"
cmd_temp="cat hg19_chrom"+str(chrom)+"clean.fa | bcftools consensus "+file_name_temp+" >Deni_chr"+str(chrom)+".fa"
CheckAndRun()

#Print message if fasta file was successfully generated
if os.path.isfile("Deni_chr"+str(chrom)+".fa") and os.stat("Deni_chr"+str(chrom)+".fa").st_size > 0:
    print("Deni_chr"+str(chrom)+".fa"+ " was successfully created!!\n")
    os.system(str("rm Deni_chr"+str(chrom)+".vcf.gz"))
    os.system(str("rm Deni_chr" + str(chrom) + ".vcf.gz.tbi"))
    os.system(str("rm Deni_chr"+str(chrom)+"_SNPs.recode.vcf"))
    os.system(str("rm Deni_chr"+str(chrom)+"_SNPs.recode.vcf.gz"))
    os.system(str("rm Deni_chr"+str(chrom)+"_SNPs.recode.vcf.gz.tbi"))
else:
    print("ERROR: Deni_chr"+str(chrom)+".fa"+" was NOT created\nSorry about that :(\nMoving on...")

####################################
###Download VCF files for French####
####################################

#create commands
french_cmd=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004468.hg19_1000g."+str(chrom)+".mod.vcf.gz> French_chr"+str(chrom)+".vcf.gz")
french_cmd_tbi=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004468.hg19_1000g."+str(chrom)+".mod.vcf.gz.tbi> French_chr"+str(chrom)+".vcf.gz.tbi")

print("Preparing to download French files with the following commands...")
print(french_cmd)
print(french_cmd_tbi)

#call the command to download the file
os.system(french_cmd)
os.system(french_cmd_tbi)

###Remove indels from French file
file_name_temp=str("French_chr"+str(chrom)+".vcf.gz")
cmd_temp=("vcftools --gzvcf "+file_name_temp+" --remove-indels --recode --recode-INFO-all --out French_chr"+str(chrom)+"_SNPs")
CheckAndRun()

#Rezip and index cleaned/"recoded" (i.e. no indels) vcf file
file_name_temp="French_chr"+str(chrom)+"_SNPs.recode.vcf"
cmd_temp=("bgzip "+file_name_temp+"; tabix -p vcf "+file_name_temp+".gz")
CheckAndRun()

#Extract fasta file
file_name_temp="French_chr"+str(chrom)+"_SNPs.recode.vcf.gz"
cmd_temp="cat hg19_chrom"+str(chrom)+"clean.fa | bcftools consensus "+file_name_temp+" >French_chr"+str(chrom)+".fa"
CheckAndRun()

#Print message if fasta file was successfully generated
if os.path.isfile("French_chr"+str(chrom)+".fa") and os.stat("French_chr"+str(chrom)+".fa").st_size > 0:
    print("French_chr"+str(chrom)+".fa"+ " was successfully created!!\n")
    os.system(str("rm French_chr"+str(chrom)+".vcf.gz"))
    os.system(str("rm French_chr" + str(chrom) + ".vcf.gz.tbi"))
    os.system(str("rm French_chr"+str(chrom)+"_SNPs.recode.vcf"))
    os.system(str("rm French_chr"+str(chrom)+"_SNPs.recode.vcf.gz"))
    os.system(str("rm French_chr"+str(chrom)+"_SNPs.recode.vcf.gz.tbi"))
else:
    print("ERROR: French_chr"+str(chrom)+".fa"+" was NOT created\nSorry about that :(\nMoving on...")

####################################
###Download VCF files for Yoroba####
####################################

#create commands
yoroba_cmd=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004475.hg19_1000g."+str(chrom)+".mod.vcf.gz> Yoroba_chr"+str(chrom)+".vcf.gz")
yoroba_cmd_tbi=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004475.hg19_1000g."+str(chrom)+".mod.vcf.gz.tbi> Yoroba_chr"+str(chrom)+".vcf.gz.tbi")

print("Preparing to download Yoroba files with the following commands...")
print(yoroba_cmd)
print(yoroba_cmd_tbi)

#call the command to download the file
os.system(yoroba_cmd)
os.system(yoroba_cmd_tbi)

###Remove indels from Yoroba file
file_name_temp=str("Yoroba_chr"+str(chrom)+".vcf.gz")
cmd_temp=("vcftools --gzvcf "+file_name_temp+" --remove-indels --recode --recode-INFO-all --out Yoroba_chr"+str(chrom)+"_SNPs")
CheckAndRun()

#Rezip and index cleaned/"recoded" (i.e. no indels) vcf file
file_name_temp="Yoroba_chr"+str(chrom)+"_SNPs.recode.vcf"
cmd_temp=("bgzip "+file_name_temp+"; tabix -p vcf "+file_name_temp+".gz")
CheckAndRun()

#Extract fasta file
file_name_temp="Yoroba_chr"+str(chrom)+"_SNPs.recode.vcf.gz"
cmd_temp="cat hg19_chrom"+str(chrom)+"clean.fa | bcftools consensus "+file_name_temp+" >Yoroba_chr"+str(chrom)+".fa"
CheckAndRun()

#Print message if fasta file was successfully generated
if os.path.isfile("Yoroba_chr"+str(chrom)+".fa") and os.stat("Yoroba_chr"+str(chrom)+".fa").st_size > 0:
    print("Yoroba_chr"+str(chrom)+".fa"+ " was successfully created!!\n")
    os.system(str("rm Yoroba_chr"+str(chrom)+".vcf.gz"))
    os.system(str("rm Yoroba_chr" + str(chrom) + ".vcf.gz.tbi"))
    os.system(str("rm Yoroba_chr"+str(chrom)+"_SNPs.recode.vcf"))
    os.system(str("rm Yoroba_chr"+str(chrom)+"_SNPs.recode.vcf.gz"))
    os.system(str("rm Yoroba_chr"+str(chrom)+"_SNPs.recode.vcf.gz.tbi"))
else:
    print("ERROR: Yoroba_chr"+str(chrom)+".fa"+" was NOT created\nSorry about that :(\nMoving on...")


####################################
### Download VCF files for San  ####
####################################

#create commands
san_cmd=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004473.hg19_1000g."+str(chrom)+".mod.vcf.gz> San_chr"+str(chrom)+".vcf.gz")
san_cmd_tbi=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004473.hg19_1000g."+str(chrom)+".mod.vcf.gz.tbi> San_chr"+str(chrom)+".vcf.gz.tbi")

print("Preparing to download San files with the following commands...")
print(san_cmd)
print(san_cmd_tbi)

#call the command to download the file
os.system(san_cmd)
os.system(san_cmd_tbi)

###Remove indels from San file
file_name_temp=str("San_chr"+str(chrom)+".vcf.gz")
cmd_temp=("vcftools --gzvcf "+file_name_temp+" --remove-indels --recode --recode-INFO-all --out San_chr"+str(chrom)+"_SNPs")
CheckAndRun()

#Rezip and index cleaned/"recoded" (i.e. no indels) vcf file
file_name_temp="San_chr"+str(chrom)+"_SNPs.recode.vcf"
cmd_temp=("bgzip "+file_name_temp+"; tabix -p vcf "+file_name_temp+".gz")
CheckAndRun()

#Extract fasta file
file_name_temp="San_chr"+str(chrom)+"_SNPs.recode.vcf.gz"
cmd_temp="cat hg19_chrom"+str(chrom)+"clean.fa | bcftools consensus "+file_name_temp+" >San_chr"+str(chrom)+".fa"
CheckAndRun()

#Print message if fasta file was successfully generated
if os.path.isfile("San_chr"+str(chrom)+".fa") and os.stat("San_chr"+str(chrom)+".fa").st_size > 0:
    print("San_chr"+str(chrom)+".fa"+ " was successfully created!!\n")
    os.system(str("rm San_chr"+str(chrom)+".vcf.gz"))
    os.system(str("rm San_chr" + str(chrom) + ".vcf.gz.tbi"))
    os.system(str("rm San_chr"+str(chrom)+"_SNPs.recode.vcf"))
    os.system(str("rm San_chr"+str(chrom)+"_SNPs.recode.vcf.gz"))
    os.system(str("rm San_chr"+str(chrom)+"_SNPs.recode.vcf.gz.tbi"))
else:
    print("ERROR: San_chr"+str(chrom)+".fa"+" was NOT created\nSorry about that :(\nMoving on...")

####################################
### Download VCF files for Han  ####
####################################

#create commands
han_cmd=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004469.hg19_1000g."+str(chrom)+".mod.vcf.gz> Han_chr"+str(chrom)+".vcf.gz")
han_cmd_tbi=("curl http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004469.hg19_1000g."+str(chrom)+".mod.vcf.gz.tbi> Han_chr"+str(chrom)+".vcf.gz.tbi")

print("Preparing to download Han files with the following commands...")
print(han_cmd)
print(han_cmd_tbi)

#call the command to download the file
os.system(han_cmd)
os.system(han_cmd_tbi)

###Remove indels from Han file
file_name_temp=str("Han_chr"+str(chrom)+".vcf.gz")
cmd_temp=("vcftools --gzvcf "+file_name_temp+" --remove-indels --recode --recode-INFO-all --out Han_chr"+str(chrom)+"_SNPs")
CheckAndRun()

#Rezip and index cleaned/"recoded" (i.e. no indels) vcf file
file_name_temp="Han_chr"+str(chrom)+"_SNPs.recode.vcf"
cmd_temp=("bgzip "+file_name_temp+"; tabix -p vcf "+file_name_temp+".gz")
CheckAndRun()

#Extract fasta file
file_name_temp="Han_chr"+str(chrom)+"_SNPs.recode.vcf.gz"
cmd_temp="cat hg19_chrom"+str(chrom)+"clean.fa | bcftools consensus "+file_name_temp+" >Han_chr"+str(chrom)+".fa"
CheckAndRun()

#Print message if fasta file was successfully generated
if os.path.isfile("Han_chr"+str(chrom)+".fa") and os.stat("Han_chr"+str(chrom)+".fa").st_size > 0:
    print("Han_chr"+str(chrom)+".fa"+ " was successfully created!!\n")
    os.system(str("rm Han_chr"+str(chrom)+".vcf.gz"))
    os.system(str("rm Han_chr" + str(chrom) + ".vcf.gz.tbi"))
    os.system(str("rm Han_chr"+str(chrom)+"_SNPs.recode.vcf"))
    os.system(str("rm Han_chr"+str(chrom)+"_SNPs.recode.vcf.gz"))
    os.system(str("rm Han_chr"+str(chrom)+"_SNPs.recode.vcf.gz.tbi"))
else:
    print("ERROR: Han_chr"+str(chrom)+".fa"+" was NOT created\nSorry about that :(\nMoving on...")

#################################################################################################################################

#####           Now check if each fasta file is there and concatenate together into one big file if they are                #####

#################################################################################################################################

if (os.path.isfile("Vindija_chr"+str(chrom)+".fa") and os.stat("Vindija_chr"+str(chrom)+".fa").st_size > 0
and os.path.isfile("Altai_chr"+str(chrom)+".fa") and os.stat("Altai_chr"+str(chrom)+".fa").st_size > 0
and os.path.isfile("Deni_chr"+str(chrom)+".fa") and os.stat("Deni_chr"+str(chrom)+".fa").st_size > 0
and os.path.isfile("French_chr"+str(chrom)+".fa") and os.stat("French_chr"+str(chrom)+".fa").st_size > 0
and os.path.isfile("Yoroba_chr"+str(chrom)+".fa") and os.stat("Yoroba_chr"+str(chrom)+".fa").st_size > 0
and os.path.isfile("San_chr"+str(chrom)+".fa") and os.stat("San_chr"+str(chrom)+".fa").st_size > 0
and os.path.isfile("Han_chr"+str(chrom)+".fa") and os.stat("Han_chr"+str(chrom)+".fa").st_size > 0):
    print("All fasta files exist\nCreating concatenated file...")
    #find and replace all fasta files
    os.system("sed 's/" + str(chrom) + "/Vindija_chr" + str(chrom) + "/g' Vindija_chr" + str(chrom) + ".fa >Vindija_chr" + str(chrom) + "clean.fa")
    os.system("sed 's/" + str(chrom) + "/Altai_chr" + str(chrom) + "/g' Altai_chr" + str(chrom) + ".fa >Altai_chr" + str(chrom) + "clean.fa")
    os.system("sed 's/" + str(chrom) + "/Deni_chr" + str(chrom) + "/g' Deni_chr" + str(chrom) + ".fa >Deni_chr" + str(chrom) + "clean.fa")
    os.system("sed 's/" + str(chrom) + "/French_chr" + str(chrom) + "/g' French_chr" + str(chrom) + ".fa >French_chr" + str(chrom) + "clean.fa")
    os.system("sed 's/" + str(chrom) + "/Yoroba_chr" + str(chrom) + "/g' Yoroba_chr" + str(chrom) + ".fa >Yoroba_chr" + str(chrom) + "clean.fa")
    os.system("sed 's/" + str(chrom) + "/San_chr" + str(chrom) + "/g' San_chr" + str(chrom) + ".fa >San_chr" + str(chrom) + "clean.fa")
    os.system("sed 's/" + str(chrom) + "/Han_chr" + str(chrom) + "/g' Han_chr" + str(chrom) + ".fa >Han_chr" + str(chrom) + "clean.fa")

    #Cat all files
    os.system("cat Vindija_chr"+str(chrom)+"clean.fa Altai_chr"+str(chrom)+"clean.fa Deni_chr"+str(chrom)+"clean.fa French_chr"+str(chrom)+"clean.fa Yoroba_chr"+str(chrom)+"clean.fa San_chr"+str(chrom)+"clean.fa Han_chr"+str(chrom)+"clean.fa>ALL_SPECIES_"+str(chrom)+".fa")

else:
    print("ERROR: Not all fasta files exist. Exiting...")



if os.path.isfile("ALL_SPECIES_"+str(chrom)+".fa") and os.stat("ALL_SPECIES_"+str(chrom)+".fa").st_size > 0:
    print("ALL_SPECIES_"+str(chrom)+".fa successfully created!\n...removing temporary fasta files and exiting!")
    os.system("rm Vindija_chr" + str(chrom) + "*.fa")
    os.system("rm Altai_chr" + str(chrom) + "*.fa")
    os.system("rm Deni_chr" + str(chrom) + "*.fa")
    os.system("rm French_chr" + str(chrom) + "*.fa")
    os.system("rm Yoroba_chr" + str(chrom) + "*.fa")
    os.system("rm San_chr" + str(chrom) + "*.fa")
    os.system("rm Han_chr" + str(chrom) + "*.fa")
    #Also remove log files
    os.system("rm *log")

    #Finally, copy output fasta file to a folder outside of /scratch
    os.system("cp ALL_SPECIES_"+str(chrom)+".fa /home/esforsythe/Documents/Download_and_processVCF/Neander180205_FINAL_ALN/")

else:
    print("ERROR: ALL_SPECIES_"+str(chrom)+".fa with NOT created!")


print("Exiting analysis of chromosome: "+str(chrom))

exit()
