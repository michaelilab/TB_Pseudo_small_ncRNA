#SCRIPT TO calculate Normalized U profile
rm(list=ls(all=TRUE)) #REMOVE ALL Variables

####### Input area #######
##path to input LAPTOP
#extDataDir <- "C:/Users/tirza/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/Input"
##path to output
#ProjectDir <- "C:/Users/tirza/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/"

#biu

##goodpairs
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/Input"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/"

#R1
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/R1"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/R1"
#R2
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/R2"

#Ld rRNA
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LD_rRNA/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LD_rRNA/GoodPairs_Results"
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LD_rRNA/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LD_rRNA/R2_Results"

##goodpairs
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_21Feb2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_21Feb2021/GoodPairs_Results"

##R2
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_21Feb2021/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_21Feb2021/R2_Results"

##goodpairs
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_21Feb2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_21Feb2021/GoodPairs_Results"

##R2
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_21Feb2021/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_21Feb2021/R2_Results"


##goodpairs
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_21Feb2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_21Feb2021/GoodPairs_Results"

##R2
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_21Feb2021/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_21Feb2021/R2_Results"

#C:\Users\User\OneDrive - Bar Ilan University\Bioinfo Computer\Shula\HydraPsi\Sequencing_21March2021\LD_rRNA
##goodpairs
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LD_rRNA/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LD_rRNA/GoodPairs_Results"

##R2
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LD_rRNA/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LD_rRNA/R2_Results"

#C:\Users\User\OneDrive - Bar Ilan University\Bioinfo Computer\Shula\HydraPsi\Sequencing_21March2021\LM_rRNA
##goodpairs
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LM_rRNA/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LM_rRNA/GoodPairs_Results"

##R2
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LM_rRNA/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Sequencing_21March2021/LM_rRNA/R2_Results"

#July 15 2021
##goodpairs
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_15July2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_15July2021/GoodPairs_Results"

#Aug09 2021
extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_09Aug2021/GoodPairs"
ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_09Aug2021/GoodPairs_Results"

#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_09Aug2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_09Aug2021/GoodPairs_Results"
#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_09Aug2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/LM_rRNA_09Aug2021/GoodPairs_Results"

#sept 14, 2021

#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_14Sept2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_14Sept2021/GoodPairs_Results"

extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_14Sept2021/R2"
ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Hs_rRNA_14Sept2021/R2_Results"

#sept 14, 2021

#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_14Sept2021/GoodPairs"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_14Sept2021/GoodPairs_Results"

#extDataDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_14Sept2021/R2"
#ProjectDir <- "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_14Sept2021/R2_Results"

#JAn 09 2021

extDataDir <- "C:/Users/tirza/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_09Jan2022/GoodPairs"
ProjectDir <- "C:/Users/tirza/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/TB_rRNA_09Jan2022/GoodPairs_Results"

###########################################################
annot_file="C:/Users/tirza/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/TB_rRNA_annot.txt"
#annot_file="C:/Users/user/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/LD_rRNA_annot.txt"
#annot_file="C:/Users/user/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/LM_rRNA_annot.txt"
#annot_file="C:/Users/user/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/LM_rRNA_annot.txt"

#annot_file="C:/Users/user/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/Hs_rRNA_annot.txt"
annot<-read.table(annot_file, sep="\t", header = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {install.packages("ggplot2")}
library(ggplot2)   

###### Treatment area ######
setwd(extDataDir)

project <-sub('.*\\/', '', extDataDir)  #FIND THE LAST "/" in the path and return the right string
print (extDataDir)
print (project)

#################CREATE FOLDER FOR ANALYSIS
dir.create (ProjectDir, showWarnings = FALSE)
ResultsDir <-paste0(ProjectDir,"/PreAnalysis") # Sub output folder
dir.create (ResultsDir, showWarnings = FALSE)  #CREATE PreAnalysis folder
###########################################

####### MAIN BLOCK #######
#LIST FOLDERS TO TREAT to check if OK
DataDirs <- list.dirs(extDataDir, full.names = TRUE, recursive = FALSE)
DataDirs					#FULL path format
#INDICATE pattern for COUNTand COVERAGE  files, the list of analyzed RNAs is the same for all folders
#count_files <- "^UCount5prime_" #indicate just a pattern to use
coverage_files <- "^coverage_" #indicate just a pattern to use

#mylibs=gsub(pattern = ".sorted.init",replacement = "",x=list.files(pattern = "init"))
count_files<- ".sorted.init$"
#coverage_files <- ".sorted.genomecov$"

#z=DataDirs[1] #for tests

ReadsStat<- data.frame(matrix(ncol = 0, nrow = 0))

for (z in DataDirs) {
# READ RNA LIST, create the list of RNA sequences to treat
Sample_files <- list.files(z, pattern=count_files, recursive=F, full.names = TRUE)
#RNA_list<-basename(Sample_files)
#type <- gsub(".csv$","",gsub("^UCount5prime_","", as.list(RNA_list)))
#type<-
#type #PRINT RNA list to check

ReadsStatSample<- data.frame(matrix(ncol = 1, nrow = 0))
colnames(ReadsStatSample)<-c("RNA")

#y=type[1] #for tests
y="rRNA_good_pairs"#only type
#y="rRNA"
#y="rRNA_R2"
#for (y in type) {
  Count_file    <- list.files(z, pattern=paste0("*",y,count_files), recursive=F, full.names = TRUE) #".*" one ore more any characters
  Coverage_file <- list.files(z, pattern=paste0(coverage_files,"*"), recursive=F, full.names = TRUE)  
  statRNA <-NA

  #reading file from csv with names of colums and rows
  sample=read.csv(Count_file, sep="\t", header = FALSE)
  sample <- sample[,c(2:3)] #keep only col 2 and 3
  sample <- na.omit(sample)
  colnames(sample) <-c("position","counts")

  #READING COVERAGE FILE$
  coverage=read.csv(Coverage_file, sep="\t", header = FALSE)
  coverage<-coverage[,c(1:3)] #skip last column with coverage value
  colnames(coverage) <-c("RNA_name","position","base")
  coverage$base<-gsub("t","T",coverage$base)  #replace eventual t in the coverage file sequence


  FullRNA<-merge(coverage,sample, by="position", all=TRUE) #merge coverage and counts
  FullRNA<-FullRNA[-1,] #SKIP First row
  FullRNA[is.na(FullRNA)] <- 0 #replace missing values by zero

  statRNA[1]<-sum(FullRNA$counts) #NUMBER OF COUNTS for RNA
  statRNA[2] <- nrow(FullRNA[which(FullRNA$counts == 0),])  #ZERO COVERAGE positions


  #NORMALIZATION by window of 11 with non-T signals (special treatment for 5 first and 5 last positions)

  FullRNA$NormMedian <- NA
  FullRNA$Nwindow <- NA

  for (x in 6:(nrow(FullRNA)-5)) {
    region <-FullRNA[c((x-5):(x+5)),]
    NonT_region<-region[which(region$base!="T"),]
    FullRNA$NormMedian[x]<-median(NonT_region[,4], na.rm = T)
    FullRNA$Nwindow[x] <- nrow(NonT_region)
  }

    #TREAT FIRST AND LAST ROWs
    regionS <- FullRNA[c(1:6),]
    FullRNA$NormMedian[1:5]<-median(regionS[which(regionS$base!="T"),][,4], na.rm = T)

    regionL <- FullRNA[c(nrow(FullRNA)-4):nrow(FullRNA),]
    FullRNA$NormMedian[(nrow(FullRNA)-4):(nrow(FullRNA))] <- median(regionL[which(regionL$base!="T"),][,4], na.rm = T)

    FullRNA$NormUcount <- FullRNA$counts/FullRNA$NormMedian

    #drop non T data
    Uprofile <-FullRNA[which(FullRNA$base=="T"),]
    Uprofile$NormUcount[is.infinite(Uprofile$NormUcount)] <- 0 #REMOVE inf
    Uprofile$NormUcount[is.nan(Uprofile$NormUcount)] <- 0 #REMOVE NAN

    statRNA[3]<-sum(Uprofile$counts) #NUMBER OF COUNTS for RNA
    statRNA[4] <- nrow(Uprofile[which(Uprofile$counts == 0),])  #ZERO COVERAGE positions

    dataStatRNA<-data.frame(matrix(ncol = 1, nrow = 4))
    colnames(dataStatRNA)<-c("stat")
    row.names(dataStatRNA) <-c(paste0(y,"_total_reads"),paste0(y,"_pos_Zero_cov"),paste0(y,"_totalU_reads"), paste0(y,"_posU_Zero_cov") )
    dataStatRNA$stat <- statRNA

    ReadsStatSample<-rbind(ReadsStatSample,dataStatRNA )


    #SAVE DATA TO DISK in extDataDir

    output_file <- paste0(z,"/Norm_Uprofile_",y,".csv")
    write.csv(Uprofile,file=output_file, quote =F)
    output_file <- paste0(z,"/FullRNA_profile_",y,".csv")
    FullRNA1<- merge(FullRNA,annot, by="position", all=TRUE) #merge coverage and count
    write.csv(FullRNA1,file=output_file, quote =F)
    #write.csv(FullRNA,file=output_file, quote =F)


    #CREATE GRAPHICS
    #DISTRIBUTION OF different values
    #SAVE IN RESULTS DIR

    SampleDir <- paste0(ResultsDir,"/", basename(z))
    dir.create (SampleDir, showWarnings = TRUE) #CREATE Sample FOLDER IN Preanalysis


    df <- FullRNA[,c(4:7)]

    msize.x=1 #parameters used in inch
    msize.y=1
    pdf.width=2.5
    pdf.height=2.5

    p<-ggplot(df, aes(x=Nwindow)) + 
      geom_histogram(binwidth=1,color="black", fill="white")+
      labs(title=y)

    pdf(paste0(SampleDir,"/NwindowDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
    print(p)
    dev.off()

    p<-ggplot(df, aes(x=NormMedian)) + 
      geom_histogram(binwidth=1,color="black", fill="lightblue")+
      labs(title=y)

    pdf(paste0(SampleDir,"/NormMedianDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
    print(p)
    dev.off()

    msize.x=1 #parameters used in inch
    msize.y=1
    pdf.width=4
    pdf.height=4
    p<-ggplot(df, aes(x=NormUcount)) + 
    #geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666")+
      labs(title=y)
    pdf(paste0(SampleDir,"/NormUcountDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
    print(p)
    dev.off()

    df <- Uprofile[,c(4:7)]
p<-ggplot(df, aes(x=Nwindow)) + 
  geom_histogram(binwidth=1,color="black", fill="white")+
  labs(title=y)
pdf(paste0(SampleDir,"/U_NwindowDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

p<-ggplot(df, aes(x=NormMedian)) + 
  geom_histogram(binwidth=1,color="black", fill="lightblue")+
  labs(title=y)
pdf(paste0(SampleDir,"/U_NormMedianDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

msize.x=1 #parameters used in inch
msize.y=1
pdf.width=4
pdf.height=4
p<-ggplot(df, aes(x=NormUcount)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title=y)
pdf(paste0(SampleDir,"/U_NormUcountDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

#BLOCK FOR rRNA ONLY

#if (y=="sc_rRNA_18S"| y=="sc_rRNA_25S" | y=="sc_rRNA_5.8S" | y=="sc_rRNA_5S") {

df<-FullRNA
  msize.x=1 #parameters used in inch
  msize.y=1
  pdf.width=4
  pdf.height=4
  p<-ggplot(df, aes(x=NormUcount,colour=base, fill=base)) + 
    #geom_histogram(aes(y=..density..))+
    geom_density(alpha=.3)+
    scale_fill_manual(values=c("red", "blue", "green", "yellow","black" ))+
    labs(title=y)
  p
  pdf(paste0(SampleDir,"/AllBasesDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
  print(p)
  dev.off()
  
  p<-ggplot(df, aes(x=NormUcount,colour=base, fill=base)) + 
    #geom_histogram(aes(y=..density..))+
    geom_density(alpha=.3)+
    scale_fill_manual(values=c("red", "blue", "green", "yellow","black" ))+
    coord_cartesian(xlim=c(0,20))+
    labs(title=y)
  p
  pdf(paste0(SampleDir,"/ZoomAllBasesDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
  print(p)
  dev.off()
#}



  colnames(ReadsStatSample) <- c(paste0(basename(z)))
  output_file <- paste0(SampleDir,"/Stat_",basename(z),".csv")
  write.csv(ReadsStatSample,file=output_file, quote =F)

  ReadsStat <-merge(ReadsStat,ReadsStatSample, by="row.names", all=T)
  row.names(ReadsStat)<-ReadsStat$Row.names
  ReadsStat$Row.names <-NULL


}


output_file <- paste0(ResultsDir,"/GlobalStat.csv")
write.csv(ReadsStat,file=output_file, quote =F)



#END


