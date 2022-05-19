##########################################################################
#this script should be run from the directory containing the count files
#it calculates the scores A,B and C and writes them into a table
##########################################################################


directory <- "RiboMethSeq/RiboMethSeq_01March2022/TB_rRNA/Input"
setwd(directory)


ResultsDir<- "RiboMethSeq/RiboMethSeq_01March2022/TB_rRNA/Output"


annot_file="DB/TB/TB_rRNA_annot_Nm.txt"
annot<-read.table(annot_file, sep="\t", header = TRUE)

### dependencies
if (!require(seqinr)) install.packages("seqinr")
library(seqinr)
### parameters

win_size=6
rRNA_length=11171

### references 

myfasta<-read.fasta("DB/TB_rRNA_chr2.fa",as.string = FALSE,set.attributes = FALSE)
### read all of the outputs from pre-processing
mylibs=gsub(pattern = ".sorted.init",replacement = "",x=list.files(pattern = "init"))

for (thisLib in mylibs){
  #read the table as it came from the pre-processing , must include the counts at the third column
  #or just change the code from [,3] to whatever column you need
  pre_myinit<-read.table(paste0(thisLib,".sorted.init"), header = FALSE)[,3]
  pre_myends<-read.table(paste0(thisLib,".sorted.3p"), header = FALSE)[,3]
  #shifting reads 1 bp 
  myinit<-c(0,pre_myinit,rep(0,rRNA_length-length(pre_myinit)-1))
  myends<-c(pre_myends[2:length(pre_myends)],rep(0,rRNA_length-length(pre_myends)+1))
  mycov<-myinit+myends
  mylength=length(myinit)
  Sc<-rep(NA,mylength)
  Sa<-rep(NA,mylength)
  Sb<-rep(NA,mylength)
  
  #sum the weights for each position - avoid overweighted results
  W<-(1+(1-0.1*(win_size)))*win_size/2
  
  
  #calculate the scores
  
  for (i in (win_size+1):(mylength-win_size-1)){
    # A score
    #  calculate average and sd for each side 
    M_l<-mean(mycov[(i-win_size/2):(i-1)])
    S_l<-sd(mycov[(i-win_size/2):(i-1)])
    
    M_r<-mean(mycov[(i+1):(i+win_size/2)])
    S_r<-sd(mycov[(i+1):(i+win_size/2)])
    
    Sa[i]<-max(0,1-(2*mycov[i]+1)/(0.5*abs(M_l-S_l)+mycov[i]+0.5*abs(M_r-S_r)+1))
    
    # B +C scores
    
    
    S1<-0
    for (j in 1:win_size){
      S1<-S1+(1-0.1*(j-1))*mycov[i-j]
    }
    
    S1<-S1/W
    
    S2<-0
    for (j in 1:win_size){
      S2<-S2+(1-0.1*(j-1))*mycov[i+j]
    }
    
    S2<-S2/W
    
    Sc[i]<-max(0,1-2*mycov[i]/(S1+S2))
    
    Sb[i]<-abs((mycov[i]-0.5*(S1+S2))/(mycov[i]+1))
  }
  
  
 # extracting fasta base by base to be printed within the final table
#doesnt work- need to fix line 
  bp<-(unlist(unname(myfasta[1]),rep(NA,rRNA_length-length(unlist(unname(myfasta[1])))))) #move base pair by one in order to know that is the bp
  mydf<-cbind.data.frame(myinit,myends,mycov,Sa,Sb,Sc,bp)
  #colnames(mydf)<- c("5p","3p","cov","Sa","Sb","Sc","bp")
  #write.csv(x =mydf ,file=paste0(ResultsDir,thisLib,".csv"))
  
  mydf<-cbind.data.frame(mydf,annot)
  colnames(mydf)<- c("5p","3p","cov","Sa","Sb","Sc","bp","rRNA","Loc","Count","modified",	"snoRNA",	"BP")
  write.csv(x =mydf ,file=paste0(ResultsDir,thisLib,".csv"))
}
