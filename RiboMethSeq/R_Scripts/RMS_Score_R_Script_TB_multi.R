##########################################################################
#this script should be run from the directory containing the count files
#it calculates the scores A,B and C and writes them into a table
##########################################################################

directory <- "Input"
setwd(directory)

ResultsDir<- "Output"

myfasta<-read.fasta("DB/TB_small_RNAs_DB.fa",as.string = FALSE,set.attributes = FALSE)

### dependencies
if (!require(seqinr)) install.packages("seqinr")
library(seqinr)
### parameters
win_size=6
mylibs=gsub(pattern = ".sorted.init",replacement = "",x=list.files(pattern = "init"))

#for testing
#thisLib="NOPI-PRSmtet-2min01032022_S7_vs_smallRNA_good_pairs"

for (thisLib in mylibs){
  
  init_file<-read.table(paste0(thisLib,".sorted.init"), header = FALSE)  
  ends_file<-read.table(paste0(thisLib,".sorted.3p"), header = FALSE)
 
  gene_list<-unique(init_file[1])
  #gene= gene_list$V1[1]
  for (gene in gene_list$V1) {
    print(gene)
    pre_myinit<-init_file[init_file$V1==gene,3]
    pre_myends<-ends_file[ends_file$V1==gene,3]
    #shifting reads 1 bp 
    myinit<-c(0,pre_myinit,rep(0,length(pre_myinit)-length(pre_myinit)))
    mydiff<-length(pre_myinit)-length(pre_myends)+2
    print(mydiff)
    myends<-c(pre_myends[2:length(pre_myends)],rep(0,length(pre_myends)-length(pre_myends)+mydiff))
    mycov<-myinit+myends
    mylength=length(myinit)
    g<-rep(gene,mylength)
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
    bp<-(unlist(unname(myfasta[gene]),rep(NA,length(myinit)-length(unlist(unname(myfasta[gene][[1]])))))) #move base pair by one in order to know that is the bp
   
    bp<-c(NA,bp)
    mydf<-cbind.data.frame(g,myinit,myends,mycov,Sa,Sb,Sc,bp)
    
    #mydf<-cbind.data.frame(mydf,annot)
    #colnames(mydf)<- c("5p","3p","cov","Sa","Sb","Sc","bp","rRNA","Loc","Count","modified",	"snoRNA",	"BP")
    colnames(mydf)<- c("Gene","5p","3p","cov","Sa","Sb","Sc","bp")
    write.csv(x =mydf ,file=paste0(ResultsDir,thisLib,"_",gene,".csv"))
  }  #for each gene
}#for each sample
  
