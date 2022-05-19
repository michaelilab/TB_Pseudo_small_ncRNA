#SCRIPT TO Extract Psi values with list works with Norm_Uprofile_ files
rm(list=ls(all=TRUE))
#INDICATE pattern for COUNT files, the list of analyzed RNAs read for folder
count_files <- "^Norm_Uprofile_" # Output of the previous script
# INDICATE PATTERN FOR REFERENCE POSITIONS FILE
pos_list_pattern <- "modList"  #actual file: Hs_rRNA_modList.csv in exttDataDir //!\\
# INDICATE RNA_List pattern
#RNA_list_pattern <-"RNA_list"

####### Input area #######
extDataDir <- "Path/to/input"	//!\\
ProjectDir <- "Path/to/output" //!\\

AnnotDir<-"C:/Users/tirza/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/Input"


####### Treatment area ######
setwd(extDataDir)

#Select the Folder Name and use as Project Name
project <-sub('.*\\/', '', extDataDir) #FIND THE LAST "/" in the path and return the right string
project

#################CREATE FOLDER FOR ANALYSIS
dir.create (ProjectDir, showWarnings = FALSE) #CREATE PROJECT FOLDER IN D:/R_folder
#ResultsDir <-paste0(ProjectDir,"/HydraPsiSeq")
ResultsDir <-paste0(ProjectDir,"/rRNA")
#ResultsDir <-paste0(ProjectDir,"/smallRNA")
dir.create (ResultsDir, showWarnings = FALSE)  #CREATE HydraPsiSeq folder
###########################################

#FIND REFERENCE positions FILE

ModList_File<-"DB/TB_rRNA_modList_pseudo_sites_05July2015.txt"

#LIST FOLDERS
DataDirs <- list.dirs(extDataDir, full.names = TRUE, recursive = FALSE)
DataDirs					#FULL path format

ModList<- read.csv(ModList_File, sep="\t", header = TRUE)
ModList$new<-NULL

setwd(ResultsDir)
#CREATE EMPTY FILES TO MERGE DATA by SAMPLE
bilanA<-setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("position", "RNA", "Known_Modif","class", "IVT","snoRNA"))
bilanC<-setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("position", "RNA", "Known_Modif", "class", "IVT","snoRNA"))
bilanM<-setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("position", "RNA", "Known_Modif", "class", "IVT","snoRNA"))

RNA_list<-basename(DataDirs)
#z=DataDirs[1]
for (z in DataDirs) {
  results <- NA

  bilanPsi<-as.data.frame(matrix(,ncol=18,nrow=0)) # with 18 columns, may also define rows

  # RNA reference table creation
  Sample_files <- list.files(z, pattern=count_files, recursive=F, full.names = TRUE)
  RNA_list<-basename(Sample_files)
  type <- as.list (RNA_list)
 
  for (x in Sample_files ){
      y=basename(x)
      print(y)
      RNA_name <-gsub(".csv$","",gsub("^Norm_Uprofile_","",y))  #remove ^Norm_Uprofile_ at the beginning and .csv$ at the end
      print(RNA_name)
      #RNA_file <- paste0(z,"/",y)
      RNA_file <- paste0(z,"/",RNA_list)
      #CREATE SampleDir
      SampleDir <- paste0(ResultsDir,"/", basename(z))
      dir.create (SampleDir, showWarnings = FALSE)

      #Select known mods in RNA 
      #ModListRNA <- ModList[grep(RNA_name,ModList$RNA), ]
      #ModListRNAPsi<- ModListRNA[grep("Psi",ModListRNA$Known_Modif),]

    
      #reading file from csv with names of colums and rows
      sample<-read.csv(x, sep=",", header = TRUE, stringsAsFactors = FALSE)
      print(x)
      sample <- na.omit(sample)
      #colnames(sample)[3]<-c("RNA") #rename column number 3
      names(sample)[names(sample) == 'RNA_name'] <- 'RNA'  #rename column with old name 'RNA_name' to "RNA'
      print (head(basename(z))) 

      #CALCULATIONS OF TREATMENT Table
      rnames <- sample$position
      row.names(sample) <- rnames		# rownames reorganisation
  

      #results <- merge(sample,ModListRNAPsi[,c(2:6)], by=c("position"), all=TRUE)  
      results<-sample
      results$NormUcount [is.na(results$NormUcount)] <- 1	# NA replaced by 1 to avoid crashes
      n1<-nrow(results)

      results$ratio <- NA
  
      for (i in 2:n1) {results$ratio[i] = (results$NormUcount[i]/results$NormUcount[i-1])}		# Ratio of counts from left to right

      results$around <- NA
      for (i in 3:n1) {
          j <- i-2
          k <- i-1 
          l <- i+1 
          m <- i+2
          results$around[i] = (1 - results$ratio[i]/(mean(results$ratio[j:k])+ mean (results$ratio[l:m])))}		# Ratio normalization

      results$deviation <- NA
      for (i in 4:n1) {
              j <- i-3 
              k <- i-1 
              l <- i+1 
              m <- i+3
              results$deviation[i] = (1 - results$around[i]/(sd(results$around[j:k])+ sd(results$around[l:m])))
            }	# Ratio standard deviation (not needed)

      results$ratio2 <- NA	
      for (i in 2:n1) {results$ratio2[i] = 1/results$ratio[i]}		# Ratio from right to left
              
      # Ratio2 normalization
      results$around2 <- NA
      for (i in 3:n1) {
                    j <- i-2
                    k <- i-1 
                    l <- i+1 
                    m <- i+2
                    results$around2[i] = (1 - results$ratio2[i]/(mean(results$ratio2[j:k])+ mean (results$ratio2[l:m])))
                  }		# Ratio2 normalization

      results$cumul <- NA
      for (i in 3:n1) {results$cumul[i] = (results$around[i]+results$around2[i+1])/2}		# Mean of the two normalization

      results$mean <- NA
      for (i in 3:n1) {results$mean[i] = results$cumul[i]}		# Score MEAN calculation

      results$scoreA <- NA
      for (i in 3:n1) {
                      j = i-2
                      k = i-1
                      l = i+1
                      m = i+2
                      results$scoreA[i] = 1-(2*results$NormUcount[i]+1)/(0.5*abs(mean(results$NormUcount[j:k])-sd(results$NormUcount[j:k]))+results$NormUcount[i]+0.5*abs(mean(results$NormUcount[l:m])-sd(results$NormUcount[l:m])))
                    }			# Score A calculation

      results$scoreB <- NA
      for (i in 3:n1) {
                      e = i-2
                      f = i-1
                      j = i+1
                      k = i+2
                      results$scoreB[i] = abs(results$NormUcount[i]-0.5*((results$NormUcount[e]*0.9+results$NormUcount[f]*1)/(0.9+1)+(results$NormUcount[j]*1+results$NormUcount[k]*0.9)/(0.9+1)))/(results$NormUcount[i]+1)
                    }		# Score B calculation

      results$scoreC <- NA
      for (i in 3:n1) {
                        e = i-2
                        f = i-1
                        j = i+1
                        k = i+2
                        results$scoreC[i] = 1-results$NormUcount[i]/(0.5*((results$NormUcount[e]*0.9+results$NormUcount[f]*1)/(0.9+1)+(results$NormUcount[j]*1+results$NormUcount[k]*0.9)/(0.9+1)))
                    }		# Score C calculation
      #print(results)
      #WRITE TREATMENT FILE
      #write.table (results, file=paste0(SampleDir,"/Treatment_",basename(z),"_", RNA_name,".csv"), sep = ",", row.names = FALSE, quote = FALSE)
      #write.table (results, file=paste0(ResultsDir,"/Treatment_",basename(z),".csv"), sep = ",", row.names = FALSE, quote = FALSE,append =TRUE)
      #write.table(results, file=paste0(ResultsDir,"/Treatment_",basename(z),".csv"), sep = ",",row.names = FALSE, quote = FALSE, append=TRUE, col.names=!file.exists(file=paste0(ResultsDir,"/Treatment_",basename(z),".csv")))
      write.table(results, file=paste0(SampleDir,"/Treatment_",basename(z),".csv"), sep = ",",row.names = FALSE, quote = FALSE, append=TRUE, col.names=!file.exists(file=paste0(ResultsDir,"/Treatment_",basename(z),".csv")))      
      #print(RNA_name)
  }#for x
}#for z

####################SELECT ONLY MODIFIED position with list

best_scorePsi <- merge(ModListRNAPsi, results, by =c("position","RNA","Known_Modif", "class", "IVT", "snoRNA"))
best_scorePsi$X<-NULL


bilanPsi <- rbind(best_scorePsi,bilanPsi)  #compilation Psi modifs by RNA

#} #end cycle y in type

if (length(row.names(bilanPsi)) != 0) {
  write.table (bilanPsi, file = paste0 (ResultsDir,"/", basename(z), "/SummaryCandidatesPsi_", basename(z),".csv"), sep = " ", row.names = FALSE , quote = F)
                                      }
#Extract Scores Mean, ScoreA and ScoreC

bilanScoreC <- bilanPsi[,c("position","RNA","Known_Modif","class", "IVT", "snoRNA","scoreC")]
colnames(bilanScoreC)[7]<- paste0("C_",basename(z))
bilanScoreA <- bilanPsi[,c("position","RNA","Known_Modif","class", "IVT", "snoRNA","scoreA")]
colnames(bilanScoreA)[7]<- paste0("A_", basename(z))
bilanScoreMean <- bilanPsi[,c("position","RNA","Known_Modif","class", "IVT", "snoRNA","mean")]
colnames(bilanScoreMean)[7]<- paste0("M_", basename(z))

#MERGE in TABLES by Score Mean, ScoreA and ScoreC
bilanA <-merge(bilanA, bilanScoreA, by =c("position","RNA", "Known_Modif","class", "IVT", "snoRNA"), all=TRUE) 
bilanC <-merge(bilanC, bilanScoreC, by =c("position","RNA", "Known_Modif","class", "IVT", "snoRNA"), all=TRUE) 
bilanM <-merge(bilanM, bilanScoreMean, by =c("position","RNA", "Known_Modif","class", "IVT", "snoRNA"), all=TRUE) 



} #end cycle z in Datadirs

bilanA <- bilanA[order( bilanA$RNA, bilanA$position ),]
bilanC <- bilanC[order( bilanC$RNA, bilanC$position ),]
bilanM <- bilanM[order( bilanM$RNA, bilanM$position ),]


if (length(row.names(bilanA)) != 0) {
  write.table (bilanA, file = paste0 (ResultsDir, "/SummaryCandidates_ScoreA_", project,".csv"), sep = ",", row.names = FALSE , quote = F)
                                     }

if (length(row.names(bilanC)) != 0) {
  write.table (bilanC, file = paste0 (ResultsDir, "/SummaryCandidates_PsiScore", project,".csv"), sep = ",", row.names = FALSE , quote = F)
}
  
if (length(row.names(bilanM)) != 0) {
    write.table (bilanM, file = paste0 (ResultsDir, "/SummaryCandidates_ScoreMean_", project,".csv"), sep = ",", row.names = FALSE , quote = F)
} 


#Variable bilanC to use for heatmap
#PREPARATION FOR HEATMAP SCRIPT

#FORMAT FILE FinalTable_PsiScore_MixProjects_New_pipeline.csv (data, name, Mewname, snoRNA)

PsiScore <- bilanC [, grep("Sample", colnames(bilanC))]  #data part only
PsiLevel <- (PsiScore-bilanC$IVT)/(1-bilanC$IVT)  # Calculate PsiLevel using IVT data

#ADD columns name, Newname = name + position, snORNA
PsiScore$name <- paste0(bilanC$RNA,"_",bilanC$Known_Modif)
PsiScore$Newname <- paste0(PsiScore$name,'_',bilanC$position)
PsiScore$snoRNA <- bilanC$snoRNA

PsiLevel$name <- paste0(bilanC$RNA,"_",bilanC$Known_Modif)
PsiLevel$Newname <- paste0(PsiLevel$name,'_',bilanC$position)
PsiLevel$snoRNA <- bilanC$snoRNA


if (length(row.names(PsiScore)) != 0) {
  write.table (PsiScore, file = paste0 (ResultsDir, "/FinalTable_PsiScore_", project,".csv"), sep = ",", row.names = FALSE , quote = F)
} 


if (length(row.names(PsiLevel)) != 0) {
  write.table (PsiLevel, file = paste0 (ResultsDir, "/FinalTable_PsiLevel_", project,".csv"), sep = ",", row.names = FALSE , quote = F)
} 



# END
