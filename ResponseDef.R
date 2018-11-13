#Nov. 13th, 2018
#2018-11-13 11:14:17 EST
#Scipher Medicine/MAD Team

################################################################################
#The script reads the Corrona metadata file from the Scipher Medicine Dropbox 
#source, reformats the data from long to wide and derives some additional 
#disease activity quantities, including endpoint CDAI, endpoint das28_crp, 
#ACR50 and MCID.
#The Corrona matadata file features variables whose values change from 0, to 3
#and/or 6 months,
################################################################################
#reading the latest Corrona metadata file from source.
CorronaMeta <- read.csv("/Users/manuel/Dropbox (DZZOM)/Scipher Team Folder/
                        _DEVELOPMENT (DX)/
                        Corrona/Observational Samples/
                        Scipher_data_export_11_08_2018_185subjects.csv", 
                        header = TRUE, na.strings = c("NA", ""))
#dim(CorronaMeta)
#[1] 518  44
#length(unique(CorronaMeta$subject_id))
#[1] 185
#~There are values for 0, 3 and 6 months.  For 6 months, values only 
#~for 148 patients.


###~~~Data reformating: from long to wide

#checking which variable changes values at 0 and 3 months.
#retrieving the 1st subset of variables whose value remains the same at 0 and 
#3 months
ContantVar1 <- colnames(CorronaMeta)[sapply(
  1:dim(CorronaMeta)[2], function(i) sapply(
    colnames(CorronaMeta), function(x) 
      sum(as.numeric((na.omit(
        CorronaMeta[(CorronaMeta$visit ==0),][,c(x)]) != na.omit(
          CorronaMeta[(CorronaMeta$visit == 3),][,c(x)])))))[[i]] == 0)]

#temporarily subsetting the CorronaMeta df for records at 0 and 6 months
CorronaMetaTmp <- CorronaMeta[(
  CorronaMeta$subject_id %in% 
    CorronaMeta[(CorronaMeta$visit == 6),]$subject_id),]

#retrieving the 2nd subset of variables whose value remains the same at 0 
#and 6 months
ContantVar2 <- colnames(CorronaMetaTmp)[sapply(
  1:dim(CorronaMetaTmp)[2], function(j) sapply(
    colnames(CorronaMetaTmp), function(x) sum(
      as.numeric((na.omit(
        CorronaMetaTmp[(CorronaMetaTmp$visit == 0),][,c(x)]) != na.omit(
          CorronaMetaTmp[(CorronaMetaTmp$visit == 6),][,c(x)])))))[[j]] == 0)]

#instantiating a list on constant variables
ContantVar <- intersect(ContantVar1, ContantVar2)

#subsetting the Corrona metadata for these non variant attributes 
CorronaMetaCore <- CorronaMeta[,ContantVar]

#removing the eular_dascrp attribute 
CorronaMetaCore <- CorronaMetaCore[,!(names(CorronaMetaCore) %in%
                                        c("eular_dascrp"))]
#setting the CorronaMetaCore as non-redundant data set
CorronaMetaCore <- unique(CorronaMetaCore)

#subsetting the Corrona metadata for variant attributes 
CorronaVar <- CorronaMeta[,names(CorronaMeta) %in% c("subject_id") |
                            !names(CorronaMeta) %in% colnames(CorronaMetaCore)]

#subsetting for visit at baseline
CorronaVar_0m <- CorronaVar[(CorronaVar$visit == 0),]
#adding a 0m suffix
colnames(CorronaVar_0m)[-c(1,2)] <-
  paste(colnames(CorronaVar_0m)[-c(1,2)], "_0m", sep = "")
#removing the visit variable
CorronaVar_0m$visit <- NULL

#subsetting for visit at 3 months
CorronaVar_3m <- CorronaVar[(CorronaVar$visit == 3),]
#adding a 3m suffix
colnames(CorronaVar_3m)[-c(1,2)] <- 
  paste(colnames(CorronaVar_3m)[-c(1,2)], "_3m", sep = "")
#removing the visit variable
CorronaVar_3m$visit <- NULL

#subsetting for visit at 6 months
CorronaVar_6m <- CorronaVar[(CorronaVar$visit == 6),]
#adding a 6m suffix
colnames(CorronaVar_6m)[-c(1,2)] <- 
  paste(colnames(CorronaVar_6m)[-c(1,2)], "_6m", sep = "")
#removing the visit variable
CorronaVar_6m$visit <- NULL

#instantiating a list object with the 4 dataframes
metaList <-list(CorronaMetaCore, CorronaVar_0m, CorronaVar_3m, CorronaVar_6m)
#joining
CorronaMetaWide <- join_all(metaList, by = 'subject_id', type = 'full')

#dim(CorronaMetaWide)
#185 101

###~~~End of data reformating
################################################################################
###~~~Response Definitions calculation

####CDAI: Clinical  Disease  Activity  Index  
CorronaMetaWide <- mutate(
  CorronaMetaWide, endpoint_cdai_3m = case_when(cdai_3m <= 10 ~ "Responder", 
                                                cdai_3m > 10 ~ "Non Responder"))
CorronaMetaWide <- mutate(
  CorronaMetaWide, endpoint_cdai_6m = case_when(cdai_6m <= 10 ~ "Responder", 
                                                cdai_6m > 10 ~ "Non Responder"))

#####ENDPOINT_DAS28CRP
CorronaMetaWide <- mutate(
  CorronaMetaWide, 
  endpoint_das28crp_3m = case_when(das28crp_3m <= 2.9 ~ "Responder", 
                                   das28crp_3m > 2.9 ~ "Non Responder"))
CorronaMetaWide <- mutate(
  CorronaMetaWide,
  endpoint_das28crp_6m = case_when(das28crp_6m <= 2.9 ~ "Responder", 
                                   das28crp_6m > 2.9 ~ "Non Responder"))

#####ENDPOINT_ACR50
CorronaMetaWide <- mutate(
  CorronaMetaWide, 
  acr50_3m = case_when(
    CorronaMetaWide$swollen_jts_28_3m <= .5 * CorronaMetaWide$swollen_jts_28_0m 
    | CorronaMetaWide$tender_jts_28_3m <= .5 * CorronaMetaWide$tender_jts_28_0m 
    & sum(as.numeric(CorronaMetaWide$di_0m <= .5 * CorronaMetaWide$di_3m),
      as.numeric(CorronaMetaWide$pt_pain_0m <= .5 * CorronaMetaWide$pt_pain_3m), 
      as.numeric(CorronaMetaWide$pt_global_assess_0m <= .5 * 
                   CorronaMetaWide$pt_global_assess_3m), as.numeric(
                     CorronaMetaWide$md_global_assess_0m <= .5 * 
                       CorronaMetaWide$md_global_assess_3m), 
      as.numeric(CorronaMetaWide$usresultsCRP_0m <= .5 * 
                   CorronaMetaWide$usresultsCRP_3m)) < 3 ~ "Responder", 
    CorronaMetaWide$swollen_jts_28_3m > .5 * CorronaMetaWide$swollen_jts_28_0m 
    | CorronaMetaWide$tender_jts_28_3m > .5 * CorronaMetaWide$tender_jts_28_0m &
      sum(as.numeric(CorronaMetaWide$di_0m <= .5 * CorronaMetaWide$di_3m), 
          as.numeric(CorronaMetaWide$pt_pain_0m <= .5 * 
                       CorronaMetaWide$pt_pain_3m), 
          as.numeric(CorronaMetaWide$pt_global_assess_0m <= .5 * 
                       CorronaMetaWide$pt_global_assess_3m), 
          as.numeric(CorronaMetaWide$md_global_assess_0m <= .5 * 
                       CorronaMetaWide$md_global_assess_3m), 
          as.numeric(CorronaMetaWide$usresultsCRP_0m <= .5 * 
                       CorronaMetaWide$usresultsCRP_3m)) >= 3  ~ "Non Responder"))

CorronaMetaWide <- mutate(CorronaMetaWide, acr50_6m = case_when(CorronaMetaWide$swollen_jts_28_6m <= .5 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_6m <= .5 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .5 * CorronaMetaWide$di_6m), as.numeric(CorronaMetaWide$pt_pain_0m <= .5 * CorronaMetaWide$pt_pain_6m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .5 * CorronaMetaWide$pt_global_assess_6m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .5 * CorronaMetaWide$md_global_assess_6m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .5 * CorronaMetaWide$usresultsCRP_6m)) < 3 ~ "Responder", CorronaMetaWide$swollen_jts_28_6m > .5 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_6m > .5 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .5 * CorronaMetaWide$di_6m), as.numeric(CorronaMetaWide$pt_pain_0m <= .5 * CorronaMetaWide$pt_pain_6m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .5 * CorronaMetaWide$pt_global_assess_6m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .5 * CorronaMetaWide$md_global_assess_6m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .5 * CorronaMetaWide$usresultsCRP_6m)) >= 3  ~ "Non Responder"))

#####ENDPOINT_ACR20
CorronaMetaWide <- mutate(CorronaMetaWide, acr20_3m = case_when(CorronaMetaWide$swollen_jts_28_3m <= .2 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_3m <= .2 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .2 * CorronaMetaWide$di_3m), as.numeric(CorronaMetaWide$pt_pain_0m <= .2 * CorronaMetaWide$pt_pain_3m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .2 * CorronaMetaWide$pt_global_assess_3m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .2 * CorronaMetaWide$md_global_assess_3m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .2 * CorronaMetaWide$usresultsCRP_3m)) < 3 ~ "Responder", CorronaMetaWide$swollen_jts_28_3m > .2 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_3m > .2 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .2 * CorronaMetaWide$di_3m), as.numeric(CorronaMetaWide$pt_pain_0m <= .2 * CorronaMetaWide$pt_pain_3m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .2 * CorronaMetaWide$pt_global_assess_3m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .2 * CorronaMetaWide$md_global_assess_3m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .2 * CorronaMetaWide$usresultsCRP_3m)) >= 3  ~ "Non Responder"))

CorronaMetaWide <- mutate(CorronaMetaWide, acr20_6m = case_when(CorronaMetaWide$swollen_jts_28_6m <= .2 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_6m <= .2 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .2 * CorronaMetaWide$di_6m), as.numeric(CorronaMetaWide$pt_pain_0m <= .2 * CorronaMetaWide$pt_pain_6m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .2 * CorronaMetaWide$pt_global_assess_6m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .5 * CorronaMetaWide$md_global_assess_6m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .2 * CorronaMetaWide$usresultsCRP_6m)) < 3 ~ "Responder", CorronaMetaWide$swollen_jts_28_6m > .2 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_6m > .2 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .2 * CorronaMetaWide$di_6m), as.numeric(CorronaMetaWide$pt_pain_0m <= .2 * CorronaMetaWide$pt_pain_6m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .2 * CorronaMetaWide$pt_global_assess_6m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .2 * CorronaMetaWide$md_global_assess_6m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .2 * CorronaMetaWide$usresultsCRP_6m)) >= 3  ~ "Non Responder"))

#####ENDPOINT_ACR70
CorronaMetaWide <- mutate(CorronaMetaWide, acr70_3m = case_when(CorronaMetaWide$swollen_jts_28_3m <= .7 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_3m <= .7 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .7 * CorronaMetaWide$di_3m), as.numeric(CorronaMetaWide$pt_pain_0m <= .7 * CorronaMetaWide$pt_pain_3m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .7 * CorronaMetaWide$pt_global_assess_3m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .7 * CorronaMetaWide$md_global_assess_3m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .7 * CorronaMetaWide$usresultsCRP_3m)) < 3 ~ "Responder", CorronaMetaWide$swollen_jts_28_3m > .7 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_3m > .7 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .7 * CorronaMetaWide$di_3m), as.numeric(CorronaMetaWide$pt_pain_0m <= .7 * CorronaMetaWide$pt_pain_3m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .7 * CorronaMetaWide$pt_global_assess_3m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .7 * CorronaMetaWide$md_global_assess_3m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .7 * CorronaMetaWide$usresultsCRP_3m)) >= 3  ~ "Non Responder"))

CorronaMetaWide <- mutate(CorronaMetaWide, acr70_6m = case_when(CorronaMetaWide$swollen_jts_28_6m <= .7 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_6m <= .7 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .7 * CorronaMetaWide$di_6m), as.numeric(CorronaMetaWide$pt_pain_0m <= .7 * CorronaMetaWide$pt_pain_6m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .7 * CorronaMetaWide$pt_global_assess_6m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .5 * CorronaMetaWide$md_global_assess_6m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .7 * CorronaMetaWide$usresultsCRP_6m)) < 3 ~ "Responder", CorronaMetaWide$swollen_jts_28_6m > .7 * CorronaMetaWide$swollen_jts_28_0m | CorronaMetaWide$tender_jts_28_6m > .7 * CorronaMetaWide$tender_jts_28_0m & sum(as.numeric(CorronaMetaWide$di_0m <= .7 * CorronaMetaWide$di_6m), as.numeric(CorronaMetaWide$pt_pain_0m <= .7 * CorronaMetaWide$pt_pain_6m), as.numeric(CorronaMetaWide$pt_global_assess_0m <= .7 * CorronaMetaWide$pt_global_assess_6m), as.numeric(CorronaMetaWide$md_global_assess_0m <= .7 * CorronaMetaWide$md_global_assess_6m), as.numeric(CorronaMetaWide$usresultsCRP_0m <= .7 * CorronaMetaWide$usresultsCRP_6m)) >= 3  ~ "Non Responder"))

#####MCID
CorronaMetaWide <- mutate(
  CorronaMetaWide, mcid_3m = case_when((cdai_0m <= 22 & cdai_0m - cdai_3m >= 6)
                                       | (cdai_0m > 22 & cdai_0m - cdai_3m >= 12)
 ~ "Responder",  
 (cdai_0m <= 22 & cdai_0m - cdai_3m < 6) | 
   (cdai_0m > 22 & cdai_0m - cdai_3m < 12) ~ "Non Responder"))

CorronaMetaWide <- mutate(
  CorronaMetaWide, mcid_6m = case_when((cdai_0m <= 22 & cdai_0m - cdai_6m >= 6)
                      | (cdai_0m > 22 & cdai_0m - cdai_6m >= 12) ~ "Responder",
                      (cdai_0m <= 22 & cdai_0m - cdai_6m < 6) | 
                    (cdai_0m > 22 & cdai_0m - cdai_6m < 12) ~ "Non Responder"))

###~~~End of Response Definitions calculation
################################################################################

#reading Q2 STAR/RSEM quant norm gene expression matrix
rsemDfnorm <- read.table("Q2/rsem_gnorm_matrix.txt", header = T, sep = "\t")
row.names(rsemDfnorm) <- rsemDfnorm$X
#transposing and editing the sample identifier
rsemDfnormt <- as.data.frame(t(rsemDfnorm[,-1]))

row.names(rsemDfnormt) <- gsub("_L1\\.\\w+|H\\d*\\.", "", 
                               row.names(rsemDfnormt))
#reading Corrona Sample Manifest
path = "/Users/manuel/Dropbox (DZZOM)/Scipher Team Folder/_DEVELOPMENT (DX)/Corrona/Observational Samples/"
sampManif <- "N0140014-CERTAIN-SDM.xlsx"
CorronaSampleManifDat <- read.xlsx(paste(path, sampManif, sep = ""), 
                                   sheet = "Sheet1", 
                                   startRow = 4, rows = c(1:203))
row.names(CorronaSampleManifDat) <- CorronaSampleManifDat$SUBJECT
row.names(CorronaMetaWide) <- CorronaMetaWide$subject_id

CorronaMetaWide <- merge(CorronaMetaWide, 
                         CorronaSampleManifDat[,c("SUBJECT", "SAMPLE.ID")], 
                         by = "row.names")

row.names(CorronaMetaWide) <- CorronaMetaWide$SAMPLE.ID
RAPrismDf_133 <- merge(CorronaMetaWide, rsemDfnormt, by = "row.names")
#removing out the Row.names columns
RAPrismDf_133 <- RAPrismDf_133[,!(colnames(RAPrismDf_133) 
                                  %in% c("Row.names","Row.names.1", "SUBJECT"))]

write.table(RAPrismDf_133, "RAPrism133PatientsMasterFileNorm.txt", sep = "\t", 
            quote = F, row.names = F)
write.table(CorronaMetaWide, "RAPrism185PatientsMetadataWide.txt", sep = "\t", 
            quote = F, row.names = F)


