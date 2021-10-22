library(GEOquery)
library(dplyr)
library("org.Hs.eg.db")
library(rqdatatable)
library(Boruta)
library(randomForest)
library(svm)
library(caTools)
library(e1071)
library(caret)


my_id <-"GSE72540"
gse <- getGEO(my_id)

##SUMMARY OF DATA###
length(gse)
gse<-gse[[1]]
print(gse) #31 samples, 62976 features
pData(gse)[1:2,]  ## print the sample information 
pData(gse)$data_processing[1] #Custom cdf used. 

fData(gse)[1,] ## print the gene annotation
exprs(gse)[1,] ## print the expression data

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse)) #log 2 scaled

#Populated empty cells without gene names by Agilent Ids
for (i in 1:62976){
  if(fData(gse)[i,"GENE_SYMBOL"]=="")
  {
    
    fData(gse)[i,"GENE_SYMBOL"]=as.character(fData(gse)[i,"ID"])
  }
  
}

##Change row and column names of expression dataframe##

exp_data<-exprs(gse)

rownames(exp_data)<-fData(gse)["GENE_SYMBOL"][[1]]
colnames(exp_data)<-pData(gse)[,39]

#create two class comparison dataset 
#wheal vs healthy: w_vs_h
#nonlesional vs healthy: nl_vs_h
#wheal vs nonlesional: w_vs_nl



healthy<-exp_data[, grep(pattern="^skin$", colnames(exp_data))]
wheal<-exp_data[, grep(pattern="wheal", colnames(exp_data))]
nonlesional<-exp_data[, grep(pattern="non-lesional", colnames(exp_data))]

w_vs_h<-cbind(wheal,healthy)
colnames(w_vs_h)[which(colnames(w_vs_h) == "skin")] <- "healthy"
nl_vs_h<-cbind(nonlesional,healthy)
w_vs_nl<-cbind(wheal,nonlesional)

w_vs_h_transpose<-t(w_vs_h)
colnames(w_vs_h_transpose)<-make.names(colnames(w_vs_h_transpose), unique=TRUE)
w_vs_h_group<-rep(c(0,1),times=c(10,8))
w_vs_h_boruta<-as.data.frame(cbind(w_vs_h_group,w_vs_h_transpose))
w_vs_h_boruta$w_vs_h_group<-ifelse(test=w_vs_h_boruta$w_vs_h_group==0,yes="wheal",no="healthy")
w_vs_h_boruta$w_vs_h_group<-as.factor(w_vs_h_boruta$w_vs_h_group)


nl_vs_h_transpose<-t(nl_vs_h)
colnames(nl_vs_h_transpose)<-make.names(colnames(nl_vs_h_transpose), unique=TRUE)
nl_vs_h_group<-rep(c(0,1),times=c(13,8))
nl_vs_h_boruta<-as.data.frame(cbind(nl_vs_h_group,nl_vs_h_transpose))
nl_vs_h_boruta$nl_vs_h_group<-ifelse(test=nl_vs_h_boruta$nl_vs_h_group==0,yes="non-lesional skin",no="healthy")
nl_vs_h_boruta$nl_vs_h_group<-as.factor(nl_vs_h_boruta$nl_vs_h_group)


w_vs_nl_transpose<-t(w_vs_nl)
colnames(w_vs_nl_transpose)<-make.names(colnames(w_vs_nl_transpose), unique=TRUE)
w_vs_nl_group<-rep(c(0,1),times=c(10,13))
w_vs_nl_boruta<-as.data.frame(cbind(w_vs_nl_group,w_vs_nl_transpose))
w_vs_nl_boruta$w_vs_nl_group<-ifelse(test=w_vs_nl_boruta$w_vs_nl_group==0,yes="wheal",no="non-lesional skin")
w_vs_nl_boruta$w_vs_nl_group<-as.factor(w_vs_nl_boruta$w_vs_nl_group)

##BORUTA FEATURE SELECTION##

#feature selection
set.seed(111)
options(expressions = 5e5)
#system("R--max-ppsize=500000")


boruta_w_vs_h<-Boruta(w_vs_h_group ~ .,data=w_vs_h_boruta,doTrace=2,maxRuns=2000) #change max runs according to paper

saveRDS(boruta_w_vs_h, file = "../data/boruta_w_vs_h.rds")

boruta_nl_vs_h<-Boruta(nl_vs_h_group ~ .,data=nl_vs_h_boruta,doTrace=2,maxRuns=2000) #change max runs according to paper

saveRDS(boruta_nl_vs_h, file = "../data/boruta_nl_vs_h.rds")

boruta_w_vs_nl<-Boruta(w_vs_nl_group ~ .,data=w_vs_nl_boruta,doTrace=2,maxRuns=2000) #change max runs according to paper

saveRDS(boruta_w_vs_nl, file = "../data/boruta_w_vs_nl.rds")


##Wheal vs healthy analysis##
boruta_w_vs_h<-readRDS("../data/boruta_w_vs_h.rds")
#Boruta performed 1999 iterations in 6.922073 mins.
#159 attributes confirmed important: ADAMTS9.AS2, AKR1B15, AOC4, ASB3, ATF7IP and 154 more;
#62784 attributes confirmed unimportant: A1BG, A1BG.1, A1BG.AS1, A1CF, A1CF.1 and 62779 more;
#33 tentative attributes left: ACR, ARL6IP4.1, C15orf63.8, CCT2.8, CDH6.1 and 28 more;


bor_w_vs_h<-TentativeRoughFix(boruta_w_vs_h) #categorize tentative feature into important/not important
print(bor_w_vs_h)
#Get confirmed attribtutes
getConfirmedFormula(bor_w_vs_h)
#Boruta performed 1999 iterations in 6.922073 mins.
#Tentatives roughfixed over the last 1999 iterations.
#165 attributes confirmed important: ADAMTS9.AS2, AKR1B15, AOC4, ASB3, ATF7IP and 160 more;
#62811 attributes confirmed unimportant: A1BG, A1BG.1, A1BG.AS1, A1CF, A1CF.1 and 62806 more;

# w_vs_h_group ~ CCL2.1 + TCP1 + HOXA9 + CCL2.2 + SORBS2 + ROR1.3 + 
#   MYC + R3HDM2 + JMJD6 + RIOK2 + SORBS2.1 + XLOC_l2_004771 + 
#   WDSUB1 + X6248 + PTP4A1 + LOC100652772 + XLOC_004598 + X7400 + 
#   ATF7IP + R3HDM2.4 + SIGIRR + RBM8A + LOC400456 + LOC100505601 + 
#   HDAC11 + ERO1L.1 + PFN4 + ERO1L.2 + ERO1L.3 + UBE2A + KRT16P3 + 
#   HOXD9 + RNF115 + IL1B.2 + X12599 + MRPL17 + GUCY1A3.2 + ZNF562 + 
#   UMPS + LOC100009676 + RG9MTD1.1 + XLOC_005355 + RHOBTB2 + 
#   CCL2.4 + SH3GLB1 + RND3.1 + SORBS2.2 + XLOC_006883 + MRPL14 + 
#   WDSUB1.4 + RG9MTD1.2 + XLOC_004598.1 + RIOK2.3 + RND3.2 + 
#   COX4I2 + TAGLN.1 + ASB3 + RND3.3 + X22388 + HIF3A.2 + NUDT17 + 
#   THG1L.3 + CAV1.4 + R3HDM2.7 + MFAP3.1 + RBM8A.2 + ADAMTS9.AS2 + 
#   PPP4R2 + DCUN1D5 + IL1B.3 + CMTM5 + RAPGEF3.1 + RBM8A.3 + 
#   PDGFB + COMMD7 + IFNGR1 + MRPL52 + SNRPG.1 + NFKB1.5 + TNRC6A + 
#   PECAM1.2 + ZBTB4.8 + IL6.1 + KDM4B.1 + VDAC1.1 + RG9MTD1.8 + 
#   KATNA1.6 + C15orf63.4 + MYC.4 + LBH.5 + MIR143HG.1 + AOC4 + 
#   CCL2.6 + FPR1 + FAM162B + CAV1.8 + OSR2 + TDP1 + S100A12 + 
#   RERGL + RPF2 + X43166 + PPP1R15B + C1orf115 + CROCC + LLPH.1 + 
#   RNF149.6 + CYGB + RBM8A.7 + S100A12.1 + MMRN2.1 + RND3.4 + 
#   NCOA7 + NUP88 + RND3.5 + TCF4 + UTP11L.5 + OGG1.1 + RNF38.1 + 
#   ECT2 + AKR1B15 + RIOK2.8 + X49868 + CEBPB.7 + CSRNP3.1 + 
#   RBM8A.8 + NKD2.1 + STARD9.3 + MGAT2 + VMP1 + PIAS1 + WDSUB1.8 + 
#   PARD3.2 + XLOC_008667 + ATF7IP.8 + CASP10.3 + CCT2.8 + NME4 + 
#   LRRC59 + SNORA41 + PAR.SN + LBH.9 + RNF149.8 + RBM8A.10 + 
#   WDR3.1 + ZNF503.2 + FAM98A + WBP11.1 + MYC.9 + MAML3.2 + 
#   FAM126A.2 + R3HDM2.9 + ERO1L.9 + GFOD1.8 + XLOC_004598.4 + 
#   RG9MTD1.9 + CCDC130.7 + XLOC_l2_004640.4 + DNAJC5.1 + CCL2.9 + 
#   HSPA8.2 + TRIAP1 + COX4NB.9 + RASL12 + FFAR3.1

stat_bor_w_vs_h<-attStats(bor_w_vs_h)  #gives summary stats of each attribute with norm value showing how one features faired against its shadow attribute

geneID_w_vs_h<-rownames(stat_bor_w_vs_h)
normhits_score_w_vs_h<-stat_bor_w_vs_h$normHits
decision_w_vs_h<-stat_bor_w_vs_h$decision
GSE72540_boruta_allfeatures_w_vs_h<-data.frame(geneID_w_vs_h,normhits_score_w_vs_h,decision_w_vs_h)
write.csv(GSE72540_boruta_allfeatures_w_vs_h,"../data/GSE72540_boruta_allfeatures_w_vs_h.csv")

borutaVars_w_vs_h <- getSelectedAttributes(bor_w_vs_h)
boruta.formula_w_vs_h <- formula(paste("w_vs_h_group ~ ", 
                                paste(borutaVars_w_vs_h, collapse = " + ")))


#With K-fold cross validation##

w_vs_h_traincontrol<-trainControl(method="cv",
                               number=5,
                               search="random",
                               savePredictions=T
)

w_vs_h_modelfit<-train(boruta.formula_w_vs_h ,data=na.exclude(w_vs_h_boruta),
                    method="rf",
                    trControl=w_vs_h_traincontrol,tuneLength=10,ntree=1000)

saveRDS(w_vs_h_modelfit,"../models/GSE72540_w_vs_h.rds")
w_vs_h_modelfit<-readRDS("../models/GSE72540_w_vs_h.rds")
w_vs_h_modelfit$bestTune #mtry #1    7

png(file=paste("../images/ml_images/w_vs_h_boruta_varImportance.png"))
plot(varImp(w_vs_h_modelfit,scale=F,type=2
                   
),main="Wheal vs Healthy Variable Importance: RF 5 fold CV",top=20)
dev.off()
plot(varImp(w_vs_h_modelfit,scale=F,type=2,
            
),main="Wheal vs Healthy Variable Importance: RF 5 fold CV",top=20) # Mean Decrease Gini importance

##r<-varImp(w_vs_h_modelfit,scale=F,type=2,useModel = FALSE)
w_vs_h_sub<-subset(w_vs_h_modelfit$pred,w_vs_h_modelfit$pred$mtry==w_vs_h_modelfit$bestTune$mtry)
caret::confusionMatrix(table(w_vs_h_sub$pred,w_vs_h_sub$obs))
# Confusion Matrix and Statistics
# 
# 
# healthy wheal
# healthy       8     0
# wheal         0    10
# 
# Accuracy : 1          
# 95% CI : (0.8147, 1)
# No Information Rate : 0.5556     
# P-Value [Acc > NIR] : 2.542e-05  
# 
# Kappa : 1          
# 
# Mcnemar's Test P-Value : NA         
#                                      
#             Sensitivity : 1.0000     
#             Specificity : 1.0000     
#          Pos Pred Value : 1.0000     
#          Neg Pred Value : 1.0000     
#              Prevalence : 0.4444     
#          Detection Rate : 0.4444     
#    Detection Prevalence : 0.4444     
#       Balanced Accuracy : 1.0000     
#                                      
#        'Positive' Class : healthy    
                                     


##Wheal vs non-lesional analysis##
boruta_w_vs_nl<-readRDS("../data/boruta_w_vs_nl.rds")
# Boruta performed 1999 iterations in 7.162584 mins.
# 202 attributes confirmed important: ADAM33, ADRM1, ANKAR, ANKRD36.6, APBB2 and 197 more;
# 62761 attributes confirmed unimportant: A1BG, A1BG.1, A1BG.AS1, A1CF, A1CF.1 and 62756
# more;
# 13 tentative attributes left: ABCF1, AXIN2, FLJ13197, FOSL1, GTF2E2 and 8 more;


bor_w_vs_nl<-TentativeRoughFix(boruta_w_vs_nl) #categorize tentative feature into important/not important
print(bor_w_vs_nl)
#Get confirmed attribtutes
getConfirmedFormula(bor_w_vs_nl)
# Boruta performed 1999 iterations in 7.162584 mins.
# Tentatives roughfixed over the last 1999 iterations.
# 204 attributes confirmed important: ADAM33, ADRM1, ANKAR, ANKRD36.6, APBB2 and 199 more;
# 62772 attributes confirmed unimportant: A1BG, A1BG.1, A1BG.AS1, A1CF, A1CF.1 and 62767
# more;
# w_vs_nl_group ~ CCL2 + CCL2.1 + CXCL1 + XLOC_005473 + GTPBP4 + 
#   X1808 + MTHFR + C2CD4B + TMEM185B + CCL2.2 + TSTA3 + LOC728392 + 
#   PHC3 + XLOC_l2_015760 + XLOC_011765 + C14orf93 + JMJD6 + 
#   CCL2.3 + RPL23AP32.1 + ZER1 + MTHFD2.1 + STAT3.1 + ORMDL1.3 + 
#   RBM8A + LOC645722 + XLOC_007996 + MOB3B + RPP38.1 + ADAM33 + 
#   LOC100131564 + CXCL2 + X10471 + DDX21 + TUG1.3 + LILRB3.1 + 
#   EIF4A1.1 + CNOT6L.1 + RPL23AP32.2 + PMM2 + MRPL17 + XLOC_l2_004640.1 + 
#   SNORA31 + RPL23AP32.3 + KIAA0141 + MFSD8 + TUBA1C + GPIHBP1 + 
#   CXCL2.1 + XLOC_003482 + ATP1B2 + SERTAD1 + RND3.1 + XLOC_011313 + 
#   PPWD1.2 + C12orf44 + XLOC_000683 + KIAA1267.1 + RPL23AP32.4 + 
#   C1orf74 + C22orf13 + SYNJ2BP + LCMT2.2 + LOC441204.6 + ZNF251 + 
#   TMEM185B.4 + RPL23AP32.5 + LOC100506778 + RND3.2 + MYC.2 + 
#   SLC3A2 + TRIB1 + X22388 + SERPINA3 + BHLHE40 + ICAM1.4 + 
#   NR3C1.1 + AZI1 + X24215 + BBX.1 + TMEM185B.5 + XLOC_000683.1 + 
#   C1orf74.2 + DKFZp686M1136 + C1orf74.3 + CEBPD + PTGS2.2 + 
#   SNORD77 + LOC153684 + ZNF259 + HSPA8.1 + X28969 + CCT7.4 + 
#   PTGS2.3 + MAPK1.4 + APBB2 + GATSL3 + RBM8A.3 + DUSP5 + C8orf44 + 
#   ZNF302.1 + CH25H + EPHA2 + LOC401127 + SNRPG.1 + X32960 + 
#   PTBP2 + FLJ45248 + STAT3.6 + LOC100652757 + PHF20L1.4 + CCL2.5 + 
#   ULK2 + XLOC_007839 + RPL23AP32.7 + HN1.5 + IL6.1 + MYC.3 + 
#   PRELID2 + POLK + LOC100506965.1 + TIA1 + C1orf74.4 + MYC.4 + 
#   LOC100505648.1 + RBM8A.4 + TTC28.AS1.2 + CCL2.6 + ZBED5 + 
#   FPR1 + X39823 + C17orf96 + NET1 + RBM8A.5 + XLOC_011763 + 
#   NRAS.1 + UTP14A + MICAL3.1 + DDX19B + S100A12 + RPF2 + NME3 + 
#   IL20 + WTAP.3 + E2F6.1 + RBM8A.7 + ISG20L2 + RND3.4 + VPS36 + 
#   CCL2.7 + OSM.1 + ZNF33B + RICTOR + PTGS2.7 + XLOC_008666 + 
#   RND3.6 + UTP11L.5 + C19orf68.1 + MYC.7 + EIF4A1.3 + OGG1.1 + 
#   XLOC_001558.4 + ZNF846 + PTGS2.8 + FUZ + ADRM1 + LOC100271722.5 + 
#   UBE2J2.1 + ICAM1.9 + C1orf74.6 + VWF.10 + MAK16.7 + PTGS2.9 + 
#   SOCS1 + C17orf90 + MGEA5.1 + CCL2.8 + NOB1.9 + PRMT5 + LRRC59 + 
#   ANKRD36.6 + SNORA41 + MAML3.1 + LDB2.1 + RND3.8 + GAR1 + 
#   LOC285147 + ZNF503.2 + X57915 + ZFP36 + C7orf40 + MYC.9 + 
#   RPL23AP32.9 + POLR3D + NAMPT.2 + ANKAR + PRMT5.1 + TRNT1.2 + 
#   ZNF626 + XLOC_l2_004640.4 + CCL2.9 + TXNIP + UTP11L.10 + 
#   HSPA8.2 + MYC.10


stat_bor_w_vs_nl<-attStats(bor_w_vs_nl)  #gives summary stats of each attribute with norm value showing how one features faired against its shadow attribute

geneID_w_vs_nl<-rownames(stat_bor_w_vs_nl)
normhits_score_w_vs_nl<-stat_bor_w_vs_nl$normHits
decision_w_vs_nl<-stat_bor_w_vs_nl$decision
GSE72540_boruta_allfeatures_w_vs_nl<-data.frame(geneID_w_vs_nl,normhits_score_w_vs_nl,decision_w_vs_nl)
write.csv(GSE72540_boruta_allfeatures_w_vs_nl,"../data/GSE72540_boruta_allfeatures_w_vs_nl.csv")

borutaVars_w_vs_nl <- getSelectedAttributes(bor_w_vs_nl)
boruta.formula_w_vs_nl <- formula(paste("w_vs_nl_group ~ ", 
                                       paste(borutaVars_w_vs_nl, collapse = " + ")))


#With K-fold cross validation##

w_vs_nl_traincontrol<-trainControl(method="cv",
                                  number=5,
                                  search="random",
                                  savePredictions=T
)

w_vs_nl_modelfit<-train(boruta.formula_w_vs_nl ,data=na.exclude(w_vs_nl_boruta),
                       method="rf",
                       trControl=w_vs_nl_traincontrol,tuneLength=10,ntree=1000)

saveRDS(w_vs_nl_modelfit,"../models/GSE72540_w_vs_nl.rds")
w_vs_nl_modelfit<-readRDS("../models/GSE72540_w_vs_nl.rds")
w_vs_nl_modelfit$bestTune #mtry #1   37

png(file=paste("../images/ml_images/w_vs_nl_boruta_varImportance.png"))
plot(varImp(w_vs_nl_modelfit,scale=F,type=2
            
),main="Wheal vs Healthy Variable Importance: RF 5 fold CV",top=20)
dev.off()
plot(varImp(w_vs_nl_modelfit,scale=F,type=2,
            
),main="Wheal vs Healthy Variable Importance: RF 5 fold CV",top=20) # Mean Decrease Gini importance

##r<-varImp(w_vs_nl_modelfit,scale=F,type=2,useModel = FALSE)
w_vs_nl_sub<-subset(w_vs_nl_modelfit$pred,w_vs_nl_modelfit$pred$mtry==w_vs_nl_modelfit$bestTune$mtry)
caret::confusionMatrix(table(w_vs_nl_sub$pred,w_vs_nl_sub$obs))

# Confusion Matrix and Statistics
# 
# 
# non-lesional skin wheal
# non-lesional skin                13     0
# wheal                             0    10
# 
# Accuracy : 1                
# 95% CI : (0.8518, 1)      
# No Information Rate : 0.5652           
# P-Value [Acc > NIR] : 2e-06            
# 
# Kappa : 1                
# 
# Mcnemar's Test P-Value : NA               
#                                            
#             Sensitivity : 1.0000           
#             Specificity : 1.0000           
#          Pos Pred Value : 1.0000           
#          Neg Pred Value : 1.0000           
#              Prevalence : 0.5652           
#          Detection Rate : 0.5652           
#    Detection Prevalence : 0.5652           
#       Balanced Accuracy : 1.0000           
#                                            
#        'Positive' Class : non-lesional skin


##Non lesional vs healthy analysis##
boruta_nl_vs_h<-readRDS("../data/boruta_nl_vs_h.rds")
# Boruta performed 1999 iterations in 5.836892 mins.
# 46 attributes confirmed important: ALPK1.1, ANKRD10.1, ATP6V1C2, AXIN2, BDKRB2 and 41
# more;
# 62927 attributes confirmed unimportant: A1BG, A1BG.1, A1BG.AS1, A1CF, A1CF.1 and 62922
# more;
# 3 tentative attributes left: CTSK, DNAJC30.1, TDRD6;


bor_nl_vs_h<-TentativeRoughFix(boruta_nl_vs_h) #categorize tentative feature into important/not important
print(bor_nl_vs_h)
#Get confirmed attribtutes
getConfirmedFormula(bor_nl_vs_h)
# Boruta performed 1999 iterations in 5.836892 mins.
# Tentatives roughfixed over the last 1999 iterations.
# 46 attributes confirmed important: ALPK1.1, ANKRD10.1, ATP6V1C2, AXIN2, BDKRB2 and 41
# more;
# 62930 attributes confirmed unimportant: A1BG, A1BG.1, A1BG.AS1, A1CF, A1CF.1 and 62925
# more;
# nl_vs_h_group ~ CCDC69 + ATP6V1C2 + PTPRK + SNORA61 + AXIN2 + 
#   OR5B17 + XLOC_001099 + ETAA1 + LOC731779 + X10342 + PIAS4 + 
#   EPM2AIP1 + X14216 + TUG1.4 + SPATA18 + BDKRB2 + XLOC_002356 + 
#   ZZZ3.1 + NUDT17 + DENND2D + XLOC_000698 + NRXN2 + FPGT + 
#   IFT43 + SPHK1 + CCAR1.5 + SKIV2L2.6 + THEM4.1 + RNF149.6 + 
#   TUG1.7 + LOC100507473 + XLOC_l2_012023.3 + RAD54B.9 + ANKRD10.1 + 
#   UNC93A + TNS1.1 + ING5.1 + HNRNPH1.1 + ORMDL1.8 + NMD3 + 
#   ALPK1.1 + SMC6.7 + ITPR2.1 + RNF149.8 + PLAT.8 + RNF149.9

stat_bor_nl_vs_h<-attStats(bor_nl_vs_h)  #gives summary stats of each attribute with norm value showing how one features faired against its shadow attribute

geneID_nl_vs_h<-rownames(stat_bor_nl_vs_h)
normhits_score_nl_vs_h<-stat_bor_nl_vs_h$normHits
decision_nl_vs_h<-stat_bor_nl_vs_h$decision
GSE72540_boruta_allfeatures_nl_vs_h<-data.frame(geneID_nl_vs_h,normhits_score_nl_vs_h,decision_nl_vs_h)
write.csv(GSE72540_boruta_allfeatures_nl_vs_h,"../data/GSE72540_boruta_allfeatures_nl_vs_h.csv")

borutaVars_nl_vs_h <- getSelectedAttributes(bor_nl_vs_h)
boruta.formula_nl_vs_h <- formula(paste("nl_vs_h_group ~ ", 
                                        paste(borutaVars_nl_vs_h, collapse = " + ")))


#With K-fold cross validation##

nl_vs_h_traincontrol<-trainControl(method="cv",
                                   number=5,
                                   search="random",
                                   savePredictions=T
)

nl_vs_h_modelfit<-train(boruta.formula_nl_vs_h ,data=na.exclude(nl_vs_h_boruta),
                        method="rf",
                        trControl=nl_vs_h_traincontrol,tuneLength=10,ntree=1000)

saveRDS(nl_vs_h_modelfit,"../models/GSE72540_nl_vs_h.rds")
nl_vs_h_modelfit<-readRDS("../models/GSE72540_nl_vs_h.rds")
nl_vs_h_modelfit$bestTune #mtry #1    4

png(file=paste("../images/ml_images/nl_vs_h_boruta_varImportance.png"))
plot(varImp(nl_vs_h_modelfit,scale=F,type=2
            
),main="Wheal vs Healthy Variable Importance: RF 5 fold CV",top=20)
dev.off()
plot(varImp(nl_vs_h_modelfit,scale=F,type=2,
            
),main="Wheal vs Healthy Variable Importance: RF 5 fold CV",top=20) # Mean Decrease Gini importance

##r<-varImp(nl_vs_h_modelfit,scale=F,type=2,useModel = FALSE)
nl_vs_h_sub<-subset(nl_vs_h_modelfit$pred,nl_vs_h_modelfit$pred$mtry==nl_vs_h_modelfit$bestTune$mtry)
caret::confusionMatrix(table(nl_vs_h_sub$pred,nl_vs_h_sub$obs))

# Confusion Matrix and Statistics
# 
# 
# healthy non-lesional skin
# healthy                 7                 0
# non-lesional skin       1                13
# 
# Accuracy : 0.9524          
# 95% CI : (0.7618, 0.9988)
# No Information Rate : 0.619           
# P-Value [Acc > NIR] : 0.0005888       
# 
# Kappa : 0.8966          
# 
# Mcnemar's Test P-Value : 1.0000000       
#                                           
#             Sensitivity : 0.8750          
#             Specificity : 1.0000          
#          Pos Pred Value : 1.0000          
#          Neg Pred Value : 0.9286          
#              Prevalence : 0.3810          
#          Detection Rate : 0.3333          
#    Detection Prevalence : 0.3333          
#       Balanced Accuracy : 0.9375          
#                                           
#        'Positive' Class : healthy       



Imp_genes<- data.frame(Wheal_vs_nonlesional=borutaVars_w_vs_nl)  
Wheal_vs_Healthy <- borutaVars_w_vs_h
Nonlesional_vs_Healthy<-borutaVars_nl_vs_h
#METHOD 1
Imp_genes$Wheal_vs_Healthy <- c(Wheal_vs_Healthy, rep(NA, nrow(Imp_genes)-length(Wheal_vs_Healthy)))  #keep as integer
Imp_genes$Nonlesional_vs_Healthy <- c(Nonlesional_vs_Healthy, rep(NA, nrow(Imp_genes)-length(Nonlesional_vs_Healthy))) 
write.csv(Imp_genes,file="../data/significant_genes.rds",row.names = FALSE)
