library(GEOquery)
library(dplyr)
library("org.Hs.eg.db")
library(rqdatatable)
library(Boruta)
library(randomForest)
library(caTools)
library(e1071)
library(caret)
library(writexl)
require(stats)
require(factoextra)
library(gplots)
library(readxl)
library(grDevices)
require(scales)
library(ggsci)

##br_rf_public<-read.csv(file="../data/significant_genes.csv")
w_vs_nl_modelfit<-readRDS("../models/GSE72540_w_vs_nl.rds")
w_vs_nl_imp<-varImp(w_vs_nl_modelfit,scale=F,type=2)
w_vs_nl_scores<-w_vs_nl_imp$importance
w_vs_nl_scores$genes<-rownames(w_vs_nl_scores)
write_xlsx(w_vs_nl_scores,"../data/w_vs_nl_scores.xlsx")
sorted_scores<-read_excel("../data/w_vs_nl_scores_sorted.xlsx")
##grepl(".",rownames(w_vs_nl_imp$importance),fixed=TRUE)


#CSU gene expression(rlog)
csu_gex<-read.csv(file="../data/CSU_2_genes.rlog.csv") #27914 rows   17 columns( gene name + 16 patients)
colnames(csu_gex)<-gsub("X","",colnames(csu_gex))


#Top 25 genes

png(file=paste("../images/ml_images/w_vs_nl_boruta_varImportance.png"))
plot(varImp(w_vs_nl_modelfit,scale=F,type=2
),main="Wheal vs non-lesional Variable Importance: RF 5 fold CV",top=35) #10 more for replicates
dev.off()

first_25<-sorted_scores[1:25,]
colnames(first_25)[2]<-"Genes"

final_geneset_exp<-merge(csu_gex,first_25,by="Genes") #"UTP11L" "X24215" gene in first 25 from GSE72540 but not in CSU dataset
final_geneset_exp$Genes<-as.character(x=final_geneset_exp$Genes)
final_geneset_exp<-final_geneset_exp[,1:17]
rownames(final_geneset_exp)<-final_geneset_exp$Genes
final_geneset_exp<-final_geneset_exp[,-1]
f_gse_t<-t(final_geneset_exp) #16 patients 23 genes
f_gse_t<-na.omit(f_gse_t)

#save data for for 23 Gex of 16 samples
write.csv(f_gse_t,file="../data/hrc_data_23genes_CSU.csv")
saveRDS(f_gse_t,file="../data/hrc_data_23genes_CSU.rds")
#standardize dataset 
f_gse_t_scaled<-scale(f_gse_t)

gse_dist<-dist(x=f_gse_t_scaled,method="euclidean")



#Hierarchical clustering
gse_hc<-hclust(gse_dist,method="ward.D")
png(file=paste("../images/ml_images/w_vs_nl_CSU_hierarchical_23genes.png"))
plot(x=gse_hc)
rect.hclust(gse_hc,h=15)
dev.off()

colors()
palette()
show_col(palette(rainbow(6)))
show_col(pal_jco(palette = c("default"))(10))
show_col(pal_jco("default", alpha = 1)(10))

png(file=paste("../images/ml_images/w_vs_nl_CSU_hierarchical_23genes_color.png"))
fviz_dend(x = gse_hc, cex = 0.8, lwd = 0.8, k = 2,
         
          k_colors = c("#0073C2FF","#CD534CFF"),
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE)
dev.off()

fviz_dend(x = gse_hc, cex = 0.8, lwd = 0.8, k = 2,
          
          k_colors = c("#0073C2FF","#CD534CFF"),
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE)

png(file=paste("../images/ml_images/w_vs_nl_CSU_hierarchical_23genes_circular.png"))
fviz_dend(x = gse_hc, cex = 0.8, lwd = 0.8, k = 2,
          
          k_colors = c("#0073C2FF","#CD534CFF"),
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE,
          type="circular"
          
          )
dev.off()
#heatmap

palette <- colorRampPalette(c("green", "black", "red"))(n = 200)
png(file=paste("../images/ml_images/w_vs_nl_CSU_heatmap_23genes.png"))
heatmap.2(f_gse_t_scaled, trace="none",
          col=palette,
          key=FALSE,
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"))
dev.off()



#To 50 genes

##plot(varImp(w_vs_nl_modelfit,scale=F,type=2
##),main="Wheal vs non-lesional Variable Importance: RF 5 fold CV",top=50)
first_50<-sorted_scores[1:50,]
colnames(first_50)[2]<-"Genes"

final_geneset_exp_50<-merge(csu_gex,first_50,by="Genes") #41 genes out of 50 founds in CSU
#"UTP11L"       "X24215"       "X39823"       "X10471"      
#"LOC100131564" "X32960"       "ZNF259"       "X1808"       
#"C22orf13"  not in CSU but in boruta/RF top 50 genes


final_geneset_exp_50$Genes<-as.character(x=final_geneset_exp_50$Genes)
final_geneset_exp_50<-final_geneset_exp_50[,1:17]
rownames(final_geneset_exp_50)<-final_geneset_exp_50$Genes
final_geneset_exp_50<-final_geneset_exp_50[,-1]
f_gse_t_50<-t(final_geneset_exp_50) #16 patients 23 genes
f_gse_t_50<-na.omit(f_gse_t_50)

#save data for for 23 Gex of 16 samples
write.csv(f_gse_t_50,file="../data/hrc_data_41genes_CSU.csv")
#standardize dataset 
f_gse_t_50_scaled<-scale(f_gse_t_50)

gse_dist_50<-dist(x=f_gse_t_50_scaled,method="euclidean")

#Hierarchical clustering
gse_hc_50<-hclust(gse_dist_50,method="ward.D")
plot(x=gse_hc_50)
rect.hclust(gse_hc_50,h=15)


#heatmap 41 genes
png(file=paste("../images/ml_images/w_vs_nl_CSU_heatmap_41genes.png"))
heatmap.2(f_gse_t_50_scaled, trace="none",
          col=palette,
          key=FALSE,
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"))
dev.off()

heatmap.2(f_gse_t_50_scaled, trace="none",
          col=palette,
          key=FALSE,
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"))





##top 20


#first_25<-sorted_scores[1:20,]
#colnames(first_25)[2]<-"Genes"

final_geneset_exp<-merge(csu_gex,first_25,by="Genes") #"UTP11L" "X24215" gene in first 25 from GSE72540 but not in CSU dataset
final_geneset_exp$Genes<-as.character(x=final_geneset_exp$Genes)
final_geneset_exp<-final_geneset_exp[,1:17]
rownames(final_geneset_exp)<-final_geneset_exp$Genes
final_geneset_exp<-final_geneset_exp[,-1]
f_gse_t<-t(final_geneset_exp) #16 patients 23 genes
f_gse_t<-na.omit(f_gse_t)

#save data for for 23 Gex of 16 samples
write.csv(f_gse_t,file="../data/hrc_data_23genes_CSU.csv")
saveRDS(f_gse_t,file="../data/hrc_data_23genes_CSU.rds")
#standardize dataset 
f_gse_t_scaled<-scale(f_gse_t)

gse_dist<-dist(x=f_gse_t_scaled,method="euclidean")



#Hierarchical clustering
gse_hc<-hclust(gse_dist,method="ward.D")
png(file=paste("../images/ml_images/w_vs_nl_CSU_hierarchical_23genes.png"))
plot(x=gse_hc)
rect.hclust(gse_hc,h=15)
dev.off()

colors()
palette()
show_col(palette(rainbow(6)))
show_col(pal_jco(palette = c("default"))(10))
show_col(pal_jco("default", alpha = 1)(10))

png(file=paste("../images/ml_images/w_vs_nl_CSU_hierarchical_23genes_color.png"))
fviz_dend(x = gse_hc, cex = 0.8, lwd = 0.8, k = 2,
          
          k_colors = c("#0073C2FF","#CD534CFF"),
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE)
dev.off()

fviz_dend(x = gse_hc, cex = 0.8, lwd = 0.8, k = 2,
          
          k_colors = c("#0073C2FF","#CD534CFF"),
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE)


