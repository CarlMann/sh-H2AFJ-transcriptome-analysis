###########################################################################
###      Install Bioconductor and required packages if necessary   ########
##########################################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("beadarray","limma", "sva", "qvalue", "illuminaHumanv4.db", "GO.db", "pheatmap", "gplots",
## "dplyr", "data.table", "calibrate")) 


###########################################################################
###      Read Illumina Bead Chip Array Data and Normalize     #############
##########################################################################

# download GSE62701_RAW.tar from the Supplementary files of GSE62701 in GEO:
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62701
# unpacking the .tar file should give 41 idat.gz files GSM1531691 to GSM1531731
# containing the raw Illumina HT12v4 Bead Chip array data. Transfer these files to
# a folder labelled "GSE62701" in your R working directory.

idatFiles <- list.files(path = file.path(getwd(),"GSE62701"), pattern = ".idat.gz", full.names=TRUE) 
sapply(idatFiles, gunzip) # decompress the .gz compressed files

# read idat files with beadarray package
library("beadarray")
library("illuminaHumanv4.db")
idatFiles <- list.files(path = file.path(getwd(),"GSE62701"), pattern = ".idat", full.names=TRUE)
idata <- readIdatFiles(idatFiles) # idata is an ExpressionSetIllumina object (beadarray package)
dim(idata) # 48107 features, 41 samples, 1 channel (green)
table(fData(idata)$Status) # display different classes of probe in the data set
idatan <- normaliseIllumina(idata, method="neqc", transform="none") # normalizes and log2 transforms the expression levels
# The neqc normalisation is described in Shi, Oshlack, Smyth. Nucl. Acid Res. 2010 (PMID: 20929874) and documented in the limma package. 
# Output from this method has control probes removed. 47224 features, 41 samples, 1 channel (green). 
table(fData(idatan)$Status) # confirm that all control probes have been removed
pData(idatan) #indicates the name of the sample arrays
# Important: ensure that the sample arrays are ordered by increasing GSM number: GSM1531691-GSM1531731 
# Change column names to those of the samples. First 3 numbers+letter correspond to the Illumina Bead Chip Array number. NT= sh-NoTarget control,
# sh2= sh2-H2AFJ, sh3=sh3-H2AFJ, Pro=Proliferating cells, Sen= Etoposide-induced senescence. The final number represents the biological replicate.
colnames(idatan) <- c("122_K_NTSen5", "122_J_sh2Pro3", "122_H_sh2Sen1","122_G_NTSen1","122_F_NTPro2","122_E_sh2Pro2","122_D_NTPro1",
                      "122_C_sh2Pro1","116_L_sh2Pro5", "116_K_NTPro4", "116_J_sh2Sen4","116_H_sh2Sen3","116_G_NTSen8","116_F_sh2Pro4",
                      "116_E_NTPro3","116_D_NTSen7","116_C_sh2Sen2","116_B_NTSen6","046_K_Rsh2Pro3","046_J_Rsh2Sen5","046_I_sh3Pro1",
                      "046_H_sh3Sen4", "046_G_sh3Sen2", "046_F_NTSen4","046_E_Rsh2Sen1","046_D_sh3Pro4", "046_C_Rsh2Pro1", "046_B_Rsh2Sen3",
                      "046_A_Rsh2Pro4", "048_L_sh3Sen3", "048_K_sh3Pro2","048_J_NTSen3", "048_I_sh3Sen1","048_H_sh3Pro3","048_G_Rsh2Sen4",
                      "048_F_Rsh2Pro5","048_E_Rsh2Sen2","048_D_sh3Sen5","048_C_Rsh2Pro2","048_B_sh3Pro5","048_A_NTSen2")
pData(idatan) # verify correct correspondances. The 3 numbers_capital letter in front of each sample name should be identical to the last 3 numbers
# and letter of the corresponding sectionNames (array names). If not, go back and order idatFiles by GSM number.                 

# filter out probe sequences that have been annotated "Bad", "No match", or no quality score available (beadarray package).
ids <- as.character(featureNames(idatan))
qual <- unlist(mget(ids, illuminaHumanv4PROBEQUALITY, ifnotfound=NA))
table(qual) # observe distribution of probe quality scores
rem <- qual == "No match" | qual == "Bad" | is.na(qual) 
filtdata <- idatan[!rem,]  # filtered data (removed probe sequences with "Bad", "No match", or no quality score available)
dim(filtdata) # 34469 features, 41 samples, 1 channel (green)

groupSamples <- c("122_D_NTPro1", "122_F_NTPro2", "116_E_NTPro3", "116_K_NTPro4", "122_G_NTSen1", "048_A_NTSen2", "048_J_NTSen3", "046_F_NTSen4",
                   "122_K_NTSen5", "116_B_NTSen6", "116_D_NTSen7", "116_G_NTSen8", "122_C_sh2Pro1", "122_E_sh2Pro2", "122_J_sh2Pro3", 
                   "116_F_sh2Pro4", "116_L_sh2Pro5", "122_H_sh2Sen1", "116_C_sh2Sen2", "116_H_sh2Sen3","116_J_sh2Sen4", "046_I_sh3Pro1", 
                   "048_K_sh3Pro2", "048_H_sh3Pro3", "046_D_sh3Pro4", "048_B_sh3Pro5", "048_I_sh3Sen1", "046_G_sh3Sen2", "048_L_sh3Sen3", 
                   "046_H_sh3Sen4", "048_D_sh3Sen5", "046_C_Rsh2Pro1", "048_C_Rsh2Pro2", "046_K_Rsh2Pro3", "046_A_Rsh2Pro4", "048_F_Rsh2Pro5",
                   "046_E_Rsh2Sen1","048_E_Rsh2Sen2", "046_B_Rsh2Sen3", "048_G_Rsh2Sen4", "046_J_Rsh2Sen5")

filtdata <- filtdata[,groupSamples] # regroup samples according to condition
sampleNames(filtdata) # verify regrouping
boxplot(filtdata) # verify good normalization


###########################################################################
###  SVA and Batch Correction   #############
##########################################################################
library("sva")
library("limma")
library("qvalue")
# visual evaluation of sample clustering
plotMDS(exprs(filtdata), cex=0.75, col=c(rep('sandybrown', 4), rep('darkgreen', 8), rep('purple', 5), rep('red', 4), rep('violet', 5),
                                          rep('blue', 5), rep('black', 5), rep('orange', 5)))
# Note batch effect by date for the shNT-Sen samples (darkgreen) that will be corrected by sva

# Make sample description Table:
Array_Name <- c("GSM1531697_9463775122_D","GSM1531695_9463775122_F","GSM1531705_9463775116_E","GSM1531700_9463775116_K","GSM1531694_9463775122_G",
                "GSM1531731_9440308048_A","GSM1531722_9440308048_J","GSM1531714_9440308046_F","GSM1531691_9463775122_K","GSM1531708_9463775116_B",
                "GSM1531706_9463775116_D","GSM1531703_9463775116_G","GSM1531698_9463775122_C","GSM1531696_9463775122_E","GSM1531692_9463775122_J",
                "GSM1531704_9463775116_F","GSM1531699_9463775116_L","GSM1531693_9463775122_H","GSM1531707_9463775116_C","GSM1531702_9463775116_H",
                "GSM1531701_9463775116_J","GSM1531711_9440308046_I","GSM1531721_9440308048_K","GSM1531724_9440308048_H","GSM1531716_9440308046_D",
                "GSM1531730_9440308048_B","GSM1531723_9440308048_I","GSM1531713_9440308046_G","GSM1531720_9440308048_L","GSM1531712_9440308046_H",
                "GSM1531728_9440308048_D","GSM1531717_9440308046_C","GSM1531729_9440308048_C","GSM1531709_9440308046_K","GSM1531719_9440308046_A",
                "GSM1531726_9440308048_F","GSM1531715_9440308046_E","GSM1531727_9440308048_E","GSM1531718_9440308046_B","GSM1531725_9440308048_G",
                "GSM1531710_9440308046_J")
Sample <- c("NTPro1","NTPro2","NTPro3","NTPro4","NTSen1","NTSen2","NTSen3","NTSen4","NTSen5","NTSen6","NTSen7","NTSen8","sh2Pro1","sh2Pro2",
            "sh2Pro3","sh2Pro4","sh2Pro5","sh2Sen1","sh2Sen2","sh2Sen3","sh2Sen4","sh3Pro1","sh3Pro2","sh3Pro3","sh3Pro4","sh3Pro5","sh3Sen1",
            "sh3Sen2","sh3Sen3","sh3Sen4","sh3Sen5","Rsh2Pro1","Rsh2Pro2","Rsh2Pro3","Rsh2Pro4","Rsh2Pro5","Rsh2Sen1","Rsh2Sen2","Rsh2Sen2",
            "Rsh2Sen4","Rsh2Sen5")
Batch <- factor(c(1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))
Covariate1 <- c("shNT_Pro","shNT_Pro","shNT_Pro","shNT_Pro","shNT_Sen","shNT_Sen","shNT_Sen","shNT_Sen","shNT_Sen","shNT_Sen","shNT_Sen","shNT_Sen",
                "sh2_Pro","sh2_Pro","sh2_Pro","sh2_Pro","sh2_Pro","sh2_Sen","sh2_Sen","sh2_Sen","sh2_Sen","sh3_Pro","sh3_Pro","sh3_Pro","sh3_Pro",
                "sh3_Pro","sh3_Sen","sh3_Sen","sh3_Sen","sh3_Sen","sh3_Sen","Rsh2_Pro","Rsh2_Pro","Rsh2_Pro","Rsh2_Pro","Rsh2_Pro","Rsh2_Sen","Rsh2_Sen",
                "Rsh2_Sen","Rsh2_Sen","Rsh2_Sen")
Covariate1 <- as.factor(Covariate1)
Covariate2 <- c("Pro","Pro","Pro","Pro","Sen","Sen","Sen","Sen","Sen","Sen","Sen","Sen","Pro","Pro","Pro","Pro","Pro","Sen","Sen","Sen","Sen","Pro","Pro",
                "Pro","Pro","Pro","Sen","Sen","Sen","Sen","Sen","Pro","Pro","Pro","Pro","Pro","Sen","Sen","Sen","Sen","Sen")
Covariate2 <- as.factor(Covariate2)
sampleTable <- data.frame(Array_Name, Sample, Batch, Covariate1, Covariate2)
sampleTable # verify the sample table

# SVA analysis
fac <- factor(sampleTable$Covariate1, levels=c("shNT_Pro","shNT_Sen","sh2_Pro", "sh2_Sen", "sh3_Pro", "sh3_Sen","Rsh2_Pro","Rsh2_Sen"))
design <- model.matrix(~fac) # model matrix assigning the 8 distinct classes of sample
colnames(design) <- levels(fac)
mod0B <- model.matrix(~as.factor(Batch),data=sampleTable) # assigning the 2 different sample batch dates to the data
fdata <- exprs(filtdata)
svobj <- sva(fdata, mod = design, mod0B) # sva object taking into account the variables of interest (8 sample classes in design)
                                        # and the 2 batch dates as adjustment variables (mod0B). 5 surrogate variables are proposed.
modSv <- cbind(design,svobj$sv) # design model matrix combined with 5 surrogate variables.
fit <- lmFit(fdata,modSv) # use limma package to fit the linear model with surrogate variables included.
mod.sv <- coefficients(fit)[,9:13] %*% t(svobj$sv) # calculate sva-estimated batch effects for each gene (Illumina probe)
expsva <- fdata - mod.sv # gene expression values after substraction of the sva estimated non-specific batch effects.
svadat <- filtdata # transfer metadata to new EspressionSetIllumina
exprs(svadat) <- expsva # transfer the sva-corrected expression data to the new EspressionSetIllumina

annot <- select(illuminaHumanv4.db, keys=rownames(fData(svadat)), columns="SYMBOL", keytype="PROBEID") # extract gene symbols for 
# each Illumina probe from the illuminaHumanv4.db. Warning message indicated some probe IDs map to multiple Entrez IDs or gene symbols.
dim(annot)# 36911 indicates that 36911-34469 = 2442 probe IDs map to several Entrez IDs or gene symbols (exact number may vary with your versions).
# We will arbitrarily remove the duplicated annotations to obtain the 1 to 1 correspondance for the majority of probes.
annodup <- duplicated(annot$PROBEID)
sum(annodup) # verify 2442 duplicated probe IDs
annot <- annot[!annodup,]
dim(annot) # verify that we know have same number of annotations as for feature data in svadat = 34469
fData(svadat)$ENTREZID <- annot$ENTREZID # add Entrez ID column to the svadat feature data
fData(svadat)$SYMBOL <- annot$SYMBOL # add gene symbol column to the svadat feature data

# Good correction of batch effects by the SVA algorithm as shown by plotMDS:
plotMDS(expsva, cex=0.75, col=c(rep('sandybrown', 4), rep('darkgreen', 8), rep('purple', 5), rep('red', 4), rep('violet', 5),
                                rep('blue', 5), rep('black', 5), rep('orange', 5)))

#########################################################################################################################################################
###  Analysis using the limma package  to identify genes that are differentially expressed in senescence in BOTH sh2+sh3 versus shNT-control  ############
#########################################################################################################################################################

sensh <- svadat[,c(5:12,18:21,27:31)]  # subset the expression data 
pData(sensh) # verify the samples
write.table(exprs(sensh),file="senshsva.txt",quote=F,sep="\t") # write txt file of sva-corrected sensh expression data for GSEA analysis
design2 <- model.matrix(~ 0+factor(c(rep(1,8), rep(2,9))))
colnames(design2) <- c("NTSen","shSen")
fitSensh <- lmFit(exprs(sensh),design2)
contrast.matrix <- makeContrasts(shSen-NTSen, levels=design2)
fit2 <- contrasts.fit(fitSensh, contrast.matrix)
fit2 <- eBayes(fit2)
senshVP <- topTable(fit2, number=Inf, genelist=fData(svadat)$SYMBOL) # table of logFC and p-values for volcano plot
senshVP$qValue <- qvalue(senshVP$P.Value)$qvalues  # add column of q-values
sum(senshVP$qValue<.05) # 284 DE genes using qvalue <.05 (FDR<.05)
senshDE <- senshVP[which(senshVP$qValue<.05),] # extract the table of 285 DE genes with qvalue <.05 (FDR<.05)
senshDE <- senshDE[order(senshDE$logFC),] # order by increasing logFC
dup <- duplicated(senshDE$ID) # are any genes duplicated (more than one Illumina probe sequence)?
sum(dup) # 36 duplicates corresponding to genes with 2 or more Illumina probe sequences for the same gene
senshDE <- senshDE[!dup,] # remove duplicate gene entries with lowest FC
dim(senshDE) # 248 genes differentially expressed in sh-H2AFJ senescent cells versus sh-NoTarget (FDR <.05)
sum(senshDE$logFC<0) # 165 genes down-regulated in senescence due to knock-down of H2AFJ (FDR <.05)
DownsenshDE <- senshDE[senshDE$logFC<0,]
sum(senshDE$logFC>0) #  83 genes up-regulated in senescence due to knock-down of H2AFJ (FDR <.05)
UpsenshDE <- senshDE[senshDE$logFC>0,]
write.table(DownsenshDE ,file="DownsenshDE.txt",quote=F,sep="\t") # write table for publication
write.table(UpsenshDE ,file="UpsenshDE.txt",quote=F,sep="\t") # write table for publication

# Making a Volcano Plot of genes differentially regulated by H2AFJ knock-down with FDR<.05 for publication
h2afj <- senshVP[2,] # values for H2AFJ
sigcol <- senshVP$qValue<.05 #index vector for assigning red color to DE genes with FDR<.05
cols <- ifelse(senshVP$ID=='H2AFJ', "blue", ifelse(sigcol, "red", "black"))
par(mar=c(5.1, 5.1, 2.1, 2.1)) # set graph margins to increase space on left and decrease space top and right relative to default
plot(senshVP$logFC, -log10(senshVP$P.Val),cex=1, pch=16, xlim=c(-3,3), xlab="log2 (Fold Change)", ylab="-log10 (p-value)",
     cex.lab=2, col=cols, cex.axis=2) # make volcano plot: log2(Fold Change) on x axis versus -log10 (p-value) on y axis
text((senshVP[2,2]+1),-log10(senshVP[2,5]), labels='H2AFJ', col="blue", cex=2) # add 'H2AFJ' label next to its data point


###############################################################################################################
###  Gene Ontology Enrichment for Genes Differentially Affected by H2AFJ knock-down in Senescence  ############
###############################################################################################################
# Use goana function in limma for GO enrichment analysis

library(GO.db)
eIDs<- select(illuminaHumanv4.db, keys=rownames(senshDE), columns="ENTREZID", keytype="PROBEID") # select Entrez IDs for the differentially expressed genes
EntrezIDs <- na.omit(eIDs$ENTREZID) # 249 Entrez IDs (one NA)
EntrezIDs <- as.vector(EntrezIDs) # goana function requires Entrez IDs as a vector
senshgoana <- goana(EntrezIDs, species = "Hs", FDR=0.01)
topGOsensh <- topGO(senshgoana, number=50) # top 50 gene ontology classes

# Recover the list of genes that are differentially expressed and in the GO categories, and append the 
# number of genes and the list of genes to topGOsensh:

GoMapping <- as.list(org.Hs.egGO2ALLEGS)
myGOs <- rownames(topGOsensh)

# create function "topGOgenes" to extract DE genes belonging to the GO term:
topGOgenes <- function(topGO) {
        myGOs <- rownames(topGO)
        myResult <- c()
        for ( i in 1:nrow(topGO) ) {
                currentIndex <- which(names(GoMapping)==myGOs[i])
                currentIDList <- GoMapping[[currentIndex]]
                currentIDList.DE <- intersect(currentIDList, EntrezIDs)
                
                currentIDNb <- length(currentIDList.DE)
                currentIDNames <- select(illuminaHumanv4.db, keys=currentIDList.DE, columns="SYMBOL", keytype="ENTREZID")
                
                myResult <- rbind( myResult, c( currentIDNb, paste(currentIDNames$SYMBOL,collapse=", ") ) )
                
        }
        colnames(myResult) <- c( "DE", "Symbol" )
        myResult <- cbind( topGO, myResult )
}


topGOsensh <- topGOgenes(topGOsensh) # apply the function to topGOsensh
topGOsensh # visualize the gene ontology classes that are most enriched in genes differentially affected by H2A.J knock-down

write.table(topGOsensh,file="topGOsensh.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication


#########################################################################################################################################################
###  Analysis using the limma package  to identify genes that are differentially expressed in senescence in sh3 versus shNT-control  ############
#########################################################################################################################################################

sensh3 <- svadat[,c(5:12,27:31)]  # subset the expression data 
pData(svadat[,c(5:12,27:31)]) # verify the samples
write.table(exprs(sensh3),file="sensh3sva.txt",quote=F,sep="\t") # write txt file of sva-corrected sensh expression data for GSEA analysis
design13 <- model.matrix(~ 0+factor(c(rep(1,8), rep(2,5))))
colnames(design13) <- c("NTSen","sh3Sen")
fitSensh3 <- lmFit(exprs(sensh3),design13)
contrast.matrix <- makeContrasts(sh3Sen-NTSen, levels=design13)
fit13 <- contrasts.fit(fitSensh3, contrast.matrix)
fit13 <- eBayes(fit13)
sensh3VP <- topTable(fit13, number=Inf, genelist=fData(svadat)$SYMBOL) # table of logFC and p-values for volcano plot
sensh3VP$qValue <- qvalue(sensh3VP$P.Value)$qvalues  # add column of q-values
sum(sensh3VP$qValue<.05) # 2736 DE genes using qvalue <.05 (FDR<.05)
sensh3DE <- sensh3VP[which(sensh3VP$qValue<.05),] # extract the table of 2736 DE genes with qvalue <.05 (FDR<.05)
sensh3DE <- sensh3DE[order(sensh3DE$logFC),] # order by increasing logFC
dup3 <- duplicated(sensh3DE$ID) # are any genes duplicated (more than one Illumina probe sequence)?
sum(dup3) # 431 duplicates corresponding to genes with 2 or more Illumina probe sequences for the same gene
sensh3DE <- sensh3DE[!dup3,] # remove duplicate gene entries with lowest FC
dim(sensh3DE) # 2305 genes differentially expressed in sh-H2AFJ senescent cells versus sh-NoTarget (FDR <.05)
sum(sensh3DE$logFC<0) # 1283 genes down-regulated in senescence due to knock-down of H2AFJ (FDR <.05)
Downsensh3DE <- sensh3DE[sensh3DE$logFC<0,]
sum(sensh3DE$logFC>0) #  1022 genes up-regulated in senescence due to knock-down of H2AFJ (FDR <.05)
Upsensh3DE <- sensh3DE[sensh3DE$logFC>0,]
write.table(Downsensh3DE ,file="Downsensh3DE.txt",quote=F,sep="\t") # write table for publication
write.table(Upsensh3DE ,file="Upsensh3DE.txt",quote=F,sep="\t") # write table for publication

# Making a Volcano Plot of genes differentially regulated by H2AFJ knock-down with FDR<.05 
which(sensh3VP$ID=="H2AFJ") # 2 probes at rows 21 and 210
sensh3VP[c(21,210),] # probe at row 21 shows much higher avg expression and will thus be used as best representative for H2AFJ
sigcol <- sensh3VP$qValue<.05 #index vector for assigning red color to DE genes with FDR<.05
cols <- ifelse(sensh3VP$ID=='H2AFJ', "blue", ifelse(sigcol, "red", "black"))
par(mar=c(5.1, 5.1, 2.1, 2.1)) # set graph margins to increase space on left and decrease space top and right relative to default
plot(sensh3VP$logFC, -log10(sensh3VP$P.Val),cex=1, pch=16, xlim=c(-3,3), xlab="log2 (Fold Change)", ylab="-log10 (p-value)",
     cex.lab=2, col=cols, cex.axis=2) # make volcano plot: log2(Fold Change) on x axis versus -log10 (p-value) on y axis
text((sensh3VP[21,2]+0.4),-log10(sensh3VP[21,5]), labels='H2AFJ', col="blue", cex=1.0) # add 'H2AFJ' label next to its data point


###############################################################################################################
###  Gene Ontology Enrichment for Genes Differentially Affected by sh3-H2AFJ knock-down in Senescence  ############
###############################################################################################################
# Use goana function in limma for GO enrichment analysis

eIDs<- select(illuminaHumanv4.db, keys=rownames(sensh3DE), columns="ENTREZID", keytype="PROBEID") # select Entrez IDs for the differentially expressed genes
EntrezIDs <- na.omit(eIDs$ENTREZID) # 2397 Entrez IDs (one NA)
EntrezIDs <- as.vector(EntrezIDs) # goana function requires Entrez IDs as a vector
sensh3goana <- goana(EntrezIDs, species = "Hs", FDR=0.01)
topGOsensh3 <- topGO(sensh3goana, number=50) # top 50 gene ontology classes

# Recover the list of genes that are differentially expressed and in the GO categories, and append the 
# number of genes and the list of genes to topGOsensh:
myGOs <- rownames(topGOsensh3)
topGOsensh3 <- topGOgenes(topGOsensh3)
topGOsensh3 # visualize the gene ontology classes that are most enriched in genes differentially affected by H2A.J knock-down
write.table(topGOsensh3,file="topGOsensh3.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication



#########################################################################################################################################################
###  Analysis using the limma package  to identify genes that are differentially expressed in senescence in sh2 versus shNT-control  ############
#########################################################################################################################################################

sensh2 <- svadat[,c(5:12,18:21)]  # subset the expression data 
pData(svadat[,c(5:12,18:21)]) # verify the samples
write.table(exprs(sensh2),file="sensh2sva.txt",quote=F,sep="\t") # write txt file of sva-corrected sensh expression data for GSEA analysis
design12 <- model.matrix(~ 0+factor(c(rep(1,8), rep(2,4))))
colnames(design12) <- c("NTSen","sh2Sen")
fitSensh2 <- lmFit(exprs(sensh2),design12)
contrast.matrix <- makeContrasts(sh2Sen-NTSen, levels=design12)
fit12 <- contrasts.fit(fitSensh2, contrast.matrix)
fit12 <- eBayes(fit12)
sensh2VP <- topTable(fit12, number=Inf, genelist=fData(svadat)$SYMBOL) # table of logFC and p-values for volcano plot
sensh2VP$qValue <- qvalue(sensh2VP$P.Value)$qvalues  # add column of q-values
sum(sensh2VP$qValue<.05) # 1407 DE genes using qvalue <.05 (FDR<.05)
sensh2DE <- sensh2VP[which(sensh2VP$qValue<.05),] # extract the table of 1407 DE genes with qvalue <.05 (FDR<.05)
sensh2DE <- sensh2DE[order(sensh2DE$logFC),] # order by increasing logFC
dup2 <- duplicated(sensh2DE$ID) # are any genes duplicated (more than one Illumina probe sequence)?
sum(dup2) # 214 duplicates corresponding to genes with 2 or more Illumina probe sequences for the same gene
sensh2DE <- sensh2DE[!dup2,] # remove duplicate gene entries with lowest FC
dim(sensh2DE) # 1193 genes differentially expressed in sh-H2AFJ senescent cells versus sh-NoTarget (FDR <.05)
sum(sensh2DE$logFC<0) # 566 genes down-regulated in senescence due to knock-down of H2AFJ (FDR <.05)
Downsensh2DE <- sensh2DE[sensh2DE$logFC<0,]
sum(sensh2DE$logFC>0) #  627 genes up-regulated in senescence due to knock-down of H2AFJ (FDR <.05)
Upsensh2DE <- sensh2DE[sensh2DE$logFC>0,]
write.table(Downsensh2DE ,file="Downsensh2DE.txt",quote=F,sep="\t") # write table for publication
write.table(Upsensh2DE ,file="Upsensh2DE.txt",quote=F,sep="\t") # write table for publication

# Making a Volcano Plot of genes differentially regulated by H2AFJ knock-down with FDR<.05 
which(sensh2VP$ID=="H2AFJ") # 2 probes at rows 1 and 18163
sensh2VP[c(1,18163),] # probe at row 1 shows much higher avg expression and will thus be used as best representative for H2AFJ
sigcol <- sensh2VP$qValue<.05 #index vector for assigning red color to DE genes with FDR<.05
cols <- ifelse(sensh2VP$ID=='H2AFJ', "blue", ifelse(sigcol, "red", "black"))
par(mar=c(5.1, 5.1, 2.1, 2.1)) # set graph margins to increase space on left and decrease space top and right relative to default
plot(sensh3VP$logFC, -log10(sensh3VP$P.Val),cex=1, pch=16, xlim=c(-3,3), xlab="log2 (Fold Change)", ylab="-log10 (p-value)",
     cex.lab=2, col=cols, cex.axis=2) # make volcano plot: log2(Fold Change) on x axis versus -log10 (p-value) on y axis
text((sensh3VP[1,2]+0.4),-log10(sensh3VP[1,5]), labels='H2AFJ', col="blue", cex=1.0) # add 'H2AFJ' label next to its data point


###############################################################################################################
###  Gene Ontology Enrichment for Genes Differentially Affected by sh3-H2AFJ knock-down in Senescence  ############
###############################################################################################################
# Use goana function in limma for GO enrichment analysis

eIDs<- select(illuminaHumanv4.db, keys=rownames(sensh2DE), columns="ENTREZID", keytype="PROBEID") # select Entrez IDs for the differentially expressed genes
EntrezIDs <- na.omit(eIDs$ENTREZID) # 1238 Entrez IDs (one NA)
EntrezIDs <- as.vector(EntrezIDs) # goana function requires Entrez IDs as a vector
sensh2goana <- goana(EntrezIDs, species = "Hs", FDR=0.01)
topGOsensh2 <- topGO(sensh2goana, number=50) # top 50 gene ontology classes

# Recover the list of genes that are differentially expressed and in the GO categories, and append the 
# number of genes and the list of genes to topGOsensh:
myGOs <- rownames(topGOsensh2)
topGOsensh2 <- topGOgenes( topGOsensh2)
topGOsensh2 # visualize the gene ontology classes that are most enriched in genes differentially affected by H2A.J knock-down
write.table(topGOsensh2,file="topGOsensh2.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication


# compare sensh2 and sensh3: use data.table and dplyr functions to merge appropriate matrix data
        
intersect(Downsensh2DE$ID, Downsensh3DE$ID) # 114 down-reg genes shared between the 2 sets
intersect(sensh2DE$ID,sensh3DE$ID) #318 genes shared between the 2 sets
intersect(rownames(topGOsensh2), rownames(topGOsensh3)) #21/50 shared GO terms

library("dplyr")
library("data.table")
library("calibrate")

# plot DownSensh2DElogFC versus DownSensh3DElogFC:
Downsensh2DElogFC <- Downsensh2DE[,1:2] # make matrix with just gene IDs and sh2-log2FC (566 genes)
colnames(Downsensh2DElogFC) <- c("ID","sh2logFC") # change column name 
Downsensh3DElogFC <- Downsensh3DE[,1:2] # 1283 genes
colnames(Downsensh3DElogFC) <- c("ID","sh3-logFC")
Downsensh2DT <- data.table(Downsensh2DElogFC, keep.rownames=TRUE) # convert to data.table to use dplyr join function
Downsensh3DT <- data.table(Downsensh3DElogFC, keep.rownames=TRUE)
shComp <- inner_join(Downsensh2DT, Downsensh3DT, by="rn") # merge data.tables using rownamnes as common factor. 
# Only 93 common Illumina Probe IDs
shComp$ID.y <- NULL # remove redundant column of gene symbols
colnames(shComp) <- c('ProbeID',"Gene","sh2logFC","sh3logFC") # rename columns
shComp <- arrange(shComp, sh3logFC) # order data.table by increasing sh3logFC
shComptext <- filter(shComp,sh3logFC<sh2logFC) 
plot(shComp[[3]],shComp[[4]], xlab="sh2logFC", ylab="sh3logFC", abline(0,1)) #plot sh2logC vs sh3logFC with abline showing equality
textxy(shComptext[[3]],shComptext[[4]],shComptext[[2]])

# plot Sensh2DElogFC versus Sensh3DElogFC:
sensh2DElogFC <- sensh2DE[,1:2] # make matrix with just gene IDs and sh2-log2FC
colnames(sensh2DElogFC) <- c("ID","sh2logFC") # change column name
sensh3DElogFC <- sensh3DE[,1:2] # make matrix with just gene IDs and sh3-log2FC
colnames(sensh3DElogFC) <- c("ID","sh3logFC") # change column name
sensh2DT <- data.table(sensh2DElogFC, keep.rownames=TRUE) # convert to data.table to use dplyr join function
sensh3DT <- data.table(sensh3DElogFC, keep.rownames=TRUE) # convert to data.table to use dplyr join function
shAllComp <- inner_join(sensh2DT, sensh3DT, by="rn") # merge data.tables using rownamnes as common factor. 236 common Illumina Probe IDs
shAllComp$ID.y <- NULL # remove redundant column of gene symbols
colnames(shAllComp) <- c('ProbeID',"Gene","sh2logFC","sh3logFC") # rename columns
shAllComp <- arrange(shAllComp, sh3logFC) # order data.table by increasing sh3logFC
plot(shAllComp[[3]],shAllComp[[4]], xlab="sh2logFC", ylab="sh3logFC", abline(0,1))
textxy(shAllComp[[3]],shAllComp[[4]],shAllComp[[2]])

# Venn Diagrams of sh2DE versus sh3DE and sh2topGO versus sh3topGO
# install.packages('VennDiagram') install if necessary
library('VennDiagram')
grid.newpage()
venn.plot1 <- draw.pairwise.venn(area1=2305,area2=1193,cross.area=318, category=c('sh3DE','sh2DE'), fill=c('yellow','blue'),
                                 cex=2, cat.cex=2, cat.pos=c(0,0), rotation.degree=180)
grid.newpage()
venn.plot2 <- draw.pairwise.venn(area1=50,area2=50,cross.area=21, category=c('sh2topGO','sh3topGO'), fill=c('blue','yellow'),
                                 cex=2, cat.cex=2, cat.pos=0)

#################################################################################################################################################
###  Differential Gene Expression Analysis for NTPro vs NTSen in order to identify genes that are up-regulated in senescence  ############
################################################################################################################################################

senPro <- svadat[,c(1:12)]
pData(svadat[,c(1:12)]) # verify the samples
design3 <- model.matrix(~ 0+factor(c(rep(1,4), rep(2,8))))
colnames(design3) <- c("NTPro","NTSen")
fitSenPro <- lmFit(exprs(senPro),design3)
contrast.matrix2 <- makeContrasts(NTSen-NTPro, levels=design3)
fit3 <- contrasts.fit(fitSenPro, contrast.matrix2)
fit4 <- eBayes(fit3)
senProDE <- topTable(fit4, p.value=.01, number=Inf, genelist=fData(svadat)$SYMBOL) #find DE genes with adjusted p-values<.01 (increase stringency because many genes affected)
dim(senProDE) #  19999 DE probes with adjusted p-values<.01 ! (Benjamini-Hochberg FDR correction). 
sum(duplicated(senProDE$ID)) # 6557 duplicated IDs. Thus, 19999-6557= 13442 unique genes differentially expressed
uniSenProDE <- senProDE[!duplicated(senProDE$ID),] # 13442 unique genes differentially expressed
sum(uniSenProDE$logFC<0) # 6887 genes down-regulated in senescence (B-H corrected p value < .01)
sum(uniSenProDE$logFC>0) #  6555 genes up-regulated in senescence (B-H corrected p value < .01)

sortSenPro <- senProDE[order(-senProDE$logFC),] # order by decreasing logFC
write.table(sortSenPro,file="senPro.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication
UpSenPro <- sortSenPro[which(sortSenPro$logFC>2),] # 449 Illumina probes up-regulated > 4-fold in senescence
sum(duplicated(UpSenPro$ID)) # 76 duplicated genes (more than 1 Illumina probe/gene). 
UpSenPro <- UpSenPro[!duplicated(UpSenPro$ID),] # remove duplicated gene probes with lowest FC. 373 genes up-regulated >4-fold in senescence.
UpSenPro2 <- sortSenPro[which(sortSenPro$logFC>1),] # 2268 Illumina probes up-regulated > 2-fold in senescence
sum(duplicated(UpSenPro2$ID)) # 447 duplicated genes (more than 1 Illumina probe/gene).
UpSenPro2 <- UpSenPro2[!duplicated(UpSenPro2$ID),] # remove duplicated gene probes with lowest FC. 1821 genes up-regulated >2-fold in senescence.
write.table(UpSenPro2,file="UpSenPro2.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication
DownSenPro <- sortSenPro[which(sortSenPro$logFC < -2.0),] #  327 Illumina probes down-regulated > 4-fold in senescence
DownSenPro <- DownSenPro[order(DownSenPro$logFC),] # order by decreasing negative logFC
sum(duplicated(DownSenPro$ID)) # 52 duplicated genes (more than 1 Illumina probe/gene)
DownSenPro <- DownSenPro[!duplicated(DownSenPro$ID),] # remove duplicated gene probes with lowest negative FC. 275 genes down-regulated >4-fold in senescence.
DownSenPro2 <- sortSenPro[which(sortSenPro$logFC < -1.0),] #  2030 Illumina probes down-regulated > 2-fold in senescence
DownSenPro2 <- DownSenPro2[order(DownSenPro2$logFC),] # order by decreasing negative logFC
sum(duplicated(DownSenPro2$ID)) # 380 duplicated genes (more than 1 Illumina probe/gene)
DownSenPro2 <- DownSenPro2[!duplicated(DownSenPro2$ID),]  # 1650 Illumina probes down-regulated > 2-fold in senescence after removing duplicate probes.1,821 1,821 
write.table(DownSenPro2,file="DownSenPro2.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication

# Define set of SASP genes proposed by the Campisi lab based mainly on antibody arrary data (Freund et al. 2010. Trends Mol. Med. PMID 20444648):

SASP <- c("AREG", "CCL2", "CCL3", "CCL5", "CCL8", "CCL26", "CCR3", "CSF2", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CCL20",
          "HGF", "ICAM1", "IGFBP2", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL1A", "IL1B", "IL6", "IL6ST", "IL7", "IL13", "MIF",  
          "MMP1", "MMP3", "MMP10", "MMP12", "MMP13", "MMP14", "OSM", "TIMP2", "TNFRSF10C", "TNFRSF11B", "TNFRSF18")

mySASP <- intersect(UpSenPro2$ID, SASP)
mySASP # 20 SASP genes up-regulated at least 2x in etoposide-induced senescence

# Making a Volcano Plot of genes differentially expressed in senescence versus proliferation with FDR<.01 for publication
senProVP <- topTable(fit4, p.value=1, number=Inf, genelist=fData(svadat)$SYMBOL) # expression data for all genes for Volcano plot
senProVP <- senProVP[!duplicated(senProVP$ID),] # 19,623 unique genes after eliminating duplicate probes
mySASPids  <- senProVP$ID %in% mySASP # index vector identifying the ids of the 20 genes in mySASP
sigcol2 <- senProVP$adj.P.Val<.01 #index vector for assigning red color to DE genes with FDR<.05
cols <- ifelse(sigcol2, "red", "black")
par(mar=c(5.1, 5.1, 2.1, 2.1)) # set graph margins to increase space on left and decrease space top and right relative to default
plot(senProVP$logFC, -log10(senProVP$P.Val),cex=1, pch=16, xlim=c(-10,10), ylim=c(0,20), xlab="log2 (Fold Change)", ylab="-log10 (p-value)",
     cex.lab=2, col=cols, cex.axis=2) # make volcano plot: log2(Fold Change) on x axis versus -log10 (p-value) on y axis
points (senProVP[mySASPids,2], -log10(senProVP[mySASPids,5]),cex=1, pch=16, col="blue") # superpose blue points representing the mySASP genes


## Finding secreted, inflammatory, immune genes regulated by H2A.J
## 1) Identify genes encoding secreted gene products with GO:

x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- "GO:0005576" # GO term for extracellular gene products
SecEIDs <- mappedLkeys(x) # Entrez IDs for all genes in GO term of extracellular region

# 2) select Entrez IDs for genes whose expression is decreased by H2AJ knock-down in senescence.
detach(package:calibrate) # detach packages that interfere with use of select for AnnotationDb
detach(package:MASS) # detach packages that interfere with use of select for AnnotationDb
detach(package:dplyr) # detach packages that interfere with use of select for AnnotationDb
eIDdownshSen <- select(illuminaHumanv4.db, keys=rownames(DownsenshDE), columns=c("ENTREZID","SYMBOL"), keytype="PROBEID")
EntrezIDdownshSen <- na.omit(eIDdownshSen$ENTREZID) # 166 genes with Entrez IDs
eIDdownshSen <- as.vector(EntrezIDdownshSen)
SASPJdep <- intersect(SecEIDs, EntrezIDdownshSen) # 63 genes dependent on H2A.J associated with extracellular gene product GO term
        
eIDUpSenPro2 <- select(illuminaHumanv4.db, keys=rownames(UpSenPro2), columns=c("ENTREZID", "SYMBOL"), keytype="PROBEID") # select Entrez IDs for genes induced >2X in senescence.
EntrezIDUpSenPro2 <- na.omit(eIDUpSenPro2$ENTREZID) # 1915 Entrez IDs
EntrezIDUpSenPro2 <- as.vector(EntrezIDUpSenPro2)
mySASP0 <- intersect(SecEIDs, EntrezIDUpSenPro2)
indids <- eIDUpSenPro2$ENTREZID %in% mySASP0
myStab <- eIDUpSenPro2[indids,]
mySASPgenes <- myStab$SYMBOL  # 545 genes in  gene ontology "secreted" whose expression is increased at least 2-fold in senescence 
mySASPJDown <- intersect(DownsenshDE$ID, mySASPgenes) # 32 "SASP genes" whose expression is decreased by H2A.J knock-down
mySASPJDown <-mySASPJDown[c(1:2, 5:13, 15:18, 21, 25, 27:32)] # remove 6 genes that aren't obviously related to secretion, immune response, senescence
manualSASP <- c("TM4SF4", "OAS2", "P2RX6", "KCNK12", "SLC2A6", "PLA2G4C") # add 6 genes from manual check of intersect(DownsenshDE$ID, UpSenPro2$ID)
                # these genes are on cell surface or implicated in inflammation.Note that HGF was just below 2X induction in senescence and was left out.
mySASPJDown <-c(mySASPJDown, manualSASP) # final list of 29 "SASP genes" whose expression is decreased by H2A.J knock-down according to the Bead Chip analysis. 
                                  

# Use goana function in limma to search for GO enrichment for genes up-regulated in senescence

eIDs<- select(illuminaHumanv4.db, keys=rownames(UpSenPro), columns="ENTREZID", keytype="PROBEID") # select Entrez IDs for the differentially expressed genes
EntrezIDs <- na.omit(eIDs$ENTREZID) #393 Entrez IDs (one NA)
EntrezIDs <- as.vector(EntrezIDs) # goana function requires Entrez IDs as a vector
UpSenProgoana <- goana(EntrezIDs, species = "Hs")
topGOUpSenPro <- topGO(UpSenProgoana, number=50) # top 50 gene ontology classes

# Recover the list of genes that are differentially expressed and in the GO categories, and append the 
# number of genes and the list of genes to topGOUpSenPro:
myGOs <- rownames(topGOUpSenPro)
topGOUpSenPro <- topGOgenes(topGOUpSenPro) # apply the topGOgenes function to topGOUpSenPro
topGOUpSenPro

# visualize the gene ontology classes that are most enriched in genes up-regulated in senescent cells
write.table(topGOUpSenPro,file="topGO_UpSenPro.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication

# Use goana function in limma to search for GO enrichment for genes down-regulated in senescence

eIDs<- select(illuminaHumanv4.db, keys=rownames(DownSenPro), columns="ENTREZID", keytype="PROBEID") # select Entrez IDs for the differentially expressed genes
EntrezIDs <- na.omit(eIDs$ENTREZID) #281 Entrez IDs (one NA)
EntrezIDs <- as.vector(EntrezIDs) # goana function requires Entrez IDs as a vector
DownSenProgoana <- goana(EntrezIDs, species = "Hs")
topGODownSenPro <- topGO(DownSenProgoana, number=50) # top 50 gene ontology classes
topGODownSenPro # visualize the gene ontology classes that are most enriched in genes down-regulated in senescent cells


# Recover the list of genes that are differentially expressed and in the GO categories, and append the 
# number of genes and the list of genes to topGODownSenPro:
myGOs <- rownames(topGODownSenPro)
topGODownSenPro <- topGOgenes(topGODownSenPro) # apply the topGOgenes function to topGODownSenPro
topGODownSenPro # visualize the gene ontology classes that are most enriched in genes down-regulated in senescent cells
write.table(topGODownSenPro,file="topGO_DownSenPro.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication


#################################################################################################################################################
#        Differential Gene Expression Analysis for sh(2+3)-Sen vs Rsh2Sen in order to compare effects of sh-H2AFJ knock-down
#        in senescence with sh2-H2AFJ cells that also express an sh2-resistant H2AFJ cDNA (rescue experiment)                       
################################################################################################################################################

Rsh2Sen <- svadat[,c(18:21,27:31,37:41)]
pData(svadat[,c(18:21,27:31,37:41)]) # verify the samples
design4 <- model.matrix(~ 0+factor(c(rep(1,9), rep(2,5))))
colnames(design4) <- c("shSen","Rsh2Sen")
fitRsh2Sen <- lmFit(exprs(Rsh2Sen), design4)
contrast.matrix3 <- makeContrasts(Rsh2Sen-shSen, levels=design4)
fit5 <- contrasts.fit(fitRsh2Sen, contrast.matrix3)
fit6 <- eBayes(fit5)
Rsh2Sen <- topTable(fit6, p.value=1, number=Inf, genelist=fData(svadat)$SYMBOL) # tabulate all gene expression values for heat map below
Rsh2SenDE <- topTable(fit6, p.value=.01, number=Inf, genelist=fData(svadat)$SYMBOL) #find DE genes with adjusted p-values<.01
dim(Rsh2SenDE) # 1372 DE probes with adjusted p-values<.01 (Benjamini & Hochberg correction: FDR=.01)
sum(duplicated(Rsh2SenDE$ID)) # 215 duplicated IDs. Thus, 2436-386= 2050 unique genes differentially expressed
uniRsh2SenDE <- Rsh2SenDE[!duplicated(Rsh2SenDE$ID),] # 1157 unique genes differentially expressed
UpRsh2SenDE <- uniRsh2SenDE[which(uniRsh2SenDE$logFC>1),] # 120 genes up-regulated > 2-fold in senescence by ectopic expression of H2A.J
UpRsh2SenDE <- UpRsh2SenDE[order(-UpRsh2SenDE$logFC),] # order genes by decreasing logFC
write.table(UpRsh2SenDE, file="UpRsh2SenDE.txt",quote=F,sep="\t") # write table for publication
DownRsh2SenDE <- Rsh2SenDE[which(uniRsh2SenDE$logFC < -1.0),] #  61 Illumina probes down-regulated > 2-fold in senescence
DownRsh2SenDE <- DownRsh2SenDE[order(DownRsh2SenDE$logFC),] # order by decreasing negative logFC
DownRsh2SenDE <- DownRsh2SenDE[!duplicated(DownRsh2SenDE$ID),]  # 59 genes down-regulated > 2-fold in senescence after removing duplicate probes.
write.table(DownRsh2SenDE,file="DownRsh2SenDE.txt",quote=F,sep="\t") # write table for publication

# Use goana function in limma to search for GO enrichment for genes up-regulated in senescent cells that ectopically over-express H2A.J
eIDs<- select(illuminaHumanv4.db, keys=rownames(UpRsh2SenDE), columns="ENTREZID", keytype="PROBEID") # select Entrez IDs for the differentially expressed genes
EntrezIDs <- na.omit(eIDs$ENTREZID) #248 Entrez IDs (one NA)
EntrezIDs <- as.vector(EntrezIDs) # goana function requires Entrez IDs as a vector
UpRsh2SenDEgoana <- goana(EntrezIDs, species = "Hs")
topGOUpRsh2SenDE <- topGO(UpRsh2SenDEgoana, number=50) # top 50 gene ontology classes
topGOUpRsh2SenDE # visualize the gene ontology classes that are most enriched in genes up-regulated in senescent cells

# Recover the list of genes that are differentially expressed and in the GO categories, and append the 
# number of genes and the list of genes to topGOUpRsh2SenDE:
myGOs <- rownames(topGOUpRsh2SenDE)
topGOUpRsh2SenDE <- topGOgenes(topGOUpRsh2SenDE)
write.table(topGOUpRsh2SenDE,file="topGOUpRsh2SenDE.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication


#################################################################################################################################################
###     Differential Gene Expression Analysis for shNTPro vs sh(2+3)Pro in order to identify genes that are differentially expressed
###     in proliferating cells upon knock-down of H2AFJ               
################################################################################################################################################
proshcomb <- svadat[,c(1:4,13:17,22:26)]
pData(svadat[,c(1:4,13:17,22:26)]) # verify the samples
design5 <- model.matrix(~ 0+factor(c(rep(1,4), rep(2,10))))
colnames(design5) <- c("NTPro","shcombPro")
fitProsh <- lmFit(exprs(proshcomb),design5)
contrast.matrix4 <- makeContrasts(shcombPro-NTPro, levels=design5)
fit7 <- contrasts.fit(fitProsh, contrast.matrix4)
fit8 <- eBayes(fit7)
proshcombDE <- topTable(fit8, p.value=.05, number=Inf, genelist=fData(svadat)$SYMBOL) #find DE genes with adjusted p-values<.05
dim(proshcombDE) # Only 34 DE probes with adjusted p-values<.05 (Benjamini & Hochberg correction: FDR<.05)!
proshcombDE <- proshcombDE[order(proshcombDE$logFC),] # order by increasing logFC. H2AFJ down 2**2 as expected. KBTBD8 down 2**3.6. 
sum(duplicated(proshcombDE$ID)) # 2 genes represented by 2 Illumina probes
proshcombDE <- proshcombDE[!duplicated(proshcombDE$ID),] # 32 unique genes differentially expressed in proliferation (one probe without annotation)
sum(proshcombDE$logFC<0) # only 12 genes down-regulated in proliferation after H2AFJ knock-down
sum(proshcombDE$logFC>0) # only 20 genes up-regulated in proliferation after H2AFJ knock-down
downProshJ <- proshcombDE[which(proshcombDE$logFC < 0),]  # table of 12 genes down-regulated in proliferation after H2AFJ knock-down
upProshJ <- proshcombDE[which(proshcombDE$logFC > 0),] # table of 20 genes down-regulated in proliferation after H2AFJ knock-down
write.table(downProshJ,file="downProshJ.txt",quote=F,sep="\t") # write table for publication
write.table(upProshJ,file="upProshJ.txt",quote=F,sep="\t") # write table for publication

########################################################################################################################################################
# Differential Gene Expression Analysis of proliferating cells expressing shNT or sh-H2AFJ versus sh2-H2AFJ cells expressing an sh2-resistant H2AFJ cDNA              
########################################################################################################################################################

allproRsh2 <- svadat[,c(1:4,13:17,22:26,32:36)]
write.table(exprs(allproRsh2),file="allproRsh2.txt",quote=F,sep="\t") # write table for GSEA
pData(svadat[,c(1:4,13:17,22:26,32:36)]) # verify the samples
design6 <- model.matrix(~ 0+factor(c(rep(1,14), rep(2,5))))
colnames(design6) <- c("shPro","Rsh2Pro")
fitallProRsh2 <- lmFit(exprs(allproRsh2),design6)
contrast.matrix5 <- makeContrasts(Rsh2Pro-shPro, levels=design6)
fit9 <- contrasts.fit(fitallProRsh2, contrast.matrix5)
fit10 <- eBayes(fit9)
Rsh2Pro <- topTable(fit10, p.value=1, number=Inf, genelist=fData(svadat)$SYMBOL) # tabulate all gene expression values for heat map below
Rsh2ProDE <- topTable(fit10, p.value=.01, number=Inf, genelist=fData(svadat)$SYMBOL) #find DE genes with adjusted p-values<.01
dim(Rsh2ProDE) # 3384 DE probes with adjusted p-values<.01  (Benjamini & Hochberg correction: FDR=.01)

UpRsh2Pro2DE <- Rsh2ProDE[which(Rsh2ProDE$logFC>1),] # 427 Illumina probes up-regulated > 2-fold by ectopic expression of H2A.J in proliferating cells
sum(duplicated(UpRsh2Pro2DE$ID))  # 81 duplicated probes
UpRsh2Pro2DE <- UpRsh2Pro2DE[!duplicated(UpRsh2Pro2DE$ID),] # 346 unique genes up-regulated > 2-fold by ectopic expression of H2A.J in proliferating cells
UpRsh2Pro2DE  <- UpRsh2Pro2DE[order(-UpRsh2Pro2DE$logFC),] # order by decreasing logFC
write.table(UpRsh2Pro2DE,file="UpRsh2ProDE.txt",quote=F,sep="\t") # write table for publication
DownRsh2ProDE <- Rsh2ProDE[which(Rsh2ProDE$logFC < -1),] # 110 Illumina probes down-regulated > 2-fold by ectopic expression of H2A.J in proliferating cells
sum(duplicated(DownRsh2ProDE$ID)) # 18 duplicated genes (more than 1 Illumina probe/gene)
DownRsh2ProDE <- DownRsh2ProDE[!duplicated(DownRsh2ProDE$ID),]  # 92 genes down-regulated > 2-fold in senescence after removing duplicate probes.
DownRsh2ProDE <- DownRsh2ProDE[order(DownRsh2ProDE$logFC),] # order by decreasing negative logFC
write.table(DownRsh2ProDE,file="DownRsh2ProDE.txt",quote=F,sep="\t") # write table for publication

# Use goana function in limma to search for GO enrichment for genes over-regulated in proliferating cells that ectopically over-express H2A.J
eIDs<- select(illuminaHumanv4.db, keys=rownames(UpRsh2Pro2DE), columns="ENTREZID", keytype="PROBEID") # select Entrez IDs for the differentially expressed genes
EntrezIDs <- na.omit(eIDs$ENTREZID) #248 Entrez IDs (one NA)
EntrezIDs <- as.vector(EntrezIDs) # goana function requires Entrez IDs as a vector
UpRsh2ProDEgoana <- goana(EntrezIDs, species = "Hs")
topGOUpRsh2ProDE <- topGO(UpRsh2ProDEgoana, number=50) # top 50 gene ontology classes
topGOUpRsh2ProDE # visualize the gene ontology classes that are most enriched in genes up-regulated in senescent cells

# Recover the list of genes that are differentially expressed and in the GO categories, and append the 
# number of genes and the list of genes to topGOUpRsh2ProDE:
myGOs <- rownames(topGOUpRsh2ProDE)
topGOUpRsh2ProDE <- topGOgenes(topGOUpRsh2ProDE)
write.table(topGOUpRsh2ProDE,file="topGOUpRsh2ProDE.txt",quote=F,sep="\t",col.names=NA,row.names=TRUE) # write table for publication

########################################################################################################################################################
# Analysis of overlaps in up-regulated transcription by ectopic H2A.J over-expression in proliferating and senescent cells
########################################################################################################################################################

dim(UpRsh2Pro2DE)  # 346 genes
dim(UpRsh2SenDE) # 120 genes
intersect(UpRsh2Pro2DE$ID, UpRsh2SenDE$ID) # 55/120 genes in common up-regulated in proliferation and senescence by ectopic H2A.J over-expression


########################################################################################################################################################
# Analysis of overlaps in GO terms
########################################################################################################################################################
intersect(rownames(topGOUpSenPro), rownames(topGOsensh))
intersect(rownames(topGOUpRsh2SenDE), rownames(topGOsensh))
intersect(rownames(topGOUpRsh2SenDE), rownames(topGOUpSenPro))
intersect(rownames(topGOUpRsh2ProDE), rownames(topGOsensh))
intersect(rownames(topGOUpRsh2ProDE), rownames(topGOUpSenPro))

########################################################################################################################################################
# Heat maps using mySASP gene set and all genes down-regulated by H2A.J in senescence
########################################################################################################################################################

# all genes down-regulated by H2A.J in senescent cells (DownsenshDE):
DownSenJind <- which(senProDE$ID %in% na.omit(DownsenshDE$ID)) # get indices DownSenJ genes from senProDE dataframe. Use na.omit to exclude
# one entry with NA gene ID. 238 IDs
DownSenJ <- senProDE[DownSenJind,] # get senProDE data for the DownSenJ genes 
sortDownSenJ <- DownSenJ[order(-DownSenJ$logFC),] # order by decreasing logFC
sDSJind <- (!duplicated(sortDownSenJ$ID)) # index for eliminating duplicates with lowest FC
sortDownSenJ <- sortDownSenJ[sDSJind,] # eliminate duplicates with lowest FC. 140 genes. Note that 24 genes that are down-regulated by H2A.J knockdown
# are not differentially regulated in senescent versus proliferating cells (no entry in senProDE)
sDownSenJ <- sortDownSenJ[which(sortDownSenJ$logFC>0),] # 104 genes that are up-regulated in senescence of the 165 genes that are down-regulated by H2A.J depletion
sDownSenJ2 <- sortDownSenJ[which(sortDownSenJ$logFC>0.99),] # 71 of these genes are up-regulated in senescence at least 2-fold

shDownJ <- senshVP[row.names(sDownSenJ2),] #from senshVP, extract probe values for DownSenJ2 genes
Rsh2SenDownJ <- Rsh2Sen[row.names(sDownSenJ2),] #from Rsh2Sen, extract probe values for DownSenJ2 genes
Rsh2ProDownJ <- Rsh2Pro[row.names(sDownSenJ2),] #from Rsh2Pro, extract probe values for DownSenJ2 genes
myDownJHMM <- cbind(sDownSenJ2$logFC,shDownJ$logFC,Rsh2SenDownJ$logFC,Rsh2ProDownJ$logFC) # matrix with logFC for mySASP, corresponding sh(2+3)H2AFJ in senescence
# and with Rescue by expression of an sh-resistant H2AFJ cDNA
colnames(myDownJHMM) <- c('SenvsPro', 'shH2AFJ-Sen', 'Rescue-sh2-Sen', 'Rescue-sh2-Pro')
row.names(myDownJHMM) <- sDownSenJ2$ID # recover gene names for the row names

# Heatmap2 with symmetry feature: blue for values <0 and red for values >0:
library("gplots")
heatmap.2(myDownJHMM, dendrogram="none", Rowv=FALSE, col = bluered(256), Colv=FALSE, scale="none", key=TRUE, density.info="none",
          trace="none", cexRow=1, cexCol=1, symm=F,symkey=T,symbreaks=T, labRow=rownames(myDownJHMM), labCol=colnames(myDownJHMM),
          lhei = c(1,8), lwid = c(0.5,4))

# sessionInfo()

# R version 3.2.1 (2015-06-18)
# Platform: x86_64-apple-darwin10.8.0 (64-bit)
# Running under: OS X 10.6.8 (Snow Leopard)
# 
# locale:
#         [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8
# 
# attached base packages:
#         [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#         [1] gplots_3.0.1              pheatmap_1.0.8            VennDiagram_1.6.17        futile.logger_1.4.1       data.table_1.9.6          GO.db_3.1.2               qvalue_2.0.0             
# [8] limma_3.24.15             sva_3.14.0                genefilter_1.50.0         mgcv_1.8-12               nlme_3.1-128              illuminaHumanv4.db_1.26.0 org.Hs.eg.db_3.1.2       
# [15] RSQLite_1.0.0             DBI_0.4-1                 AnnotationDbi_1.30.1      GenomeInfoDb_1.4.3        IRanges_2.2.9             S4Vectors_0.6.6           beadarray_2.18.0         
# [22] ggplot2_2.1.0             Biobase_2.28.0            BiocGenerics_0.14.0      
# 
# loaded via a namespace (and not attached):
#         [1] gtools_3.5.0         reshape2_1.4.1       splines_3.2.1        lattice_0.20-33      colorspace_1.2-6     chron_2.3-47         BeadDataPackR_1.20.0 survival_2.39-4     
# [9] XML_3.98-1.2         RColorBrewer_1.1-2   calibrate_1.7.2      lambda.r_1.1.7       plyr_1.8.3           stringr_1.0.0        munsell_0.4.3        gtable_0.2.0        
# [17] caTools_1.17.1       labeling_0.3         illuminaio_0.10.0    Rcpp_0.12.5          KernSmooth_2.23-15   xtable_1.8-2         openssl_0.9.2        scales_0.4.0        
# [25] gdata_2.17.0         base64_2.0           annotate_1.46.1      XVector_0.8.0        stringi_1.1.1        dplyr_0.4.3          GenomicRanges_1.20.8 bitops_1.0-6        
# [33] tools_3.2.1          magrittr_1.5         lazyeval_0.1.10      futile.options_1.0.0 MASS_7.3-45          Matrix_1.2-6         assertthat_0.1       R6_2.1.2    
