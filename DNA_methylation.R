---
title: "**DNA methylation pipeline analysis in R code**"
author:  "**Francesca Catalogne**"
output:
  prettydoc::html_pretty:
    theme: cayman
---



### Step 1

*Load raw data with minfi and create an object called RGset storing the RGChannelSet object.* 

**This block will load minfi, a R library that provides the tools needed to analyze the Ilumina Methylation arrays. In this case we will be utilizing Illumina 450K array.  Analysis and setting of the working directory is then done. Loading of raw data with minfi and creation of an object called RGset storing the RGChannelSet data that are the two color channels each sample is measured on within the single array. **

```{r message=FALSE, warning=FALSE}
rm(list=ls())
suppressMessages(library(minfi))
setwd("/Users/fcatalogne/Desktop/DRD/Report/")
baseDir <- ("./Input_data")
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
RGset
```


### Step 2  

**Creation of the dataframes Red and Green to store the red and green fluorescences respectively within the array of 450,000 CpG positions.**

```{r message=FALSE, warning=FALSE}
Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))
```


### Step 3  

**Retrieving of red and green fluorescence for the address *18744490*:**  

```{r message=FALSE, warning=FALSE}
probe_red <- Red[rownames(Red)=="18744490",]
probe_green <- Green[rownames(Green)=="18744490",]
probe_red
probe_green
```

|       Sample      |Red fluor| Green fluor | Type | 
|-------------------|---------|-------------|------|
|X5775278051_R01C02 |  3723   |    10539    |  II  |
|X5775278078_R04C01 |  824    |    1649     |  II  |
|X5775278078_R02C01 |  4109	  |    9977	    |  II  |
|X5775278078_R05C02 |  844    |    1551     |  II  |
|X5775278078_R04C02 |  525	  |    474	    |  II  |
|X5930514034_R01C02 |  3761   |    9184	    |  II  |
|X5930514035_R06C02 |  4422	  |    13333	  |  II  |
|X5930514054_R06C01 |  576    |    622      |  II  |


#### *Optional*

**Checking in the manifest file if the address corresponds to a *Type I* or a *Type II* probe: **   
**Since each CpG is associated with two measurements, either "methylated" or "un-methylated", there needs to be a distinction between the two that is recognizable, this is through Type I or Type II design. **  
**CpGs measured with Type I design use a single color, with two different probes in the same color channel providing the methylated and the un-methylated measurements. **   
**CpGs of Type II are measured using a single probe, and two different colors provide the methylated and the un-methylated measurements, therefore on this type of array, there is not a one-to-one correspondence between probes and CpG positions.**

```{r message=FALSE, warning=FALSE}
load("./Illumina450Manifest_clean.RData")
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="18744490",]
```

*Within the column titled "type", this specific address corresponds to a Type II probe, therefore no color channel is specified.*   


### Step 4

**Creation of the object MSet.raw from the preprocessRaw function that is transforming the Red/Green channel for an Illumina methylation array into an actual methylation signal, without normalization.** 

```{r, echo=FALSE, warning=FALSE}
MSet.raw <- preprocessRaw(RGset)
MSet.raw
```


### Step 5

**Performing the following quality checks on the channel:**  
*QCplot*
```{r message=FALSE, warning=FALSE}
qc <- getQC(MSet.raw)
plotQC(qc)
```

**This Quality Control plot contains the ratio of unmethylated to methylated medians, taking its binary log value, we see that the values are localized above the line, indicating high signals, therefore high quality results.**

**Strip plots are made for each control probe type specified.**
```{r message=FALSE, warning=FALSE}
controlStripPlot(RGset, controls="NEGATIVE")
```

  *For what concerns the intensity, it is ideal to be below 1000, which evidently is the case, as the data falls below 10 since it is log base two. It should be noted that the green and red labels are swapped unintentionally in this schematic.*

**Next is the calculation of the detection pValues; for each sample, we calculate the probes that have a detection p-value higher than the threshold *0.05*.**

```{r message=FALSE, warning=FALSE}
detP <- detectionP(RGset)
failed <- detP>0.05
summary(failed)
```

|Sample           |Failed  |
|:-----------------:|:------:|
|5775278051_R01C02| 237    |
|5775278078_R01C01| 580    |
|5775278078_R02C01| 264    |
|5775278078_R02C02| 596    |
|5775278078_R04C02| 593    |
|5930514034_R01C02| 91     |
|5930514035_R06C02| 115    |
|5930514054_R06C01| 457    |


### Step 6

**Calculation of the raw beta and M values with their plot of their densities of mean methylation values using the divided wild type (WT) and down syndrome (DS) sample** 
**It is recommended to subset the beta and M values matrices in order to retain DS or WT subjects and apply the function mean to the 2 subsets.** 

**Subset the DS and WT and assign  MSet to each:**
```{r message=FALSE, warning=FALSE}
MSet.raw
csv <- read.csv("./Input_data/Samplesheet_report_2021.csv")
WT <- csv[csv$Group=='WT', 'Basename']
DS <- csv[csv$Group=='DS', 'Basename']
wtSet <- MSet.raw[,colnames(MSet.raw) %in% WT]
dsSet <- MSet.raw[,colnames(MSet.raw) %in% DS]
```

**Beta and M are calculated for each DS and WT subsets:**

```{r message=FALSE, warning=FALSE}
wtBeta <- getBeta(wtSet)
wtM <- getM(wtSet)

dsBeta <- getBeta(dsSet)
dsM <- getM(dsSet)
```

**Beta and M mean is calculated for each DS and WT, stripping NA values with na.rm=T:**

```{r message=FALSE, warning=FALSE}
mean_wtBeta <- apply(wtBeta,MARGIN=1,mean,na.rm=T)
mean_dsBeta <- apply(dsBeta,MARGIN=1,mean,na.rm=T)
mean_wtM <- apply(wtM,MARGIN=1,mean,na.rm=T)
mean_dsM <- apply(dsM,MARGIN=1,mean,na.rm=T)
```

**Calculation of the density distribution on Beta values:**

```{r message=FALSE, warning=FALSE}
d_mean_wtBeta <- density(mean_wtBeta,na.rm=T)
d_mean_dsBeta <- density(mean_dsBeta, na.rm=T)
```


**Apply the same steps on M values:**

```{r message=FALSE, warning=FALSE}
d_mean_dsM <- density(mean_dsM)
d_mean_wtM <- density(mean_wtM)
```

**A plot is created using an overlay of "mfrow"**

```{r message=FALSE, warning=FALSE}
par(mfrow=c(1,2))
par("mar")
plot(d_mean_wtBeta,main="Density of Beta Values",col="orange")
lines(d_mean_dsBeta,main="Density of Beta Values",col="purple")
plot(d_mean_dsM,main="Density of M Values",col="red")
lines(d_mean_wtM,main="Density of M Values",col="blue")
```



### Step 7

**Normalize the data using the function preprocessQuantile and compare raw data and normalized data.**

**This function implements stratified quantile normalization preprocessing. The normalization procedure is applied to the Meth and Unmeth intensities separately. The distribution of type I and type II signals is forced to be the same by first quantile normalizing the type II probes across samples and then interpolating a reference distribution to which we normalize the type I probes.**

**Subset the manifest into dfI and dfII in order to obtain raw beta probe values next by getting only Type I or Type II probes:**

```{r message=FALSE, warning=FALSE}
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)
```

**Names of the probes stored in row names and beta values from the beta matrix that are also found in the manifest from each probe:**

```{r message=FALSE, warning=FALSE}
beta <- getBeta(MSet.raw)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
```


**Mean of the unprocessed data:**

```{r message=FALSE, warning=FALSE}
mean_beta_I <- apply(beta_I,1,mean)
mean_beta_II <- apply(beta_II,1,mean)
```


**Density of the mean of the unprocessed data:**
```{r message=FALSE, warning=FALSE}
d_mean_beta_I <- density(mean_beta_I,na.rm=T)
d_mean_beta_II <- density(mean_beta_II,na.rm=T)
```

**For the Raw data, components included are the:matrix of beta values, the d_mean_beta_I and d_mean_beta_II objects**

**Standard deviations:** 

```{r message=FALSE, warning=FALSE}
sd_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_beta_II <- apply(beta_II,1,sd,na.rm=T)
```

**Density of the standard deviation:**

```{r message=FALSE, warning=FALSE}
d_sd_beta_I <- density(sd_beta_I,)
d_sd_beta_II <- density(sd_beta_II)
```

**Normalization is applied using the RGSet input which is the raw intensities of the red and green channels, and normalization function 'preprocessQuantile' that will stratify the microarray by CpG region such as island, shore, sea:**

```{r message=FALSE, warning=FALSE}
preprocessQUANT_results <- preprocessQuantile(RGset) 
```

**Repeat the process for the normalized data:  **   
**Assign the normalized function to beta probes:**

```{r message=FALSE, warning=FALSE}
beta_preprocessQUANT <- getBeta(preprocessQUANT_results)

```

**Assign those found in the array and the corresponding row names in the beta probe values:**

```{r message=FALSE, warning=FALSE}
beta_preprocessQUANT_I <- beta_preprocessQUANT[rownames(beta_preprocessQUANT) %in% dfI$IlmnID,]
beta_preprocessQUANT_II <- beta_preprocessQUANT[rownames(beta_preprocessQUANT) %in% dfII$IlmnID,]
```

**Mean of the normalized values: **

```{r message=FALSE, warning=FALSE}
mean_beta_preprocessQUANT_I <- apply(beta_preprocessQUANT_I,1,mean)
mean_beta_preprocessQUANT_II <- apply(beta_preprocessQUANT_II,1,mean)
```

**Density of the mean of the normalized values:**

```{r message=FALSE, warning=FALSE}
d_mean_beta_preprocessQUANT_I <- density(mean_beta_preprocessQUANT_I,na.rm=T)
d_mean_beta_preprocessQUANT_II <- density(mean_beta_preprocessQUANT_II,na.rm=T)
```

**Standard deviation of the values:**

```{r message=FALSE, warning=FALSE}
sd_beta_preprocessQUANT_I <- apply(beta_preprocessQUANT_I,1,sd)
sd_beta_preprocessQUANT_II <- apply(beta_preprocessQUANT_II,1,sd)
```

**Density of the standard deviation:**

```{r message=FALSE, warning=FALSE}
d_sd_beta_preprocessQUANT_I <- density(sd_beta_preprocessQUANT_I,na.rm=T)
d_sd_beta_preprocessQUANT_II <- density(sd_beta_preprocessQUANT_II,na.rm=T)
```

**A plot with 6 panels in which, for both raw and normalized data, is shown along with the density plots of beta mean values according to the chemistry of the probes.**

**The density plot of beta standard deviation values according to the chemistry of the probes and the boxplot of beta values:**

```{r message=FALSE, warning=FALSE}
par(mfrow=c(2,3))
plot(d_mean_beta_I,col="blue",main="raw beta")
lines(d_mean_beta_II,col="red")

plot(d_sd_beta_I,col="blue",main="raw sd")
lines(d_sd_beta_II,col="red")

boxplot(beta)
plot(d_mean_beta_preprocessQUANT_I,col="blue",main="preprocessQuantile beta")
lines(d_mean_beta_preprocessQUANT_II,col="red")
plot(d_sd_beta_preprocessQUANT_I,col="blue",main="preprocessQuantile sd")
lines(d_sd_beta_preprocessQUANT_II,col="red")
boxplot(beta_preprocessQUANT)
```
 
 
*The main objective of this comparison is between the raw data and the normalized data. From what is graphed, data coming from the beta values and their density mean show the the raw data has a greater variation at higher densities, with a lower type II probe density range towards lower methylation levels and at higher methylation levels. However where density is lowest and where the methylation levels are centered, it seems the point of intersection from the raw data occurs at lower methylation levels, where as after normalization it is more centered. The values of type I for high density at high methylation levels in normalized data seemed to stretch wider showing greater levels of high density at higher methylation levels, whereas for the type II probe, this peak at high methlyation levels seems to narrow, thus overall creating a closer overlay between the two probes at higher methylation levels within the normalized data. At low methylation levels, type II increased in density while the peak of type I falls at 4 for normalized data, where it was 5 for raw data. In regards to the box plots, it is evident that the data falls closer along the same line with a decreased number of outliers or standard deviation thus indicating higher quality and more reliable data.*


### Step 8

**Performance of a PCA on the beta matrix generated in step 7.**

```{r message=FALSE, warning=FALSE}
pca_results <- prcomp(t(beta_preprocessQUANT),scale=T)
print(summary(pca_results))
plot(pca_results)

```

*Plotting the PCA shows that variance for each component is virtually distributed equally amongst PC 1 to 7.  The first PC shows the highest proportion of variance, while the 8th PC has a variance equal to 0. *

```{r message=FALSE, warning=FALSE}
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2)
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),pos=1)
```

```{r message=FALSE, warning=FALSE}
pheno <- read.csv("./Samplesheet_report_2021.csv",header=T, stringsAsFactors=T)
str(pheno)
pheno$Group
```

```{r message=FALSE, warning=FALSE}
levels(pheno$Group)
palette(c("blue","red"))
plot(pca_results$x[,1], pca_results$x[,2],cex=1.5,pch=20,col=pheno$Group,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.4,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=16)

```

*The plots shows the distribution of PC1 to PC2. There is an outlier for W present at an abnormally high PC for both. DS seems to be equally distributed amongst the two components, while W seems to be associated with higher PC1 values and a greater range of variance for PC1 with a peak at 500 and minimum at -500.*  



### Step 9

**Using the matrix of normalized beta values generated in step 7, it is possible to identify differentially methylated probes between group DS and group WT using the function t-test. This will be done only for the the first 2000 probes which yields only one chromosome; chr 1 due to memory overloads.**

``` {r message=FALSE, warning=FALSE}
pheno <- read.csv("./Input_data/Samplesheet_report_2021.csv",header=T, stringsAsFactors=T)
My_ttest_function <- function(x) {
  t_test <- t.test(x~ pheno$Group)
  return(t_test$p.value)
} 

first20k_beta_preprocessQUANT <- beta_preprocessQUANT[1:2000,]
pValues_ttest_first20k <- apply(first20k_beta_preprocessQUANT,1, My_ttest_function)
length(pValues_ttest_first20k)

final_ttest_first20k <- data.frame(first20k_beta_preprocessQUANT, pValues_ttest_first20k)
head(final_ttest_first20k)
summary(final_ttest_first20k)

```


### Step 10

**Application of multiple test corrections with a significance threshold of 0.05.**

```{r message=FALSE, warning=FALSE}
raw_pValues <- final_ttest_first20k$pValues_ttest
corrected_pValues_BH <- p.adjust(raw_pValues,"BH")
corrected_pValues_Bonf <- p.adjust(raw_pValues,"bonferroni")
final_ttest_first20k_corrected <- data.frame(final_ttest_first20k, corrected_pValues_BH, corrected_pValues_Bonf)
summary(final_ttest_first20k_corrected)
boxplot(final_ttest_first20k_corrected[,9:11])
dim(final_ttest_first20k_corrected[final_ttest_first20k_corrected$pValues_ttest<=0.05,])
dim(final_ttest_first20k_corrected[final_ttest_first20k_corrected$corrected_pValues_BH<=0.05,])
dim(final_ttest_first20k_corrected[final_ttest_first20k_corrected$corrected_pValues_Bonf<=0.05,])
```
**How many probes do you identify as differentially methylated considering nominal pValues** 
*162 probes remain considering the nomial pValues which are indeed lower before the correction.*  

**How many after Bonferroni and BH correction?**
*There are no signingicant probes due to the threshold of 0.5 after the correction for both corerections Bonferrini and BH and we must also take into consideration the quantity of samples we are using, which is very small at only 2000.*


### Step 11

**Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis**

```{r message=FALSE, warning=FALSE}
beta_first20k <- final_ttest_first20k_corrected[,1:8]

beta_first20k_groupDS <- beta_first20k[,pheno$Group=="DS"]
mean_beta_first20k_groupDS <- apply(beta_first20k_groupDS,1,mean)
beta_first20k_groupWT <- beta_first20k[,pheno$Group=="WT"]
mean_beta_first20k_groupWT <- apply(beta_first20k_groupWT,1,mean)
 
delta_first20k <- mean_beta_first20k_groupWT-mean_beta_first20k_groupDS
head(delta_first20k)

toVolcPlot <- data.frame(delta_first20k, -log10(final_ttest_first20k_corrected$pValues_ttest))
head(toVolcPlot)

plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5, xlim= c(-0.2,0.3), ylim = c(0,4))
abline(a=-log10(0.01),b=0,col="red")

```

*In the plot is it clear that the data is distributed fairly and has reached the threshold line in red. The difference between the two groups against the pValue is made and those probes that are on the upper area are of interest because they show the largest difference between the two groups with the highest significance due to the lowest pValue. *


**Manhattan plot:**

```{r message=FALSE, warning=FALSE}
library(qqman)
final_ttest_first20k_corrected <- data.frame(rownames(final_ttest_first20k_corrected),final_ttest_first20k_corrected)
colnames(final_ttest_first20k_corrected)[1] <- "IlmnID"
final_ttest_first20k_corrected_annotated <- merge(final_ttest_first20k_corrected, Illumina450Manifest_clean,by="IlmnID")

str(final_ttest_first20k_corrected_annotated)

input_Manhattan <- final_ttest_first20k_corrected_annotated[colnames(final_ttest_first20k_corrected_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_ttest_first20k")]
dim(input_Manhattan)
str(input_Manhattan)
levels(input_Manhattan$CHR)

input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
levels(input_Manhattan$CHR)
input_Manhattan$CHR = as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest_first20k")
```
**It is important to note that along this chromosome 1, there are less genes and therefore less probes the closer we are to the telomer and this is visible with the plot as you can see the relative number of methylated sites is highly sparse towards the left side.**

### Step 12

**Production of a heatmap of the top 100 differentially methylated probes**

```{r message=FALSE, warning=FALSE}
library(gplots)
input_heatmap=as.matrix(final_ttest_first20k[1:100,1:8])
pheno$Group
colorbar <- c("green","green","orange","orange","orange","green","green","orange")
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
```


*Single*

```{r message=FALSE, warning=FALSE}
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
```

*Average:*

```{r message=FALSE, warning=FALSE}
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
```

### *Optional:*

**As DS is caused by the trisomy of chromosome 21, we can plot the density of the methylation values of the probes mapping only on chromosome 21.**

*First the manifest is accessed and specifically chromosome 21 is pulled and assigned, then drop levels is applied to remove unused levels*

```{r message=FALSE, warning=FALSE}
chr21 <- Illumina450Manifest_clean[Illumina450Manifest_clean$CHR=="21",]
chr21 <- droplevels(chr21)
dim(chr21)
str(chr21)
```

*This is taking the probes from previous normalization procedure.*

```{r message=FALSE, warning=FALSE}
betaWT_21 <- wtBeta[rownames(wtBeta) %in% chr21$IlmnID,]
betaDS_21 <- dsBeta[rownames(dsBeta) %in% chr21$IlmnID,]
mWT_21 <- wtM[rownames(wtM) %in% chr21$IlmnID,]
mDS_21 <- dsM[rownames(dsM) %in% chr21$IlmnID,]

mean_betaWT_21 <- apply(betaWT_21,1,mean)
mean_betaDS_21 <- apply(betaDS_21,1,mean)
mean_mWT_21 <- apply(mWT_21,1,mean)
mean_mDS_21 <- apply(mDS_21,1,mean)

d_mean_betaWT_21 <- density(mean_betaWT_21)
d_mean_betaDS_21 <- density(mean_betaDS_21)
d_mean_mWT_21 <- density(mean_mWT_21)
d_mean_mDS_21 <- density(mean_mDS_21)

par(mfrow = c(1,2))
plot(d_mean_betaWT_21,col="orange", main="Beta Values chr21")
lines(d_mean_betaDS_21, col = "blue")
plot(d_mean_mWT_21, col="orange", main="M Values chr21")
lines(d_mean_mDS_21, col="blue")
```


*The WT and DS sample overlap eachother almost identicallly in this instance.The number of differentially methylated probes is 4243.*





