###downloading and selecting environmental data
setwd("...")
#install.packages("colorRamps")
library(geodata)
library(raster)
library(sp)
library(rgdal)
library(gtools)
library(colorRamps)
library(data.table)

# download / prepare current/past environmental data
past <- getData("worldclim",var="bio",res=5) ###saved in "wc5"
coords<-read.table('locations_coordinates_allkokanee.txt', header=TRUE)
onlycoords <- data.frame(x=coords$x,y=coords$y)
points <- SpatialPoints(onlycoords, proj4string = past@crs) 
past_env <- raster::extract(past,points) 
bio_past<-cbind.data.frame(coords, past_env)
bio_past[,11:21]<-bio_past[,11:21]/10

###scale past environmental variables
bio_past_factors<-bio_past[c(1:4)]
bio_past_varonly<-bio_past[,5:29]
bio_past_varonly<-bio_past_varonly[,-4] ###remove tds

## Standardization of the variables
bio_past_varonly_scale<-scale(bio_past_varonly)
## Recovering scaling coefficients
scale_env <- attr(bio_past_varonly_scale, 'scaled:scale')
center_env <- attr(bio_past_varonly_scale, 'scaled:center')
bio_past_varonly_scale <- as.data.frame(bio_past_varonly_scale)
bio_past_scale<-cbind(bio_past_factors,bio_past_varonly_scale)


###save env data
write.table(bio_past, "env_variables_past.txt", sep='\t', row.names=FALSE)
write.table(bio_past_scale, "env_variables_past_scale.txt", sep='\t', row.names=FALSE)


# download / prep past climate data
climRasts <- list.files(path="./wc5",
                        full.names=T,
                        pattern=".bil")
#get files in correct order
climRasts <- mixedsort(climRasts)

# stack, subset, and clip to species range polygon
climRasts <- stack(climRasts)
names(climRasts) <- paste("bio", 01:19, sep="")
climRasts$bio1<-climRasts$bio1/10
climRasts$bio2<-climRasts$bio2/10
climRasts$bio3<-climRasts$bio3/10
climRasts$bio5<-climRasts$bio5/10
climRasts$bio6<-climRasts$bio6/10
climRasts$bio7<-climRasts$bio7/10
climRasts$bio8<-climRasts$bio8/10
climRasts$bio9<-climRasts$bio9/10
climRasts$bio10<-climRasts$bio10/10
climRasts$bio11<-climRasts$bio11/10

#climRasts <- climRasts[[which(names(climRasts) %in% colnames(envGF))]]
extent<-c(-143,-111,47,61)
#climRasts <- raster::mask(stack(climRasts), extent)
climRasts <- crop(climRasts, extent)
climRasts <- raster::trim(climRasts)
plot(climRasts, col=rgb.tables(1000))

mask <- climRasts[[1]]>-100 


###### genotypic data #############

###preparation input filed for RDA in linux/unix
#vcftools --vcf allkokanee_nomiss30_norep_nakokanee_ld_maf05.vcf --maf 0.05 --012 --out allkokanee_nomiss30_norep_nakokanee_ld_maf05
#replace -1 with NA, already done with 
#sed 's/-1/NA/g' allkokanee_nomiss30_norep_nakokanee_ld_maf05.012 > allkokanee_nomiss30_norep_nakokanee_ld_maf05_na.012
library(data.table)
library(tidyverse)
library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
#library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
#library(qvalue)
library(robustbase)
#library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)
setwd("...")
###import genotypic data
geno.012_kok <- fread("allkokanee_nomiss30_norep_nakokanee_ld_maf05_na.012")[,-1] #load genotype matrix
geno.012_kok.pos <- fread("allkokanee_nomiss30_norep_nakokanee_ld_maf05.012.pos", header=FALSE) %>% #load SNPs info
  mutate(., locus=paste(V1,V2,sep='_')) #create a new column for SNP info name (CHR + position)
geno.012_kok.indv <- fread("allkokanee_nomiss30_norep_nakokanee_ld_maf05.012.indv", header=FALSE) #load individuals info
#evaluate % of missing
sum(is.na(geno.012_kok))/(dim(geno.012_kok)[1]*dim(geno.012_kok)[2]) #[1] 0.09359583
#impute missing with the most common geno
set.seed(1)
geno.012_kok.imp <- apply(geno.012_kok, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(geno.012_kok.imp))

geno.012_kok_df<-as.data.frame(geno.012_kok.imp)
geno.012_kok.pos_df<-as.data.frame(geno.012_kok.pos)
geno.012_kok.indv_df<-as.data.frame(geno.012_kok.indv)

#Set rownames and colnames to the geno matrix
dimnames(geno.012_kok_df) <- list(geno.012_kok.indv_df$V1, geno.012_kok.pos_df$locus)
#check the geno matrix
geno.012_kok_df[1:12,1:9]

allfreq_geno<-geno.012_kok_df

###remove all objects that are not necessary and only take a lot of space
rm(list=ls(pattern="^geno"))

##import popmap
popmap_noam<-read.table("popmap_4RDA_ecotypes.txt", header= TRUE)
popmap_noam_pop <- popmap_noam %>% select(ind, pop)



######RDA based on individual genotypes ###
###inferring population structure
### In this case from the same dataset
## Running a PCA on all markers
pca_geno <- rda(allfreq_geno, scale=T)
screeplot(pca_geno, type = "barplot", npcs=10, main="PCA Eigenvalues")
## population structure table
PCs_geno <- scores(pca_geno, choices=c(1:3), display="sites", scaling=0)
PopStruct_geno <- data.frame(Ind = row.names(allfreq_geno), PCs_geno)
colnames(PopStruct_geno) <- c("Ind", "PC1", "PC2", "PC3")
rownames(PopStruct_geno) <- NULL
plot(PopStruct_geno$PC1, PopStruct_geno$PC2)
PopStruct_geno$location<-popmap_noam$pop
PopStruct_geno$ecospawning<-popmap_noam$spawning
PopStruct_geno$ecodepth<-popmap_noam$depth

bio_past_scale<-read.table("env_variables_past_scale.txt", header = TRUE)
variables<-merge(PopStruct_geno,bio_past_scale, by = "location")
write.table(variables, "env_variables_past_scale_pca_ecotype.txt ", sep="\t", row.names= FALSE)
variables<-read.table("env_variables_past_scale_pca_ecotype.txt", header=TRUE)

### plot PCA pretty
new_bg<-c("#001219", "#5a189a", "#B785FF", "#00509d", "#4A9FFF", "#AAD8FF", "#005E29", "#00C05A", "#83D477", "#b9e769", "#e9d8a6", "#efea5a", "#f1c453", "#f29e4c", "#f3722c", "#f94144", "#d00000", "#9d0208", "#6a040f", "#990066", "#d94a8c", "#fad2e1")

pc12<-variables %>%
  mutate(location= factor(location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Bobtail", "Arctic", "Puntzi", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%>%
  ggplot( aes(PC1, PC2, color = location, fill = location)) + 
  geom_point(shape=21, alpha=0.7) + scale_color_manual(values=new_bg) + scale_fill_manual(values=new_bg) + theme_bw() 

pc13<-variables %>%
  mutate(location= factor(location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Bobtail", "Arctic", "Puntzi", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%>%
  ggplot( aes(PC1, PC3, color = location, fill = location)) + 
  geom_point(shape=21, alpha=0.7) + scale_color_manual(values=new_bg) + scale_fill_manual(values=new_bg)+ theme_bw() 
# library(plotly)
# plot_ly(x=variables$PC1, y=variables$PC2, z=variables$PC3, type="scatter3d", mode="markers", color= variables$location)
#library(patchwork)
combined_pcas<-pc12 + pc13 & theme(legend.position = "bottom")
combined_pcas_legend<- combined_pcas + plot_layout(guides="collect")
combined_pcas_legend + plot_annotation(title="Sequence variation", tag_levels = 'A')


## Null model
RDA0 <- rda(allfreq_geno ~ 1,  variables) 
## Full model
RDA_6var_nopop <- rda(allfreq_geno ~  surface_area + pH + bio5 + bio6 + bio15 + bio16, variables)

RDA_6var_pop <- rda(allfreq_geno ~  surface_area + pH + bio5 + bio6 + bio15 + bio16 + Condition(PC1 + PC2 + PC3), variables)
RDA_4var_pop <- rda(allfreq_geno ~ bio5 + bio6 + bio15 + bio16 + Condition(PC1 + PC2 + PC3), variables)

RsquareAdj(RDA_6var_nopop)
# $r.squared
# [1] 0.0957709
# 
# $adj.r.squared
# [1] 0.07076917

RsquareAdj(RDA_6var_pop)
# $r.squared
# [1] 0.07289394
# 
# $adj.r.squared
# [1] 0.04956888

RsquareAdj(RDA_4var_pop)
# $r.squared
# [1] 0.05092257
# 
# $adj.r.squared
# [1] 0.03514195


## Stepwise procedure with ordiR2step function
mod1 <- ordiR2step(RDA0, RDA_6var_nopop, Pin = 0.01, R2permutations = 1000, R2scope = T)
# Step: R2.adj= 0 
# Call: allfreq_geno ~ 1 
# 
# R2.adjusted
# <All variables> 0.070769172
# + bio05         0.020813380
# + bio06         0.015198532
# + pH            0.013717024
# + bio15         0.013423403
# + bio16         0.010858343
# + surface_area  0.008154532
# <none>          0.000000000
# 
# Df    AIC    F Pr(>F)   
# + bio05  1 2750.5 5.74  0.002 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step: R2.adj= 0.02081338 
# Call: allfreq_geno ~ bio05 
# 
# R2.adjusted
# <All variables>  0.07076917
# + bio06          0.03297964
# + bio15          0.03230131
# + pH             0.03166568
# + bio16          0.03011155
# + surface_area   0.02835078
# <none>           0.02081338
# 
# Df    AIC     F Pr(>F)   
# + bio06  1 2748.7 3.793  0.002 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step: R2.adj= 0.03297964 
# Call: allfreq_geno ~ bio05 + bio06 
# 
# R2.adjusted
# <All variables>  0.07076917
# + bio15          0.04506641
# + pH             0.04258495
# + bio16          0.04159249
# + surface_area   0.04061492
# <none>           0.03297964
# 
# Df    AIC      F Pr(>F)   
# + bio15  1 2746.9 3.7972  0.002 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step: R2.adj= 0.04506641 
# Call: allfreq_geno ~ bio05 + bio06 + bio15 
# 
# R2.adjusted
# <All variables>  0.07076917
# + pH             0.05478873
# + bio16          0.05328317
# + surface_area   0.05312728
# <none>           0.04506641
# 
# Df    AIC      F Pr(>F)   
# + pH  1 2745.6 3.2629  0.002 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step: R2.adj= 0.05478873 
# Call: allfreq_geno ~ bio05 + bio06 + bio15 + pH 
# 
# R2.adjusted
# <All variables>  0.07076917
# + bio16          0.06299901
# + surface_area   0.06258895
# <none>           0.05478873
# 
# Df    AIC      F Pr(>F)   
# + bio16  1 2744.6 2.9189  0.002 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step: R2.adj= 0.06299901 
# Call: allfreq_geno ~ bio05 + bio06 + bio15 + pH + bio16 
# 
# R2.adjusted
# <All variables>  0.07076917
# + surface_area   0.07076917
# <none>           0.06299901
# 
# Df    AIC      F Pr(>F)   
# + surface_area  1 2743.7 2.8229  0.002 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step: R2.adj= 0.07076917 
# Call: allfreq_geno ~ bio05 + bio06 + bio15 + pH + bio16 + surface_area 

screeplot(RDA_6var_nopop, main="Eigenvalues of constrained axes")
screeplot(RDA_6var_pop, main="Eigenvalues of constrained axes") #2 axes
screeplot(RDA_4var_pop, main="Eigenvalues of constrained axes") #2 axes

levels(variables$location)<-c("Anderson","Arctic","Arrow_Hill","Arrow_Mosq","Bobtail","Bonaparte",   
                              "Christina","Cluculz","Cowichan","Dunn","EastBarriere","Kalamalka",   
                              "Kootenay","La_Hache","Nicola","Okanagan","Puntzi","Shawningan",  
                              "Sockeye","Tchesinkut","Thutade", "Wood")
bg<-c("#efea5a","#2c699a","#9d0208","#6a040f","#048ba8","#83e377","#d00000","#013a63","#d94a8c","#b9e769","#e9d8a6","#f94144","#990066","#16db93","#f1c453","#f29e4c","#0db39e","#fad2e1","#001219","#815ac0","#3c096c","#f3722c")
names(bg)<-levels(variables$location)

#plot1
pdf(file="plot_RDA_6var_nopop.pdf")
plot(RDA_6var_nopop, type="n", scaling=3)
points(RDA_6var_nopop, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_6var_nopop, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables$location])
text(RDA_6var_nopop, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

#plot2
pdf(file="plot_RDA_6var_pop.pdf")
plot(RDA_6var_pop, type="n", scaling=3)
points(RDA_6var_pop, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_6var_pop, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables$location])
text(RDA_6var_pop, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

#plot3
pdf(file="plot_RDA_4var_pop.pdf")
plot(RDA_4var_pop, type="n", scaling=3)
points(RDA_4var_pop, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_4var_pop, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables$location])
text(RDA_4var_pop, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

###variance partitioning 
## full model
pRDAfull<-rda(allfreq_geno ~ ecodepth + ecospawning + x + y + PC1 + PC2 + PC3 + surface_area + pH + bio5 + bio6 + bio15 + bio16, variables)
### climate model
pRDAclim<-rda(allfreq_geno ~  surface_area + pH + bio5 + bio6 + bio15 + bio16 + Condition(ecodepth + ecospawning + x + y + PC1 + PC2 + PC3), variables)
### pop structure model
pRDAstruc<-rda(allfreq_geno ~ PC1 + PC2 + PC3 + Condition(ecodepth + ecospawning + x + y + surface_area + pH + bio5 + bio6 + bio15 + bio16), variables)
### geography model
pRDAgeog<-rda(allfreq_geno ~  x + y + Condition(ecodepth + ecospawning + PC1 +PC2 +PC3 + surface_area + pH + bio5 + bio6 + bio15 + bio16), variables)
## ecotype model
pRDAecot<-rda(allfreq_geno ~  ecodepth + ecospawning + Condition(x + y +  PC1 + PC2 + PC3 + surface_area + pH + bio5 + bio6 + bio15 + bio16), variables)


RsquareAdj(pRDAfull)
RsquareAdj(pRDAclim)
RsquareAdj(pRDAstruc)
RsquareAdj(pRDAgeog)
RsquareAdj(pRDAecot)
# RsquareAdj(pRDAfull)
# $r.squared
# [1] 0.1845951
# 
# $adj.r.squared
# [1] 0.1341176
# 
# > RsquareAdj(pRDAclim)
# $r.squared
# [1] 0.06291756
# 
# $adj.r.squared
# [1] 0.04090427
# 
# > RsquareAdj(pRDAstruc)
# $r.squared
# [1] 0.03992637
# 
# $adj.r.squared
# [1] 0.02960532
# 
# > RsquareAdj(pRDAgeog)
# $r.squared
# [1] 0.02247363
# 
# $adj.r.squared
# [1] 0.01547101
# 
# > RsquareAdj(pRDAecot)
# $r.squared
# [1] 0.01686154
# 
# $adj.r.squared
# [1] 0.009567733

# anova(pRDAfull)
# anova(pRDAclim)
# anova(pRDAstruc)
# anova(pRDAgeog)


###Identifying candidate genes based on 6 variables (including area and ph)
### method 1 - outlier based on SD
### correlation each SNP
load.rda <- scores(RDA_6var_pop, choices=c(1:3), display="species")
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}


cand1 <- outliers(load.rda[,1],3) # 
cand2 <- outliers(load.rda[,2],3) # 
cand3 <- outliers(load.rda[,3],3) # 

ncand <- length(cand1) + length(cand2) + length(cand3) 
ncand #63473
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2)  <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

### 
foo <- matrix(nrow=ncand, ncol=6)

colnames(foo) <- c( "surface_area" , "pH" , "bio5", "bio6", "bio15",	"bio16")
var_env<-variables[c(7,8,14,15,24,25)] 
#var_env_num<-as.numeric(unlist(var_env))
# for (i in 1:length(cand$snp)) {
#   nam <- cand[i,2]
#   snp.gen <- allfreq_geno[,nam]
#   foo[i,] <- apply(var_env,2,function(x) cor(x,snp.gen))
# }


for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  # print(nam)
  snp.gen <- allfreq_geno[,nam]
  #  print(snp.gen)
  foo[i,] <- apply(var_env,2,function(x) cor(x,snp.gen))
  #  print(foo[i,])
}
head(foo)

cand_corr <- cbind.data.frame(cand,foo)  
head(cand_corr)
write.table(cand_corr, "outliers_RDA_6var_pop.txt", sep = "\t",row.names = FALSE)

for (i in 1:length(cand_corr$snp)) {
  bar <- cand_corr[i,]
  cand_corr[i,10] <- names(which.max(abs(bar[4:9]))) # gives the variable
  cand_corr[i,11] <- max(abs(bar[4:9]))              # gives the correlation
}

colnames(cand_corr)[10] <- "predictor"
colnames(cand_corr)[11] <- "correlation"

cand_corr_uni <- cand_corr[!duplicated(cand_corr$snp),]
write.table(cand_corr_uni, "outliers_RDA_6var_pop_unique.txt", sep = "\t",row.names = FALSE)
table(cand_corr_uni$predictor) 
# 
# bio05        bio06        bio15        bio16           pH surface_area 
# 11547        14624         7868         3302        14197         5930 

cand_corr_strong <- cand_corr_uni[cand_corr_uni$correlation > 0.5,]
write.table(cand_corr_strong, "outliers_RDA_6var_pop_unique_strong.txt", sep = "\t",row.names = FALSE)
table(cand_corr_strong$predictor) 
# bio05        bio06        bio15        bio16           pH surface_area 
# 521          339          347           23          182           32 


sel <- cand_corr$snp
env <- cand_corr$predictor
#"surface_area" , "pH" , "bio5", "bio6", "bio15",	"bio16"
env[env=="surface_area"] <- '#1f78b4'
env[env=="pH"] <- '#a6cee3'
env[env=="bio5"] <- '#6a3d9a'
env[env=="bio6"] <- '#e31a1c'
env[env=="bio15"] <- '#33a02c'
env[env=="bio16"] <- '#ffff33'


# color by predictor:
col.pred <- rownames(RDAclim$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("NC",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg2 <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c', '#33a02c','#ffff33')

# axes 1 & 2
pdf(file="plot_RDA_6var_pop_onlyoutliers.pdf")

plot(RDA_6var_pop, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(RDA_6var_pop, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(RDA_6var_pop, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(RDA_6var_pop, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("surface_area" , "pH" , "bio5", "bio6", "bio15",	"bio16"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg2)

dev.off()


###### Adaptive landscape
## Adaptively enriched RDA
cand_corr_uni<-read.table("outliers_RDA_6var_pop_unique.txt", header = TRUE)
cand_corr_strong <- cand_corr_uni[cand_corr_uni$correlation > 0.5,]
#allfreq_geno_out<-allfreq_geno[,cand_corr_uni$snp]
RDA_outliers <- rda(allfreq_geno[,cand_corr_uni$snp] ~ surface_area  + pH + bio5 + bio6 + bio15 + bio16, variables) 
#RDA_outliers <- rda(allfreq_geno_out ~ surface_area  + pH + bio5 + bio6 + bio15 + bio16, variables) 

RDA0_out <- rda(allfreq_geno[,cand_corr_uni$snp] ~ 1,  variables) 
mod1_out <- ordiR2step(RDA0_out, RDA_outliers, Pin = 0.01, R2permutations = 1000, R2scope = T)

# RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=RDA1, y=RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  #xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))


pdf("plot_rda_outliers.pdf")
plot(RDA_outliers, type="n", scaling=3)
points(RDA_outliers, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_outliers, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables$location])
text(RDA_outliers, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

## Adaptively enriched RDA with strong candidates
RDA_outliers_strong <- rda(allfreq_geno[,cand_corr_strong$snp] ~ surface_area  + pH + bio5 + bio6 + bio15 + bio16, variables)
RDA0_out_strong <- rda(allfreq_geno[,cand_corr_strong$snp] ~ 1,  variables) 
mod1_out_strong <- ordiR2step(RDA0_out_strong, RDA_outliers_strong, Pin = 0.01, R2permutations = 1000, R2scope = T)


TAB_loci_strong <- as.data.frame(scores(RDA_outliers_strong, choices=c(1:2), display="species", scaling="none"))
TAB_var_strong <- as.data.frame(scores(RDA_outliers_strong, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_strong, aes(x=RDA1*10, y=RDA2*10), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  # geom_segment(data = TAB_var_strong, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_strong, aes(x=RDA1, y=RDA2, label = row.names(TAB_var_strong)), size = 2.5, family = "Times") +
  #xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space - strong outliers") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

pdf("plot_rda_outliers_strong.pdf")
plot(RDA_outliers_strong, type="n", scaling=3)
points(RDA_outliers_strong, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_outliers_strong, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables$location])
text(RDA_outliers_strong, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

###select SNPs associated only with 4 env variables (t and prec)
head(cand_corr_strong)
cand_corr_uni_4var<-cand_corr_uni %>% filter(predictor != "surface_area") %>% filter(predictor != "pH") 
cand_corr_strong_4var <- cand_corr_uni_4var[cand_corr_uni_4var$correlation > 0.5,]
RDA_outliers_4var <- rda(allfreq_geno[,cand_corr_uni_4var$snp] ~  bio5 + bio6 + bio15 + bio16, variables) 

## Adaptively enriched RDA only temperature and precipitation - strong outliers only
RDA_outliers_strong_4var <- rda(allfreq_geno[,cand_corr_strong_4var$snp] ~  bio5 + bio6 + bio15 + bio16, variables) 
TAB_loci_strong_4var <- as.data.frame(scores(RDA_outliers_strong_4var, choices=c(1:2), display="species", scaling="none"))
TAB_var_strong_4var <- as.data.frame(scores(RDA_outliers_strong_4var, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_strong_4var, aes(x=RDA1*10, y=RDA2*10), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  # geom_segment(data = TAB_var_strong, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_strong_4var, aes(x=RDA1, y=RDA2, label = row.names(TAB_var_strong_4var)), size = 2.5, family = "Times") +
  #xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space - strong outliers for temp and prec") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

pdf("plot_rda_outliers_strong_4var.pdf")
plot(RDA_outliers_strong_4var, type="n", scaling=3)
points(RDA_outliers_strong_4var, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_outliers_strong_4var, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables$location])
text(RDA_outliers_strong_4var, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

###adaptive index
source("./RDA-landscape-genomics-main/src/adaptive_index.R")
library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)

res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers_strong_4var, K = 2, env_pres = climRasts, method = "loadings", scale_env = scale_env, center_env = center_env)
## Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}
## Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
ggplot(data = TAB_RDA) + 
  #geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  #geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-143, -111), ylim = c(47, 61), expand = FALSE) + 
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

ggplot() + 
  #geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(data = TAB_RDA,aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  #geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-143, -111), ylim = c(47, 61), expand = FALSE) + 
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11), strip.background =element_rect(fill="white")) +
  geom_point(variables, mapping=aes(x,y), shape=21,fill="lightblue",size=2)



### Genomic offset
### prepare inoput files for genomic offset calculations in GradientForest 
###save subset of dataset - strong outliers
allfreq_geno_out<-allfreq_geno[,cand_corr_uni$snp]
allfreq_geno_out[1:12,1:9]
write.table(allfreq_geno_out, "allfreq_geno_outliers.txt", sep='\t')
dim(allfreq_geno_out)


allfreq_geno_strong<-allfreq_geno[,cand_corr_strong$snp]
allfreq_geno_strong[1:12,1:9]
write.table(allfreq_geno_strong, "allfreq_geno_strong_outliers.txt", sep='\t')
dim(allfreq_geno_strong)

allfreq_geno_strong_4var<-allfreq_geno[,cand_corr_strong_4var$snp]
allfreq_geno_strong_4var[1:12,1:9]
write.table(allfreq_geno_strong_4var, "allfreq_geno_strong_outliers_4var.txt", sep='\t')
dim(allfreq_geno_strong_4var)

###sample 1000 random SNPs
set.seed(1)
allfreq_geno_1k<-cbind(sample(allfreq_geno[-1], 1000))
allfreq_geno_1k[1:12,1:9]
write.table(allfreq_geno_1k, "allfreq_geno_1k_random.txt", sep='\t')

####import popmap for pop allele frequency calculation
popmap_noam<-read.table("popmap_allkokanee_nojapan.txt", header= TRUE)
popmap_noam_pop <- popmap_noam %>% select(ind, pop)
###estimated population allele frequencies
allfreq_geno_out_AF <- aggregate(allfreq_geno_out, by = list(popmap_noam_pop$pop), function(x) mean(x, na.rm = T)/2)
allfreq_geno_strong_AF <- aggregate(allfreq_geno_strong, by = list(popmap_noam_pop$pop), function(x) mean(x, na.rm = T)/2)
allfreq_geno_1k_AF <- aggregate(allfreq_geno_1k, by = list(popmap_noam_pop$pop), function(x) mean(x, na.rm = T)/2)
allfreq_geno_strong_4var_AF <- aggregate(allfreq_geno_strong_4var, by = list(popmap_noam_pop$pop), function(x) mean(x, na.rm = T)/2)

write.table(allfreq_geno_out_AF, "allfreq_geno_out_AF.txt", sep='\t')
write.table(allfreq_geno_strong_AF, "allfreq_geno_strong_AF.txt", sep='\t')
write.table(allfreq_geno_1k_AF, "allfreq_geno_1k_AF.txt", sep='\t')
write.table(allfreq_geno_strong_4var_AF, "allfreq_geno_strong_4var_AF", sep='\t')
