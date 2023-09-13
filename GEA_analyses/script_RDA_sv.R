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


###f###prepare inout files for RDA (vcf --> 012)
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
library(dplyr)
setwd("...")
###import genotypic data
sv.012_kok <- fread("structural variation/delly_allkokanee_gea_sv_filtered_pass_genotype_miss80_nobnd_na.012")[,-1] #load genotype matrix
sv.012_kok.pos <- fread("structural variation/delly_allkokanee_gea_sv_filtered_pass_genotype_miss80_nobnd.012.pos", header=FALSE) %>% #load SNPs info
  mutate(., locus=paste(V1,V2,sep='_')) #create a new column for SNP info name (CHR + position)
sv.012_kok.indv <- fread("structural variation/delly_allkokanee_gea_sv_filtered_pass_genotype_miss80_nobnd.012.indv", header=FALSE) #load individuals info
#evaluate % of missing
sum(is.na(sv.012_kok))/(dim(sv.012_kok)[1]*dim(sv.012_kok)[2]) #[1] 0.09359583
#impute missing with the most common geno
set.seed(1)
sv.012_kok.imp <- apply(sv.012_kok, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(sv.012_kok.imp))

sv.012_kok_df<-as.data.frame(sv.012_kok.imp)
sv.012_kok.pos_df<-as.data.frame(sv.012_kok.pos)
sv.012_kok.indv_df<-as.data.frame(sv.012_kok.indv)

#Set rownames and colnames to the geno matrix
dimnames(sv.012_kok_df) <- list(sv.012_kok.indv_df$V1, sv.012_kok.pos_df$locus)
#check the geno matrix
sv.012_kok_df[1:12,1:9]

allfreq_sv<-sv.012_kok_df

###remove all objects that are not necessary and only take a lot of space
#rm(list=ls(pattern="^sv"))

# ##import popmap
# popmap_noam<-read.table("popmap_4RDA_ecotypes_sv.txt", header= TRUE)
#popmap_noam_pop <- popmap_noam %>% select(ind, pop)
####import popmap for pop allele frequency calculation
popmap_noam_sv<-read.table("popmap_allkokanee_nojapan_sv.txt", header= TRUE)
popmap_noam_pop_sv <- popmap_noam_sv %>% dplyr::select(ind, pop)



###estimated population allele frequencies
allfreq_sv_AF <- aggregate(allfreq_sv, by = list(popmap_noam_pop_sv$pop), function(x) mean(x, na.rm = T)/2)

allfreq_sv_AF[1:12,1:9]
write.table(allfreq_sv_AF, "geno.012_kok_sv_af.txt", sep='\t')
write.table(allfreq_sv_AF[,cand_corr_sv_uni$snp], "geno.012_kok_sv_af_outliers.txt", sep='\t')



###  to calulate allele frequency for each ecotype in the  okanagan separately
popmap_eco_sv<-read.table("popmap_allkokanee_nojapan_okav_eco_sv.txt", header= TRUE)
popmap_eco_pop_sv <- popmap_eco_sv %>% dplyr::select(ind, pop)

allfreq_sv_eco_AF <- aggregate(allfreq_sv, by = list(popmap_eco_pop_sv$pop), function(x) mean(x, na.rm = T)/2)

allfreq_sv_eco_AF[1:12,1:9]
write.table(allfreq_sv_eco_AF, "geno.012_kok_sv_eco_af.txt", sep='\t')
cand_corr_sv_uni<-read.table("outliers_RDA_6var_pop_unique_sv.txt", header=TRUE)
write.table(allfreq_sv_eco_AF[,cand_corr_sv_uni$snp], "geno.012_kok_sv_eco_af_outliers.txt", sep='\t')

######RDA based on individual genotypes ###
###inferring population structure
### In this case from the same dataset
## Running a PCA on neutral genetic markers
variables_sv<-fread("env_variables_past_scale_pca_ecotype_sv.txt", header = TRUE)

pca_sv <- rda(allfreq_sv, scale=T)
screeplot(pca_sv, type = "barplot", npcs=10, main="PCA Eigenvalues")
## Neutral population structure table
PCs_sv <- scores(pca_sv, choices=c(1:3), display="sites", scaling=0)
PopStruct_sv <- data.frame(Ind = row.names(allfreq_sv), PCs_sv)
colnames(PopStruct_sv) <- c("Ind", "PC1", "PC2", "PC3")
rownames(PopStruct_sv) <- NULL
plot(PopStruct_sv$PC1, PopStruct_sv$PC2)
PopStruct_sv$location<-popmap_noam_sv$pop
PopStruct_sv$ecospawning<-popmap_noam_sv$spawning
PopStruct_sv$ecodepth<-popmap_noam_sv$depth

p1<-ggplot(PopStruct_sv, aes(PC1,PC2, color=location)) + geom_point()
p2<-ggplot(PopStruct_sv, aes(PC1,PC3, color=location)) + geom_point()
p3<-ggplot(PopStruct_sv, aes(PC2,PC3, color=location)) + geom_point()
ggplotly(p1)
ggplotly(p2)
ggplotly(p3)

# variables_sv
# bg<-c("#efea5a","#2c699a","#9d0208","#6a040f","#048ba8","#83e377","#d00000","#013a63","#d94a8c","#b9e769","#e9d8a6","#f94144","#990066","#16db93","#f1c453","#f29e4c","#0db39e","#fad2e1","#001219","#815ac0","#3c096c","#f3722c")
# new_bg<-c("#001219", "#3c096c", "#B785FF", "#013a63", "#4A9FFF", "#AAD8FF", "#005E29", "#00C05A", "#83D477", "#b9e769", "#e9d8a6", "#efea5a", "#f1c453", "#f29e4c", "#f3722c", "#f94144", "#d00000", "#9d0208", "#6a040f", "#990066", "#d94a8c", "#fad2e1")
# p1<-ggplot(variables_sv, aes(PC1, PC2, color = Population)) + 
#   geom_point() + scale_color_manual(values=bg) + theme_bw()

#install.packages("plotly")

new_bg<-c("#001219", "#5a189a", "#B785FF", "#00509d", "#4A9FFF", "#AAD8FF", "#005E29", "#00C05A", "#83D477", "#b9e769", "#e9d8a6", "#efea5a", "#f1c453", "#f29e4c", "#f3722c", "#f94144", "#d00000", "#9d0208", "#6a040f", "#990066", "#d94a8c", "#fad2e1")

pc12<-variables_sv %>%
  mutate(location= factor(location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Bobtail", "Arctic", "Puntzi", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%>%
  ggplot( aes(PC1, PC2, color = location, fill=location)) + 
  geom_point(shape=21, alpha=0.7) + scale_color_manual(values=new_bg) + scale_fill_manual(values=new_bg) + theme_bw() 

pc13<-variables_sv %>%
  mutate(location= factor(location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Bobtail", "Arctic", "Puntzi", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%>%
  ggplot( aes(PC1, PC3, color = location, fill=location)) + 
  geom_point(shape=21, alpha=0.7) + scale_color_manual(values=new_bg) + scale_fill_manual(values=new_bg) + theme_bw() 
# library(plotly)
# plot_ly(x=variables$PC1, y=variables$PC2, z=variables$PC3, type="scatter3d", mode="markers", color= variables$location)
library(patchwork)
combined_pcas<-pc12 + pc13 & theme(legend.position = "bottom")
combined_pcas_legend<- combined_pcas + plot_layout(guides="collect")
combined_pcas_legend + plot_annotation(title="Structural variation", tag_levels = 'A')






library(plotly)
ggplotly(p1)

variables_sv %>%
  mutate(location= factor(location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Arctic", "Puntzi", "Bobtail", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%>%
  ggplot( aes(PC1, PC2, color = location)) + 
  geom_point() + scale_color_manual(values=new_bg) + theme_bw() 

variables_sv %>%
  mutate(location= factor(location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Arctic", "Puntzi", "Bobtail", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%>%
  ggplot( aes(PC1, PC3, color = location)) + 
  geom_point() + scale_color_manual(values=new_bg) + theme_bw() 

variables_sv %>%
  mutate(location= factor(location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Arctic", "Puntzi", "Bobtail", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%>%
  ggplot( aes(PC2, PC3, color = location)) + 
  geom_point() + scale_color_manual(values=new_bg) + theme_bw() 


plot_ly(x=variables_sv$PC1, y=variables_sv$PC2, z=variables_sv$PC3, type="scatter3d", mode="markers", color= variables_sv$location)

data %>%
  arrange(val) %>%
  mutate(name = factor(name, levels=c("north", "north-east", "east", "south-east", "south", "south-west", "west", "north-west"))) %>%
  ggplot( aes(x=name, y=val)) +
  geom_segment( aes(xend=name, yend=0)) +
  geom_point( size=4, color="orange") +
  theme_bw() +
  
ggplot(variables_sv, aes(PC1, PC2, color = Population)) + 
  geom_point() 


# bio_past_scale<-read.table("env_variables_past_scale.txt", header = TRUE)
# variables_sv<-merge(PopStruct_sv,bio_past_scale, by = "location")
# write.table(variables_sv, "env_variables_past_scale_pca_ecotype_sv.txt ", sep="\t", row.names= FALSE)
variables_sv<-fread("env_variables_past_scale_pca_ecotype_sv.txt", header = TRUE)
#variables_sv<-read.table("env_variables_past_scale_pca_ecotype_sv.txt", header=TRUE)
## Null model
RDA0_sv <- rda(allfreq_sv ~ 1,  variables_sv) 
## Full model
RDA_6var_nopop_sv <- rda(allfreq_sv ~  surface_area + pH + bio5 + bio6 + bio15 + bio16, variables_sv)

RDA_6var_pop_sv <- rda(allfreq_sv ~  surface_area + pH + bio5 + bio6 + bio15 + bio16 + Condition(PC1 + PC2 + PC3), variables_sv)
RDA_4var_pop_sv <- rda(allfreq_sv ~ bio5 + bio6 + bio15 + bio16 + Condition(PC1 + PC2 + PC3), variables_sv)

RsquareAdj(RDA_6var_nopop_sv)
# $r.squared
# [1] 0.136594
# 
# $adj.r.squared
# [1] 0.112499

RsquareAdj(RDA_6var_pop_sv)
# $r.squared
# [1] 0.08998682
# 
# $adj.r.squared
# [1] 0.06832477

RsquareAdj(RDA_4var_pop_sv)
# $r.squared
# [1] 0.0631016
# 
# $adj.r.squared
# [1] 0.04833627



## Stepwise procedure with ordiR2step function
mod1_sv <- ordiR2step(RDA0_sv, RDA_6var_nopop_sv, Pin = 0.01, R2permutations = 1000, R2scope = T)
screeplot(RDA_6var_nopop_sv, main="Eigenvalues of constrained axes")
screeplot(RDA_6var_pop_sv, main="Eigenvalues of constrained axes") #2 axes
screeplot(RDA_4var_pop_sv, main="Eigenvalues of constrained axes") #2 axes

levels(variables_sv$location)<-c("Anderson","Arctic","Arrow_Hill","Arrow_Mosq","Bobtail","Bonaparte",   
                              "Christina","Cluculz","Cowichan","Dunn","EastBarriere","Kalamalka",   
                              "Kootenay","La_Hache","Nicola","Okanagan","Puntzi","Shawningan",  
                              "Sockeye","Tchesinkut","Thutade", "Wood")
bg<-c("#efea5a","#2c699a","#9d0208","#6a040f","#048ba8","#83e377","#d00000","#013a63","#d94a8c","#b9e769","#e9d8a6","#f94144","#990066","#16db93","#f1c453","#f29e4c","#0db39e","#fad2e1","#001219","#815ac0","#3c096c","#f3722c")
names(bg)<-levels(variables_sv$location)

#plot1
pdf(file="plot_RDA_6var_nopop_sv.pdf")
plot(RDA_6var_nopop_sv, type="n", scaling=3)
points(RDA_6var_nopop_sv, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_6var_nopop_sv, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables_sv$location])
text(RDA_6var_nopop_sv, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables_sv$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

#plot2
pdf(file="plot_RDA_6var_pop_sv.pdf")
plot(RDA_6var_pop_sv, type="n", scaling=3)
points(RDA_6var_pop_sv, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_6var_pop_sv, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables_sv$location])
text(RDA_6var_pop_sv, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables_sv$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

#plot3
pdf(file="plot_RDA_4var_pop_sv.pdf")
plot(RDA_4var_pop_sv, type="n", scaling=3)
points(RDA_4var_pop_sv, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_4var_pop_sv, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables_sv$location])
text(RDA_4var_pop_sv, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables_sv$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

###variance partitioning 
## full model
pRDAfull_sv<-rda(allfreq_sv ~ ecodepth + ecospawning + x + y + PC1 + PC2 + PC3 + surface_area + pH + bio5 + bio6 + bio15 + bio16, variables_sv)
### climate model
pRDAclim_sv<-rda(allfreq_sv ~  surface_area + pH + bio5 + bio6 + bio15 + bio16 + Condition(ecodepth + ecospawning + x + y + PC1 + PC2 + PC3), variables_sv)
### pop structure model
pRDAstruc_sv<-rda(allfreq_sv ~ PC1 + PC2 + PC3 + Condition(ecodepth + ecospawning + x + y + surface_area + pH + bio5 + bio6 + bio15 + bio16), variables_sv)
### geography model
pRDAgeog_sv<-rda(allfreq_sv ~  x + y + Condition(ecodepth + ecospawning + PC1 +PC2 +PC3 + surface_area + pH + bio5 + bio6 + bio15 + bio16), variables_sv)
## ecotype model
pRDAecot_sv<-rda(allfreq_sv ~  ecodepth + ecospawning + Condition(x + y +  PC1 + PC2 + PC3 + surface_area + pH + bio5 + bio6 + bio15 + bio16), variables_sv)


RsquareAdj(pRDAfull_sv)
RsquareAdj(pRDAclim_sv)
RsquareAdj(pRDAstruc_sv)
RsquareAdj(pRDAgeog_sv)
RsquareAdj(pRDAecot_sv)
# > RsquareAdj(pRDAfull_sv)
# $r.squared
# [1] 0.2528753
# 
# $adj.r.squared
# [1] 0.20618
# 
# > RsquareAdj(pRDAclim_sv)
# $r.squared
# [1] 0.07311407
# 
# $adj.r.squared
# [1] 0.05324902
# 
# > RsquareAdj(pRDAstruc_sv)
# $r.squared
# [1] 0.04948814
# 
# $adj.r.squared
# [1] 0.04054701
# 
# > RsquareAdj(pRDAgeog_sv)
# $r.squared
# [1] 0.0264955
# 
# $adj.r.squared
# [1] 0.02032317
# 
# > RsquareAdj(pRDAecot_sv)
# $r.squared
# [1] 0.02505478
# 
# $adj.r.squared
# [1] 0.01880698
###Identifying candidate genes based on 6 variables (including area and ph)
### method 1 - outlier based on SD
### correlation each SNP
load.rda_sv <- scores(RDA_6var_pop_sv, choices=c(1:3), display="species")
hist(load.rda_sv[,1], main="Loadings on RDA1")
hist(load.rda_sv[,2], main="Loadings on RDA2")
hist(load.rda_sv[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}


cand1_sv <- outliers(load.rda_sv[,1],3) # 
cand2_sv <- outliers(load.rda_sv[,2],3) # 
cand3_sv <- outliers(load.rda_sv[,3],3) # 

ncand_sv <- length(cand1_sv) + length(cand2_sv) + length(cand3_sv) 
ncand_sv #329
cand1_sv <- cbind.data.frame(rep(1,times=length(cand1_sv)), names(cand1_sv), unname(cand1_sv))
cand2_sv <- cbind.data.frame(rep(2,times=length(cand2_sv)), names(cand2_sv), unname(cand2_sv))
cand3_sv <- cbind.data.frame(rep(3,times=length(cand3_sv)), names(cand3_sv), unname(cand3_sv))

colnames(cand1_sv) <- colnames(cand2_sv)  <- colnames(cand3_sv) <- c("axis","snp","loading")

cand_sv <- rbind(cand1_sv, cand2_sv, cand3_sv)
cand_sv$snp <- as.character(cand_sv$snp)

### 
foo_sv <- matrix(nrow=ncand_sv, ncol=6)

colnames(foo_sv) <- c( "surface_area" , "pH" , "bio5", "bio6", "bio15",	"bio16")
library(dplyr)
var_env_sv<-variables_sv %>% dplyr::select(surface_area , pH , bio5, bio6, bio15,	bio16)
#var_env_num<-as.numeric(unlist(var_env))
# for (i in 1:length(cand$snp)) {
#   nam <- cand[i,2]
#   snp.gen <- allfreq_geno[,nam]
#   foo[i,] <- apply(var_env,2,function(x) cor(x,snp.gen))
# }


for (i in 1:length(cand_sv$snp)) {
  nam <- cand_sv[i,2]
  # print(nam)
  snp.gen <- allfreq_sv[,nam]
  #  print(snp.gen)
  foo_sv[i,] <- apply(var_env_sv,2,function(x) cor(x,snp.gen))
  #  print(foo[i,])
}
head(foo_sv)

cand_corr_sv <- cbind.data.frame(cand_sv,foo_sv)  
head(cand_corr_sv)
write.table(cand_corr_sv, "outliers_RDA_6var_pop_sv.txt", sep = "\t",row.names = FALSE)

for (i in 1:length(cand_corr_sv$snp)) {
  bar <- cand_corr_sv[i,]
  cand_corr_sv[i,10] <- names(which.max(abs(bar[4:9]))) # gives the variable
  cand_corr_sv[i,11] <- max(abs(bar[4:9]))              # gives the correlation
}

colnames(cand_corr_sv)[10] <- "predictor"
colnames(cand_corr_sv)[11] <- "correlation"

cand_corr_sv_uni <- cand_corr_sv[!duplicated(cand_corr_sv$snp),]
write.table(cand_corr_sv_uni, "outliers_RDA_6var_pop_unique_sv.txt", sep = "\t",row.names = FALSE)
table(cand_corr_sv_uni$predictor) 
# 
# bio15        bio16         bio5         bio6           pH surface_area 
# 29           22           47           76           94           54 

cand_corr_strong_sv <- cand_corr_sv_uni[cand_corr_sv_uni$correlation > 0.5,]
write.table(cand_corr_strong_sv, "outliers_RDA_6var_pop_unique_strong_sv.txt", sep = "\t",row.names = FALSE)
table(cand_corr_strong_sv$predictor) 
# bio5 bio6 
# 2    3 
plot(density(cand_corr_sv_uni$correlation))


sel_sv <- cand_corr_sv$snp
env_sv <- cand_corr_sv$predictor
#"surface_area" , "pH" , "bio5", "bio6", "bio15",	"bio16"
env_sv[env_sv=="surface_area"] <- '#1f78b4'
env_sv[env_sv=="pH"] <- '#a6cee3'
env_sv[env_sv=="bio5"] <- '#6a3d9a'
env_sv[env_sv=="bio6"] <- '#e31a1c'
env_sv[env_sv=="bio15"] <- '#33a02c'
env_sv[env_sv=="bio16"] <- '#ffff33'


# color by predictor:
col.pred <- rownames(RDA_6var_pop_sv$CCA$v) # pull the SNP names

for (i in 1:length(sel_sv)) {           # color code candidate SNPs
  foo <- match(sel_sv[i],col.pred)
  col.pred[foo] <- env_sv[i]
}

col.pred[grep("NC",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg2 <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c', '#33a02c','#ffff33')

# axes 1 & 2
pdf(file="plot_RDA_6var_pop_onlyoutliers_sv.pdf")

plot(RDA_6var_pop_sv, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(RDA_6var_pop_sv, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(RDA_6var_pop_sv, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(RDA_6var_pop_sv, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("surface_area" , "pH" , "bio5", "bio6", "bio15",	"bio16"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg2)

dev.off()


###### Adaptive landscape
## Adaptively enriched RDA
#cand_corr_sv_uni<-read.table("outliers_RDA_6var_pop_unique_sv.txt", header = TRUE)
#cand_corr_strong_sv <- cand_corr_uni_sv[cand_corr_uni_sv$correlation > 0.5,]
cand_corr_sv_uni<-fread("outliers_RDA_6var_pop_unique_sv.txt", header=TRUE)
RDA_outliers_sv <- rda(allfreq_sv[,cand_corr_sv_uni$snp] ~ surface_area  + pH + bio5 + bio6 + bio15 + bio16, variables_sv) 
# RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers_sv, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers_sv, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*5, y=RDA2*5), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=RDA1, y=RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  #xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))


pdf("plot_rda_outliers_sv.pdf")
plot(RDA_outliers_sv, type="n", scaling=3)
points(RDA_outliers_sv, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables_sv$location])
points(RDA_outliers_sv, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
text(RDA_outliers_sv, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables_sv$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

## Adaptively enriched RDA with strong candidates
RDA_outliers_strong_sv <- rda(allfreq_sv[,cand_corr_strong_sv$snp] ~ surface_area  + pH + bio5 + bio6 + bio15 + bio16, variables_sv) 
TAB_loci_strong_sv <- as.data.frame(scores(RDA_outliers_strong_sv, choices=c(1:2), display="species", scaling="none"))
TAB_var_strong_sv <- as.data.frame(scores(RDA_outliers_strong_sv, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_strong_sv, aes(x=RDA1*10, y=RDA2*10), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  # geom_segment(data = TAB_var_strong, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_strong_sv, aes(x=RDA1, y=RDA2, label = row.names(TAB_var_strong_sv)), size = 2.5, family = "Times") +
  #xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space - strong outliers") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

pdf("plot_rda_outliers_strong_sv.pdf")
plot(RDA_outliers_strong_sv, type="n", scaling=3)
points(RDA_outliers_strong_sv, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_outliers_strong_sv, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables_sv$location])
text(RDA_outliers_strong_sv, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables_sv$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()


###### 
###select SNPs associated only with 4 env variables (t and prec)
head(cand_corr_strong_sv)
cand_corr_uni_4var_sv<-cand_corr_sv_uni %>% filter(predictor != "surface_area") %>% filter(predictor != "pH") 
#cand_corr_strong_4var <- cand_corr_uni_4var[cand_corr_uni_4var$correlation > 0.5,]

## Adaptively enriched RDA only temperature and precipitation
RDA_outliers_4var_sv <- rda(allfreq_sv[,cand_corr_uni_4var_sv$snp] ~  bio5 + bio6 + bio15 + bio16, variables_sv) 
TAB_loci_4var_sv <- as.data.frame(scores(RDA_outliers_4var_sv, choices=c(1:2), display="species", scaling="none"))
TAB_var_4var_sv <- as.data.frame(scores(RDA_outliers_4var_sv, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_4var_sv, aes(x=RDA1*10, y=RDA2*10), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var_4var_sv, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_4var_sv, aes(x=RDA1, y=RDA2, label = row.names(TAB_var_4var_sv)), size = 2.5, family = "Times") +
  #xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space - SVs outliers for temp and prec") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

pdf("plot_rda_outliers_4var_sv.pdf")
plot(RDA_outliers_4var_sv, type="n", scaling=3)
points(RDA_outliers_4var_sv, display="species", pch=20, cex=0.7, col="gray32", scaling=3)  
points(RDA_outliers_4var_sv, display="sites", pch=21, cex=1.3, col="lightgrey", scaling=3, bg=bg[variables_sv$location])
text(RDA_outliers_4var_sv, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(variables_sv$location), bty="n", col="lightgrey", pch=21, cex=1, pt.bg=bg)
dev.off()

###adaptive index
source("../RDA-landscape-genomics-main/src/adaptive_index.R")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")
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
#library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)

res_RDA_proj_current_sv <- adaptive_index(RDA = RDA_outliers_4var_sv, K = 2, env_pres = climRasts, method = "loadings", scale_env = scale_env, center_env = center_env)
## Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_proj_current_sv$RDA1, res_RDA_proj_current_sv$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}
## Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
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
  geom_point(variables_sv, mapping=aes(x,y), shape=21,fill="lightblue",size=2)

