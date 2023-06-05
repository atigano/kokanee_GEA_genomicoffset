setwd("~/Dropbox/postdoc UBC/analyses_allkokanee")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
library(dplyr)
library(tidyverse)
detach("package:tidyverse", unload=TRUE)
setwd("~/Dropbox/postdoc UBC/analyses_allkokanee//")
#install.packages(c("geodata","raster","sp","rgdal"))
library(geodata)
library(raster)
library(sp)
library(rgdal)
library(data.table)
library(gtools)
library(colorRamps)

#libraries <- c("rgdal", "raster", "gtools", "adespatial", "ade4", "adegraphics", "spdep", "maptools", "adegenet", "qvalue", "pcadapt", "gdm", "gradientForest", "vegan", "rangeBuilder", "rgeos", "assigner", "ggfortify")
#sapply(libraries, function(x) { suppressMessages(require(x, character.only=T)) } )


env_scale<-read.table("env_variables_pastfuture.txt", header = TRUE)
setnames(env_scale, new = c('bio1','bio2','bio3','bio4','bio5', 'bio6', 'bio7', 'bio8', 'bio9'), 
         old = c('bio01','bio02','bio03','bio04','bio05', 'bio06', 'bio07', 'bio08', 'bio09'))
env_past<-env_scale %>% dplyr::filter(dataset=="past") %>% dplyr::select(-x, -y, -degree_days,-location,-region, -dataset, -tds)

###analysis with strong sequence outliers
geno_outliers_strong<-read.table("allfreq_geno_strong_AF_adj.txt", header =TRUE)
snps_out_strong<-geno_outliers_strong%>%dplyr::select(contains("NC"))

maxLevel<-log2(0.368*nrow(env_past)/2)
gf_out_strong<-gradientForest(data=cbind(env_past, snps_out_strong),
                              predictor.vars=colnames(env_past),
                              response.vars=colnames(snps_out_strong),
                              ntree=500,
                              maxLevel=maxLevel,
                              trace = T,
                              corr.threshold = 0.50)
types<-c("Overall.Importance",
         "Split.Density",
         "Cumulative.Importance",
         "Performance")
sel<-1
plot(gf_out_strong,plot.type=types[sel])

  


by.importance_out_strong<-names(importance(gf_out_strong))
      
plot(gf_out_strong,plot.type = "C",
     show.species = F,
     common.scale=T,
     imp.vars=by.importance_out_strong,
     lwd=2,
     col="blue",
     cex.lab=1.5,
     line.ylab=0.9,
     par.args=list(mgp=c(1.5,0.5,0),
                   mar=c(2.5,1,1,1),
                   omi=c(0,0.3,0,0)))

by.importance_ran<-names(importance(gf_ran))

plot(gf_ran,plot.type = "C",
     show.species = F,
     common.scale=T,
     imp.vars=by.importance_ran,
     lwd=2,
     col="blue",
     cex.lab=1.5,
     line.ylab=0.9,
     par.args=list(mgp=c(1.5,0.5,0),
                   mar=c(2.5,1,1,1),
                   omi=c(0,0.3,0,0)))

###sv analysis - only outliers

geno_outliers_sv<-read.table("geno.012_kok_sv_af_outliers.txt", header =TRUE)
snps_out_sv<-geno_outliers_sv%>%dplyr::select(contains("NC"))

maxLevel<-log2(0.368*nrow(env_past)/2)
gf_out_sv<-gradientForest(data=cbind(env_past, snps_out_sv),
                              predictor.vars=colnames(env_past),
                              response.vars=colnames(snps_out_sv),
                              ntree=500,
                              maxLevel=maxLevel,
                              trace = T,
                              corr.threshold = 0.50)
types<-c("Overall.Importance",
         "Split.Density",
         "Cumulative.Importance",
         "Performance")
sel<-1
plot(gf_out_sv,plot.type=types[sel])


### strong outliers but with ecotypes split in okanagan lakes
env_scale_oka_eco<-read.table("env_variables_pastfuture_oka_eco.txt", header = TRUE)
setnames(env_scale_oka_eco, new = c('bio1','bio2','bio3','bio4','bio5', 'bio6', 'bio7', 'bio8', 'bio9'), 
         old = c('bio01','bio02','bio03','bio04','bio05', 'bio06', 'bio07', 'bio08', 'bio09'))
env_past_oka_eco<-env_scale_oka_eco %>% dplyr::filter(dataset=="past") %>% dplyr::select(-x, -y, -degree_days,-location,-region, -dataset, -tds)

geno_outliers_oka<-read.table("allfreq_geno_strong_AF_oka_eco.txt", header =TRUE)
snps_out_oka<-geno_outliers_oka%>%dplyr::select(contains("NC"))
maxLevel<-log2(0.368*nrow(env_past_oka_eco)/2)
gf_out_oka_eco<-gradientForest(data=cbind(env_past_oka_eco, snps_out_oka),
                          predictor.vars=colnames(env_past_oka_eco),
                          response.vars=colnames(snps_out_oka),
                          ntree=500,
                          maxLevel=maxLevel,
                          trace = T,
                          corr.threshold = 0.50)


geno_outliers_oka_sv<-read.table("geno.012_kok_sv_eco_af_outliers.txt", header =TRUE)
sv_out_oka<-geno_outliers_oka_sv%>%dplyr::select(contains("NC"))
maxLevel<-log2(0.368*nrow(env_past_oka_eco)/2)
gf_out_oka_eco_sv<-gradientForest(data=cbind(env_past_oka_eco, sv_out_oka),
                               predictor.vars=colnames(env_past_oka_eco),
                               response.vars=colnames(sv_out_oka),
                               ntree=500,
                               maxLevel=maxLevel,
                               trace = T,
                               corr.threshold = 0.50)


###simulate 


### prepare map data
### current
rastFiles <- list.files(path="./wc5",
                        full.names=T,
                        pattern=".bil")
#get files in correct order
rastFiles <- mixedsort(rastFiles)

# stack, subset, and clip to species range polygon
climRasts <- stack(rastFiles)
names(climRasts) <- paste("bio", 1:19, sep="")
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
#### future projections 
###2040-2060
### worst scenario
cmip6Rasts <- stack("./wc2.1_5m/wc2.1_5m_bioc_UKESM1-0-LL_ssp585_2041-2060.tif")
names(cmip6Rasts) <- paste("bio", 1:19, sep="")
#cmip6Rasts <- cmip6Rasts[[which(names(cmip6Rasts) %in% colnames(envGF))]]
cmip6Names <- names(cmip6Rasts)
cmip6Rasts <- cmip6Rasts*mask
names(cmip6Rasts) <- cmip6Names
plot(cmip6Rasts, col=rgb.tables(1000))
deltarastw<-cmip6Rasts-climRasts
plot(deltarastw, col=rgb.tables(1000))

sum(is.na(climRasts[[1]][]))
sum(is.na(cmip6Rasts[[1]][]))
climRasts[[6]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[5]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[4]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[3]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[2]][which(is.na(cmip6Rasts[[1]][]))] <- NA
climRasts[[1]][which(is.na(cmip6Rasts[[1]][]))] <- NA
cmip6Rasts <- crop(cmip6Rasts, extent)

### best scenario
cmip6Rastsb <- stack("./wc2.1_5m/wc2.1_5m_bioc_UKESM1-0-LL_ssp245_2041-2060.tif")
names(cmip6Rastsb) <- paste("bio", 1:19, sep="")
#cmip6Rasts <- cmip6Rasts[[which(names(cmip6Rasts) %in% colnames(envGF))]]
cmip6Namesb <- names(cmip6Rastsb)
cmip6Rastsb <- cmip6Rastsb*mask
names(cmip6Rastsb) <- cmip6Namesb
plot(cmip6Rastsb, col=rgb.tables(1000))
deltarastb<-cmip6Rastsb-climRasts
plot(deltarastb, col=rgb.tables(1000))

sum(is.na(climRasts[[1]][]))
sum(is.na(cmip6Rastsb[[1]][]))
climRasts[[6]][which(is.na(cmip6Rastsb[[1]][]))] <- NA
climRasts[[5]][which(is.na(cmip6Rastsb[[1]][]))] <- NA
climRasts[[4]][which(is.na(cmip6Rastsb[[1]][]))] <- NA
climRasts[[3]][which(is.na(cmip6Rastsb[[1]][]))] <- NA
climRasts[[2]][which(is.na(cmip6Rastsb[[1]][]))] <- NA
climRasts[[1]][which(is.na(cmip6Rastsb[[1]][]))] <- NA
cmip6Rastsb <- crop(cmip6Rastsb, extent)


# GF functions will not take rasters, so need to extract data from
# the rasters (cell ID as well)
env_trns <- raster::extract(climRasts, 1:ncell(climRasts[[1]]))
env_trns <- data.frame(cell=1:ncell(climRasts[[1]]), env_trns)
env_trns <- na.omit(env_trns)

# same for future
# worst scenario
env_trns_futurew <- raster::extract(cmip6Rasts, 1:ncell(cmip6Rasts[[1]]))
env_trns_futurew <- data.frame(cell=1:ncell(cmip6Rasts[[1]]), env_trns_futurew)
env_trns_futurew <- na.omit(env_trns_futurew)

# worst scenario
env_trns_futureb <- raster::extract(cmip6Rastsb, 1:ncell(cmip6Rastsb[[1]]))
env_trns_futureb <- data.frame(cell=1:ncell(cmip6Rastsb[[1]]), env_trns_futureb)
env_trns_futureb <- na.omit(env_trns_futureb)


# check congruence
dim(env_trns)
dim(env_trns_futurew)
dim(env_trns_futureb)
#transform env using gf models
predCand<-predict(gf_out_strong, env_trns[,-1])
#plot(predCand,col=rgb.tables(1000)) ### this is not working


###worst scenario
transCandw<-predict(gf_out_strong, env_trns_futurew[,-1])
genoffset_outw<-sqrt(rowSums((transCandw-predCand)^2))
genoffsetw<-mask
genoffsetw[env_trns_futurew$cell] <-genoffset_outw
plot(genoffsetw,zlim = c(0,0.5),col=rgb.tables(1000))

####sequence variation
###best scenario
transCandb<-predict(gf_out_strong, env_trns_futureb[,-1])
genoffset_outb<-sqrt(rowSums((transCandb-predCand)^2))
genoffsetb<-mask
genoffsetb[env_trns_futureb$cell] <-genoffset_outb
plot(genoffsetb,zlim = c(0,0.5),col=rgb.tables(1000))
###best scenario only bio5
imp.vars<-c("bio5")
predCand_temp<-predict(gf_out_strong, env_trns[imp.vars])
transCand_tempb<-predict(gf_out_strong, env_trns_futureb[imp.vars])
genoffset_out_tempb<-sqrt(rowSums((transCand_tempb-predCand_temp)^2))
genoffset_tempb<-mask
genoffset_tempb[env_trns_futureb$cell] <-genoffset_out_tempb
plot(genoffset_tempb,zlim = c(0,0.038), col=rgb.tables(1000))

###worst scenario
#imp.vars<-c("bio5")
predCand_temp<-predict(gf_out_strong, as.data.frame(env_trns[imp.vars]))
transCand_tempw<-predict(gf_out_strong, env_trns_futurew[imp.vars])
genoffset_out_tempw<-sqrt(rowSums((transCand_tempw-predCand_temp)^2))
genoffset_tempw<-mask
genoffset_tempw[env_trns_futurew$cell] <-genoffset_out_tempw
plot(genoffset_tempw,zlim = c(0,0.038), col=rgb.tables(1000))


coords<-env_scale %>% dplyr::filter(dataset=="past") %>% dplyr::select(x, y)
points<-sp::SpatialPoints(coords, proj4string = genoffset_tempw@crs)
genomicofsset_values_b<-raster::extract(genoffset_tempb, coords)
genomicofsset_values_w<-raster::extract(genoffset_tempw, coords)

genomicofsset_b<-cbind(coords, genomicofsset_values_b)
genomicofsset_w<-cbind(coords, genomicofsset_values_w)

names(genomicofsset_b)[3] ="gen_off"
names(genomicofsset_w)[3] ="gen_off"

genomicofsset_b<-cbind(env_scale[1:22,1], genomicofsset_b)
names(genomicofsset_b)[1] ="location"
genomicofsset_w<-cbind(env_scale[1:22,1], genomicofsset_w)
names(genomicofsset_w)[1] ="location"

max(genomicofsset_b$gen_off) ### tchesinkut
max(genomicofsset_w$gen_off) ### bonaparte

###### genomic offset structural variation
###best scenario
#transform env using gf models
predCand_sv<-predict(gf_out_sv, env_trns[,-1])
#plot(predCand,col=rgb.tables(1000)) ### this is not working

transCandb_sv<-predict(gf_out_sv, env_trns_futureb[,-1])
genoffset_outb_sv<-sqrt(rowSums((transCandb_sv-predCand_sv)^2))
genoffsetb_sv<-mask
genoffsetb_sv[env_trns_futureb$cell] <-genoffset_outb_sv
plot(genoffsetb_sv,zlim = c(0,0.25),col=rgb.tables(1000))

###best scenario only bio5
imp.vars<-c("bio5")
predCand_temp_sv<-predict(gf_out_sv, env_trns[imp.vars])
transCand_tempb_sv<-predict(gf_out_sv, env_trns_futureb[imp.vars])
genoffset_out_tempb_sv<-sqrt(rowSums((transCand_tempb_sv-predCand_temp_sv)^2))
genoffset_tempb_sv<-mask
genoffset_tempb_sv[env_trns_futureb$cell] <-genoffset_out_tempb_sv
plot(genoffset_tempb_sv,zlim = c(0,0.012), col=rgb.tables(1000))

###worst scenario
#imp.vars<-c("bio5")
predCand_temp_sv<-predict(gf_out_sv, as.data.frame(env_trns[imp.vars]))
transCand_tempw_sv<-predict(gf_out_sv, env_trns_futurew[imp.vars])
genoffset_out_tempw_sv<-sqrt(rowSums((transCand_tempw_sv-predCand_temp_sv)^2))
genoffset_tempw_sv<-mask
genoffset_tempw_sv[env_trns_futurew$cell] <-genoffset_out_tempw_sv
plot(genoffset_tempw_sv,zlim = c(0,0.012), col=rgb.tables(1000))


coords<-env_scale %>% dplyr::filter(dataset=="past") %>% dplyr::select(x, y)
points<-sp::SpatialPoints(coords, proj4string = genoffset_tempb@crs)
genomicofsset_values_b_sv<-raster::extract(genoffset_tempb_sv, coords)
genomicofsset_values_w_sv<-raster::extract(genoffset_tempw_sv, coords)

genomicofsset_b_sv<-cbind(coords, genomicofsset_values_b_sv)
genomicofsset_w_sv<-cbind(coords, genomicofsset_values_w_sv)

names(genomicofsset_b_sv)[3] ="gen_off"
names(genomicofsset_w_sv)[3] ="gen_off"

genomicofsset_b_sv<-cbind(env_scale[1:22,1], genomicofsset_b_sv)
names(genomicofsset_b_sv)[1] ="location"
genomicofsset_w_sv<-cbind(env_scale[1:22,1], genomicofsset_w_sv)
names(genomicofsset_w_sv)[1] ="location"


deltarastb_df<-as.data.frame(deltarastb, xy=TRUE)
deltarastw_df<-as.data.frame(deltarastw, xy=TRUE)


###genomic offset okanagan ecotypes separated
#transform env using gf models
predCand_oka<-predict(gf_out_oka_eco, env_trns[,-1])
#plot(predCandoka,col=rgb.tables(1000)) ### this is not working


###worst scenario
transCandw_oka<-predict(gf_out_oka_eco, env_trns_futurew[,-1])
genoffset_outw_oka<-sqrt(rowSums((transCandw_oka-predCand_oka)^2))
genoffsetw_oka<-mask
genoffsetw_oka[env_trns_futurew$cell] <-genoffset_outw_oka
plot(genoffsetw_oka,zlim = c(0,0.6),col=rgb.tables(1000))

###best scenario
transCandb_oka<-predict(gf_out_oka_eco, env_trns_futureb[,-1])
genoffset_outb_oka<-sqrt(rowSums((transCandb_oka-predCand_oka)^2))
genoffsetb_oka<-mask
genoffsetb_oka[env_trns_futureb$cell] <-genoffset_outb_oka
plot(genoffsetb,zlim = c(0,0.25),col=rgb.tables(1000))
###best scenario only bio5
imp.vars<-c("bio5")
predCand_temp_oka<-predict(gf_out_oka_eco, env_trns[imp.vars])
transCand_tempb_oka<-predict(gf_out_oka_eco, env_trns_futureb[imp.vars])
genoffset_out_tempb_oka<-sqrt(rowSums((transCand_tempb_oka-predCand_temp_oka)^2))
genoffset_tempb_oka<-mask
genoffset_tempb_oka[env_trns_futureb$cell] <-genoffset_out_tempb_oka
plot(genoffset_tempb_oka,zlim = c(0,0.04), col=rgb.tables(1000))
###worst scenario
#imp.vars<-c("bio5")
predCand_temp_oka<-predict(gf_out_oka_eco, as.data.frame(env_trns[imp.vars]))
transCand_tempw_oka<-predict(gf_out_oka_eco, env_trns_futurew[imp.vars])
genoffset_out_tempw_oka<-sqrt(rowSums((transCand_tempw_oka-predCand_temp_oka)^2))
genoffset_tempw_oka<-mask
genoffset_tempw_oka[env_trns_futurew$cell] <-genoffset_out_tempw_oka
plot(genoffset_tempw_oka,zlim = c(0,0.04), col=rgb.tables(1000))

coords_eco<-env_scale_oka_eco %>% dplyr::filter(dataset=="past") %>% dplyr::select(x, y)
point_okaeco<-sp::SpatialPoints(coords_eco, proj4string = genoffset_tempb_oka@crs)
genomicofsset_values_b_oka<-raster::extract(genoffset_tempb_oka, coords_eco)
genomicofsset_values_w_oka<-raster::extract(genoffset_tempw_oka, coords_eco)

genomicofsset_b_oka<-cbind(coords_eco, genomicofsset_values_b_oka)
genomicofsset_w_oka<-cbind(coords_eco, genomicofsset_values_w_oka)

names(genomicofsset_b_oka)[3] ="gen_off"
names(genomicofsset_w_oka)[3] ="gen_off"

genomicofsset_b_oka<-cbind(env_scale_oka_eco[1:25,1], genomicofsset_b_oka)
names(genomicofsset_b_oka)[1] ="location"
genomicofsset_w_oka<-cbind(env_scale_oka_eco[1:25,1], genomicofsset_w_oka)
names(genomicofsset_w_oka)[1] ="location"

genomicofsset_b_oka
genomicofsset_w_oka

plot(deltarastb[["bio5"]],zlim = c(-2,13), col=rgb.tables(1000))
plot(deltarastw[["bio5"]],zlim = c(-2,13), col=rgb.tables(1000))


### okanagan ecortpes split structural outliers
###genomic offset okanagan ecotypes separated
#transform env using gf models
predCand_oka_sv<-predict(gf_out_oka_eco_sv, env_trns[,-1])
#plot(predCandoka,col=rgb.tables(1000)) ### this is not working


###worst scenario
transCandw_oka_sv<-predict(gf_out_oka_eco_sv, env_trns_futurew[,-1])
genoffset_outw_oka_sv<-sqrt(rowSums((transCandw_oka_sv-predCand_oka_sv)^2))
genoffsetw_oka_sv<-mask
genoffsetw_oka_sv[env_trns_futurew$cell] <-genoffset_outw_oka_sv
plot(genoffsetw_oka_sv,zlim = c(0,0.30),col=rgb.tables(1000))

###best scenario
transCandb_oka_sv<-predict(gf_out_oka_eco_sv, env_trns_futureb[,-1])
genoffset_outb_oka_sv<-sqrt(rowSums((transCandb_oka_sv-predCand_oka_sv)^2))
genoffsetb_oka_sv<-mask
genoffsetb_oka_sv[env_trns_futureb$cell] <-genoffset_outb_oka_sv
plot(genoffsetb_oka_sv,zlim = c(0,0.3),col=rgb.tables(1000))
###best scenario only bio5
imp.vars<-c("bio5")
predCand_temp_oka_sv<-predict(gf_out_oka_eco_sv, env_trns[imp.vars])
transCand_tempb_oka_sv<-predict(gf_out_oka_eco_sv, env_trns_futureb[imp.vars])
genoffset_out_tempb_oka_sv<-sqrt(rowSums((transCand_tempb_oka_sv-predCand_temp_oka_sv)^2))
genoffset_tempb_oka_sv<-mask
genoffset_tempb_oka_sv[env_trns_futureb$cell] <-genoffset_out_tempb_oka_sv
plot(genoffset_tempb_oka_sv,zlim = c(0,0.01), col=rgb.tables(1000))
###worst scenario
#imp.vars<-c("bio5")
predCand_temp_oka_sv<-predict(gf_out_oka_eco_sv, as.data.frame(env_trns[imp.vars]))
transCand_tempw_oka_sv<-predict(gf_out_oka_eco_sv, env_trns_futurew[imp.vars])
genoffset_out_tempw_oka_sv<-sqrt(rowSums((transCand_tempw_oka_sv-predCand_temp_oka_sv)^2))
genoffset_tempw_oka_sv<-mask
genoffset_tempw_oka_sv[env_trns_futurew$cell] <-genoffset_out_tempw_oka_sv
plot(genoffset_tempw_oka_sv,zlim = c(0,0.012), col=rgb.tables(1000))

coords_eco<-env_scale_oka_eco %>% dplyr::filter(dataset=="past") %>% dplyr::select(x, y)
point_okaeco<-sp::SpatialPoints(coords_eco, proj4string = genoffset_tempb_oka_sv@crs)
genomicofsset_values_b_oka_sv<-raster::extract(genoffset_tempb_oka_sv, coords_eco)
genomicofsset_values_w_oka_sv<-raster::extract(genoffset_tempw_oka_sv, coords_eco)

genomicofsset_b_oka_sv<-cbind(coords_eco, genomicofsset_values_b_oka_sv)
genomicofsset_w_oka_sv<-cbind(coords_eco, genomicofsset_values_w_oka_sv)

names(genomicofsset_b_oka_sv)[3] ="gen_off"
names(genomicofsset_w_oka_sv)[3] ="gen_off"

genomicofsset_b_oka_sv<-cbind(env_scale_oka_eco[1:25,1], genomicofsset_b_oka_sv)
names(genomicofsset_b_oka_sv)[1] ="location"
genomicofsset_w_oka_sv<-cbind(env_scale_oka_eco[1:25,1], genomicofsset_w_oka_sv)
names(genomicofsset_w_oka_sv)[1] ="location"

genomicofsset_b_oka_sv
genomicofsset_w_oka_sv

plot(deltarastb[["bio5"]],zlim = c(-2,13), col=rgb.tables(1000))
plot(deltarastw[["bio5"]],zlim = c(-2,13), col=rgb.tables(1000))

###extract deltabioclim fro each pop
deltab_pop<-raster::extract(deltarastb, coords)
deltaw_pop<-raster::extract(deltarastw, coords)
deltab_pop_df<-as.data.frame(deltab_pop)
deltaw_pop_df<-as.data.frame(deltaw_pop)

min(deltab_pop_df$bio1)
max(deltab_pop_df$bio1)
min(deltab_pop_df$bio5)
max(deltab_pop_df$bio5)

min(deltaw_pop_df$bio1)
max(deltaw_pop_df$bio1)


library(ggplot2)
library(scales)
#install.packages("sf")
#install.packages("ggnewscale")
library(ggnewscale)
library(sf)
#install.packages("patchwork")
library(patchwork)
###sequence plots
best_snp<-ggplot() +
  geom_raster(data=deltarastb_df, aes(x=x,y=y, fill=bio5)) + coord_fixed() +
  coord_sf(xlim = c(-143,-111), ylim=c(47,61), expand=FALSE) +
  #scale_fill_viridis_c(option = "D",limits=c(-3,13),na.value = NA,name=paste("\u0394T (\u00b0C)\n")) +
  #scale_fill_gradient2(low = "blue", mid='white',high = "red", midpoint=0) +
  #scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), values=rescale(c(-1,0,1))) +
  scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), limits=c(-13,13),na.value=NA, name=paste("bio5\n\u0394T (\u00b0C)")) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid=element_blank(), plot.background = element_blank(), 
        panel.background = element_blank()) +
  new_scale_fill() +
  geom_point(data=genomicofsset_b, mapping=aes(x,y, fill=gen_off, size=gen_off), shape=21, col="black") +
  scale_fill_viridis_c(option = "H", limits =c(0.007,0.032), guide = "legend", name ="Genomic offset") +
  scale_size_continuous(limits = c(0.007,0.032), name ="Genomic offset", range=c(0.1,5)) +
  ggtitle("Sequence variation - best scenario (RCP2.6)") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

worst_snp<-ggplot() +
  geom_raster(data=deltarastw_df, aes(x=x,y=y, fill=bio5)) + coord_fixed() +
  coord_sf(xlim = c(-143,-111), ylim=c(47,61), expand=FALSE) +
  #scale_fill_viridis_c(option = "D",limits=c(-3,13),na.value = NA,name=paste("\u0394T (\u00b0C)\n")) +
  #scale_fill_gradient2(low = "blue", mid='white',high = "red", midpoint=0) +
  #scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), values=rescale(c(-1,0,1))) +
  scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), limits=c(-13,13),na.value=NA, name=paste("bio5\n\u0394T (\u00b0C)")) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid=element_blank(), plot.background = element_blank(), 
        panel.background = element_blank()) +
  new_scale_fill() +
  geom_point(data=genomicofsset_w, mapping=aes(x,y, fill=gen_off, size=gen_off), shape=21, col="black") +
  scale_fill_viridis_c(option = "H", limits =c(0.007,0.032), guide = "legend", name ="Genomic offset") +
  scale_size_continuous(limits = c(0.007,0.032), name ="Genomic offset", range=c(0.1,5)) +
  ggtitle("Sequence variation - worst scenario (RCP8.5)") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

combined_snp <- best_snp/worst_snp & theme(legend.position = "right")
combined_snp_legend<-combined_snp + plot_layout(guides = "collect")
combined_snp_legend
#### SVs plots
best_sv<-ggplot() +
  geom_raster(data=deltarastb_df, aes(x=x,y=y, fill=bio5)) + coord_fixed() +
  coord_sf(xlim = c(-143,-111), ylim=c(47,61), expand=FALSE) +
  #scale_fill_viridis_c(option = "D",limits=c(-3,13),na.value = NA,name=paste("\u0394T (\u00b0C)\n")) +
  #scale_fill_gradient2(low = "blue", mid='white',high = "red", midpoint=0) +
  #scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), values=rescale(c(-1,0,1))) +
  scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), limits=c(-13,13),na.value=NA, name=paste("bio5\n\u0394T (\u00b0C)")) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid=element_blank(), plot.background = element_blank(), 
        panel.background = element_blank(), text = element_text(size = 15)) +
  new_scale_fill() +
  geom_point(data=genomicofsset_b_sv, mapping=aes(x,y, fill=gen_off, size=gen_off), shape=21, col="black") +
  scale_fill_viridis_c(option = "H", limits =c(0.0017,0.0086), guide = "legend", name ="Genomic offset") +
  scale_size_continuous(limits = c(0.0017,0.0086), name ="Genomic offset", range=c(0.1,5)) +
  ggtitle("Structural variation - best scenario (RCP2.6)") +
  theme(plot.title = element_text(hjust = 0.5))

worst_sv<-ggplot() +
  geom_raster(data=deltarastw_df, aes(x=x,y=y, fill=bio5)) + coord_fixed() +
  coord_sf(xlim = c(-143,-111), ylim=c(47,61), expand=FALSE) +
  #scale_fill_viridis_c(option = "D",limits=c(-3,13),na.value = NA,name=paste("\u0394T (\u00b0C)\n")) +
  #scale_fill_gradient2(low = "blue", mid='white',high = "red", midpoint=0) +
  #scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), values=rescale(c(-1,0,1))) +
  scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), limits=c(-13,13),na.value=NA, name=paste("bio5\n\u0394T (\u00b0C)")) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid=element_blank(), plot.background = element_blank(), 
        panel.background = element_blank()) +
  new_scale_fill() +
  geom_point(data=genomicofsset_w_sv, mapping=aes(x,y, fill=gen_off, size=gen_off), shape=21, col="black") +
  scale_fill_viridis_c(option = "H", limits =c(0.0017,0.0086), guide = "legend", name ="Genomic offset") +
  scale_size_continuous(limits = c(0.0017,0.0086), name ="Genomic offset", range=c(0.1,5)) +
  ggtitle("Structural variation - worst scenario (RCP8.5)") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) 

combined_sv <- best_sv/worst_sv & theme(legend.position = "right")
combined_sv_legend<-combined_sv + plot_layout(guides = "collect")

combined_sv_legend
combined_snp_legend

####correlation with latitude
ggplot(genomicofsset_b, aes(y,gen_off)) + geom_point()
ggplot(genomicofsset_w, aes(y,gen_off)) + geom_point()

ggplot(genomicofsset_b_sv, aes(y,gen_off)) + geom_point()
ggplot(genomicofsset_w_sv, aes(y,gen_off)) + geom_point()

lm_snp_b_lat<-lm(y ~ gen_off, data=genomicofsset_b)
lm_snp_w_lat<-lm(y ~ gen_off, data=genomicofsset_w)
lm_sv_b_lat<-lm(y ~ gen_off, data=genomicofsset_b_sv)
lm_sv_w_lat<-lm(y ~ gen_off, data=genomicofsset_w_sv)

summary(lm_snp_b_lat)
summary(lm_snp_w_lat)
summary(lm_sv_b_lat)
summary(lm_sv_b_lat)
summary(lm_snp_b_lat)

# Call:
#   lm(formula = y ~ gen_off, data = genomicofsset_b)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.2872 -1.6768 -0.6348  1.6909  4.9761 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    59.84       2.78  21.522 2.65e-15 ***
#   gen_off      -432.77     145.88  -2.967  0.00763 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.461 on 20 degrees of freedom
# Multiple R-squared:  0.3056,	Adjusted R-squared:  0.2709 
# F-statistic: 8.801 on 1 and 20 DF,  p-value: 0.007627
# 
# > summary(lm_snp_w_lat)
# 
# Call:
#   lm(formula = y ~ gen_off, data = genomicofsset_w)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.3196 -1.3105 -0.1236  1.5928  3.8165 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   64.988      3.573  18.187 6.57e-14 ***
#   gen_off     -539.813    144.235  -3.743  0.00128 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.265 on 20 degrees of freedom
# Multiple R-squared:  0.4119,	Adjusted R-squared:  0.3825 
# F-statistic: 14.01 on 1 and 20 DF,  p-value: 0.001283
# 
# > summary(lm_sv_b_lat)
# 
# Call:
#   lm(formula = y ~ gen_off, data = genomicofsset_b_sv)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.8600 -1.6918 -0.4077  1.1328  4.6868 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    59.523      2.549  23.352  5.5e-16 ***
#   gen_off     -1351.883    433.417  -3.119  0.00541 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.422 on 20 degrees of freedom
# Multiple R-squared:  0.3273,	Adjusted R-squared:  0.2936 
# F-statistic: 9.729 on 1 and 20 DF,  p-value: 0.005405
# 
# > summary(lm_sv_b_lat)
# 
# Call:
#   lm(formula = y ~ gen_off, data = genomicofsset_b_sv)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.8600 -1.6918 -0.4077  1.1328  4.6868 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    59.523      2.549  23.352  5.5e-16 ***
#   gen_off     -1351.883    433.417  -3.119  0.00541 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.422 on 20 degrees of freedom
# Multiple R-squared:  0.3273,	Adjusted R-squared:  0.2936 
# F-statistic: 9.729 on 1 and 20 DF,  p-value: 0.005405

deltarastw_df
deltat_b<-raster::extract(deltarastb, coords)
deltat_w<-raster::extract(deltarastw, coords)

genomicofsset_b<-cbind(genomicofsset_b,deltat_b)
genomicofsset_w<-cbind(genomicofsset_w,deltat_w)
genomicofsset_b_sv<-cbind(genomicofsset_b_sv,deltat_b)
genomicofsset_w_sv<-cbind(genomicofsset_w_sv,deltat_w)

lm_snp_b<-lm(bio5 ~ gen_off, data=genomicofsset_b)
lm_snp_w<-lm(bio5 ~ gen_off, data=genomicofsset_w)
lm_sv_b<-lm(bio5 ~ gen_off, data=genomicofsset_b_sv)
lm_sv_w<-lm(bio5 ~ gen_off, data=genomicofsset_w_sv)
summary(lm_snp_b)
summary(lm_snp_w)
summary(lm_sv_b)
summary(lm_sv_b)

# > summary(lm_snp_b)
# 
# Call:
#   lm(formula = bio5 ~ gen_off, data = genomicofsset_b)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.9454 -0.4271  0.1757  0.7292  1.2360 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.313      1.042   2.220  0.03813 *  
#   gen_off      251.878     54.660   4.608  0.00017 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9221 on 20 degrees of freedom
# Multiple R-squared:  0.515,	Adjusted R-squared:  0.4907 
# F-statistic: 21.23 on 1 and 20 DF,  p-value: 0.0001701
# 
# > summary(lm_snp_w)
# 
# Call:
#   lm(formula = bio5 ~ gen_off, data = genomicofsset_w)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.6054 -0.4073  0.2532  0.5228  0.8802 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.1823     1.3024   -0.14     0.89    
# gen_off     381.1200    52.5712    7.25 5.16e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8255 on 20 degrees of freedom
# Multiple R-squared:  0.7244,	Adjusted R-squared:  0.7106 
# F-statistic: 52.56 on 1 and 20 DF,  p-value: 5.156e-07
# 
# > summary(lm_sv_b)
# 
# Call:
#   lm(formula = bio5 ~ gen_off, data = genomicofsset_b_sv)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.90599 -0.08382  0.16815  0.52521  0.84566 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.1818     0.8468   2.576    0.018 *  
#   gen_off     841.3543   143.9889   5.843 1.02e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8047 on 20 degrees of freedom
# Multiple R-squared:  0.6306,	Adjusted R-squared:  0.6121 
# F-statistic: 34.14 on 1 and 20 DF,  p-value: 1.023e-05
# 
# > summary(lm_sv_b)
# 
# Call:
#   lm(formula = bio5 ~ gen_off, data = genomicofsset_b_sv)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.90599 -0.08382  0.16815  0.52521  0.84566 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.1818     0.8468   2.576    0.018 *  
#   gen_off     841.3543   143.9889   5.843 1.02e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8047 on 20 degrees of freedom
# Multiple R-squared:  0.6306,	Adjusted R-squared:  0.6121 
# F-statistic: 34.14 on 1 and 20 DF,  p-value: 1.023e-05
# 



ggplot(genomicofsset_b, aes(bio5,gen_off)) + geom_point()
ggplot(genomicofsset_w, aes(bio5,gen_off)) + geom_point()

ggplot(genomicofsset_b_sv, aes(bio5,gen_off)) + geom_point()
ggplot(genomicofsset_w_sv, aes(bio5,gen_off)) + geom_point()

genomicofsset_b$model<-"best"
genomicofsset_b_sv$model<-"best"
genomicofsset_w$model<-"worst"
genomicofsset_w_sv$model<-"worst"

genomicofsset_b$variant<-"snp"
genomicofsset_b_sv$variant<-"sv"
genomicofsset_w$variant<-"snp"
genomicofsset_w_sv$variant<-"sv"


all_genoff_deltaclim<-rbind(genomicofsset_b, genomicofsset_w, genomicofsset_b_sv, genomicofsset_w_sv)
write.table(all_genoff_deltaclim, "kokanee_gea_genomicoffset.txt", sep = "\t")
all_genoff_deltaclim<-read.table("kokanee_gea_genomicoffset.txt", header=TRUE)
###plot relatrionship climate change - genomic offset  
go_deltat<-ggplot(all_genoff_deltaclim,mapping=aes(bio5,gen_off,colour=model, shape=variant)) +
    geom_point(size=2) + geom_smooth(aes(linetype = variant),method="lm",se=FALSE) +
    ylab("Genomic offset") + xlab(paste("bio5 - \u0394T (\u00b0C)")) +
    theme_bw() +
    scale_shape_manual(name = "Variant", 
                       values=c(snp = 16, sv = 17), 
                       labels=c(snp = "SNPs", sv = "SVs"), 
                       guide = guide_legend(order = 1))+
    scale_linetype_manual(name = "Variant", 
                          values=c(snp = "solid", sv = "dotted"), 
                          labels=c(snp = "SNPs", sv = "SVs"), 
                          guide = guide_legend(order = 1)) +
    scale_color_manual(name = "Climate scenario", 
                       labels=c(best = "Best (RCP2.6)", worst = "Worst (RCP8.5)"), 
                      values=c(best = "#219ebc", worst = "#fb8500"),
                       guide = guide_legend(order = 2)) 
 
go_lat<-ggplot(all_genoff_deltaclim,mapping=aes(y,gen_off,colour=model, shape=variant)) +
  geom_point(size=2) + geom_smooth(aes(linetype = variant),method="lm",se=FALSE) +
  ylab("Genomic offset") + xlab(paste("Latitude")) +
  theme_bw() +
  scale_shape_manual(name = "Variant", 
                     values=c(snp = 16, sv = 17), 
                     labels=c(snp = "SNPs", sv = "SVs"), 
                     guide = guide_legend(order = 1))+
  scale_linetype_manual(name = "Variant", 
                        values=c(snp = "solid", sv = "dotted"), 
                        labels=c(snp = "SNPs", sv = "SVs"), 
                        guide = guide_legend(order = 1)) +
  scale_color_manual(name = "Climate scenario", 
                     labels=c(best = "Best (RCP2.6)", worst = "Worst (RCP8.5)"), 
                     values=c(best = "#219ebc", worst = "#fb8500"),
                     guide = guide_legend(order = 2)) 

combined_go_corr <- go_lat+go_deltat & theme(legend.position = "right")
combined_go_corr_legend<-combined_go_corr + plot_layout(guides = "collect")
combined_go_corr_legend


#### add relationship gen offset and diversity
#install.packages("MuMIn")
library(data.table)
library(tidyverse)
library(lme4)
library(predictmeans)
library(lmerTest)
library(MuMIn)
library(data.table)
#het <-fread("env_variable_sel_scale_pop_coords_het.txt", header= TRUE)

het<-fread("env_variable_all_pop_coords_het.txt", header = TRUE)
het_svs<-fread("env_variable_all_pop_coords_het_svs.txt", header = TRUE)
# Load dplyr
library(dplyr)

# Group by mean using dplyr
mean_het <- het %>% group_by(pop) %>% 
  summarise(mean_het=mean(het_prop),
            .groups = 'drop')
mean_het

# Convert tibble to df
mean_het_df <- mean_het %>% as.data.frame()
mean_het_df

# Group by mean using dplyr
mean_het_svs <- het_svs %>% group_by(pop) %>% 
  summarise(mean_het_svs=mean(het_prop),
            .groups = 'drop')
mean_het_svs

# Convert tibble to df
mean_het_svs_df <- mean_het_svs %>% as.data.frame()
mean_het_svs_df
data.table::setnames(mean_het_svs,'mean_het_svs', 'mean_het')

mean_het_df4<-rbind(mean_het_df,mean_het_df,mean_het_svs_df,mean_het_svs_df)

all_genoff_deltaclim_het<-cbind(all_genoff_deltaclim,mean_het_df4)

go_het<-ggplot(all_genoff_deltaclim_het,mapping=aes(mean_het,gen_off,colour=model, shape=variant)) +
  geom_point(size=2) + #geom_smooth(aes(linetype = variant),method="lm",se=FALSE) +
  ylab("Genomic offset") + xlab(paste("Mean proportion of heterozygous sites")) +
  theme_bw() +
  scale_shape_manual(name = "Variant", 
                     values=c(snp = 16, sv = 17), 
                     labels=c(snp = "SNPs", sv = "SVs"), 
                     guide = guide_legend(order = 1))+
  scale_linetype_manual(name = "Variant",
                        values=c(snp = "solid", sv = "dotted"),
                        labels=c(snp = "SNPs", sv = "SVs"),
                        guide = guide_legend(order = 1)) +
  scale_color_manual(name = "Climate scenario", 
                     labels=c(best = "Best (RCP2.6)", worst = "Worst (RCP8.5)"), 
                     values=c(best = "#219ebc", worst = "#fb8500"),
                     guide = guide_legend(order = 2)) 

combined_go_corr <- go_lat+go_deltat + go_het & theme(legend.position = "right")
combined_go_corr_legend<-combined_go_corr + plot_layout(guides = "collect")
combined_go_corr_legend




genoff<-fread("kokanee_gea_genomicoffset.txt", header=TRUE)


###
plot(genomicofsset_b$gen_off, genomicofsset_b_sv$gen_off)
plot(genomicofsset_w$gen_off, genomicofsset_w_sv$gen_off)

go_b_lm<-lm(genomicofsset_b$gen_off~ genomicofsset_b_sv$gen_off)
go_w_lm<-lm(genomicofsset_w$gen_off~ genomicofsset_w_sv$gen_off)

summary(go_b_lm)
# Call:
#   lm(formula = genomicofsset_b$gen_off ~ genomicofsset_b_sv$gen_off)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0013784 -0.0004953 -0.0001163  0.0004402  0.0016012 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.0017600  0.0008749   2.012   0.0579 .  
# genomicofsset_b_sv$gen_off 2.9443229  0.1487681  19.791 1.32e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0008315 on 20 degrees of freedom
# Multiple R-squared:  0.9514,	Adjusted R-squared:  0.949 
# F-statistic: 391.7 on 1 and 20 DF,  p-value: 1.319e-14

summary(go_w_lm)
# Call:
#   lm(formula = genomicofsset_w$gen_off ~ genomicofsset_w_sv$gen_off)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -1.831e-03 -4.679e-04  2.452e-05  4.621e-04  1.714e-03 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.001852   0.001339   1.384    0.182    
# genomicofsset_w_sv$gen_off 3.002568   0.175354  17.123 2.05e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0008873 on 20 degrees of freedom
# Multiple R-squared:  0.9361,	Adjusted R-squared:  0.9329 
# F-statistic: 293.2 on 1 and 20 DF,  p-value: 2.047e-13


install.packages("RColorBrewer")
library(RColorBrewer)
climRasts_df<-as.data.frame(climRasts, xy=TRUE)


library(ggrepel)
library(ggplot2)
library(ggnewscale)
all_genoff_deltaclim<-read.table("kokanee_gea_genomicoffset.txt", header=TRUE)
location_coord<-all_genoff_deltaclim[1:22,2:4]
climRasts_df<-as.data.frame(climRasts, xy=TRUE)
location_coord$location = factor(location_coord$location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Bobtail", "Arctic", "Puntzi", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan"))

new_bg<-c("#001219", "#5a189a", "#B785FF", "#00509d", "#4A9FFF", "#AAD8FF", "#005E29", "#00C05A", "#83D477", "#b9e769", "#e9d8a6", "#efea5a", "#f1c453", "#f29e4c", "#f3722c", "#f94144", "#d00000", "#9d0208", "#6a040f", "#990066", "#d94a8c", "#fad2e1")


ggplot() +
  geom_raster(data=climRasts_df, aes(x=x,y=y, fill=bio5), alpha=0.5) + coord_fixed() +
  coord_sf(xlim = c(-143,-111), ylim=c(47,61), expand=FALSE) +
  #scale_fill_viridis_c(option = "D",limits=c(-3,13),na.value = NA,name=paste("\u0394T (\u00b0C)\n")) +
  #scale_fill_gradient2(low = "blue", mid='white',high = "red", midpoint=0) +
  #scale_fill_gradientn(colours=c("blue","cyan","white", "yellow","red"), values=rescale(c(-1,0,1))) +
  scale_fill_gradientn(colours=c("blue","cyan", "yellow","red"), limits=c(0,31),na.value=NA, name=paste("bio5 (\u00b0C)")) +
  #scale_fill_viridis_c(option="heat",na.value=NA, name=paste("bio5 (\u00b0C)"), direction=-1)+
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid=element_blank(), plot.background = element_blank(), 
        panel.background = element_blank(), legend.position = "none", text = element_text(size = 15)) +
  borders("world", size=0.1) +
  new_scale_fill() +
  geom_point(data=location_coord, mapping=aes(x,y, fill=location), shape=21, size=3, col="black") + 
  geom_text_repel(data=location_coord, mapping=aes(x,y,label=location), max.overlaps = Inf, size=5) + 
  scale_fill_manual(values=new_bg) 




location_coord$location = factor(location_coord$location, levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Bobtail", "Arctic", "Puntzi", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan"))

new_bg<-c("#001219", "#3c096c", "#B785FF", "#013a63", "#4A9FFF", "#AAD8FF", "#005E29", "#00C05A", "#83D477", "#b9e769", "#e9d8a6", "#efea5a", "#f1c453", "#f29e4c", "#f3722c", "#f94144", "#d00000", "#9d0208", "#6a040f", "#990066", "#d94a8c", "#fad2e1")
levels=c("Sockeye", "Thutade", "Tchesinkut", "Cluculz", "Arctic", "Puntzi", "Bobtail", "La_Hache", "Bonaparte", "Dunn", "EastBarriere", "Anderson", "Nicola", "Okanagan", "Wood", "Kalamalka", "Christina", "Arrow_Hill", "Arrow_Mosq", "Kootenay", "Cowichan", "Shawningan")))%
+
  new_scale_fill() +
   +
  scale_fill_viridis_c(option = "H", limits =c(0.0018,0.01), guide = "legend", name ="Genomic offset") +
  scale_size_continuous(limits = c(0.0018,0.01), name ="Genomic offset", range=c(0.1,5)) +
  ggtitle("Structural variation - worst scenario(RCP8.5)") +
  theme(plot.title = element_text(hjust = 0.5)) 

geno_outliers_oka_new<-as.data.frame(t(geno_outliers_oka[-1,-1]))
geno_outliers_oka[1:10,1:10]
rownames(geno_outliers_oka)<-NULL
