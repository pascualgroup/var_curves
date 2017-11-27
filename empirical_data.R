# Define functions --------------------------------------------------------
source('mtn_functions.R')

cleanSurveyData <- function(df, min.var.per.isolate=40, max.var.per.isolate=60, plotit=F){
  require(ggplot2)
  require(cowplot)
  df <- bipartite::empty(df)
  if(plotit){
    p1 <- qplot(rowSums(df), main = 'Number of times a var type occurs in the host population')
    p2 <- qplot(colSums(df), binwidth=min.var.per.isolate, main='number of var types in an isolate')+geom_vline(xintercept = min.var.per.isolate, color='red')
    p <- plot_grid(p1, p2,
                   labels = c("A", "B"),
                   ncol = 1, nrow = 2)
    print(p)
  }
  df <- df[, which(colSums(df)>=min.var.per.isolate & colSums(df)<=max.var.per.isolate)]  # Remove isolates by constrains on number of var genes
  df <- bipartite::empty(df)
  # df <- df[which(rowSums(df)>1),] # Remove var types that appear only once
  # df <- bipartite::empty(df)
  cat(nrow(df),' var types in ',ncol(df),' isolates')
  return(data.matrix(df))
}

# Initialize --------------------------------------------------------------

prep.packages(c('sqldf')) # Packages for SQL
prep.packages(c("stringr","reshape2","splitstackshape","plyr","dplyr")) # Packages for data manipulation
prep.packages(c("ggplot2","RColorBrewer","grid","gtable","gridExtra","gplots","cowplot", "ggdendro","ggpubr")) # Packages for plotting
prep.packages(c("igraph","bipartite")) # Packages for network analysis
prep.packages(c("entropy","vegan")) # Packages for analysis


# Get and clean empirical data ---------------------------------------------------

mainData <- read.csv('THESIS_FINAL_S1_S2_45487_DBLa_types_1284_isolates.csv', header = T)
names(mainData)[1] <- 'OTU_ID'
anyDuplicated(mainData$OTU_ID)
anyDuplicated(colnames(mainData))
varGroupings <- read.csv('DBLa_types_and_Ups_45487.csv', header = T)
names(varGroupings)[1] <- 'OTU_ID'
mainData <- subset(mainData, OTU_ID%in%subset(varGroupings, Ups=='BC')$OTU_ID)
rownames(mainData) <- mainData[,1]
mainData <- mainData[,-1]

### Separate by surveys
isolates_all <- colnames(mainData)
table(str_sub(isolates_all,1,4))
S1_isolates_idx <- which(str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='S1' | 
                           str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='SBS1' |
                           str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='RS1' | # RS1 and R2S1 are isolates for which infectino was not detect with the first PCR but with subsequent PCRs
                           str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='R2S1')
S2_isolates_idx <- which(str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='S2' |
                           str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='SBS2')
intersect(S1_isolates_idx,S2_isolates_idx) # Any overlap in isolates?
length(isolates_all) - (length(S1_isolates_idx)+length(S2_isolates_idx)) # Are all isolates included?


S1_isolates <- str_sub(isolates_all[S1_isolates_idx],str_locate(isolates_all[S1_isolates_idx], 'MRS')[,1],str_locate(isolates_all[S1_isolates_idx], '\\.')[,1]-1)
anyDuplicated(S1_isolates) # make sure that hosts with RS1 or R2S1 are not in the S1 set.
S2_isolates <- str_sub(isolates_all[S2_isolates_idx],str_locate(isolates_all[S2_isolates_idx], 'MRS')[,1],str_locate(isolates_all[S2_isolates_idx], '\\.')[,1]-1)
anyDuplicated(S2_isolates)
Data_S1 <- mainData[,S1_isolates_idx]
Data_S2 <- mainData[,S2_isolates_idx]

# Limit to 55 vars (MOI=1)
maxVarIsolate <- 55 # 55 because we do not take grouping A
Data_S1 <- cleanSurveyData(Data_S1, min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S2 <- cleanSurveyData(Data_S2, min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)

mainDataClean <- mainData[unique(c(rownames(Data_S1),rownames(Data_S2))),
                          unique(c(colnames(Data_S1),colnames(Data_S2)))]
mainDataClean <- bipartite::empty(mainDataClean)
dim(mainDataClean)



# EMPIRICAL: Define Layers -----------------------------------------------------------
# # Remove the single individual with a chronic infection
# Data_S1 <- Data_S1[,-grep(pattern = 'MRS1011', colnames(Data_S1))]
# Data_S2 <- Data_S2[,-grep(pattern = 'MRS1011', colnames(Data_S2))]

empiricalLayer_1 <- overlapAlleleAdj(t(Data_S1))
empiricalLayer_2 <- overlapAlleleAdj(t(Data_S2))
diag(empiricalLayer_1) <- 0
diag(empiricalLayer_2) <- 0

layersEmpirical <- list(empiricalLayer_1, empiricalLayer_2)
similarityMatrixEmpirical <- overlapAlleleAdj(t(mainDataClean))

sapply(layersEmpirical, gdensity)

plotSurveyLayer(layersEmpirical[[1]], main='S1, Oct 2012', vertex.color='navy')
plotSurveyLayer(layersEmpirical[[2]], main='S2, Jun 2013', vertex.color='navy')


# EMPIRICAL: Define cutoff -----------------------------------------------------------
write.csv(similarityMatrixEmpirical, 'similarityMatrixEmpirical.csv') # save this so to be able to plot edge weight distributons fast

cutoffProbability <- 0.95

# links <- as.vector(similarityMatrixEmpirical[similarityMatrixEmpirical!=0])
# length(links)
# number_links_to_preserve <- ceiling(length(links)*0.05)

# Explore the similarity between isoaltes
cutoffValue <- quantile(as.vector(similarityMatrixEmpirical[similarityMatrixEmpirical!=0]), probs = cutoffProbability)
p_edge_weight_histogram_empirical <- ggHistogram(as.vector(similarityMatrixEmpirical[similarityMatrixEmpirical!=0]), xlab = 'Non-zero similarity values')+geom_vline(xintercept = cutoffValue, color='red')
similarityMatrixEmpirical_cutoff <- similarityMatrixEmpirical
similarityMatrixEmpirical_cutoff[similarityMatrixEmpirical_cutoff<cutoffValue] <- 0
layersEmpirical_cutoff <- layersEmpirical
for (i in 1:length(layersEmpirical_cutoff)){
  x <-layersEmpirical_cutoff[[i]]
  x[x<cutoffValue] <- 0
  layersEmpirical_cutoff[[i]] <- x
}
# plotSurveyLayer(empiricalLayer_1)

sapply(layersEmpirical_cutoff, gdensity)

# Plot layers 
pdf('~/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/Results/empirical_layers.pdf', width = 8, height = 8)
plotSurveyLayer(layersEmpirical_cutoff[[1]], main='S1, Oct 2012', vertex.color='navy')
plotSurveyLayer(layersEmpirical_cutoff[[2]], main='S2, Jun 2013', vertex.color='navy')
dev.off()


# ABM: Get data from SQL -------------------------------------------------------
iteration=1
args=c(paste('mtn_12_seasonality_empirical_irs_2_run_',iteration,sep=''),27000,30,0.9,1)
experiment <- as.character(args[1])
burnin <- as.numeric(args[2])
windowWidth <- as.numeric(args[3]) # 30 is for a month
cutoffPercentile <- as.numeric(args[4])
MOI <- as.numeric(args[5])

setwd('~/Documents/sqlite/')
# maxTime <- firstSample+windowWidth*nSamples*(aggregateResolution+samplingSpace)
maxTime <- 36000
db <- dbConnect(SQLite(), dbname = paste(experiment,'.sqlite',sep=''))
strainComposition <- dbGetQuery(db, 'SELECT strainId, geneId FROM strains')
loci <- dbGetQuery(db, 'SELECT alleleId, geneId FROM loci')

# Get host infection data from sqlite.
temporalData <- dbGetQuery(db, paste('SELECT time, hostId, strainId, active FROM sampledHostInfections WHERE time <= ',maxTime,sep=''))
temporalData <- temporalData[temporalData$time>=burnin,]
temporalData <- temporalData[!duplicated(temporalData),] #Remove duplicated infections (Ask Qixin why there are duplicated infections)
temporalData <- temporalData[!is.na(temporalData$active),] # Only include strains active in the host blood
temporalData$timeSlice <- ceiling((temporalData$time-burnin)/windowWidth) # subtract the burnin time (18,000 days)
range(temporalData$time)

# Set months
calendar <- build_calendar(num_years = 100, burnin = 27000, year_to_start = 100)
calendar$empirical_survey[calendar$month_sim=='Oct'] <- 'S1'
calendar$empirical_survey[calendar$month_sim=='Jun'] <- 'S2'

# Set the dates that correspond to the survey
temp <- subset(calendar, !is.na(empirical_survey), select=c('running_day','empirical_survey'))
names(temp)[1] <- 'time'
temporalData <- merge(temporalData, temp, all.x=T)
temporalData <- subset(temporalData, !is.na(empirical_survey))
temporalData <- with(temporalData, temporalData[order(timeSlice,hostId,strainId),])


# ABM: Build networks -------------------------------------------------

## For S1
temporalData_S1 <- subset(temporalData, empirical_survey=='S1')
timeSlices_S1 <- unique(temporalData_S1$timeSlice)
hostMOIdf <- c()
for(ts in timeSlices_S1){
  hostMOI <- binarize(xtabs(~hostId+strainId, data=subset(temporalData_S1, timeSlice==ts)))
  hostMOIdf <- rbind(hostMOIdf, data.frame(timeSlice=ts, hostId=rownames(hostMOI), MOI=rowSums(hostMOI)))
}
temporalData_S1 <- merge(temporalData_S1,hostMOIdf, by=c('hostId','timeSlice'), all.x = T)
temporalData_S1 <- subset(temporalData_S1, MOI==MOI)
# Account only for the strains which occur in the surveys
strainComposition_reduced_S1 <- subset(strainComposition, strainId%in%temporalData_S1$strainId)
strainComposition_reduced_S1 <- merge(strainComposition_reduced_S1,loci,by='geneId')
layersFull_S1 <- vector(mode = 'list', length = length(timeSlices_S1))
for (ts in timeSlices_S1){
  print(paste('[',Sys.time(), '] Building layer for S1, time slice ',ts,'...',sep=''))
  strains <- with(temporalData_S1, subset(strainId, timeSlice==ts)) # get the data for the particular slice
  strainCompositionLayer <- strainComposition_reduced_S1[strainComposition_reduced_S1$strainId%in%strains,]
  similarityMatrix <- overlapAlleleAdj(table(strainCompositionLayer$strainId, strainCompositionLayer$geneId)) # Notice that we should use geneId, like in the empirical data
  diag(similarityMatrix) <- 0
  layersFull_S1[[which(timeSlices_S1==ts)]] <- similarityMatrix
}

## For S2
temporalData_S2 <- subset(temporalData, empirical_survey=='S2')
timeSlices_S2 <- unique(temporalData_S2$timeSlice)
hostMOIdf <- c()
for(ts in timeSlices_S2){
  hostMOI <- binarize(xtabs(~hostId+strainId, data=subset(temporalData_S2, timeSlice==ts)))
  hostMOIdf <- rbind(hostMOIdf, data.frame(timeSlice=ts, hostId=rownames(hostMOI), MOI=rowSums(hostMOI)))
}
temporalData_S2 <- merge(temporalData_S2,hostMOIdf, by=c('hostId','timeSlice'), all.x = T)
temporalData_S2 <- subset(temporalData_S2, MOI==MOI)
# Account only for the strains which occur in the surveys
strainComposition_reduced_S2 <- subset(strainComposition, strainId%in%temporalData_S2$strainId)
strainComposition_reduced_S2 <- merge(strainComposition_reduced_S2,loci,by='geneId')
layersFull_S2 <- vector(mode = 'list', length = length(timeSlices_S2))
for (ts in timeSlices_S2){
  print(paste('[',Sys.time(), '] Building layer for S2, time slice ',ts,'...',sep=''))
  strains <- with(temporalData_S2, subset(strainId, timeSlice==ts)) # get the data for the particular slice
  strainCompositionLayer <- strainComposition_reduced_S2[strainComposition_reduced_S2$strainId%in%strains,]
  similarityMatrix <- overlapAlleleAdj(table(strainCompositionLayer$strainId, strainCompositionLayer$geneId)) # Notice that we should use geneId, like in the empirical data
  diag(similarityMatrix) <- 0
  layersFull_S2[[which(timeSlices_S2==ts)]] <- similarityMatrix
}

### Sub-sample networks to obtain the same number of nodes as in empirical networks
layersFull_S1_sampled <- c()
for (i in 1:length(layersFull_S1)){
  x <- layersFull_S1[[i]]
  sm <- sample(nrow(x),nrow(layersEmpirical[[1]]),F)
  layersFull_S1_sampled[[i]] <- x[sm,sm]
}  
sapply(layersFull_S1_sampled, dim)
layersFull_S2_sampled <- c()
for (i in 1:length(layersFull_S2)){
  x <- layersFull_S2[[i]]
  sm <- sample(nrow(x),nrow(layersEmpirical[[2]]),F)
  layersFull_S2_sampled[[i]] <- x[sm,sm]
}  
sapply(layersFull_S2_sampled, dim)

# ABM: Apply cutoff -------------------------------------------------------
x <- xtabs(~strainId+alleleId, strainComposition_reduced_S1)  
similarityMatrix_S1 <- overlapAlleleAdj(x)
cutoffValue_S1 <- quantile(as.vector(similarityMatrix_S1[similarityMatrix_S1!=0]), probs = cutoffProbability)
# ggHistogram(as.vector(similarityMatrix_S1[similarityMatrix_S1!=0]), xlab = 'Non-zero similarity values')+geom_vline(xintercept = cutoffValue, color='red')
similarityMatrix_S1_cutoff <- similarityMatrix_S1
similarityMatrix_S1_cutoff[similarityMatrix_S1_cutoff<cutoffValue] <- 0
layersFull_S1_cutoff <- layersFull_S1
for (i in 1:length(layersFull_S1_cutoff)){
  x <-layersFull_S1_cutoff[[i]]
  x[x<cutoffValue] <- 0
  layersFull_S1_cutoff[[i]] <- x
}
### S2
x <- xtabs(~strainId+alleleId, strainComposition_reduced_S2)  
similarityMatrix_S2 <- overlapAlleleAdj(x)
cutoffValue_S2 <- quantile(as.vector(similarityMatrix_S2[similarityMatrix_S2!=0]), probs = cutoffProbability)
# ggHistogram(as.vector(similarityMatrix_S2[similarityMatrix_S2!=0]), xlab = 'Non-zero similarity values')+geom_vline(xintercept = cutoffValue, color='red')
similarityMatrix_S2_cutoff <- similarityMatrix_S2
similarityMatrix_S2_cutoff[similarityMatrix_S2_cutoff<cutoffValue] <- 0
layersFull_S2_cutoff <- layersFull_S2
for (i in 1:length(layersFull_S2_cutoff)){
  x <-layersFull_S2_cutoff[[i]]
  x[x<cutoffValue] <- 0
  layersFull_S2_cutoff[[i]] <- x
}
# plotSurveyLayer(empiricalLayer_1)

sapply(layersFull_S1_cutoff, gdensity)
sapply(layersFull_S2_cutoff, gdensity)
sapply(layersEmpirical_cutoff, gdensity)


# ABM: Compare seasonal networks --------------------------

features_ABM_S1 <- sapply(layersFull_S1, calculateFeatures)
features_ABM_S1 <- scale(features_ABM_S1)
colnames(features_ABM_S1) <- paste('S1',timeSlices_S1,sep='_')
features_ABM_S2 <- sapply(layersFull_S2, calculateFeatures)
features_ABM_S2 <- scale(features_ABM_S2)
colnames(features_ABM_S2) <- paste('S2',timeSlices_S2,sep='_')
hc_ABM <- hclust(dist(t(cbind(features_ABM_S1,features_ABM_S2))))
plot(hc_ABM)


# Compare ABM to EMPIRICAL ------------------------------------------------
library(adegenet)

S1_networks <- layersFull_S1_sampled
S2_networks <- layersFull_S2_sampled
s_empirical <- layersEmpirical

network_features_to_include <- 1:32

features_ABM_S1 <- sapply(S1_networks, calculateFeatures)
colnames(features_ABM_S1) <- paste('S1',timeSlices_S1,sep='_')
features_ABM_S2 <- sapply(S2_networks, calculateFeatures)
colnames(features_ABM_S2) <- paste('S2',timeSlices_S2,sep='_')
x <- t(cbind(features_ABM_S1[network_features_to_include,],features_ABM_S2[network_features_to_include,]))
groups_ABM <- c(rep('S1',25),rep('S2',25))
grp <- find.clusters(x, max.n.clust=2, n.pca=10, n.clust=2)
table(groups_ABM,grp$grp)
dapc.ABM <- dapc(x, grp$grp, scale=T, center=T, n.pca=10, n.da=2)
scatter(dapc.ABM)
# Try to predict where the data is clustered
features_empirical <- sapply(s_empirical, calculateFeatures)
colnames(features_empirical) <- c('S1_E','S2_E')
pred <- predict.dapc(dapc.ABM, t(features_empirical[network_features_to_include,]), scale=T, center=T)
assignplot(dapc.ABM, new.pred = pred)



# Shuffle empirical networks ----------------------------------------------

# This procudeure allocates isolates uniformly at random to the different
# surveys. It maintains the number of isolates per survey.
allSurveyData <- list(Data_S1,Data_S2)
sapply(allSurveyData, nrow) # Number of vars per survey
sapply(allSurveyData, ncol) # Number of isolates per survey
allIsolates <- unique(unlist(sapply(allSurveyData, colnames)))

shuffled_networks_S1 <- shuffled_networks_S2 <- list()
for(r in 1:100){
  isolateSampleSize <- sapply(allSurveyData, ncol)
  isolateSet <- allIsolates
  shuffledIsolates_S1 <- sample(isolateSet, isolateSampleSize[1], replace = F)
  isolateSet <- isolateSet[!isolateSet %in% shuffledIsolates_S1]
  shuffledIsolates_S2 <- isolateSet
  
  Data_S1_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S1]
  Data_S1_shuffled <- bipartite::empty(Data_S1_shuffled)
  dim(Data_S1_shuffled)
  Data_S2_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S2]
  Data_S2_shuffled <- bipartite::empty(Data_S2_shuffled)
  dim(Data_S2_shuffled)
  
  allSurveyData_shuffled <- list(Data_S1_shuffled,Data_S2_shuffled)
  length(unique(unlist(sapply(allSurveyData_shuffled, rownames)))) # Same number of vars as in the observed data?
  
  vars <- unique(unlist(sapply(allSurveyData_shuffled,rownames)))
  isolates <- unique(unlist(sapply(allSurveyData_shuffled,colnames)))
  superMatrix <- matrix(0, nrow=length(vars), ncol=length(isolates), dimnames = list(vars, isolates))
  superMatrix[rownames(Data_S1_shuffled),colnames(Data_S1_shuffled)] <- data.matrix(Data_S1_shuffled)
  superMatrix[rownames(Data_S2_shuffled),colnames(Data_S2_shuffled)] <- data.matrix(Data_S2_shuffled)
  
  similarityMatrix_shuffled <- overlapAlleleAdj(t(superMatrix))
  identical(similarityMatrix_shuffled,similarityMatrixEmpirical)
  
  # Create the layers
  empiricalLayer_1_shuffled <- overlapAlleleAdj(t(Data_S1_shuffled))
  empiricalLayer_2_shuffled <- overlapAlleleAdj(t(Data_S2_shuffled))
  diag(empiricalLayer_1_shuffled) <- 0
  diag(empiricalLayer_2_shuffled) <- 0

  # # If applying cutoff
  # cutoffValue <- quantile(as.vector(similarityMatrixEmpirical[similarityMatrixEmpirical!=0]), probs = cutoffProbability)
  # similarityMatrix_shuffled[similarityMatrix_shuffled<cutoffValue] <- 0
  # for (i in 1:length(layersEmpirical_shuffled)){
  #   x <-layersEmpirical_shuffled[[i]]
  #   x[x<cutoffValue] <- 0
  #   layersEmpirical_shuffled[[i]] <- x
  # }
  
  shuffled_networks_S1[[r]] <- empiricalLayer_1_shuffled
  shuffled_networks_S2[[r]] <- empiricalLayer_2_shuffled
}
sapply(shuffled_networks_S1,nrow)
sapply(shuffled_networks_S2,nrow)
mean(sapply(shuffled_networks_S1,gdensity))
mean(sapply(shuffled_networks_S2,gdensity))
sapply(layersEmpirical,gdensity)

# Try to predict where the data is clustered
features_shuffled_S1 <- sapply(shuffled_networks_S1, calculateFeatures)
features_shuffled_S2 <- sapply(shuffled_networks_S2, calculateFeatures)
dim(features_shuffled_S1)
dim(features_shuffled_S2)

predictions <- c()
for (r in 1:100){
  joint <- cbind(features_shuffled_S1[,r],features_shuffled_S2[,r])
  colnames(joint) <- c('S1_shuff','S2_shuff')
  pred <- predict.dapc(dapc.ABM, t(joint[network_features_to_include,]), scale=T, center=T)
  predictions <- rbind(predictions, pred$assign)
}


assignplot(dapc.ABM, new.pred = pred)


features_shuffled_S1 <- sapply(shuffled_networks_S1, calculateFeatures)
colnames(features_shuffled_S1) <- paste('S1',1:ncol(features_shuffled_S1),sep='_')
features_shuffled_S2 <- sapply(shuffled_networks_S2, calculateFeatures)
colnames(features_shuffled_S2) <- paste('S2',1:ncol(features_shuffled_S1),sep='_')
x <- t(cbind(features_shuffled_S1[network_features_to_include,],features_shuffled_S2[network_features_to_include,]))
groups <- c(rep('S1',ncol(features_shuffled_S1)),rep('S2',ncol(features_shuffled_S1)))
grp <- find.clusters(x, max.n.clust=25, n.pca=10, n.clust=2)
table(groups, grp$grp)
dapc.shuff <- dapc(x, grp$grp, scale=T, center=T, n.pca=10, n.da=2)
scatter(dapc.shuff)
# Try to predict where the data is clustered
features_empirical <- sapply(s_empirical, calculateFeatures)
colnames(features_empirical) <- c('S1_E','S2_E')
pred <- predict.dapc(dapc.shuff, t(features_empirical[network_features_to_include,]), scale=T, center=T)
assignplot(dapc.shuff, new.pred = pred)



# Run across cutoffs ------------------------------------------------------
cutOffRange <- seq(0,0.95,by=0.05)
# Explore the similarity between isoaltes
x <- as.vector(similarityMatrixEmpirical)
cutoffValues <- quantile(as.vector(similarityMatrixEmpirical[similarityMatrixEmpirical!=0]), probs = cutOffRange)
ggHistogram(as.vector(similarityMatrixEmpirical[similarityMatrixEmpirical!=0]))+geom_vline(xintercept = cutoffValues, color='red')
plotList <- list()
for (cutoffValue in cutoffValues){
  print(cutoffValue)
  filename <- paste('Infomap_survey_empirical_MOI1_cutoff_',cutoffValue,sep='')
  similarityMatrixEmpirical_cutoff <- similarityMatrixEmpirical
  layersEmpirical_cutoff <- layersEmpirical
  similarityMatrixEmpirical_cutoff[similarityMatrixEmpirical_cutoff<cutoffValue] <- 0
  for (i in 1:length(layersEmpirical_cutoff)){
    x <-layersEmpirical_cutoff[[i]]
    x[x<cutoffValue] <- 0
    diag(x) <- 0
    layersEmpirical_cutoff[[i]] <- x
  }
  build_network_Infomap(paste(filename,'.txt',sep=''), layersEmpirical_cutoff, similarityMatrixEmpirical_cutoff)
  system(paste('./Infomap ',filename,'.txt . -2 -i multiplex --multiplex-relax-rate -1 -d -N 10 --rawdir --tree --expanded --silent',sep=''))
  modules_empirical <- infomap_readTreeFile(paste(filename,'_expanded.tree',sep=''), reorganize_modules = T, remove_buggy_instances = F,max_layers = numLayers)
  modules_empirical <- subset(modules_empirical, layer<=4)
  moduleSummary_empirical <- modules_empirical %>% group_by(layer,module) %>% summarise(numStrains=length(strain))
  moduleSummary_empirical$module <- factor(moduleSummary_empirical$module)
  plotList[[which(cutoffValues==cutoffValue)]] <- ggplot(moduleSummary_empirical, aes(x=layer, y=module, size=numStrains, color=as.factor(module)))+geom_point()+
  theme(legend.position="none")+labs(title=paste('MOI=1, cutoff=',cutOffRange[which(cutoffValues==cutoffValue)], ', mean(d)=',round(mean(sapply(layersEmpirical_cutoff,gdensity,F,T)),2),sep=''))
}


  
plotList[[20]] <- ggHistogram(as.vector(similarityMatrixEmpirical), bw = 0.02)+
  geom_vline(xintercept = cutoffValues, color='red')

cowplot::plot_grid(plotlist=plotList, ncol = 5, nrow=4)



# Synthetic random networks -----------------------------------------------

# This procedure produces random isolate-by-var networks, while maintaining only
# the size and density of the networks
Data_S1_shuffled <- matrix(rbinom(nrow(Data_S1)*ncol(Data_S1),1,gdensity(Data_S1, bipartite = T)),nrow(Data_S1),ncol(Data_S1))
dimnames(Data_S1_shuffled) <- dimnames(Data_S1)
Data_S2_shuffled <- matrix(rbinom(nrow(Data_S2)*ncol(Data_S2),1,gdensity(Data_S2, bipartite = T)),nrow(Data_S2),ncol(Data_S2))
dimnames(Data_S2_shuffled) <- dimnames(Data_S2)
Data_S3_shuffled <- matrix(rbinom(nrow(Data_S3)*ncol(Data_S3),1,gdensity(Data_S3, bipartite = T)),nrow(Data_S3),ncol(Data_S3))
dimnames(Data_S3_shuffled) <- dimnames(Data_S3)
Data_S4_shuffled <- matrix(rbinom(nrow(Data_S4)*ncol(Data_S4),1,gdensity(Data_S4, bipartite = T)),nrow(Data_S4),ncol(Data_S4))
dimnames(Data_S4_shuffled) <- dimnames(Data_S4)
allSurveyData_shuffled <- list(Data_S1_shuffled,Data_S2_shuffled,Data_S3_shuffled,Data_S4_shuffled)
length(unique(unlist(sapply(allSurveyData_shuffled, rownames)))) # Same number of vars as in the observed data?


# This procedure maintains frequency of vars (row sums)
Data_S1_shuffled <- vegan::nullmodel(Data_S1, method = 'c0')$data
Data_S2_shuffled <- vegan::nullmodel(Data_S2, method = 'c0')$data
Data_S3_shuffled <- vegan::nullmodel(Data_S3, method = 'c0')$data
Data_S4_shuffled <- vegan::nullmodel(Data_S4, method = 'c0')$data
allSurveyData_shuffled <- list(Data_S1_shuffled,Data_S2_shuffled,Data_S3_shuffled,Data_S4_shuffled)
length(unique(unlist(sapply(allSurveyData_shuffled, rownames)))) # Same number of vars as in the observed data?


vars <- unique(unlist(sapply(allSurveyData,rownames)))
isolates <- unique(unlist(sapply(allSurveyData,colnames)))
superMatrix <- matrix(0, nrow=length(vars), ncol=length(isolates), dimnames = list(vars, isolates))
superMatrix[rownames(Data_S1_shuffled),colnames(Data_S1_shuffled)] <- data.matrix(Data_S1_shuffled)
superMatrix[rownames(Data_S2_shuffled),colnames(Data_S2_shuffled)] <- data.matrix(Data_S2_shuffled)
superMatrix[rownames(Data_S3_shuffled),colnames(Data_S3_shuffled)] <- data.matrix(Data_S3_shuffled)
superMatrix[rownames(Data_S4_shuffled),colnames(Data_S4_shuffled)] <- data.matrix(Data_S4_shuffled)

similarityMatrix_shuffled <- overlapAlleleAdj(t(superMatrix))
identical(similarityMatrix_shuffled,similarityMatrixEmpirical)




# Diversity ---------------------------------------------------------------

allSurveyData <- list(Data_S1,Data_S2)
sapply(allSurveyData, nrow) # Number of vars per survey
sapply(allSurveyData, ncol) # Number of isolates per survey
allIsolates <- unique(unlist(sapply(allSurveyData, colnames)))
allVars <- unique(unlist(sapply(allSurveyData, rownames)))

## var diversity
vars_in_surveys <- matrix(0, nrow = length(allVars), ncol=4, dimnames=list(allVars,1:4))
vars_in_surveys[rownames(Data_S1),1] <- 1
vars_in_surveys[rownames(Data_S2),2] <- 1
vars_in_surveys[rownames(Data_S3),3] <- 1
vars_in_surveys[rownames(Data_S4),4] <- 1
vars_in_surveys <- bipartite::empty(vars_in_surveys)
n_vars <- nrow(vars_in_surveys)
# Proportion of vars per survey
colSums(vars_in_surveys)/n_vars
# Similarity between surveys in the var types
vegan::betadiver(t(vars_in_surveys), method = 10) 
calculatePTS(t(vars_in_surveys), return_matrix = T)
# Proportion of vars that appeared in n surveys
table(rowSums(vars_in_surveys))/n_vars

# Var persistence: How many vars pesisteed from the previous survey, and how may are new?
var_persistence_df <- data.frame(survey=1:4,persisted_prop=NA,new_prop=NA)
for (s in 2:4){
  # vars_persisted <- sum(vars_in_surveys[,2]==1 & vars_in_surveys[,1]==1)
  vars_new <- sum(vars_in_surveys[,s]-vars_in_surveys[,(s-1)]==1)
  total_vars <- sum(vars_in_surveys[,s]==1)
  vars_persisted <- total_vars-vars_new
  var_persistence_df[s,'persisted_prop'] <- vars_persisted/total_vars
  var_persistence_df[s,'new_prop'] <- vars_new/total_vars
}

var_persistence_df$N <- colSums(vars_in_surveys) # Absolute number of vars in each survey
var_persistence_df$Prop_per_survey <- colSums(vars_in_surveys)/n_vars # Proportion of vars (out of all the vars across surveys) that appear in each survey 
var_persistence_df$n_surveys <- table(rowSums(vars_in_surveys))/n_vars # Proportion of vars that appeared in n surveys 
var_persistence_df$isolates <- sapply(allSurveyData, ncol)

var_persistence_df[,c(2,3,5,6)] <- round(var_persistence_df[,c(2,3,5,6)],2)

# Calculate median MOI
var_persistence_df$MOI_median <- sapply(allSurveyData, function(x) median(colSums(x)/60))


d <- melt(var_persistence_df, id.vars='survey', measure.vars=c('N','prev_all'))
ggplot(d, aes(survey, value))+
  geom_point(size=3)+
  geom_line()+
  facet_wrap(~variable,scales='free')+
  labs(title='PTS S1-S4, upsBC, MOI=ALL')


## PTS distributions for each season
PTS_S1 <- data.frame(survey='S1', PTS=calculatePTS(t(Data_S1)))
PTS_S2 <- data.frame(survey='S2', PTS=calculatePTS(t(Data_S2)))
PTS_S3 <- data.frame(survey='S3', PTS=calculatePTS(t(Data_S3)))
PTS_S4 <- data.frame(survey='S4', PTS=calculatePTS(t(Data_S4)))
d_PTS <- rbind(PTS_S1,PTS_S2,PTS_S3,PTS_S4)
p1 <- ggplot(d_PTS, aes(survey, PTS, color=survey))+geom_violin()+labs(title='PTS S1-S4, upsBC, MOI=ALL')
p2 <- ggplot(d_PTS, aes(PTS, fill=survey))+geom_density()+facet_wrap(~survey)+labs(title='PTS S1-S4, upsBC, MOI=ALL')
p3 <- ggtexttable(var_persistence_df, rows = NULL, 
                  theme = ttheme("mOrange"))

pdf('~/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/Results/empirical_data_diversity_all_MOI.pdf', width = 11, height = 8)
gridExtra::grid.arrange(p1,p2,p3, layout_matrix = cbind(c(1,3), c(2,3)))
dev.off()


#### Plot for the grant
# var_persistence_df <- rbind(var_persistence_df,c(5,rep(NA,5)))
# var_persistence_df <- rbind(var_persistence_df,c(6,rep(NA,5)))
# var_persistence_df$prev_soe <- c(44.9,26.4,33.1,18.9,27.8,16.5)
# var_persistence_df$prev_vea <- c(38.6,27.6,24.9,27.3,26,9.4)
# var_persistence_df$prev_all <- c(41.9,27.0,29.3,23,26.9,13)

grant_plot_data <- data.frame(ID=1:8,
                              season=rep(c('Wet','Dry'),4),
                              year=c(2012,2013,2013,2014,2014,2015,2015,2016))
grant_plot_data$diversity <- NA
grant_plot_data$diversity[c(1,2,4,5)] <- var_persistence_df$N[1:4]
grant_plot_data$normalized_diversity <- NA
grant_plot_data$normalized_diversity[c(1,2,4,5)] <- var_persistence_df$N[1:4]/var_persistence_df$isolates[1:4]
grant_plot_data$MOI <- NA
grant_plot_data$MOI[c(1,2,4,5)] <- var_persistence_df$MOI_median

grant_plot_data$prevalence <- NA
grant_plot_data$prevalence[c(1,2,4,5,7,8)] <- c(41.9,27.0,29.3,23,26.9,13)


## First version
scale_factor_axis <- 500
d <- grant_plot_data[1:8,]
p1 <- ggplot(d, aes(x = ID))+
  scale_x_continuous(name='Survey', labels = d$season, breaks=1:8)
p1 <- p1 + geom_line(aes(y=diversity, colour = "Pool diversity"))+
  geom_point(aes(y=diversity, colour = "Pool diversity"), size=6)
p1 <- p1 + geom_line(aes(y=normalized_diversity*scale_factor_axis, colour = "Normalized Diversity"))+
  geom_point(aes(y=normalized_diversity*scale_factor_axis, colour = "Normalized Diversity"), size=6)
p1 <- p1 + labs(y='Pool diversity') + theme(legend.position = c(0.6, 0.9)) +
  scale_y_continuous(limits=c(10000,40000),
                     # breaks = seq(10000,40000,5000),
                     sec.axis = sec_axis(~./scale_factor_axis,
                                         name = "Normalized diversity\n[types/isolate]",
                                         breaks=seq(20,80,10))) # now adding the secondary axis, following the example in the help file ?scale_y_continuous

### Second version -- with MOI
scale_factor_axis <- 15000
p1 <- ggplot(grant_plot_data, aes(x = ID))+
  scale_x_continuous(name='Survey', labels = grant_plot_data$season, breaks=1:8)
p1 <- p1 + geom_line(aes(y=diversity, colour = "var DBLa types sampled"))+
  geom_point(aes(y=diversity, colour = "var DBLa types sampled"), size=6)
p1 <- p1 + geom_line(aes(y=MOI*scale_factor_axis, colour = "Median MOI"))+
  geom_point(aes(y=MOI*scale_factor_axis, colour = "Median MOI"), size=6)
p1 <- p1 + labs(y='var DBLa types sampled') + theme(legend.position = c(0.6, 0.9)) +
  scale_y_continuous(limits=c(10000,40000),
                     breaks = seq(10000,40000,5000),
                     sec.axis = sec_axis(~./scale_factor_axis,
                                         breaks = seq(0,2.5,0.5),
                                         name = "Median MOI"
                     )) # now adding the secondary axis, following the example in the help file ?scale_y_continuous



p2 <- ggplot(grant_plot_data, aes(x=ID, y=prevalence))+
  scale_x_continuous(name='Survey', labels = grant_plot_data$season, breaks=1:8)+
  scale_y_continuous(breaks=seq(0,100,10), limits=c(0,60))+
  geom_line()+
  geom_point(size=6)+labs(y='Prevalence [%]')


plot_grid(p1,p2, ncol=1, align = 'v')


pdf('~/Dropbox/RESEARCH/Grant_Malaria_2017/Fig_diversity_prevalence_MOI.pdf', 10,6)
plot_grid(p1,p2, ncol=1, align = 'v')
dev.off()

# ggplot(grant_plot_data, aes(ID, prevalence, fill=season))+geom_bar(stat='identity')
p_prevalence <- ggplot(grant_plot_data, aes(ID, prevalence, colour = "Prevalence"))+geom_point(size=3)+geom_line()+
  scale_x_continuous(breaks=1:8)

# var_persistence_df$persisted_prop <- var_persistence_df$persisted_prop * 100
# p <- ggplot(var_persistence_df, aes(x = survey))
# p <- p + geom_line(aes(y = N, colour = "Diversity"))+geom_point(aes(y = N, colour = "Diversity"),size=4)
# p <- p + geom_line(aes(y = prev_all*1000, colour = "Prevalence"))+geom_point(aes(y = prev_all*1000, colour = "Prevalence"),size=4)
# p <- p + scale_y_continuous(sec.axis = sec_axis(~./1000, name = "Prevalence [%]")) # now adding the secondary axis, following the example in the help file ?scale_y_continuous
# p <- p + scale_colour_manual(values = c("red", "blue")) + scale_x_continuous(breaks=1:6)
# p <- p + labs(y = "Diversity",
#               x = "Survey")
# p <- p + theme(legend.position = c(0.6, 0.9))
# ggsave('~/Dropbox/RESEARCH/Grant_Malaria_2017/empirical_diversity_prevalence.svg',p, width = 8, height = 6)

p_PTS_S1 <- ggplot(PTS_S1, aes(PTS, y=..scaled..))+geom_density(fill='purple')+xlim(c(0,1))
p_PTS_S2 <- ggplot(PTS_S2, aes(PTS, y=..scaled..))+geom_density(fill='purple')+xlim(c(0,1))
p_PTS_S3 <- ggplot(PTS_S3, aes(PTS, y=..scaled..))+geom_density(fill='purple')+xlim(c(0,1))
p_PTS_S4 <- ggplot(PTS_S4, aes(PTS, y=..scaled..))+geom_density(fill='purple')+xlim(c(0,1))
ggsave('~/Dropbox/RESEARCH/Grant_Malaria_2017/p_PTS_S1.svg',p_PTS_S1, width = 8, height = 6)
ggsave('~/Dropbox/RESEARCH/Grant_Malaria_2017/p_PTS_S2.svg',p_PTS_S2, width = 8, height = 6)
ggsave('~/Dropbox/RESEARCH/Grant_Malaria_2017/p_PTS_S3.svg',p_PTS_S3, width = 8, height = 6)
ggsave('~/Dropbox/RESEARCH/Grant_Malaria_2017/p_PTS_S4.svg',p_PTS_S4, width = 8, height = 6)

## var code diversity
# var codes are unique combinations of vars in isolates.
var_code_PTS <- calculatePTS(mainDataClean)
dim(var_code_PTS)
hist(var_code_PTS[lower.tri(var_code_PTS)])

#Between-season comparison
betweenSeasonPTS_mean <- matrix(0,4,4,dimnames = list(paste('S',1:4,sep=''),paste('S',1:4,sep='')))
for (s1 in paste('S',1:4,sep='')){
  for (s2 in paste('S',1:4,sep='')){
    x <- mainDataClean[]
    x <- var_code_PTS[str_sub(rownames(var_code_PTS),1,2)==s1,str_sub(colnames(var_code_PTS),1,2)==s2]
    betweenSeasonPTS_mean[s2,s1] <- round(mean(x[lower.tri(x)]),5)
  }
}
betweenSeasonPTS_mean[upper.tri(betweenSeasonPTS_mean)] <- 0



identical(rownames(var_code_PTS),colnames(var_code_PTS))
isolates <- rownames(var_code_PTS)
isolates <- str_sub(isolates,str_locate(isolates, 'MRS')[,1],str_locate(isolates, '\\.')[,1]-1)
isolates[duplicated(isolates)]
var_code_PTS[grep(isolates[duplicated(isolates)],rownames(var_code_PTS)),]

# calculateGeneticSimilarity_Sparse <- function(M){
#   x <- Matrix(M)
#   x <- tcrossprod(x>0)
#   x <- x/rowSums(x>0)
#   print(object.size(x))
#   print(dim(x))
#   return(x)
# }
# empiricalLayer_1 <- calculateGeneticSimilarity_Sparse(t(Data_S1))
# empiricalLayer_2 <- calculateGeneticSimilarity_Sparse(t(Data_S2))
# empiricalLayer_3 <- calculateGeneticSimilarity_Sparse(t(Data_S3))
# empiricalLayer_4 <- calculateGeneticSimilarity_Sparse(t(Data_S4))


