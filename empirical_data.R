# Initialize --------------------------------------------------------------
source('mtn_functions.R')

prep.packages(c("stringr","ggplot2","igraph","adegenet","dplyr","cowplot", "ggdendro","ggpubr"))

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

build_network_Infomap <- function(file='Infomap_survey_empirical.txt', layerList, similarityMatrix, numLayers=4){
  nodeLabel <- sort(unique(unlist(lapply(layerList,rownames))))
  nodeList <- data.frame(nodeID=1:length(nodeLabel), nodeLabel)
  # Build interlayer edges
  ile_edgelist <- list()
  for (t in 1:(numLayers-1)){
    strainCopies_t <- rownames(layerList[[t]])
    strainCopies_t1 <- rownames(layerList[[t+1]])
    m <- similarityMatrix[strainCopies_t,strainCopies_t1]
    NZ <- which(m!=0, arr.ind = T) # non-zero elements
    edges_interlayer <- data.frame(layer_s=rep(t,nrow(NZ)),
                                   node_s=rep(NA,nrow(NZ)),
                                   layer_t=rep(t+1,nrow(NZ)),
                                   node_t=rep(NA,nrow(NZ)),
                                   w=rep(NA,nrow(NZ)))
    
    edges_interlayer$node_s <- rownames(m)[NZ[,1]]
    edges_interlayer$node_t <- colnames(m)[NZ[,2]]
    for (n in 1:nrow(edges_interlayer)){
      edges_interlayer[n,'w'] <- m[as.character(edges_interlayer[n,'node_s']),as.character(edges_interlayer[n,'node_t'])]
    }
    edges_interlayer$node_s <- nodeList$nodeID[match(edges_interlayer$node_s,nodeList$nodeLabel)]
    edges_interlayer$node_t <- nodeList$nodeID[match(edges_interlayer$node_t,nodeList$nodeLabel)]
    
    ile_edgelist[[t]] <- edges_interlayer
  }
  
  inter_edges <- do.call(rbind.data.frame,ile_edgelist)
  intra_edges <- infomap_makeIntralayerEdges(layerList,nodeList)
  
  # ggHistogram(inter_edges$w)
  # ggHistogram(intra_edges$w)
  
  ## Write file for infomap
  message('Writing Infomap files')
  print(paste('Infomap file:',file))
  if (file.exists(file)){unlink(file)}
  sink(file, append = T)
  cat("# A network in a general multiplex format");cat('\n')
  cat(paste("*Vertices",nrow(nodeList)));cat('\n')
  write.table(nodeList, file, append = T,sep=' ', quote = T, row.names = F, col.names = F)
  cat("*Multiplex");cat('\n')
  cat("# layer node layer node [weight]");cat('\n')
  cat("# Intralayer edges");cat('\n')
  write.table(intra_edges, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
  cat("# Interlayer edges");cat('\n')
  write.table(inter_edges, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
  sink.reset()
}

# Get and clean empirical data ---------------------------------------------------

mainData <- read.table('~/Dropbox/Qixin_Shai_Malaria/S1-S4/S_DNAseq_all_inframe_renamed_otuTable.txt', header = T)
anyDuplicated(mainData$OTU_ID)
anyDuplicated(colnames(mainData))
varGroupings <- read.table('~/Dropbox/Qixin_Shai_Malaria/S1-S4/S_AAseq_all_inframe_renamed_centroids.fasta_classify.txt', header = T)
mainData <- subset(mainData, OTU_ID%in%subset(varGroupings, Grouping=='BC')$target)
rownames(mainData) <- mainData[,1]
mainData <- mainData[,-1]

### Separate by surveys
isolates_all <- colnames(mainData)
S1_isolates_idx <- which(str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='S1' | 
                           str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='RS1' | # RS1 and R2S1 are isolates for which infectino was not detect with the first PCR but with subsequent PCRs
                           str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='R2S1')
S2_isolates_idx <- which(str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='S2')
S3_isolates_idx <- which(str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='S3')
S4_isolates_idx <- which(str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='S4')
Pl_isolates_idx <- which(str_sub(isolates_all, 1, str_locate(isolates_all, 'MRS')[,1]-1)=='P')

# Are all isolates included?
length(isolates_all) - (length(S1_isolates_idx)+length(S2_isolates_idx)+length(S3_isolates_idx)+length(S4_isolates_idx))


S1_isolates <- str_sub(isolates_all[S1_isolates_idx],str_locate(isolates_all[S1_isolates_idx], 'MRS')[,1],str_locate(isolates_all[S1_isolates_idx], '\\.')[,1]-1)
any(duplicated(S1_isolates)) # make sure that hosts with RS1 or R2S1 are not in the S1 set.
S2_isolates <- str_sub(isolates_all[S2_isolates_idx],3,str_locate(isolates_all[S2_isolates_idx], '\\.')[,1]-1)
S3_isolates <- str_sub(isolates_all[S3_isolates_idx],3,str_locate(isolates_all[S3_isolates_idx], '\\.')[,1]-1)
S4_isolates <- str_sub(isolates_all[S4_isolates_idx],3,str_locate(isolates_all[S4_isolates_idx], '\\.')[,1]-1)
Pl_isolates <- str_sub(isolates_all[Pl_isolates_idx],3,str_locate(isolates_all[Pl_isolates_idx], '\\.')[,1]-1)

Data_S1 <- mainData[,S1_isolates_idx]
Data_S2 <- mainData[,S2_isolates_idx]
Data_S3 <- mainData[,S3_isolates_idx]
Data_S4 <- mainData[,S4_isolates_idx]

# Limit to 55 vars (MOI=1)
maxVarIsolate <- 55 # 55 because we do not take grouping A
Data_S1 <- cleanSurveyData(Data_S1, min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S2 <- cleanSurveyData(Data_S2, min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S3 <- cleanSurveyData(Data_S3, min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S4 <- cleanSurveyData(Data_S4, min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)

mainDataClean <- mainData[unique(c(rownames(Data_S1),rownames(Data_S2),rownames(Data_S3),rownames(Data_S4))),
                          unique(c(colnames(Data_S1),colnames(Data_S2),colnames(Data_S3),colnames(Data_S4)))]
mainDataClean <- bipartite::empty(mainDataClean)
dim(mainDataClean)

# mainData <- cleanSurveyData(mainData, min.var.per.isolate = 30, max.var.per.isolate = maxVarIsolate, plotit = T)
# PTS_mainData <- calculatePTS(mainData)
# ggHistogram(as.vector(PTS_mainData[lower.tri(PTS_mainData)]))
# summary(PTS_mainData[lower.tri(PTS_mainData)])


# Diversity ---------------------------------------------------------------

allSurveyData <- list(Data_S1,Data_S2,Data_S3,Data_S4)
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



# Define Layers -----------------------------------------------------------
# # Remove the single individual with a chronic infection
# Data_S1 <- Data_S1[,-grep(pattern = 'MRS1011', colnames(Data_S1))]
# Data_S2 <- Data_S2[,-grep(pattern = 'MRS1011', colnames(Data_S2))]

empiricalLayer_1 <- overlapAlleleAdj(t(Data_S1))
empiricalLayer_2 <- overlapAlleleAdj(t(Data_S2))
empiricalLayer_3 <- overlapAlleleAdj(t(Data_S3))
empiricalLayer_4 <- overlapAlleleAdj(t(Data_S4))

diag(empiricalLayer_1) <- 0
diag(empiricalLayer_2) <- 0
diag(empiricalLayer_3) <- 0
diag(empiricalLayer_4) <- 0

layersEmpirical <- list(empiricalLayer_1, empiricalLayer_2, empiricalLayer_3, empiricalLayer_4)
similarityMatrixEmpirical <- overlapAlleleAdj(t(mainDataClean))

sapply(layersEmpirical, gdensity)

pdf('~/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/Results/empirical_layers_no_cutoff.pdf', width = 8, height = 8)
plotSurveyLayer(layersEmpirical[[1]], main='S1, Oct 2012', vertex.color='navy')
plotSurveyLayer(layersEmpirical[[2]], main='S2, Jun 2013', vertex.color='navy')
plotSurveyLayer(layersEmpirical[[3]], main='S3, Jun 2014', vertex.color='navy')
plotSurveyLayer(layersEmpirical[[4]], main='S4, Oct 2014', vertex.color='navy')
dev.off()

# If not cutoff is used then set probability to 0

# Define cutoff -----------------------------------------------------------
write.csv(similarityMatrixEmpirical, 'similarityMatrixEmpirical.csv') # save this so to be able tro plt edge weight distributons fast

cutoffProbability <- 0.9

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

# Plot layers -------------------------------------------------------------

pdf('~/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/Results/empirical_layers.pdf', width = 8, height = 8)
plotSurveyLayer(layersEmpirical_cutoff[[1]], main='S1, Oct 2012', vertex.color='navy')
plotSurveyLayer(layersEmpirical_cutoff[[2]], main='S2, Jun 2013', vertex.color='navy')
plotSurveyLayer(layersEmpirical_cutoff[[3]], main='S3, Jun 2014', vertex.color='navy')
plotSurveyLayer(layersEmpirical_cutoff[[4]], main='S4, Oct 2014', vertex.color='navy')
dev.off()

# Run Infomap -- Empirical -----------------------------------------------
build_network_Infomap('Infomap_survey_empirical.txt', layersEmpirical_cutoff, similarityMatrixEmpirical_cutoff)

system('./Infomap Infomap_survey_empirical.txt . -2 -i multiplex --multiplex-relax-rate -1 -d -N 20 --rawdir --tree --expanded --silent')

modules_empirical <- infomap_readTreeFile('Infomap_survey_empirical_expanded.tree', reorganize_modules = T, remove_buggy_instances = T,max_layers = 4)
## Create a module-by-layer summary
moduleSummary_empirical <- modules_empirical %>% group_by(layer,module) %>% summarise(numStrains=length(strain))
moduleSummary_empirical <- subset(moduleSummary_empirical, layer<=4)
# moduleSummary$strain_prop <- moduleSummary$numStrains / lInfo$strains_total[moduleSummary$layer]
moduleSummary_empirical$module <- factor(moduleSummary_empirical$module)
x <- xtabs(~module+layer, data = modules_empirical)
table(rowSums(binarize(x)))/nModules # The proportion of modules that appear in n layers
prop_mods_in_layers <- data.frame(table(rowSums(binarize(x)))/nModules) # The proportion of modules that appear in n layers
p_prop_in_layers_empirical <- ggplot(prop_mods_in_layers, aes(Var1,Freq))+geom_bar(stat = 'identity')+labs(x='Number of layers')



nModules <- max(modules_empirical$module)
module_continuity_empirical <- data.frame(module=rownames(x),num_layers=rowSums(binarize(x)))
moduleSummary_empirical <- merge(moduleSummary_empirical,module_continuity_empirical, all.x=T)
p_empirical <- ggplot(moduleSummary_empirical, aes(x=layer, y=module, size=3, color=as.factor(num_layers)))+
  geom_point()+
  theme(legend.position="none")+
  scale_color_manual(values=c('#fef0d9','#fdcc8a','#fc8d59','#d7301f'))+
  labs(title=paste('Empirical, upsBC, MOI=1, cutoff=',cutoffProbability))
p_empirical

# Distribution of module sizes
total_number_strains <- length(unique(modules_empirical$strain))

moduleSize_empirical <- modules_empirical %>% group_by(module) %>% summarise(strainProp=length(strain)/total_number_strains)
p_moduleSize_empirical <- ggplot(moduleSize_empirical, aes(strainProp, y=..scaled..))+geom_density(fill = 'navy')+labs(x = 'Module size')

# # Proportiion of new moduels appearing after intervention
# moduleInLayer <- binarize(xtabs(~module+layer,modules_empirical))
# beginning <- apply(moduleInLayer,1, function(x) which(x!=0)[1]) # position of first non-zero element from each row.
# end <- apply(moduleInLayer,1, function(x) tail(which(x!=0),1)) # position of last non-zero element from each row.
# moduleIntervals <- data.frame(start=beginning,end=end)
# table(beginning)/nModules

cowplot::plot_grid(p_empirical,p_edge_weight_histogram_empirical,p_moduleSize_empirical, nrow=3)

##### Plots for presentation
p_prop_in_layers_empirical <- ggplot(prop_mods_in_layers, aes(Var1,Freq))+
  geom_bar(stat = 'identity', fill='navy',color='black')+labs(x='Number of layers')

p_empirical <- ggplot(moduleSummary_empirical, aes(x=layer, y=module))+
  geom_point(size=5, color='navy')+
  theme(legend.position="none")

p_moduleSize_empirical <- ggplot(moduleSize_empirical, aes(strainProp, y=..scaled..))+geom_density(fill = 'navy')+labs(x = 'Normalized module size', y='Frequency')

plot_grid(p_empirical,p_prop_in_layers_empirical,p_moduleSize_empirical, nrow=1)

svg('~/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/Presentations/temporal_empirical_data.svg',16,8)
plot_grid(p_empirical,p_prop_in_layers_empirical,p_moduleSize_empirical, nrow=1)
dev.off()

pdf('~/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/Presentations/temporal_empirical_data_separate.pdf',6,4)
p_empirical
p_prop_in_layers_empirical
p_moduleSize_empirical
p_edge_weight_histogram_empirical
dev.off()



# PTS within/between modules ----------------------------------------------


### Calculate PTS within modules
strainByGene <- t(mainDataClean)
moduleSize_empirical$PTS <- NA
for (m in subset(moduleSize_empirical, numStrains>2)$module){
  moduleSize_empirical[moduleSize_empirical$module==m,'PTS'] <- modulePTS(modules_empirical, moduleId = m)
}
ggplot(moduleSize_empirical, aes(x=numStrains,y=PTS))+geom_point()+geom_smooth(method='lm', se=F)
d <- moduleSize_empirical %>% group_by(numStrains) %>% summarise(numModules=length(module), meanPTS=mean(PTS))
ggplot(d, aes(x=numStrains,y=meanPTS, size=numModules))+geom_point()+geom_smooth(method='lm', se=F)
mean(moduleSize_empirical$PTS, na.rm = T)

### PTS Between modules
repertoires_vars_df <- matrix.to.coordinateList(mainDataClean)
names(repertoires_vars_df) <- c('var','strain','w')
d <- merge(modules_empirical,repertoires_vars_df, by='strain')
x <- xtabs(~module+var, d)
between_modules_PTS <- calculatePTS(x)

### Within and between
ggplot()+
  geom_density(aes(moduleSize_empirical$PTS, y=..scaled..), fill='red', alpha=0.1)+
  geom_density(aes(calculatePTS(x), y=..scaled..),fill='blue', alpha=0.1)

# Module/Repertoire distribution across isolates---------------------------

# Because we work with MOI=1 then the distribution of repertoires in isoaltes 
# will always produce maximum entropy, since an isolate is a repertoire. So it 
# remains to see how modules are distributed. The number of cases is actually
# the number of strains (or isolates); again, becaues of MOI=1.

moduleDistribution <- modules_empirical %>% group_by(layer,module) %>% summarise(numStrains=length(strain))
diversModules <- moduleDistribution %>% group_by(layer) %>% summarise(H_normalized=vegan::diversity(numStrains)/log(length(module)+1))




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


# Layer Permutations ------------------------------------------------------
library(gtools)
perms <- permutations(4,4)
permPlots <- list()
for (p in 1:nrow(perms)){
  print(p)
  layersEmpirical_perm <- layersEmpirical_cutoff[perms[p,]]
  f <- paste('perm_',p,'.txt',sep='')
  build_network_Infomap(f, layersEmpirical_perm, similarityMatrixEmpirical_cutoff)
  
  system(paste('./Infomap ',f,' . -2 -i multiplex --multiplex-relax-rate -1 -d -N 10 --rawdir --tree --map --expanded --silent',sep=''))
  
  modules_perm <- infomap_readTreeFile(paste('perm_',p,'_expanded.tree',sep=''), reorganize_modules = T, remove_buggy_instances = F,max_layers = numLayers)
  modules_perm <- subset(modules_perm, layer<=4)
  moduleSummary_perm <- modules_perm %>% group_by(layer,module) %>% summarise(numStrains=length(strain))
  moduleSummary_perm$module <- factor(moduleSummary_perm$module)
  permPlots[[p]] <- ggplot(moduleSummary_perm, aes(x=layer, y=module, size=numStrains, color=as.factor(module)))+geom_point()+
    theme(legend.position="none")+labs(title=paste(perms[p,],collapse = '-'))
  moduleSize_perm <- modules_perm %>% group_by(module) %>% summarise(numStrains=length(strain))
  print(ks.test(moduleSize_empirical$numStrains, moduleSize_perm$numStrains)$p.value)
}


cowplot::plot_grid(plotlist=permPlots, ncol = 6, nrow=4)




# Shuffle order of isolates in time ---------------------------------------

# This procudeure allocates isolates uniformly at random to the different
# surveys. It maintains the number of isolates per survey.

results_list <- list()
plots_list <- list()
ks_results <- list()
for(r in 1:20){
  isolateSampleSize <- sapply(allSurveyData, ncol)
  isolateSet <- allIsolates
  shuffledIsolates_S1 <- sample(isolateSet, isolateSampleSize[1], replace = F)
  isolateSet <- isolateSet[!isolateSet %in% shuffledIsolates_S1]
  shuffledIsolates_S2 <- sample(isolateSet, isolateSampleSize[2], replace = F)
  isolateSet <- isolateSet[!isolateSet %in% shuffledIsolates_S2]
  shuffledIsolates_S3 <- sample(isolateSet, isolateSampleSize[3], replace = F)
  isolateSet <- isolateSet[!isolateSet %in% shuffledIsolates_S3]
  shuffledIsolates_S4 <- isolateSet
  
  Data_S1_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S1]
  Data_S1_shuffled <- bipartite::empty(Data_S1_shuffled)
  dim(Data_S1_shuffled)
  Data_S2_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S2]
  Data_S2_shuffled <- bipartite::empty(Data_S2_shuffled)
  dim(Data_S2_shuffled)
  Data_S3_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S3]
  Data_S3_shuffled <- bipartite::empty(Data_S3_shuffled)
  dim(Data_S3_shuffled)
  Data_S4_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S4]
  Data_S4_shuffled <- bipartite::empty(Data_S4_shuffled)
  dim(Data_S4_shuffled)
  allSurveyData_shuffled <- list(Data_S1_shuffled,Data_S2_shuffled,Data_S3_shuffled,Data_S4_shuffled)
  length(unique(unlist(sapply(allSurveyData_shuffled, rownames)))) # Same number of vars as in the observed data?
  
  vars <- unique(unlist(sapply(allSurveyData_shuffled,rownames)))
  isolates <- unique(unlist(sapply(allSurveyData_shuffled,colnames)))
  superMatrix <- matrix(0, nrow=length(vars), ncol=length(isolates), dimnames = list(vars, isolates))
  superMatrix[rownames(Data_S1_shuffled),colnames(Data_S1_shuffled)] <- data.matrix(Data_S1_shuffled)
  superMatrix[rownames(Data_S2_shuffled),colnames(Data_S2_shuffled)] <- data.matrix(Data_S2_shuffled)
  superMatrix[rownames(Data_S3_shuffled),colnames(Data_S3_shuffled)] <- data.matrix(Data_S3_shuffled)
  superMatrix[rownames(Data_S4_shuffled),colnames(Data_S4_shuffled)] <- data.matrix(Data_S4_shuffled)
  
  similarityMatrix_shuffled <- overlapAlleleAdj(t(superMatrix))
  identical(similarityMatrix_shuffled,similarityMatrixEmpirical)
  
  # Create the layers
  empiricalLayer_1_shuffled <- overlapAlleleAdj(t(Data_S1_shuffled))
  empiricalLayer_2_shuffled <- overlapAlleleAdj(t(Data_S2_shuffled))
  empiricalLayer_3_shuffled <- overlapAlleleAdj(t(Data_S3_shuffled))
  empiricalLayer_4_shuffled <- overlapAlleleAdj(t(Data_S4_shuffled))
  diag(empiricalLayer_1_shuffled) <- 0
  diag(empiricalLayer_2_shuffled) <- 0
  diag(empiricalLayer_3_shuffled) <- 0
  diag(empiricalLayer_4_shuffled) <- 0
  layersEmpirical_shuffled <- list(empiricalLayer_1_shuffled, empiricalLayer_2_shuffled, empiricalLayer_3_shuffled, empiricalLayer_4_shuffled)
  
  # If applying cutoff
  cutoffValue <- quantile(as.vector(similarityMatrixEmpirical[similarityMatrixEmpirical!=0]), probs = cutoffProbability)
  similarityMatrix_shuffled[similarityMatrix_shuffled<cutoffValue] <- 0
  for (i in 1:length(layersEmpirical_shuffled)){
    x <-layersEmpirical_shuffled[[i]]
    x[x<cutoffValue] <- 0
    layersEmpirical_shuffled[[i]] <- x
  }
  
  ###
  build_network_Infomap('Infomap_survey_empirical_shuffled.txt', layersEmpirical_shuffled, similarityMatrix_shuffled)
  
  system('./Infomap Infomap_survey_empirical_shuffled.txt . -2 -i multiplex --multiplex-relax-rate -1 -d -N 10 --rawdir --tree --expanded --silent')
  
  modules_shuffled <- infomap_readTreeFile('Infomap_survey_empirical_shuffled_expanded.tree', reorganize_modules = T, remove_buggy_instances = F,max_layers = numLayers)
  modules_shuffled <- subset(modules_shuffled, layer<=4)
  moduleSummary_shuffled <- modules_shuffled %>% group_by(layer,module) %>% summarise(numStrains=length(strain))
  moduleSummary_shuffled$module <- factor(moduleSummary_shuffled$module)
  nModules <- max(modules_shuffled$module)
  
  moduleSize_shuffled <- modules_shuffled %>% group_by(module) %>% summarise(numStrains=length(strain))
  ks_results[[r]] <- ks.test(moduleSize_empirical$numStrains, moduleSize_shuffled$numStrains)$p.value
  
  p <- ggplot(moduleSummary_shuffled, aes(x=layer, y=module, size=numStrains, color=as.factor(module)))+geom_point()+
    theme(legend.position="none")+labs(title=r)
  plots_list[[r]] <- p
  results_list[[r]] <- moduleSummary_shuffled
  
}

cowplot::plot_grid(plotlist=plots_list[1:20], ncol = 5, nrow=4)+labs(title='')

#sapply(results_list, function(z) any(table(z$module)>2))

# Shuffle vars among isolates ---------------------------------------------
# 
# # This procedure shuffles the order of vars in time
# 
# varSampleSize <- sapply(allSurveyData, nrow)
# varSet <- unlist(sapply(allSurveyData, rownames)) # This includes duplicates in vars because some of them appear in >1 layers
# shuffledVars_S1 <- sample(varSet, varSampleSize[1], replace = F)
# varSet <- varSet[!varSet %in% shuffledVars_S1]
# shuffledVars_S2 <- sample(varSet, varSampleSize[2], replace = F)
# varSet <- varSet[!varSet %in% shuffledVars_S2]
# shuffledVars_S3 <- sample(varSet, varSampleSize[3], replace = F)
# varSet <- varSet[!varSet %in% shuffledVars_S3]
# shuffledVars_S4 <- varSet
# 
# Data_S1_shuffled <- mainDataClean[rownames(mainDataClean)%in%shuffledVars_S1,colnames(mainDataClean)%in%colnames(Data_S1)]
# Data_S1_shuffled <- bipartite::empty(Data_S1_shuffled)
# dim(Data_S1_shuffled)
# Data_S2_shuffled <- mainDataClean[rownames(mainDataClean)%in%shuffledVars_S2,colnames(mainDataClean)%in%colnames(Data_S2)]
# Data_S2_shuffled <- bipartite::empty(Data_S2_shuffled)
# dim(Data_S2_shuffled)
# Data_S3_shuffled <- mainDataClean[rownames(mainDataClean)%in%shuffledVars_S3,colnames(mainDataClean)%in%colnames(Data_S3)]
# Data_S3_shuffled <- bipartite::empty(Data_S3_shuffled)
# dim(Data_S3_shuffled)
# Data_S4_shuffled <- mainDataClean[rownames(mainDataClean)%in%shuffledVars_S4,colnames(mainDataClean)%in%colnames(Data_S4)]
# Data_S4_shuffled <- bipartite::empty(Data_S4_shuffled)
# dim(Data_S4_shuffled)
# 
# allSurveyData_shuffled <- list(Data_S1_shuffled,Data_S2_shuffled,Data_S3_shuffled,Data_S4_shuffled)
# length(unique(unlist(sapply(allSurveyData_shuffled, rownames)))) # Same number of vars as in the observed data?
# length(unique(unlist(sapply(allSurveyData_shuffled, colnames)))) # Same number of isolates as in the observed data?
# 
# vars <- unique(unlist(sapply(allSurveyData_shuffled,rownames)))
# isolates <- unique(unlist(sapply(allSurveyData_shuffled,colnames)))
# superMatrix <- matrix(0, nrow=length(vars), ncol=length(isolates), dimnames = list(vars, isolates))
# superMatrix[rownames(Data_S1_shuffled),colnames(Data_S1_shuffled)] <- data.matrix(Data_S1_shuffled)
# superMatrix[rownames(Data_S2_shuffled),colnames(Data_S2_shuffled)] <- data.matrix(Data_S2_shuffled)
# superMatrix[rownames(Data_S3_shuffled),colnames(Data_S3_shuffled)] <- data.matrix(Data_S3_shuffled)
# superMatrix[rownames(Data_S4_shuffled),colnames(Data_S4_shuffled)] <- data.matrix(Data_S4_shuffled)
# dim(superMatrix)
# similarityMatrix_shuffled <- overlapAlleleAdj(t(superMatrix))
# identical(similarityMatrix_shuffled,similarityMatrixEmpirical)

mainDataClean_shuffled <- mainDataClean[sample(rownames(mainDataClean), nrow(mainDataClean), F),]
identical(mainDataClean_shuffled,mainDataClean)
similarityMatrix_shuffled <- overlapAlleleAdj(t(mainDataClean_shuffled))
identical(similarityMatrix_shuffled,similarityMatrixEmpirical)


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

# Run Infomap for shuffled networks ---------------------------------------


# Create the layers
empiricalLayer_1_shuffled <- overlapAlleleAdj(t(Data_S1_shuffled))
empiricalLayer_2_shuffled <- overlapAlleleAdj(t(Data_S2_shuffled))
empiricalLayer_3_shuffled <- overlapAlleleAdj(t(Data_S3_shuffled))
empiricalLayer_4_shuffled <- overlapAlleleAdj(t(Data_S4_shuffled))
diag(empiricalLayer_1_shuffled) <- 0
diag(empiricalLayer_2_shuffled) <- 0
diag(empiricalLayer_3_shuffled) <- 0
diag(empiricalLayer_4_shuffled) <- 0
layersEmpirical_shuffled <- list(empiricalLayer_1_shuffled, empiricalLayer_2_shuffled, empiricalLayer_3_shuffled, empiricalLayer_4_shuffled)

# If applying cutoff
cutoffValue <- quantile(as.vector(similarityMatrix_shuffled[similarityMatrix_shuffled!=0]), probs = cutoffProbability)
ggHistogram(as.vector(similarityMatrix_shuffled[similarityMatrix_shuffled!=0]), add_vline = cutoffValue)
similarityMatrix_shuffled[similarityMatrix_shuffled<cutoffValue] <- 0
for (i in 1:length(layersEmpirical_shuffled)){
  x <-layersEmpirical_shuffled[[i]]
  x[x<cutoffValue] <- 0
  layersEmpirical_shuffled[[i]] <- x
}


###
build_network_Infomap('Infomap_survey_empirical_shuffled.txt', layersEmpirical_shuffled, similarityMatrix_shuffled)

system('./Infomap Infomap_survey_empirical_shuffled.txt . -2 -i multiplex --multiplex-relax-rate -1 -d -N 10 --rawdir --tree --expanded')

modules_shuffled <- infomap_readTreeFile('Infomap_survey_empirical_shuffled_expanded.tree', reorganize_modules = T, remove_buggy_instances = T,max_layers = 4)
moduleSummary_shuffled <- modules_shuffled %>% group_by(layer,module) %>% summarise(numStrains=length(strain))
moduleSummary_shuffled$module <- factor(moduleSummary_shuffled$module)
ggplot(moduleSummary_shuffled, aes(x=layer, y=module, size=numStrains, color=as.factor(module)))+geom_point()+
  theme(legend.position="none")+labs(title='MOI=1, cutoff=none')




# Compare networks before and after intervention --------------------------

# This section is to compare the networks without modularity to see if their structure has changed.
# It is not an explicit temporal analysis.
# Calculate structural metrics for the observed networks. Then cluster the 
# features to try to dewtect similarity between the networks. It is clear that
# the networks are clustered by seasonality. S1 and S4 are both in the end of
# the wet season, S2 and S3 are in the wnd of the dry season.
features <- sapply(layersEmpirical, calculateFeatures)
features_scaled <- scale(features)
hc_observed <- hclust(dist(t(features_scaled)), method = 'average')
p_dendrogram_empirical <- ggdendrogram(hc_observed)+labs(title='Empirical')


## Now, to test the null hypothesis that time does not play a role in the
## structure of the network, compare that observed clustering to that derived in
## 100 shuffled networks. To shuffle the networks we permute the order of isolates.

shuffled <- list()
for(r in 1:100){
  print(r)
  isolateSampleSize <- sapply(allSurveyData, ncol)
  isolateSet <- allIsolates
  shuffledIsolates_S1 <- sample(isolateSet, isolateSampleSize[1], replace = F)
  isolateSet <- isolateSet[!isolateSet %in% shuffledIsolates_S1]
  shuffledIsolates_S2 <- sample(isolateSet, isolateSampleSize[2], replace = F)
  isolateSet <- isolateSet[!isolateSet %in% shuffledIsolates_S2]
  shuffledIsolates_S3 <- sample(isolateSet, isolateSampleSize[3], replace = F)
  isolateSet <- isolateSet[!isolateSet %in% shuffledIsolates_S3]
  shuffledIsolates_S4 <- isolateSet
  
  Data_S1_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S1]
  Data_S1_shuffled <- bipartite::empty(Data_S1_shuffled)
  dim(Data_S1_shuffled)
  Data_S2_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S2]
  Data_S2_shuffled <- bipartite::empty(Data_S2_shuffled)
  dim(Data_S2_shuffled)
  Data_S3_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S3]
  Data_S3_shuffled <- bipartite::empty(Data_S3_shuffled)
  dim(Data_S3_shuffled)
  Data_S4_shuffled <- mainDataClean[,colnames(mainDataClean)%in%shuffledIsolates_S4]
  Data_S4_shuffled <- bipartite::empty(Data_S4_shuffled)
  dim(Data_S4_shuffled)
  allSurveyData_shuffled <- list(Data_S1_shuffled,Data_S2_shuffled,Data_S3_shuffled,Data_S4_shuffled)
  length(unique(unlist(sapply(allSurveyData_shuffled, rownames)))) # Same number of vars as in the observed data?
  
  vars <- unique(unlist(sapply(allSurveyData_shuffled,rownames)))
  isolates <- unique(unlist(sapply(allSurveyData_shuffled,colnames)))
  superMatrix <- matrix(0, nrow=length(vars), ncol=length(isolates), dimnames = list(vars, isolates))
  superMatrix[rownames(Data_S1_shuffled),colnames(Data_S1_shuffled)] <- data.matrix(Data_S1_shuffled)
  superMatrix[rownames(Data_S2_shuffled),colnames(Data_S2_shuffled)] <- data.matrix(Data_S2_shuffled)
  superMatrix[rownames(Data_S3_shuffled),colnames(Data_S3_shuffled)] <- data.matrix(Data_S3_shuffled)
  superMatrix[rownames(Data_S4_shuffled),colnames(Data_S4_shuffled)] <- data.matrix(Data_S4_shuffled)
  
  similarityMatrix_shuffled <- overlapAlleleAdj(t(superMatrix))
  # print(identical(similarityMatrix_shuffled,similarityMatrixEmpirical))
  
  # Create the layers
  empiricalLayer_1_shuffled <- overlapAlleleAdj(t(Data_S1_shuffled))
  empiricalLayer_2_shuffled <- overlapAlleleAdj(t(Data_S2_shuffled))
  empiricalLayer_3_shuffled <- overlapAlleleAdj(t(Data_S3_shuffled))
  empiricalLayer_4_shuffled <- overlapAlleleAdj(t(Data_S4_shuffled))
  diag(empiricalLayer_1_shuffled) <- 0
  diag(empiricalLayer_2_shuffled) <- 0
  diag(empiricalLayer_3_shuffled) <- 0
  diag(empiricalLayer_4_shuffled) <- 0
  layersEmpirical_shuffled <- list(empiricalLayer_1_shuffled, empiricalLayer_2_shuffled, empiricalLayer_3_shuffled, empiricalLayer_4_shuffled)
  
  # If applying cutoff
  # cutoffValue <- quantile(as.vector(similarityMatrix_shuffled[similarityMatrix_shuffled!=0]), probs = cutoffProbability)
  # similarityMatrix_shuffled[similarityMatrix_shuffled<cutoffValue] <- 0
  # for (i in 1:length(layersEmpirical_shuffled)){
  #   x <-layersEmpirical_shuffled[[i]]
  #   x[x<cutoffValue] <- 0
  #   layersEmpirical_shuffled[[i]] <- x
  # }
  shuffled[[r]] <- layersEmpirical_shuffled
}

# Test for signifance: In how many of the permutations the order of the clustering was the same?
counter=0
identical_clustering_plots <- clustering_plots <- list()
features_shuffled <- list()
for (r in 1:99){
  features_shuffled[[r]] <- sapply(shuffled[[r]], calculateFeatures)
  x <- scale(features_shuffled[[r]])
  x <- dist(t(x))
  hc_shuffled <- hclust(x, method = 'average')
  print(paste(r,'--',hc_shuffled$order))
  clustering_plots[[r]] <- ggdendrogram(hc_shuffled)+labs(title=r)
  if(identical(hc_shuffled$order,hc_observed$order) | identical(hc_shuffled$order,rev(hc_observed$order))){
    counter <- counter+1
    if(length(hc_shuffled$height)==length(hc_observed$height)){
      identical_clustering_plots[[counter]] <- ggdendrogram(hc_shuffled)+labs(title=r)
    }
  }  
}
counter/100
# Plot clustering results for shuffled networks similar to the observed
identical_clustering_plots[[counter+1]] <- p_dendrogram_empirical
cowplot::plot_grid(plotlist=identical_clustering_plots)
# plot all the clustering results
clustering_plots[[r+1]] <- p_dendrogram_empirical
cowplot::plot_grid(plotlist=clustering_plots)

# Should visually inspect the shuffled clustering as well!!
pdf('~/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/intermediate_results/empirical_clustering_shuffled_comparison',20,16)
cowplot::plot_grid(plotlist=clustering_plots)
dev.off()
