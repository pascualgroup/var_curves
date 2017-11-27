analyzeModularStructure <- function(modules, verbose=T){
  cat('=====================================================');cat('\n')
  modularStructureMetrics <- rep(NA,10)
  names(modularStructureMetrics) <- c('num_modules','number_layers_min','number_layers_max','layer_first','layer_last','temporal_overlap','module_size_mean','module_size_entropy','module_size_sd','mean_strain_similarity')
  
  modularStructureMetrics['num_modules'] <- length(unique(modules$module_ordered))
  
  ### 1. Module persistence and overlap
  modulesLayers <- binarize(table(modules$module_ordered,modules$layer))
  modulePersistence <- rowSums(modulesLayers)
  if (verbose) {
    cat('Modules persistence:');cat('\n')
    print(modulePersistence)
    cat('------------');cat('\n')
  }
  modularStructureMetrics['number_layers_min'] <- min(modulePersistence)
  modularStructureMetrics['number_layers_max'] <- max(modulePersistence)
  modularStructureMetrics['layer_first'] <- as.numeric(colnames(modulesLayers)[1])
  modularStructureMetrics['layer_last'] <- as.numeric(colnames(modulesLayers)[ncol(modulesLayers)])
  # The overlap in time is the similarity in the layers in which they occurred. 
  # Hence, the fact that a module did not occur throughout the life span of 
  # another is taken into account. This implicitly considers the life span of 
  # modules. In that example the overlap is <1 because module 1
  # occured throughout a shorter period of module 2's life span.
  # mod1:     +++++ 
  # mod2: +++++++++++++++++++
  moduleTemporalOverlap <- mean(round(1-vegdist(modulesLayers, method = 'jaccard'),3))
  modularStructureMetrics['temporal_overlap'] <- moduleTemporalOverlap
  
  ### 2. Number of strains in a module
  d=modules[,c('strain','module_ordered')]
  d=unique(d)
  strainsInModules <- table(d$module_ordered)
  if (verbose){
    cat('Unique strains in modules:');cat('\n')
    print(strainsInModules)
    cat('------------');cat('\n')
  }
  modularStructureMetrics['module_size_mean'] <- mean(strainsInModules)
  modularStructureMetrics['module_size_sd'] <- sd(strainsInModules)
  modularStructureMetrics['module_size_entropy'] <- entropy(strainsInModules)
  
  ### 3. Similarity in module composition
  x=table(modules$module_ordered, modules$strain)
  x[x>0] <- 1
  moduleSimilarity <- round(1-vegdist(x, method = 'jaccard'),3)
  if (verbose){
    cat('Similarity in strain structure:');cat('\n')
    print(moduleSimilarity)
  }
  modularStructureMetrics['mean_strain_similarity'] <- mean(moduleSimilarity)
  
  # heatmap(as.matrix(moduleSimilarity), Rowv = NA, Colv = NA, symm = T)
  # plt=qplot(moduleSimilarity[lower.tri(moduleSimilarity)])+labs(x='Module similarity (Jaccard)', y='Count')
  # plt=ggDefaultTheme(plt)
  
  return(modularStructureMetrics)
}


# apply a given cutoff to a list of networks
applyCutoff <- function(N,cutoff){
  newlist <- vector('list',length(N))
  for (i in 1:length(N)){
    #cat('Network ');cat(i);cat('\t')
    x <- N[[i]]
    x[x<cutoff] <- 0
    newlist[[i]] <- x
  }
  return(newlist)
}


applyCutoff2 <- function(N,percentile){
  newlist <- vector('list',length(N))
  for (i in 1:length(N)){
    #cat('Network ');cat(i);cat('\t')
    x <- N[[i]]
    cutoff=quantile(x[x>0],percentile)
    x[x<cutoff] <- 0
    newlist[[i]] <- x
  }
  return(newlist)
}


binarize <- function(mat){
  m <- mat 
  m[which(m > 0)] <- 1
  return(m)
}


build_calendar <- function(num_years = 100, burnin=27000, year_to_start=10, plotit=F){# This is the number of years from start to end, including pre-burnin
  months_in_year <- rep(c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), each=30)
  calendar <- data.frame(running_day=seq(from = 1,to = 360*num_years,by=1),
                         year_sim=rep(1:num_years, each=360),
                         month_sim=rep(months_in_year,num_years),
                         day_sim=rep(1:30,num_years))
  calendar$layer <- ceiling((calendar$running_day-burnin)/30)
  calendar$burnin <- 'No'
  calendar$burnin[1:burnin] <- 'Yes'
  calendar$empirical_survey <- NA
  calendar$empirical_IRS <- NA
  
  # Set the dates corresponding to the empirical survey dates. year_to_start is
  # the year in the simulation (post-burnin) where simulated data is starting to
  # be matched to the empirical data.
  if (year_to_start < num_years){
    calendar$empirical_survey[calendar$year_sim==burnin/360+year_to_start & calendar$month_sim=='Oct'] <- 'S1'
    calendar$empirical_survey[calendar$year_sim==burnin/360+year_to_start+1 & calendar$month_sim=='Jun'] <- 'S2'
    calendar$empirical_survey[calendar$year_sim==burnin/360+year_to_start+2 & calendar$month_sim=='Jun'] <- 'S3'
    calendar$empirical_survey[calendar$year_sim==burnin/360+year_to_start+2 & calendar$month_sim=='Oct'] <- 'S4'
    calendar$empirical_survey[calendar$year_sim==burnin/360+year_to_start+3 & calendar$month_sim=='Oct'] <- 'S5'
    calendar$empirical_survey[calendar$year_sim==burnin/360+year_to_start+4 & calendar$month_sim=='Jun'] <- 'S6'
    irs_1_start <- extract_from_calendar(calendar, burnin, year_to_start+1, 'Oct',d=15)$running_day
    irs_1_end <- extract_from_calendar(calendar, burnin, year_to_start+1, 'Dec',d=30)$running_day
    calendar$empirical_IRS[calendar$running_day%in%irs_1_start:irs_1_end] <- 'IRS'
    irs_2_start <- extract_from_calendar(calendar, burnin, year_to_start+2, 'May',d=15)$running_day
    irs_2_end <- extract_from_calendar(calendar, burnin, year_to_start+2, 'Jul',d=30)$running_day
    calendar$empirical_IRS[calendar$running_day%in%irs_2_start:irs_2_end] <- 'IRS'
    irs_3_start <- extract_from_calendar(calendar, burnin, year_to_start+2, 'Dec',d=15)$running_day
    irs_3_end <- extract_from_calendar(calendar, burnin, year_to_start+3, 'Feb',d=30)$running_day
    calendar$empirical_IRS[calendar$running_day%in%irs_3_start:irs_3_end] <- 'IRS'
  
    if (plotit){
      temp <- calendar %>% filter(!is.na(empirical_survey) | !is.na(empirical_IRS))
      p <- ggplot(temp, aes(running_day))
      p <- p+geom_point(aes(y=empirical_survey),color='blue')
      p <- p+geom_point(aes(y=empirical_IRS),color='red')
      p <- p+geom_vline(xintercept = c(irs_2_start,irs_2_end))
      print(p)
    }
  }
  
  return(calendar)
}


# Find a certain day by year (from burnin) and month.
extract_from_calendar <- function(cal, burnin=27000, y, m, d=NULL){
  if (is.null(d)){
    tmp <- (burnin/360+y)
    x <- subset(cal, year_sim==tmp & month_sim==m)
    return(x)
  } else {
    tmp <- (burnin/360+y)
    x <- subset(cal, year_sim==tmp & month_sim==m & day_sim==d)
    return(x)
  }
}

buildInfomapIntervention <- function(phase){
  if (phase=='B'){layersIntervention <- 1:layerInterventionStart}
  if (phase=='D'){layersIntervention <- (layerInterventionStart+1):layerInterventionEnd}
  if (phase=='A'){layersIntervention <- (layerInterventionEnd+1):numberOfLayers}
  theLayers <- layersFull[layersIntervention]
  nodeLabel <- sort(unique(unlist(lapply(theLayers,rownames))))
  nodeList <- data.frame(nodeID=1:length(nodeLabel), nodeLabel)
  # Build interlayer edges
  ile_edgelist <- list()
  for (t in 1:(length(theLayers)-1)){
    cat('Building interlayer edges for layers: ',t,'-->',t+1,'\n')
    strainCopies_t <- rownames(theLayers[[t]])
    strainCopies_t1 <- rownames(theLayers[[t+1]])
   
    # if intervention was too strong and there are less than 5 strains in a particular layer, it should be skipped
    if (length(strainCopies_t1)<5){
      mprint(paste('Need at least 5 strains in layer',t+1,'!!! skipping interlayer edges',t,'-->',t+1,'...'))
      next
    }
    if (length(strainCopies_t)<5){
      mprint(paste('Need at least 5 strains in layer',t,'!!! skipping interlayer edges',t,'-->',t+1,'...'))
      next
    }
    m <- similarityMatrix[splitText(strainCopies_t, after = F),splitText(strainCopies_t1, after = F)]
    rownames(m) <- strainCopies_t
    colnames(m) <- strainCopies_t1
    NZ <- which(m!=0, arr.ind = T) # non-zero elements
    edges_interlayer <- data.frame(layer_s=rep(t,nrow(NZ)),
                                   node_s=rep(NA,nrow(NZ)),
                                   layer_t=rep(t+1,nrow(NZ)),
                                   node_t=rep(NA,nrow(NZ)),
                                   w=rep(NA,nrow(NZ)))
    edges_interlayer$node_s <- rownames(m)[NZ[,1]]
    edges_interlayer$node_t <- colnames(m)[NZ[,2]]
    for (n in 1:nrow(edges_interlayer)){
      # Here use the indices instead of strain names because there are bugs in the
      # sparse matrix package which result in a value of 0 if the name in the row 
      # and column are the same (i.e., the diagonal). If a regular matrix is used
      # it doesnt really matter if using the indices or the names.
      edges_interlayer[n,'w'] <- m[which(rownames(m)==edges_interlayer[n,'node_s']),which(colnames(m)==edges_interlayer[n,'node_t'])]
      # This is to use the names of the strains
      # edges_interlayer[n,'w'] <- m[as.character(edges_interlayer[n,'node_s']),as.character(edges_interlayer[n,'node_t'])]
    }
    edges_interlayer$node_s <- nodeList$nodeID[match(edges_interlayer$node_s,nodeList$nodeLabel)]
    edges_interlayer$node_t <- nodeList$nodeID[match(edges_interlayer$node_t,nodeList$nodeLabel)]
    ile_edgelist[[t]] <- edges_interlayer
  }
  inter_edges <- do.call(rbind.data.frame,ile_edgelist)
  rm(ile_edgelist)
  intra_edges <- infomap_makeIntralayerEdges(theLayers,nodeList)
  ## Write file for infomap
  mprint('Writing Infomap files')
  file <- paste(filenameBase,'_Infomap_multilayer_',phase,'.txt',sep='')
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


calculateFeatures <- function(x){
  require(igraph) 

  giant.component <- function(g) { 
    cl <- clusters(g) 
    induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
  }

  f_01_averageLocalClusteringCoeff <- function(g,GC=F){
    if (GC) {
      gc <- giant.component(g)
      return(mean(transitivity(gc, type = 'local'), na.rm = T))
    } else {
      return(mean(transitivity(g, type = 'local'), na.rm = T))
    }
  }
  
  f_02_averageLocalClusteringCoeffWeighted <- function(g,GC=F){
    if (GC) {
      gc <- giant.component(g)
      return(mean(transitivity(gc, type = 'barrat'), na.rm = T))
    } else {
      return(mean(transitivity(g, type = 'barrat'), na.rm = T))
    }
  }
  
  
  f_03_globalClusteringCoeff <- function(g, GC=F){
    if (GC) {
      gc <- giant.component(g)
      return(transitivity(gc, type = 'global'))
    } else {
      return(transitivity(g, type = 'global'))
    }
  }
  
  f_04_gdensity <- function(g){
    return(graph.density(g))
  }
  
  f_05_proportionSingletons <- function(g){
    sum(igraph::degree(g)==0)/length(V(g))
  }
  
  f_06_proportionEndpoints <- function(g){
    sum(igraph::degree(g)==1)/length(V(g))
  }
  
  f_07_meanDegree <- function(g){
    mean(igraph::degree(g))
  }
  
  f_08_meanDegreeNotSingletons <- function(g){
    mean(igraph::degree(g)[igraph::degree(g)!=0])
  }
  
  f_09_meanStrength <- function(g){
    mean(igraph::strength(g))
  }
  
  f_10_meanStrengthNotSingletons <- function(g){
    mean(igraph::strength(g)[igraph::degree(g)!=0])
  }
  
  f_11_entropyDegreeDistribution <-  function(g,verbose=F){
    y=igraph::degree(g)
    freq=prop.table(table(y))
    if (verbose){print(freq)}
    -sum(freq * log(freq, base = 2))
  }
  
  f_12_ratioComponents <- function(g) { 
    cl <- clusters(g) 
    cl$no/length(V(g))
  }
  
 f_13_averageComponentSize <- function(g) { 
    cl <- clusters(g) 
    mean(cl$csize)
  }
  
  f_14_entropyComponentSize <-  function(g,verbose=F){
    cl <- clusters(g)
    y=cl$csize
    freq=prop.table(table(y))
    if (verbose){print(freq)}
    -sum(freq * log(freq, base = 2))
  }
  
  f_15_giantConnectedRatio <- function(g){
    cl <- clusters(g)
    max(cl$csize)/length(V(g))
  }
  
  f_16_meanEccentricity <- function(g){
    mean(eccentricity(g))
  }
  
  f_17_gdiameter <- function(g){
    return(diameter(g))
  }
  
  f_18_meanDiameterComponents <- function(g){
    cl <- clusters(g)
    d <- 0
    for (m in 1:(cl$no)){
      d <- d + diameter(induced.subgraph(g, which(cl$membership == m)))
    }
    return(d/cl$no)
  }
  
  f_19_globalEfficiency <- function(g){
    d_ij <- shortest.paths(g)
    d_ij <- d_ij[lower.tri(d_ij)]
    d_ij <- d_ij[!is.infinite(d_ij)]
    N=length(V(g))
    1/(N*(N-1))*sum(1/d_ij)
  }
  
  f_20_averageClosenessCentrality <- function(g){
    mean(igraph::closeness(g, weights = NULL))
  }
  

  motifsProportion <- function(g){ # Calculate the proportion of each of the 16 motifs out of the total motifs found
    motifs <- graph.motifs(g, size = 3)
    motifs.prop <- motifs/sum(motifs, na.rm = T)
    #names(motifs.prop) <- paste('motif',1:16,sep='')
    return(motifs.prop[-c(1,2,4,12)]) # Motifs 1,2,4 are constantly (across all cutoffs) NA and 12 is constantly 0 (I tested it)
  }
  
  
  
  # Main function starts here
  
  ## Transfomr to an igraph object
  if (!is.igraph(x)){
    g <- graph.adjacency(x, mode = 'directed', weighted = T, diag = F)
  }
  
  ## Calculate properties
  featureVector <- vector(length=32)
  # Diagnostics of transitivity
  featureVector[1] <- f_01_averageLocalClusteringCoeff(g,F)    # Clustering coefficient averaged across all nodes
  featureVector[2] <- f_02_averageLocalClusteringCoeffWeighted(g,F)    # Barrat's clustering coefficient averaged across all nodes
  featureVector[3] <- f_03_globalClusteringCoeff(g,F)
  # Diagnostics of degree/sterngth
  featureVector[4] <- f_04_gdensity(g)                    # Graph density
  featureVector[5] <- f_05_proportionSingletons(g)             # Proportion of nodes with degree 0 of all the nodes
  featureVector[6] <- f_06_proportionEndpoints(g)              # Proportion of nodes with degree 1 of all the nodes
  featureVector[7] <- f_07_meanDegree(g)                       # Average degree 
  featureVector[8] <- f_08_meanDegreeNotSingletons(g)
  featureVector[9] <- f_09_meanStrength(g)                     # Average strength
  featureVector[10] <- f_10_meanStrengthNotSingletons(g)
  featureVector[11] <- f_11_entropyDegreeDistribution(g)       # Average measurement of the heterogeneity of the network. See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  featureVector[12] <- f_12_ratioComponents(g)                 # Number of components relative to networks size
  featureVector[13] <- f_13_averageComponentSize(g)
  featureVector[14] <- f_14_entropyComponentSize(g)
  featureVector[15] <- f_15_giantConnectedRatio(g)             # Proportion of nodes in the giant component
  #Diagnostics of shortest-paths
  featureVector[16] <- f_16_meanEccentricity(g)                 # Eccentricity is the maximum shortest distance from a node to all other nodes. This is averaged across all nodes
  featureVector[17] <- f_17_gdiameter(g)                        # length of the longest geodesic
  featureVector[18] <- f_18_meanDiameterComponents(g)
  featureVector[19] <- f_19_globalEfficiency(g)                # See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  featureVector[20] <- f_20_averageClosenessCentrality(g)      # 
  # Diagnostics of motifs
  featureVector[21:32] <- motifsProportion(g)
  return(featureVector)
}




calculate_Pvalue <- function (observed,null,mode,plotit=F,min.values=100){
  if (length(null)<min.values){stop('there must be at least 100 null values. aborting')}
  if (length(null)<1000){warning('It is better to have at least 1000 null values')}
  if (mode=='h'){
    cutoffHigh <- sort(null)[ceiling(0.95*length(null))]
    print(sum(null>observed)/length(null))
    if(plotit==T) {ggHistogram(null)+ geom_vline(xintercept=observed,col='red')+geom_vline(xintercept=c(cutoffHigh),col='black', linetype='dashed')}
    return(observed > cutoffHigh)
  }
  if (mode=='l'){
    cutoffLow <- sort(null)[floor(0.05*length(null))]
    print(sum(null<observed)/length(null))
    if(plotit==T) {ggHistogram(null)+ geom_vline(xintercept=observed,col='red')+geom_vline(xintercept=c(cutoffLow),col='black', linetype='dashed')}
    return(observed < cutoffLow)
  }   
  if (mode=='2'){
    cutoffLow <- sort(null)[floor(0.025*length(null))]
    cutoffHigh <- sort(null)[ceiling(0.975*length(null))]
    cat(sum(null<observed)/length(null));cat(' (low) ');cat(sum(null>observed)/length(null));cat(' (high)')
    if(plotit==T) {ggHistogram(null)+ geom_vline(xintercept=observed,col='red')+geom_vline(xintercept=c(cutoffLow,cutoffHigh),col='black', linetype='dashed')}
    return(observed < cutoffLow | observed > cutoffHigh)
  }
}
#null=rnorm(500)
#observed=2
#calculate_Pvalue(observed,null,'2',T)


CorenessLayout <- function(g) {
  coreness <- graph.coreness(g);
  xy <- array(NA, dim=c(length(coreness), 2));
  
  shells <- sort(unique(coreness));
  for(shell in shells) {
    v <- 1 - ((shell-1) / max(shells));
    nodes_in_shell <- sum(coreness==shell);
    angles <- seq(0,360,(360/nodes_in_shell));
    angles <- angles[-length(angles)]; # remove last element
    xy[coreness==shell, 1] <- sin(angles) * v;
    xy[coreness==shell, 2] <- cos(angles) * v;
  }
  return(xy);
}

fastJaccard <- function(M){
  m_1=matrix(1, nrow(M),ncol(M))
  A <- t(M)%*%M
  B=t(m_1)%*%M-A
  C=t(B)
  J=A/(A+B+C)
  J
}

findLongestModuleSequence <- function(modulesequence){
  temp <- cumsum(c(1, diff(modulesequence) - 1))
  temp2 <- rle(temp)
  modulesequence[which(temp == with(temp2, values[which.max(lengths)]))]
}


findSegments <- function(x){
  hit <- which(is.na(x))
  n <- length(hit)
  ind <- which(hit[-1] - hit[-n] > 1)
  starts <- c(hit[1], hit[ ind+1 ])
  ends <- c(hit[ ind ], hit[n])
  cbind(starts,ends)
}


gdensity <- function(m, bipartite=F, directed=T){
  if (bipartite==T){
    return(sum(m!=0)/(nrow(m)*ncol(m)))
  }
  if (bipartite==F & directed==F){ # d=2*E/(V*(V-1)) -- this EXCLUDES THE DIAGONAL
    no_edges <- sum(m[lower.tri(m)]!=0)+sum(m[upper.tri(m)]!=0) # This is to excldue the diagonal and in undirected networks equals 2E
    no_possible_edges <- nrow(m)*(nrow(m)-1)
    return(no_edges/no_possible_edges)
  }
  if (bipartite==F & directed==T){ # d=E/(V*(V-1)) -- this EXCLUDES THE DIAGONAL
    no_edges <- sum(m[lower.tri(m)]!=0)+sum(m[upper.tri(m)]!=0) # This is to excldue the diagonal
    no_possible_edges <- nrow(m)*(nrow(m)-1)
    return(no_edges/no_possible_edges)
  }
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


ggDefaultTheme <- function(plt,plotit=T,textsize=20,...){
  plt=plt+theme_bw()+theme(strip.background=element_blank(),
                           
                           axis.line = element_line(colour = "black"),
                           # axis.ticks.y = element_blank(),
                           axis.text = element_text(size=textsize),
                           axis.title = element_text(size=textsize),
                           strip.text = element_text(size=textsize),
                           #legend.position="none",
                           #panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_rect(color='black'),
                           panel.background = element_blank(),...)
  if(plotit){plot(plt)}
  return(plt)
}


ggDefaultThemeClassic <- function(plt,plotit=T,textsize=20,...){
  plt=plt+theme_classic()+theme(strip.background=element_blank(),
                                # axis.line = element_line(colour = "black"),
                                # axis.ticks.y = element_blank(),
                                axis.text = element_text(size=textsize),
                                axis.title = element_text(size=textsize),
                                strip.text = element_text(size=textsize),...)
  if(plotit){plot(plt)}
  return(plt)
}


ggHistogram <- function(x, col='navy', bw=NULL, xlab='', title='', add_vline=NULL, tail_percentile=NULL,...) {
  require(ggplot2)
  d <- data.frame(xID=1:length(x), x)
  
  if (is.null(tail_percentile)){
    if (is.null(add_vline)){
      plt=ggplot(d, aes(x))+geom_histogram(fill=col, binwidth = bw)+labs(x=xlab, title=title)
    } else {
      plt=ggplot(d, aes(x))+geom_histogram(fill=col, binwidth = bw)+labs(x=xlab, title=title)+geom_vline(xintercept = add_vline)
    }
  }
  if (!is.null(tail_percentile)){
    print(paste('Value at ',tail_percentile,' percentile is: ',quantile(d$x, probs = tail_percentile),sep=''))
    plt=ggplot(d, aes(x))+geom_histogram(fill=col, binwidth = bw)+labs(x=xlab, title=title)+geom_vline(xintercept = quantile(d$x, probs = tail_percentile))
  }
  return(plt)
}


ggTimeSeries <- function(list_of_vars, vnames=NULL,cols=NULL,...){
  # d <- data.frame(l=1:length(x), x)
  if (class(list_of_vars)!='list'){stop('the series must be a list')}
  d <- as.data.frame(do.call('cbind', list_of_vars))
  if(!is.null(vnames)){colnames(d) <- vnames}
  d$t <- 1:nrow(d)
  d <- melt(d, id.vars='t')
  if (is.null(cols)){plt <- ggplot(d, aes(t, value, color=variable))+geom_line()
  } else
    { plt <- ggplot(d, aes(t, value, color=variable))+geom_line()+scale_color_manual(values=cols)
    }
  return(plt)
}


ggplotToBrowser <- function(p, w=35, h=16.875) {
  ggsave(filename = tf_img <- tempfile(fileext = ".svg"), plot = p, width=w, height=h, units = 'cm')
  html <- sprintf('<html><body><img src="%s"></body></html>', paste0("file:///", tf_img))
  cat(html, file = tf_html <- tempfile(fileext = ".html"))
  if(Sys.info()[1]=="Linux"){
    options(browser="google-chrome")
    browseURL(tf_html)  
  } else {
    system(sprintf("open %s", tf_html))  
  }
}


h2 <- function(x,n=10,m=10){
  if(n>nrow(x)){n <- nrow(x)}
  if(m>ncol(x)){m <- ncol(x)}
  print(x[1:n,1:m])
  invisible()
}

ht <- function(x,...){
  print(head(x,...))
  cat('------------------------------------\n')
  print(tail(x,...))
}

hmod_organizeModulesByLayer <- function(df, mod.var){
  print('Re-organizing modules by layer...')
  df=df[with(df, order(layer,get(mod.var))),]
  x=with(df, table(get(mod.var))[unique(get(mod.var))])
  df$module_ordered <- NA
  for (i in 1:nrow(df)){
    df[i,'module_ordered'] <- which(names(x)==df[i,mod.var])
  }
  return(df)
}


importMatrix <- function(file,sepchr=','){
  x <- read.table(file, header = T, sep = sepchr)
  colnames(x) <- rownames(x)
  x <- data.matrix(x)
  print(x[1:5,1:5])
  return(x)
}


infomap_makeIntralayerEdges <- function(layerList, nodeList) {
  listOfEdgesLists <- list()
  for (l in 1:length(layerList)){
    print(paste('[',Sys.time(), '] creating edge list of layer ',l,' for Infomap...',sep=''))
    
    # If there are no nodes in the layer
    if (length(layerList[[l]])==0){
      mprint(paste('No strains in layer',l,'!!! skipping intralayer edges'))
      next
    }
    if(sum(layerList[[l]])==0){
      mprint(paste('No edges in layer',l,'!!! skipping intralayer edges'))
      next
    }
    listOfEdgesLists[[l]] <- matrix.to.coordinateList2(layerList[[l]])
  }
  all(sapply(listOfEdgesLists, function(x) all(x$i%in% nodeList$nodeLabel))==T)
  all(sapply(listOfEdgesLists, function(x) all(x$j%in% nodeList$nodeLabel))==T)
  for (l in 1:length(listOfEdgesLists)){
    if(is.null(listOfEdgesLists[[l]])){
      listOfEdgesLists[[l]] <- data.frame(layer_s=NA,node_s=NA,layer_t=NA,node_t=NA)
      next
    }
    listOfEdgesLists[[l]]$layer_s <- l
    listOfEdgesLists[[l]]$node_s <- nodeList$nodeID[match(listOfEdgesLists[[l]]$i,nodeList$nodeLabel)]
    listOfEdgesLists[[l]]$layer_t <- l
    listOfEdgesLists[[l]]$node_t <- nodeList$nodeID[match(listOfEdgesLists[[l]]$j,nodeList$nodeLabel)]
    listOfEdgesLists[[l]] <- listOfEdgesLists[[l]][,c(4:7,3)]
  }
  
  edges_intralayer <- plyr::ldply(listOfEdgesLists,data.frame)
  edges_intralayer <- edges_intralayer[complete.cases(edges_intralayer),]
  return(edges_intralayer)
}


infomap_makeInterlayerEdges <- function(ileMatrix, nodeList){
  edges_interlayer <- data.frame(layer_s=as.numeric(),
                                 node_s=as.numeric(),
                                 layer_t=as.numeric(),
                                 node_t=as.numeric(),
                                 w=as.numeric())
  NZidx <- which(ileMatrix!=0,arr.ind = T) # indices of non-zero matrix cells
  pb <- txtProgressBar(min = 1, max = nrow(NZidx), style = 3)
  for (i in 1:nrow(NZidx)){
    setTxtProgressBar(pb, i)
    edges_interlayer[i, 'node_s'] <- rownames(ileMatrix)[NZidx[i,1]]
    edges_interlayer[i, 'node_t'] <- rownames(ileMatrix)[NZidx[i,1]]
    edges_interlayer[i, 'layer_s'] <- NZidx[i,2]
    edges_interlayer[i, 'layer_t'] <- NZidx[i,2]+1
    edges_interlayer[i, 'w'] <- ileMatrix[NZidx[i,1],NZidx[i,2]]
  }
  edges_interlayer$node_s <- nodeList$nodeID[match(edges_interlayer$node_s,nodeList$nodeLabel)]
  edges_interlayer$node_t <- nodeList$nodeID[match(edges_interlayer$node_t,nodeList$nodeLabel)]
  return(edges_interlayer)
}


infomap_readTreeFile <- function(file, reorganize_modules=T, max_layers=372, remove_buggy_instances=T){
  require(splitstackshape)
  lines <- readLines(file)
  #lines <- readLines('Infomap_linux/output/S11_S0_w30_300_300_0.1_expanded.tree');length(lines)
  cat(lines[1]);cat('\n')
  x=read.table(file, skip = 2, stringsAsFactors = F)
  modules <- data.frame(module=rep(NA,nrow(x)),
                        strain=rep(NA,nrow(x)),
                        layer=rep(NA,nrow(x)),
                        flow=rep(NA,nrow(x)), stringsAsFactors = F)
  
  modules$path <- x[,1]
  x.module <- x[,1]
  x.module <- cSplit(as.data.table(x.module),'x.module',':')
  x.module <- as.data.frame(x.module)
  modules$module <- x.module[,1]
  cat(nrow(x),'state nodes','in',paste(max(modules$module),'modules, organized in',length(x.module),'levels...'));cat('\n')
  modules$flow <- x[,2]
  x.strain <- x[,3]
  x.strain <- read.table(text = x.strain, sep = "|", colClasses = "character", stringsAsFactors = F, strip.white = T)
  modules$strain <- x.strain$V1
  modules$layer <- as.numeric(x$V4)
  
  # There is a bug in Infomap that assigns nodes to layers which do not exist (with higher number than the existing layers).
  
  if(remove_buggy_instances){
    buggy <- modules[modules$layer>max_layers,]
    if (nrow(buggy)>=1){
      cat('\n')
      print('---------------------------------')
      print('Some buggy instances encountered!')
      print(paste('removed ', nrow(buggy),' instances which were assigned to layers which do not exist.',sep=''))
      print(paste('total flow of removed instance: ',sum(buggy$flow)))
      print('Buggy instances written to file')
      write.table(buggy, paste(str_sub(file, 1, str_locate(file, 'output/'))[2],'buggy_instances_infomap.txt',sep=''))
      modules <- modules[modules$layer<=max_layers,]
    }
  }
  
  if(reorganize_modules){  # organize the names of modules to be consecutive
    cat('Re-organizing modules by layer...');cat('\t')
    modules=modules[with(modules, order(layer,module)),]
    x=table(modules$module)[unique(modules$module)]
    modules$module_ordered <- NA
    for (i in 1:nrow(modules)){
      modules[i,'module_ordered'] <- which(names(x)==modules[i,'module'])
    }
    modules <- modules[,-1]
    modules <- modules[,c('module_ordered','strain','layer','flow','path')]
    names(modules)[1] <- 'module'
  }
  print('Done!')
  return(modules)
}


infomap_readTreeFile_single <- function(file, layer=NULL){
  require(splitstackshape)
  lines <- readLines(file)
  #lines <- readLines('Infomap_linux/output/S11_S0_w30_300_300_0.1_expanded.tree');length(lines)
  cat(lines[1]);cat('\n')
  x=read.table(file, skip = 2, stringsAsFactors = F)
  modules <- data.frame(module=rep(NA,nrow(x)),
                        nodeLabel=rep(NA,nrow(x)),
                        flow=rep(NA,nrow(x)), stringsAsFactors = F)
  
  modules$path <- x[,1]
  x.module <- x[,1]
  x.module <- cSplit(as.data.table(x.module),'x.module',':')
  x.module <- as.data.frame(x.module)
  modules$module <- x.module[,1]
  cat(nrow(x),'state nodes','in',paste(max(modules$module),'modules, organized in',length(x.module),'levels...'));cat('\t')
  modules$flow <- x[,2]
  modules$nodeLabel <- x[,3]
  if (!is.null(layer)){modules$layer <- layer}
  # There is a bug in Infomap that assigns nodes to layers which do not exist (with higher number than the existing layers).
  print('Done!')
  return(modules)
}

is.notzero <- function(x){
  return(x[x!=0])
}




linkageDiseq <- function(m, output='D_prime', make_symmetric=T){
  # Example matrix:
  # m <- matrix(c(1,1,0,0,1,0,1,1,1,1,
  #               1,0,1,1,1,0,1,1,0,1,
  #               0,0,1,1,0,1,1,0,1,1,
  #               0,1,1,0,1,1,1,1,0,0,
  #               0,1,1,0,1,0,0,1,0,0,
  #               1,0,1,0,1,0,1,1,0,0),
  #             nrow=6, ncol=10, byrow = T, dimnames = list(1:6,letters[1:10]))
  
  # in m rows are strains cols are alleles or genes
  
  # Definitions:
  # p_i   is the probability of finding allele i in any given matrix (the number of times an allele occurs in the strain population divided by the number of strains)
  # P_ij  is a matrix in which an entry is the independent probability of encountering allele i and allele j. This is the pairwise multiplication of elements of p_i
  # J_ij  is a matrix in which an entry is the joint probability of encountering BOTH allele i and allele j IN THE SAME STRAIN. This is the number of strains where both i and j occur divided by the number of strains
  # K_ij  is a matrix whose values are (1-p_i)*(1-p_j), This is used to define Dmax. Entries are calculated as the pairwise multiplication of 1-p_i
  # L_ij  is a matrix whose values are p_i*(1-p_j). This is used to define Dmax.
  # M_ij  is a matrix whose values are p_j*(1-p_i). This is used to define Dmax.
  # D_prime is defined as: if(LD<0){Dmax <- min(pA*pB,(1-pA)*(1-pB))} else Dmax <- min(pA*(1-pB),pB*(1-pA))
  # r is: r=D/sqrt(pA(1-pA)pB(1-pB)) == D/sqrt(P_ij*K_ij)
  
  m[m>0] <- 1
  
  p_i <- colSums(m)/nrow(m)
  P_ij <- K_ij <- L_ij <- M_ij <- Dmin <- matrix(0, nrow=ncol(m), ncol=ncol(m), dimnames = list(colnames(m),colnames(m)))
  P_ij[lower.tri(P_ij)] <- combn(p_i, m = 2, FUN = prod)
  J_ij <- crossprod(m)/nrow(m)
  J_ij[upper.tri(J_ij)] <- 0; diag(J_ij) <- 0
  K_ij[lower.tri(K_ij)] <- combn(p_i, m = 2, FUN = function(x) prod(1-x))
  x=combn(p_i, m = 2)
  L_ij[lower.tri(L_ij)] <-x[1,]*(1-x[2,]) # This is pA(1-pB)
  M_ij[lower.tri(M_ij)] <-(1-x[1,])*x[2,] # This is pB(1-pA)
  
  D <- J_ij-P_ij
  
  # Return results
  if (output=='D') {
    if(make_symmetric){D[upper.tri(D)] <- t(D)[upper.tri(D)]}
    return(D)
  }
  if (output=='D_prime') {
    Dmin[,] <- 1
    Dmin1 <- pmax(-1*P_ij, -1*J_ij) # When D<0
    Dmin2 <- pmin(L_ij,M_ij) # When D>0
    Dmin[D<0] <- Dmin1[D<0]
    Dmin[D>0] <- Dmin2[D>0]
    D_prime <- D/Dmin
    D_prime[is.infinite(D_prime)] <- 0
    if(make_symmetric){D_prime[upper.tri(D_prime)] <- t(D_prime)[upper.tri(D_prime)]}
    return(D_prime)
  }
  if (output=='r') {
    r <- D/sqrt(P_ij*K_ij)
    if(make_symmetric){r[upper.tri(r)] <- t(r)[upper.tri(r)]}
    return(r)
  }
}


loadExperiment <- function(project,layersToInclude,cutoffPercentile){
  fileSignature <- paste(project,'_',layersToInclude,'_cutoffP_',cutoffPercentile,sep='')
  print('Loading layers...')
  Layers <- list()
  for (r in 1:layersToInclude){
    Layers[[r]] <- data.matrix(read.table(paste(project,'/',fileSignature,'_layer_',r,'.txt',sep='')))
    colnames(Layers[[r]]) <- rownames(Layers[[r]])
  }
  print('Loading prevalence data...')
  sFrequency <- data.matrix(read.table(paste(project,'/',fileSignature,'_sFrequency.txt',sep='')))
  colnames(sFrequency) <- 1:layersToInclude
  sFrequency_normalized <- data.matrix(read.table(paste(project,'/',fileSignature,'_sFrequency_normalized.txt',sep='')))
  colnames(sFrequency_normalized) <- 1:layersToInclude
  sOccurrence <- data.matrix(read.table(paste(project,'/',fileSignature,'_sOccurrence.txt',sep='')))
  colnames(sOccurrence) <- 1:layersToInclude
  sPrevalence <- data.matrix(read.table(paste(project,'/',fileSignature,'_sPrevalence.txt',sep='')))
  colnames(sPrevalence) <- 1:layersToInclude
  strain_allele <- data.matrix(read.table(paste(project,'/',fileSignature,'_strain_allele.txt',sep='')))
  colnames(strain_allele) <- str_replace_all(colnames(strain_allele),pattern = 'X','')
  interlayerEdges <- data.matrix(read.table(paste(project,'/',fileSignature,'_interlayerEdges.txt',sep='')))
  colnames(interlayerEdges) <- str_replace_all(colnames(interlayerEdges),pattern = 'X','')
  colnames(interlayerEdges) <- str_replace_all(colnames(interlayerEdges),pattern = '\\...','-->')
  lInfo <- read.table(paste(project,'/',fileSignature,'_layerInfo.txt',sep=''))
  
  infomapFile <- paste(project,'/',fileSignature,'_Infomap','.txt',sep='')
  
  ifelse(file.exists(infomapFile),print('Infomap file ok'),print('Infomap file does not exist'))
  nodes <- data.frame(nodeID=1:nrow(sPrevalence), nodeLabel=rownames(sPrevalence), stringsAsFactors = F)
  print('Creating igraph objects')
  ### Create a list of graphs for the layers.
  graphList <- vector('list',length(Layers))
  for (i in 1:length(Layers)){
    # cat(i);cat('\t')
    x <- Layers[[i]]
    g <- graph.adjacency(x, mode='directed', weighted = T, diag = F)
    graphList[[i]] <- g
  }
  names(graphList) <- paste('L',str_pad(1:length(Layers),width = 2,side = 'l',pad = '0'),sep='_')
  print('Done!')
  return(list(Layers=Layers,
              sFrequency=sFrequency,
              sFrequency_normalized=sFrequency_normalized,
              sOccurrence=sOccurrence,
              sPrevalence=sPrevalence,
              interlayerEdges=interlayerEdges,
              strain_allele=strain_allele,
              lInfo=lInfo,
              nodes=nodes,
              graphList=graphList,
              infomapFile=infomapFile))
}

loadExperiments_GoogleSheets <- function(workBookName='mtn_experiments',sheetID=2){
  require(googlesheets)
  GS <- gs_title(workBookName)
  col_types <- GS %>% gs_read(ws=1, col_names=T)
  col_t <- unname(as.list(col_types[1,]))
  experiments <- GS %>% gs_read(ws=sheetID, col_names=T, col_types=col_t)
  print(experiments)
  return(experiments)
}

makeCooccurMatrix <- function(m){
  # row is isolate, column is var gene
  corMat<-rcorr(m, type="pearson")
  corMat_pearsonCoeff<-corMat[[1]]
  corMat_pearsonCoeff[corMat[[3]]>0.05] <- 0 # Only take significant correlations
  return(corMat_pearsonCoeff)
}


makelayerList <- function(layerList,frequencyMatrix,ocurrenceMatrix,prevalenceMatrix,layerInfo,layersToInclude,cutoffPercentile,min.layers){
  # Possible to limit the data by setting low value for layersToInclude.
  timeSlices <- 1:layersToInclude
  layerList <- layerList[timeSlices]
  layerInfo <- layerInfo[timeSlices,]
  
  # Remove transient strains
  layerList <- removeTransientStrains(layerList,min.layers)
  
  # Set the diagonal to 0, instead of 1
  for (l in 1:length(layerList)){diag(layerList[[l]])=0}
  
  # Set cutoff
  layerList <- applyCutoff2(layerList,cutoffPercentile)
  
  # Name the layers
  names(layerList) <- paste('L',str_pad(1:length(layerList),width = 2,side = 'l',pad = '0'),sep='_')
  # match the strains
  strains <- unique(unlist(lapply(layerList,rownames)))
  frequencyMatrix <- frequencyMatrix[strains,1:layersToInclude]
  ocurrenceMatrix <- ocurrenceMatrix[strains,1:layersToInclude]
  prevalenceMatrix <- prevalenceMatrix[strains,1:layersToInclude]
  # Gather information on layers that can be gathered only after removing transient strains
  layerInfo$strains_persistent=sapply(layerList, nrow)
  layerInfo$edges=sapply(layerList, function(x) length(x[x!=0]))
  layerInfo$density=sapply(layerList, gdensity)
  layerInfo$layer <- timeSlices
  layerInfo$meanSimilarity <- sapply(layerList, function(x) mean(x[x!=0]))
  # Some checks. Anything wrong?
  ### Layers:
  print('Checking layers...  ')
  # Any rows/cols summing to zero?
  if (any(unlist(sapply(layerList, rowSums))==0)) {print('!!! rows sum to zero')}
  if (any(unlist(sapply(layerList, colSums))==0)) {print('!!! cols sum to zero')}
  # Check to see if any of the layers does not contain a network. The time window
  # should be adjusted in a way that this does not occur
  if(any(sapply(layerList, is.null)==T)) {print('!!! At least one layer contains no network')}
  if(is.null(names(layerList))){print('!!! layers are unlabeled')}
  ### Strain matrices
  print('Checking strain matrices... ')
  cat('strain occurence dim: ');print(dim(ocurrenceMatrix))
  cat('strain frequency dim: ');print(dim(frequencyMatrix))
  cat('strain prevalence dim: ');print(dim(prevalenceMatrix))
  if(any(colSums(ocurrenceMatrix)==0)){print('!!! ocurrenceMatrix: 1 or more columns sum to 0')}
  if(any(colSums(prevalenceMatrix, na.rm = T)==0)){print('!!! prevalenceMatrix: 1 or more columns sum to 0')}
  if(!all(rownames(ocurrenceMatrix)==rownames(prevalenceMatrix))){print('!!! row names of ocurrenceMatrix and prevalenceMatrix do not match')}
  if(!all(rownames(ocurrenceMatrix)==rownames(frequencyMatrix))){print('!!! row names of ocurrenceMatrix and frequencyMatrix do not match')}
  if(!all(rownames(prevalenceMatrix)==rownames(frequencyMatrix))){print('!!! row names of prevalenceMatrix and frequencyMatrix do not match')}
  
  return(list(Layers=layerList,layerInfo=layerInfo,strainFrequency=frequencyMatrix,strainOccurrence=ocurrenceMatrix,strainPrevalence=prevalenceMatrix))
}


makeInterlayerEdges <- function (inputMatrix,completeMissingStrains=T){
  
  if (completeMissingStrains){
    # Complete the input matrix: Because sampling misses some strains which do
    # occur in a layer, I complete these myself (see variable strainOccurrence).
    # However, it is impossible to calculate the prevalecne/frequency of these strains directly.
    # Therefore, I assume that in the instances (layers) where they are missing the
    # prevalence/frequency of a strain will be the mean between two observed occurrences.
    print('Completing missing strains...')
    pb <- txtProgressBar(min = 1, max = nrow(inputMatrix), style = 3)
    for (s in 1:nrow(inputMatrix)){
      setTxtProgressBar(pb, s)
      na_segments <- findSegments(inputMatrix[s,])
      #if there are no NA values in the row
      if (any(is.na(na_segments))) next
      for (i in 1:nrow(na_segments)){
        # if the first column has an NA (open segment) it is impossible to calculate an average
        if (na_segments[i,1]==1) {
          inputMatrix[s,na_segments[i,1]:na_segments[i,2]] <- inputMatrix[s,na_segments[i,2]+1]
          next
        } 
        # if there is a segment of NAs that reaches the last layer take the value of the layer before
        if (na_segments[i,2]==ncol(inputMatrix)) {  
            inputMatrix[s,na_segments[i,1]:na_segments[i,2]] <- inputMatrix[s,na_segments[i,1]-1]
            next
        }
        # otherwise proceed normally for the segment (calculate average)
        x=na_segments[i,1]-1
        y=na_segments[i,2]+1
        avg <- (inputMatrix[s,x]+inputMatrix[s,y])/2
        inputMatrix[s,na_segments[i,1]:na_segments[i,2]] <- avg
      }
    }
    if(any(is.na(inputMatrix))){stop('There are still NA values in the prevalecne matrix')}
  }
  
  # Now, calculate relative change in strain frequency/prevalence between two consecutive layers
  relativeChange <- matrix(0,nrow(inputMatrix),
                                     ncol(inputMatrix)-1,
                                     dimnames=list(rownames(inputMatrix),
                                                   colnames(inputMatrix)[-1]))
  
  print('Calculating interlayer edges...')
  pb <- txtProgressBar(min = 1, max = (ncol(inputMatrix)-1), style = 3)
  for (x in 1:(ncol(inputMatrix)-1)){
    setTxtProgressBar(pb, x)
    relativeChange[,x] <- inputMatrix[,x+1]/inputMatrix[,x]
    colnames(relativeChange)[x] <- paste(x,'-->',x+1,sep='')
  }
  relativeChange[is.na(relativeChange)] <- 0
  relativeChange[is.infinite(relativeChange)] <- 0
  relativeChange[is.nan(relativeChange)] <- 0
  return(relativeChange)
}


# Efficient function to transform a matrix to an edge list:
matrix.to.coordinateList <- function(m, default_w=1, undirected=F, verbose=F){
  if (undirected){
    m[upper.tri(m)] <- 0
  }
  NZ <- which(m!=0, arr.ind = T) # non-zero elements
  COO <- data.frame(i=rep(NA,nrow(NZ)),j=rep(NA,nrow(NZ)),w=rep(default_w,nrow(NZ)))
  COO$i <- rownames(m)[NZ[,1]]
  COO$j <- colnames(m)[NZ[,2]]
  if (verbose){print(paste('Found ',nrow(COO),' edges'))}
  for (n in 1:nrow(COO)){
    COO[n,'w'] <- m[as.character(COO[n,'i']),as.character(COO[n,'j'])]
  }
  return(COO)
}

matrix.to.coordinateList2 <- function(m, undirected=F, verbose=F, remove_negative_edges=T){
  if (undirected){
    cat('undirected ')
    m[upper.tri(m)] <- 0
    g <- graph.adjacency(m, mode = 'undirected', weighted = T)
  } else {
    g <- graph.adjacency(m, mode = 'directed', weighted = T)
  }
  COO <- get.data.frame(g)
  negative_edges <- sum(COO$weight<0)
  if (verbose){print(paste('Found ',nrow(COO),' edges (',negative_edges,' negative).',sep=''))}
  if (remove_negative_edges) {COO <- subset(COO, weight>=0)}
  names(COO) <- c('i','j','w')
  return(COO)
}




MI <- function (N) {
  S <- sum(N)
  CA = dim(N)[1]; CB = dim(N)[2]
  
  # Danon et al.'s 2005 method from http://arxiv.org/pdf/cond-mat/0505245.pdf
  Iup=0
  for (i in 1:CA){
    for (j in 1:CB){
      if (N[i,j] != 0){
        Ni.=sum(N[i,])
        N.j=sum(N[,j])
        Iup = Iup+N[i,j]*log((N[i,j]*S)/(Ni.*N.j))
      }
    }
  }
  Idown1=0;Idown2=0
  for (i in 1:CA){
    Ni.=sum(N[i,])
    Idown1=Idown1+Ni.*log(Ni./S)
  }
  for (j in 1:CB){
    N.j=sum(N[,j])
    Idown2=Idown2+N.j*log(N.j/S)
  }
  I=-2*Iup/(Idown1+Idown2)
  
  return(I)
}

mprint <- function(text){
  print(text)
  message(text)
}

# Taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#This function calculates similarity among strains. Row is each parasite, columns are gene type or allele types
# Taken from Qixin (June 2016)
overlapAlleleAdj<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(newmat)  
}
# mat <- matrix(c(1,0,1,1,0,0,1,1,1,0,1,1), ncol=6, nrow=2) # rows are strains, cols are alleles
# overlapAlleleAdj(mat)

padWithNA <- function(v, len){
  head(c(v, NA * 1:len), len)
}


plotLayer <- function(l,graph_list,moduleList,big_modules=NULL,unique_colors=F,
                      ver.col=NULL,verbose=T,cutoff_g=0,
                      edge_weight_multiply=5,...){ # requires the modules data frame!
  g=graph_list[[l]]
  layer_modules <- moduleList[moduleList$layer==l,] # the modules in that layer
  if(any(duplicated(layer_modules$strainCopy))){print('There are strains that occur in more than 1 module.')}
  if (verbose){
    cat('The following strains appear in the graph but not in the modules data frame. They will be removed:\n')
    cat(setdiff(V(g)$name,layer_modules[,'strainCopy']));cat('\n')
  }
  g <- delete_vertices(g, setdiff(V(g)$name,layer_modules[,'strainCopy']))
  
  g <- delete_edges(g, which(E(g)$weight<quantile(E(g)$weight, cutoff_g))) # remove all edges smaller than the cutoff
  #layout <-layout.kamada.kawai(g)
  
  V(g)$module <- layer_modules$module[match(V(g)$name,layer_modules$strainCopy)]
  if (verbose){
    cat(length(unique(V(g)$module)));cat(' modules are in this layer:\n')
    cat(sort(unique((V(g)$module))));cat('\n')
  }
  
  if (is.null(ver.col)){
    if (unique_colors) {
      cols = gg_color_hue(length(unique(V(g)$module)))
      names(cols) <- sort(unique((V(g)$module)))
      if (!is.null(big_modules)){cols[!names(cols)%in%big_modules] <- NA}
      V(g)$color <- cols[match(V(g)$module, names(cols))]
    } else {
      cols = gg_color_hue(max(moduleList$module))
      if (!is.null(big_modules)){cols[-big_modules] <- NA}
      V(g)$color <- cols[V(g)$module]
    }
  } else {
    V(g)$color <- ver.col
  }
    
  # plot.new()
  # par(mar=c(0,3,3,3))
  plot(g, 
       #layout=layout, 
       vertex.color=V(g)$color,
       # vertex.label=V(g)$module,
       vertex.label=NA,
       # edge.arrow.mode='-', 
       edge.arrow.width=1,
       edge.arrow.size=0.2,
       edge.curved=0.5, 
       edge.width=E(g)$weight*edge_weight_multiply,
       # asp=0, # This needed if changing margins
       ...)
}

plotSurveyLayer <- function(x, zeroDiag=T, cutoff_g=NULL,...){
  if(zeroDiag){diag(x) <- 0}
  g <- graph.adjacency(x, mode = 'directed', weighted = T, diag = F)
  if(!is.null(cutoff_g)){g <- delete_edges(g, which(E(g)$weight<quantile(E(g)$weight, cutoff_g)))} # remove all edges smaller than the cutoff
  plot(g, 
       # vertex.color='dark green',
       vertex.size=6,
       vertex.label=NA,
       # edge.arrow.mode='-', 
       edge.arrow.width=1,
       edge.arrow.size=0.2,
       edge.curved=0.5, 
       edge.color='black',
       edge.width=E(g)$weight*10,...)  
}

# These are the functions used for the malaria temporal networks
prep.packages = function(package.list) {
  loaded = package.list %in% .packages()
  if ( all(loaded) ) return(invisible())
  
  package.list = package.list[!loaded]
  installed = package.list %in% .packages(TRUE)
  if ( !all(installed) ) install.packages(package.list[!installed], repos="http://cran.rstudio.com/")
  for ( p in package.list )
  {
    print(paste("Loading package:",p))
    library(p,character.only=TRUE)
  }
}
# prep.packages(c("igraph","bipartite", #networks
#                 "ggplot2","RColorBrewer","grid","gtable","gridExtra","gplots", #graphics
#                 "stringr","googlesheets","reshape2","sqldf","RSQLite",'plyr',"foreach",'splitstackshape','data.tree', # data manipulation
#                 "betalink","TSA","Hmisc","inflection","entropy",'dynlm',"astsa","codyn")) #analysis


# If isolate A has a repertoire of nA types and isolate B has a repertoire of nB
# types, and the two isolates share a total of nAB types, and the union of the
# two sets is U_AB, we define PTS_AB = 2nAB/(nA+nB), and JS_AB  = nAB/U_AB.
calculatePTS <- function(mat, return_matrix=F) { # Row names are repertoires, col names are var genes
  newmat<-tcrossprod(mat>0) # var types shared between every two repertoires
  rs <- rowSums(mat>0)
  divMat<-matrix(rep(rs,length(rs)),length(rs),length(rs))
  divMat<-divMat+t(divMat) # This calculates the na+nb
  pts_matrix <- 2*newmat/divMat
  if (return_matrix){
    return(pts_matrix)
  } else {
    return(pts_matrix[lower.tri(pts_matrix)])
  }
}


removeNA <- function(x){
  x <- x[!is.na(x)]
}


# Remove the transient strains, defined as those which appear in <= n layers
removeTransientStrains <- function(layerlist,min.layers=3){
  strains <- unlist(lapply(layerlist, rownames))
  strains1 <- length(unique(strains))
  transientstrains <- names(which(table(strains)<=min.layers))
  layerlist <- lapply(layerlist, function(x) x[!rownames(x)%in%transientstrains,!colnames(x)%in%transientstrains])
  strains2 <- length(unique(unlist(lapply(layerlist, rownames))))
  print(paste('All:',strains1))
  print(paste('Transient removed:',length(transientstrains)))
  print(paste('left:',strains2,'that occur in at least',min.layers,'layers'))
  return(layerlist)
}


sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}


showAsDF <- function(x){
  d <- data.frame(id=1:length(x),x=x)
  print(d)
}

splitText <- function(str,after=T,splitchar='\\.'){
  if (after){
    return(sapply(strsplit(str, split=splitchar), tail, 1))
  }
  if (after==F){
    return(sapply(strsplit(str, split=splitchar), head, 1))
  }
}

stripDigits <- function(x){
  as.numeric(gsub("([0-9]+).*$", "\\1", x))
}

write_infomap <- function(x,file,nodeList,undirected = T){
  cat(paste('[',Sys.time(), '] writing Infomap file...',sep='')) 
  infomapEdgeList <- matrix.to.coordinateList2(x, undirected = undirected, verbose = T, remove_negative_edges=T)
  infomapEdgeList <- subset(infomapEdgeList, w>=0)
  infomapEdgeList$node_s <- nodes$nodeID[match(infomapEdgeList$i,nodes$nodeLabel)]
  infomapEdgeList$node_t <- nodes$nodeID[match(infomapEdgeList$j,nodes$nodeLabel)]
  if(file.exists(file)){unlink(file)}
  sink(file, append = T)
  cat(paste("*Vertices",nrow(nodes)));cat('\n')
  write.table(nodes[,c('nodeID','nodeLabel')], file, append = T,sep=' ', quote = T, row.names = F, col.names = F)
  cat(paste("*Arcs",nrow(infomapEdgeList)));cat('\n')
  write.table(infomapEdgeList[,c('node_s','node_t','w')], file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
  sink.reset()
}



