#setwd('~/Dropbox/RESEARCH/PROJECTS/Shazia')
#source('mtn_functions.R')
#library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
moi_arg <- as.numeric(args[1])
ups_arg <- as.character(args[2])
runs_arg <- as.numeric(args[3])
threshold_arg <- as.numeric(args[4])
subsample_vars_arg <- as.logical(args[5])


# Functions ---------------------------------------------------------------
makeCurve <- function(dat,threshold=0.95,subsample_vars){
  R <- colnames(dat) # The set of repertoires
  V <- rownames(dat) # The set of vars
  cumulativeSampledVars <- c()
  cumulativeNumberSampledVars <- 0
  bites <- 0
  totalVars <- 0
  threshold <- ceiling(threshold*length(V))
  while(totalVars<threshold){ # Loop until number of vars reaches the threshold
    print(paste('Bite: ',bites,' | total vars: ',totalVars,' | ',round(totalVars/length(V)*100,2),'% of vars sampled',sep=''))
    bites <- bites+1
    sampledRepertoire <- sample(R, 1) # Sample a repertoire
    sampledVars <- names(which(dat[,sampledRepertoire]==1)) # What are the var genes of this rpertoire?
    # If the sampled repertoire has a length>60 then subsample it to produce a repertoire
    # cat('# vars: ');cat(length(sampledVars));cat(' | ')
    if (subsample_vars){
      if (length(sampledVars)>60){
        # cat('subsampling... | ')
        vars2Subsample <- as.numeric(sample(names(repertoireLengthDistribution),1,prob=repertoireLengthDistribution)) # Determine the number of vars to subsample based on the length distribution of MOI=1
        sampledVars <- sample(sampledVars, vars2Subsample, F) 
      }
    }
    # Limit to groups of Ups
    # if (Ups_subset!='ABC'){
    #   sampledVars <- sampledVars[which(sampledVars%in%subset(dataUps, Ups==Ups_subset)$Type)]
    # }
    
    numberNewVars <- sum(!sampledVars%in%cumulativeSampledVars) # How many of the sampled vars were not seen before?
    cumulativeNumberSampledVars <- c(cumulativeNumberSampledVars,numberNewVars+totalVars) # Add this number to the vector
    totalVars <- tail(cumulativeNumberSampledVars,1)
    
    cumulativeSampledVars <- c(cumulativeSampledVars, sampledVars) # Add sampled vars to the cumulative list of sampled vars
    cumulativeSampledVars <- unique(cumulativeSampledVars) # Remove duplicates in the cumulative list of sampled vars
  }
  return(data.frame(bites=0:bites,propVars=cumulativeNumberSampledVars/length(V)))
}

subsetDataMOI <- function(MOI){
  if (MOI>=100){
    data_MOI <- data
  } else {
    maxVars <- MOI*60
    data_MOI <- data[,colSums(data)<=maxVars]
    data_MOI <- bipartite::empty(data_MOI)
    dim(data_MOI)
    repertoires <- colnames(data_MOI)
    vars <- rownames(data_MOI)
    print(paste('MOI<=',MOI,' | ',length(repertoires),' repertoires | ',length(vars),' vars',sep=''))
  }
  return(data_MOI)  
}

# Initialize --------------------------------------------------------------
message('Loading data...')
data <- read.csv('S1 S2 Single Infections 20-60 DBLa types ASYMP RESERVOIR.csv', row.names = 1)
data <- data.matrix(data)
repertoireLengthDistribution <- table(colSums(data))

data <- read.csv('S1 S2 All isolates 20 or more DBLa types ASYMP RESERVOIR.csv', row.names = 1)
vars <- rownames(data)
data <- data.matrix(data)
range(data)

# Get A/BC groups
dataUps <- read.csv('S1 S2 All isolates 20 or more DBLa types ASYMP RESERVOIR A NON A PRELIM.csv')
dataUps <- dataUps[,1:2]
table(dataUps$Ups)
setequal(dataUps$Type,vars)
dataUps <- dataUps[match(vars,dataUps$Type),]
all(dataUps$Type==vars)

# Run curves --------------------------------------------------------------

runCurves <- function(MOI, Ups='ABC', runs=20, threshold=0.95, subsample_vars, write2file=T){
  dataMOI <- subsetDataMOI(MOI)
  if (Ups=='A'){
    dataMOI <- dataMOI[rownames(dataMOI)%in%dataUps$Type[dataUps$Ups=='A'],]
  }
  if (Ups=='BC'){
    dataMOI <- dataMOI[rownames(dataMOI)%in%dataUps$Type[dataUps$Ups=='BC'],]
  }

  results <- list()
  for (run in 1:runs){
    print(paste('[',Sys.time(),'] run: ',run,sep=''))
    tmp <- makeCurve(dataMOI,threshold, subsample_vars = subsample_vars)
    tmp$run <- run
    results[[run]] <- tmp
  }
  df <- do.call('rbind', results)
  df$MOI <- MOI
  df$Ups <- Ups
  if(write2file){
    write.table(df, paste('curveData_',MOI,'_',Ups,'_',runs,'_',threshold,'_',subsample_vars,'.csv',sep=''),row.names = F, sep=',')
  }
  return(df)
}

#When running on midway
message('Running curves...')
curveData <- runCurves(MOI = moi_arg, Ups = ups_arg, runs=runs_arg, threshold = threshold_arg, subsample_vars = subsample_vars_arg)

# When running localy:
resultsMOI1A <- runCurves(MOI = 1, Ups = 'A', runs=20, threshold = 0.95)
resultsMOI2A <- runCurves(MOI = 2, Ups = 'A', runs=20, threshold = 0.95)
resultsMOI3A <- runCurves(MOI = 3, Ups = 'A', runs=20, threshold = 0.95)
resultsMOI200A <- runCurves(MOI = 200, Ups = 'A', runs=20, threshold = 0.95)

resultsMOI1BC <- runCurves(MOI = 1, Ups = 'BC', runs=2, threshold = 0.95)
resultsMOI2BC <- runCurves(MOI = 2, Ups = 'BC', runs=2, threshold = 0.95)
resultsMOI3BC <- runCurves(MOI = 3, Ups = 'BC', runs=2, threshold = 0.95)
resultsMOI200BC <- runCurves(MOI = 200, Ups = 'BC', runs=2, threshold = 0.95)

resultsMOI1ABC <- runCurves(MOI = 1, Ups = 'ABC', runs=2, threshold = 0.95)
resultsMOI2ABC <- runCurves(MOI = 2, Ups = 'ABC', runs=2, threshold = 0.95)
resultsMOI3ABC <- runCurves(MOI = 3, Ups = 'ABC', runs=2, threshold = 0.95)
resultsMOI200ABC <- runCurves(MOI = 200, Ups = 'ABC', runs=2, threshold = 0.95)



resultsABC <- runCurves(MOI = 200, Ups = 'ABC', runs=1, threshold = 0.95, subsample_vars = F)
resultsABC_subsample <- runCurves(MOI = 200, Ups = 'ABC', runs=1, threshold = 0.95, subsample_vars = T)
resultsA <- runCurves(MOI = 200, Ups = 'A', runs=1, threshold = 0.95, subsample_vars = F)
resultsA_subsample <- runCurves(MOI = 200, Ups = 'A', runs=1, threshold = 0.95, subsample_vars = T)
resultsBC <- runCurves(MOI = 200, Ups = 'BC', runs=1, threshold = 0.95, subsample_vars = F)
resultsBC_subsample <- runCurves(MOI = 200, Ups = 'BC', runs=1, threshold = 0.95, subsample_vars = T)

resultsMOI200BC$bound='U'
resultsMOI200BC_subsample$bound='L'
d <- rbind(resultsMOI200BC,resultsMOI200BC_subsample)
ggplot(d, aes(bites/25,propVars,group=bound,color=bound))+
  geom_point(data=subset(d, run>=1 & run<10),alpha=0.05)+
  geom_smooth(color='black')+
  theme_bw()+
  labs(x='Years', y='Proportion of vars accumulated')
