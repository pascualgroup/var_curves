# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
ups_arg <- as.character(args[1])
threshold_arg <- as.numeric(args[2])
MOI <- as.numeric(args[3])
data_file <- as.character(args[4])
 
# data_file <- 'For_Shai_1099_isolates_764_types_at_70%_seq_ID.csv'
run_arg <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Functions --------------------------------------------------------------
prep.packages = function(package.list) {
  loaded = package.list %in% .packages()
  if ( all(loaded) ) return(invisible())
  
  package.list = package.list[!loaded]
  installed = package.list %in% .packages(TRUE)
  if ( !all(installed) ) install.packages(package.list[!installed], repos="http://cran.rstudio.com/")
  for ( p in package.list )
  {
    print(paste("Loading package:",p))
    suppressMessages(library(p,character.only=TRUE))
  }
}

# This function is taken from the bipartite package
empty <- function (web, count = FALSE) {
  web[is.na(web)] <- 0
  if (NCOL(web) == 1 | NROW(web) == 1) {
    if (NCOL(web) == 1 & NROW(web) != 1) {
      nr <- sum(web > 0)
      nc <- 1
    }
    if (NROW(web) == 1 & NCOL(web) != 1) {
      nc <- sum(web > 0)
      nr <- 1
    }
    if (NROW(web) == 1 & NCOL(web) == 1) {
      nr <- 1
      nc <- 1
    }
    out <- web[1:nr, 1:nc, drop = FALSE]
    if (count) 
      attr(out, "empty") <- c(`empty rows` = NROW(web) - 
                                nr, `empty columns` = NCOL(web) - nc)
    return(out)
  }
  cempty <- which(colSums(web) == 0)
  rempty <- which(rowSums(web) == 0)
  cind <- if (length(cempty) == 0) 
    1:NCOL(web)
  else (1:NCOL(web))[-cempty]
  rind <- if (length(rempty) == 0) 
    1:NROW(web)
  else (1:NROW(web))[-rempty]
  out <- web[rind, cind, drop = FALSE]
  if (count) 
    attr(out, "empty") <- c(`empty rows` = length(rempty), 
                            `empty columns` = length(cempty))
  return(out)
}

subsetDataMOI <- function(data_input, MOI, min_repertoire_length=20){
  if (MOI>=100){
    data_MOI <- data_input
  } else {
    maxVars <- MOI*60
    data_MOI <- data_input[,colSums(data_input)>=min_repertoire_length & colSums(data_input)<=maxVars]
    data_MOI <- empty(data_MOI)
    dim(data_MOI)
    print(paste('MOI<=',MOI,' | ',ncol(data_MOI),' repertoires | ',nrow(data_MOI),' vars',sep=''))
  }
  return(data_MOI)  
}


# Libraries ---------------------------------------------------------------
library(data.table)
library(stringr)

# Load data ---------------------------------------------------------------
data <- fread(paste('Data/',data_file,sep=''))
# Get A/BC groups
dataUps <- data[,1:2]
print(table(dataUps$Ups))
# Get DBLa types
DBLa <- data[,1]
# Make a matrix
data <- data[,-c(1,2)]
data <- data.matrix(data)
rownames(data) <- DBLa$DBLa_type


# The distribution of repertoire length is based on MOI=1 and is per ups A/BC group
if (ups_arg!='ABC'){
  print(ups_arg)
  # Get the distribution of repertoire lengths
  dataMOI1 <- subsetDataMOI(data, 1)
  dataMOI1 <- dataMOI1[rownames(dataMOI1)%in%subset(dataUps, Ups==ups_arg)$DBLa_type,]
  dataMOI1 <- empty(dataMOI1)
  repertoireLengthDistribution <- table(colSums(dataMOI1))
  # Subset the data to the desired MOI
  data <- subsetDataMOI(data, MOI)
  # Now subsample the data to obtain only the desired ups group (A or BC)
  data <- data[rownames(data)%in%subset(dataUps, Ups==ups_arg)$DBLa_type,]
  data <- empty(data)
} else {
  print('ABC')
  # Get the distribution of repertoire lengths
  dataMOI1 <- subsetDataMOI(data, 1)
  repertoireLengthDistribution <- table(colSums(dataMOI1))
  # Subset the data to the desired MOI
  data <- subsetDataMOI(data, MOI)
}
print('------------------------------------------------------------')
print(paste('Ups:',ups_arg,'| MOI =',MOI,'| ',ncol(data),'repertoires |',nrow(data),'types.'))
print('------------------------------------------------------------')


vars <- rownames(data)

# Run curves --------------------------------------------------------------
# Subsampling vars is to create a lower boundry to the curve
makeCurve <- function(dat, threshold=0.95, subsample_vars){
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
    
    # When sub-sampling vars, then if the sampled repertoire has a length>60
    # then subsample it to produce a repertoire of length 60 or less. This immitates the case of
    # MOI=1
    if (subsample_vars){
      if (length(sampledVars)>60){
        # cat('subsampling... | ')
        vars2Subsample <- as.numeric(sample(names(repertoireLengthDistribution),1,prob=repertoireLengthDistribution)) # Determine the number of vars to subsample based on the length distribution of MOI=1
        sampledVars <- sample(sampledVars, vars2Subsample, F) 
      }
    }
    # If not subsampling, or in other words, if we look at isolates with MOI>1,
    # then make sure that the sample we take is not from an isolate with MOI=1.
    # This is to ensure that the upper boundary is for MOI>1.
    if (!subsample_vars){
      if (length(sampledVars)>60){
        next # If the sample was taken from an individual with MOI=1 then discard the sample.
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

pos <- str_locate(data_file,'%')[1]
resolution <- str_sub(data_file,pos-2,pos-1)

print('Making curves with subsample...')
tmp <- makeCurve(data, threshold_arg, subsample_vars=T)
tmp$run <- run_arg
tmp$Ups <- ups_arg
tmp$subsample <- T
tmp$resolution <- resolution
write.table(tmp, paste('Results/curveData_',ups_arg,'_',threshold_arg,'_subsampleT_',resolution,'_run_',run_arg,'.csv',sep=''),row.names = F, sep=',')

print('Making curves without subsample...')
tmp <- makeCurve(data, threshold_arg, subsample_vars=F)
tmp$run <- run_arg
tmp$Ups <- ups_arg
tmp$subsample <- F
tmp$resolution <- resolution
write.table(tmp, paste('Results/curveData_',ups_arg,'_',threshold_arg,'_subsampleF_',resolution,'_run_',run_arg,'.csv',sep=''),row.names = F, sep=',')
  
