args <- commandArgs(trailingOnly=TRUE)
ups_arg <- as.character(args[1])
run_arg <- as.numeric(args[2])
threshold_arg <- as.numeric(args[3])

# Initialize --------------------------------------------------------------
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


# Get A/BC groups
dataUps <- read.csv('S1 S2 All isolates 20 or more DBLa types ASYMP RESERVOIR A NON A PRELIM.csv')
dataUps <- dataUps[,1:2]
table(dataUps$Ups)
#setequal(dataUps$Type,vars)
# dataUps <- dataUps[match(vars,dataUps$Type),]
#all(dataUps$Type==vars)

if (ups_arg!='ABC'){
  print(ups_arg)
  data <- read.csv('S1 S2 Single Infections 20-60 DBLa types ASYMP RESERVOIR.csv', row.names = 1)
  data <- data.matrix(data)
  data <- data[rownames(data)%in%subset(dataUps, Ups==ups_arg)$Type,]
  data <- empty(data)
  repertoireLengthDistribution <- table(colSums(data))
  data <- read.csv('S1 S2 All isolates 20 or more DBLa types ASYMP RESERVOIR.csv', row.names = 1)
  data <- data.matrix(data)
  data <- data[rownames(data)%in%subset(dataUps, Ups==ups_arg)$Type,]
  data <- empty(data)
} else {
  print('ABC')
  data <- read.csv('S1 S2 Single Infections 20-60 DBLa types ASYMP RESERVOIR.csv', row.names = 1)
  data <- data.matrix(data)
  repertoireLengthDistribution <- table(colSums(data))
  data <- read.csv('S1 S2 All isolates 20 or more DBLa types ASYMP RESERVOIR.csv', row.names = 1)
  data <- data.matrix(data)
}

range(data)
vars <- rownames(data)



# Run curves --------------------------------------------------------------
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

message('Loading data...')

tmp <- makeCurve(data, threshold_arg, subsample_vars=T)
tmp$run <- run_arg
tmp$Ups <- ups_arg
write.table(tmp, paste('curveData_',ups_arg,'_',threshold_arg,'_subsampleT_',run_arg,'.csv',sep=''),row.names = F, sep=',')

tmp <- makeCurve(data, threshold_arg, subsample_vars=F)
tmp$run <- run_arg
tmp$Ups <- ups_arg
write.table(tmp, paste('curveData_',ups_arg,'_',threshold_arg,'_subsampleF_',run_arg,'.csv',sep=''),row.names = F, sep=',')
  
