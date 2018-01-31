# Initialize --------------------------------------------------------------
library(ggplot2)
library(data.table)

getRuns <- function(Ups, resolution, iterations){
  lower_curve <- upper_curve <- c()
  for (i in 1:iterations){
    print(i)
    lower_curve <- rbind(lower_curve, fread(paste('curveData_',Ups,'_0.95_subsampleT_',resolution,'_run_',i,'.csv',sep=''), he=T))
    upper_curve <- rbind(upper_curve, fread(paste('curveData_',Ups,'_0.95_subsampleF_',resolution,'_run_',i,'.csv',sep=''), he=T))
  }
  lower_curve$bound <- 'L'
  upper_curve$bound <- 'U'
  curveData <- rbind(lower_curve,upper_curve)
  return(curveData)
}


# Get data ----------------------------------------------------------------
setwd('/home/shai/Documents/Shazia/Results')
d <- NULL
for (r in c(70,80,90,96)){
  for (u in (c('A','BC','ABC'))){
    x <- getRuns(u,r,25)
    d <- rbind(d, x)
  }
}


# Plot --------------------------------------------------------------------
d$Ups <- factor(d$Ups, levels = c('ABC','BC','A'))
EIR=100
ggplot(d, aes(bites/EIR, propVars, color=subsample))+
  geom_smooth(se = T, method = "gam", formula = y ~ s(log(x+1)))+
  theme_bw()+
  scale_color_manual(values=c('springgreen4','brown'))+
  labs(x=paste('Years since birth (assuming EIR of ',EIR,')',sep=''), y='Proportion of DBLa types accumulated')+
  # facet_wrap(~Ups,scales='free')
  facet_grid(Ups~resolution, scales='free_x')

