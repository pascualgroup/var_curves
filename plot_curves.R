library(ggplot2)
getRuns <- function(Ups){
  lower_curve <- upper_curve <- c()
  for (i in 1:250){
    print(i)
    lower_curve <- rbind(lower_curve, read.csv(paste('curveData_',Ups,'_0.95_subsampleT_',i,'.csv',sep=''), he=T))
    upper_curve <- rbind(upper_curve, read.csv(paste('curveData_',Ups,'_0.95_subsampleF_',i,'.csv',sep=''), he=T))
  }
  lower_curve$bound <- 'L'
  upper_curve$bound <- 'U'
  curveData <- rbind(lower_curve,upper_curve)
  return(curveData)
}
curveData_A <- getRuns('A')
curveData_BC <- getRuns('BC')
curveData_ABC <- getRuns('ABC')

curveData_A$Ups <- 'A'
curveData_BC$Ups <- 'BC'
curveData_ABC$Ups <- 'ABC'

d <- rbind(curveData_A,curveData_BC,curveData_ABC)
write.csv(d,'curve_data_paper.csv')

# d <- read.csv('curve_data_paper.csv')
pdf('curves_paper.pdf', width = 14, height = 8)
ggplot(d, aes(bites/25,propVars,group=bound,color=bound))+
  # geom_point(data=subset(d, run>=1 & run<5),alpha=0.05)+
  geom_smooth(se = T, method = "gam", formula = y ~ s(log(x+1)))+
  theme_bw()+
  labs(x='Years', y='Proportion of vars accumulated')+
  facet_wrap(~Ups,scales='free')
dev.off()
# 
# resultsMOI200BC$bound='U'
# resultsMOI200BC_subsample$bound='L'
# d <- rbind(resultsMOI200BC,resultsMOI200BC_subsample)
# ggplot(d, aes(bites/25,propVars,group=bound,color=bound))+
#   geom_point(data=subset(d, run>=1 & run<10),alpha=0.05)+
#   geom_smooth(color='black')+
#   theme_bw()+
#   labs(x='Years', y='Proportion of vars accumulated')


# 
# 
# 
# d.A <- c()
# d.A <- rbind(d.A, read.csv('curveData_1_A_100_0.95.csv', he=T))
# d.A <- rbind(d.A, read.csv('curveData_2_A_100_0.95.csv', he=T))
# d.A <- rbind(d.A, read.csv('curveData_3_A_100_0.95.csv', he=T))
# d.A <- rbind(d.A, read.csv('curveData_200_A_100_0.95.csv', he=T))
# d.A$MOI <- factor(d.A$MOI)
# d.A$Ups <- factor(d.A$Ups)
# 
# d.BC <- c()
# d.BC <- rbind(d.BC, read.csv('curveData_1_BC_100_0.95.csv', he=T))
# d.BC <- rbind(d.BC, read.csv('curveData_2_BC_100_0.95.csv', he=T))
# d.BC <- rbind(d.BC, read.csv('curveData_3_BC_100_0.95.csv', he=T))
# d.BC <- rbind(d.BC, read.csv('curveData_200_BC_100_0.95.csv', he=T))
# d.BC$MOI <- factor(d.BC$MOI)
# d.BC$Ups <- factor(d.BC$Ups)
# 
# d.ABC <- c()
# d.ABC <- rbind(d.ABC, read.csv('curveData_1_ABC_100_0.95.csv', he=T))
# d.ABC <- rbind(d.ABC, read.csv('curveData_2_ABC_100_0.95.csv', he=T))
# d.ABC <- rbind(d.ABC, read.csv('curveData_3_ABC_100_0.95.csv', he=T))
# d.ABC <- rbind(d.ABC, read.csv('curveData_200_ABC_100_0.95.csv', he=T))
# d.ABC$MOI <- factor(d.ABC$MOI)
# d.ABC$Ups <- factor(d.ABC$Ups)
# 
# pdf('curves_A.pdf')
# ggplot(d.A, aes(bites/25,propVars,group=MOI,color=MOI))+
#   geom_point(data=subset(d.A, run>=1 & run<20),alpha=0.05)+
#   geom_smooth(color='black')+
#   theme_bw()+
#   labs(x='Years', y='Proportion of vars accumulated')
# dev.off()
# 
# pdf('curves_BC.pdf')
# ggplot(d.BC, aes(bites/25,propVars,group=MOI,color=MOI))+
#   geom_point(data=subset(d.BC, run>=1 & run<20),alpha=0.05)+
#   geom_smooth(color='black')+
#   theme_bw()+
#   labs(x='Years', y='Proportion of vars accumulated')
# dev.off()
# 
# pdf('curves_ABC.pdf')
# ggplot(d.ABC, aes(bites/25,propVars,group=MOI,color=MOI))+
#   geom_point(data=subset(d.ABC, run>=1 & run<10),alpha=0.05)+
#   geom_smooth(color='black')+
#   theme_bw()+
#   labs(x='Years', y='Proportion of vars accumulated')
# dev.off()
# 
# d <- rbind(d.A,d.BC,d.ABC)
# pdf('curves_ALL.pdf')
# ggplot(d, aes(bites/25,propVars,group=MOI,color=MOI))+
#   geom_point(data=subset(d, run>=1 & run<10), alpha=0.05)+
#   geom_smooth(color='black')+
#   facet_wrap(~Ups)+theme_bw()+
#   labs(x='Years', y='Proportion of vars accumulated')
# dev.off()
# 