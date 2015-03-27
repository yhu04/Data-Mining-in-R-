##################################
###### Predicting Algae Bloom#####
##################################

################
# Load Dataset #
################
algae=read.table('Analysis,txt',
                       header=F,
                       dec='.'
                       col.names=c('season','size','speed','mxPH','mn02','C1','N03','NH4',
                       'oPO4','Chla','a1','a2','a3','a4','a5','a6','a7'),
                       na.strings=c('XXXXX'))
                       
################################
# Data Visualization & Summary #
################################
summary(algae)

##plot mxPH
library(car)
par(mfrow=c(1,2))
hist(algae$mxPH,prob=T,xlab='',main='Histogram of maximum pH value',ylim=0:1)
lines(density(algae$mxPH,na.rm=T))
rug(jitter(algae$mxPH))
qq.plot(algae$mxPH,main='Normal QQ plot of maximum pH value')
par(mfrwo=c(1,1))

##plot oPO4
boxplot(algae$oPO4,ylab='Orthophosphate(oPO4)')
rug(jitter(algae$oPO4),side=2)
abline(h=mean(algae$oPO4,na.rm=T),lty=2)

##plot NH4
plot(algae$NH4,xlab='')
abline(h=mean(algae$NH4,na.rm=T),lty=1)
abline(h=mean(algae$NH4,na.rm=T)+sd(algae$NH4,na.rm=T),lty=2)
abline(h=median(algae$NH4,na.rm=T),lty=3)
identify(algae$NH4)

## plot size
library(lattice)
bwplot(size~a1,data=algae,ylab='River Size',xlab='Algal A1')
library(Hmisc)
bwplot(size~a1,data=algae,panel=panel.bpplot,probs=seq(0.01,0.49,by=0.01),datadensity=TRUE,ylab='River Size',xlab='Algal A1')

## plot mnO2
minO2=equal.count(na.omit(algae$mnO2),number=4,overlap=1/5)
stripplot(season~a3|minO2,data=algae[!is.na(algae$mnO2),])

##################
# Missing Values #
##################
#remove missing values
algae[!complete.cases(algae),]
nrow(algae[!complete.case(algae),])
algae=na.omit(algae)
apply(algae,1,function(x) sum(is.na(x)))
data(algae)
manyNAs(algae,0.2)
algae=algae[-manyNAs(algae),]

#fill in most requent values
algae[48,'mxPH']=mean(algae$mxPH,na.rm=T)
algae[is.na(algae$Chla),'Chla']=median(algae$Chla,na.rm=T)
data(algae)
algae=algae[-manyNAs(algae),]
algae=centralImputation(algae)

#fill in missing values by correlations
symnum(cor(algae[,4:18],use='complete.obs'))
data(algae)
algae=algae[-manyNAs(algae),]
fit=lm(PO4~oPO4,data=algae)
fillPO4=function(oP){
	if (is.na(oP))
	   return(NA)
	else return(42.897+1.293*oP)
}
algae[is.na(algae$PO4),'PO4']=sapply(algae[is.na(algae$PO4),'oPO4'],fillPO4)
algae$season=factor(algae$season,levels=c('spring','summer','autumn','winter'))
histogram(~mxPH|season,data=algae)
histogram(~mxPH|size,data=algae)
stripplot(size~mxPH)|speed,data=algae,jitter=T)

#fill in missing values with similarity 
data(algae)
algae=algae[-manyNAs(algae),]
clean.algae=knnImputation(algae,k=10)

#####################
# Prediction Models #
#####################

## Multivariable linear regression---without missing values 
lm.a1=lm(a1~.,data=clean.algae[,1:12])
summary(lm.a1)
anova(lm.a1)
lm2.a1=update(lm.a1,.~.-season)
anova(lm.a1,lm2.a1)
final.lm=step(lm.a1)
summary(final.lm)

## Regression Trees 
library(rpart)
data(algae)
algae=algae[-manyNAs(algae),]
rt.a1=rpart(a1~.,data=algae[,1:12])
prettyTree(rt.a1)
printcp(rt.a1)
rt2.a1=prune(rt.a1,cp=0.08)
first.tree=rpart(a1~.,data=algae[,1:12])
snip.rpart(first.tree,c(4,7))

##################################
# Model evaluation and selection #
##################################
lm.prediction.a1=predict(final.lm,clean.algae)
rt.prediction.a1=predict(rt.a1,algae)

## Mean absolute error 
mae.a1.lm=mean(abs(lm.prediction.a1-algae[,'a1']))
mae.a1.rt=mean(abs(rt.prediction.a2-algae[,'a1']))

## Mean squared error 
mse.a1.lm=mean((lm.prediction.a1-algae[,'a1'])^2)
mse.a1.rt=mean((rt.prediction.a2-algae[,'a1'])^2)

## Normalized 
nmse.a1.lm=mean((lm.prediction.a1-algae[,'a1'])^2)/mean((mean(algae[,'a1'])-algae[,'a1'])^2)
nmse.a1.rt=mean((rt.prediction.a2-algae[,'a1'])^2)/mean((mean(algae[,'a1'])-algae[,'a1'])^2)

## visualization 
old.par=par(mfrow=c(1,2))
plot(lm.prediction.a1,algae[,'a1'],main='Linear Model',xlab='Predictions',ylab='True Values')
abline=(0,1,lty=2)
plot(rt.prediction.a1,algae[,'a1'],main='Regression Tree',xlab='Predictions',ylab='True Values')
abline(0,1,lty=2)
par(old.par)

sensible.lm.predictions.a1=ifelse(lm.prediction.a1<0,0,lm.prediction.a1)

## cross-validation
cv.rpart=function(form,train,test,...){
	m=rpartXse(form,train,...)
	p=predict(m,test)
	c=(nmse=mse/mean((mean(resp(form,train))-resp(from,test)))^2)
}
cv.lm=function(form,train,test,...){
	m=lm(form,train,...)
	p=predict(m,test)
	p=ifelse(p<0,0,p)
	mse=mean((p-resp(form,test))^2)
	c=(nmse=mse/mean((mean(resp(form,train))-resp(from,test)))^2)	
}

res=experimentalComparison(
         c(dataset(a1~.,clean.algae[,1:12],'a1')),
         c(variants('cv.lm'),variants('cv.rpart',se=c(0,0.5,1))),
         cvSettings(3,10,1234)
        )
        
summary(res)       
plot(res)

DSs=sapply(names(clean.algae[12:18],function(x,names.attrs){
	                 f=as.formula(paste(x,'~.'))
	                 dataset(f,clean.algae[,c(names.attrs,x)],x)},
	                 names(clean.algae))) 
	                 
##############################	                 
# Prediction on test dataset #	                 
##############################	            
bestModelsNames=sapply(bestScores(res.all),function(x) x['nmse','system'])
learners=c(rf='randomForest',rpart='rpartXse')
funcs=leaners[sapply(strsplit(bestModesNames,'\\.'),function(x) x[2])]
parSetts=lapply(bestModelsNames,function(x) getVariant(x,res.all)@pars)
bestModels=list()
for (a in 1:7){
	form=as.formula(paste(names(clean.algae)[11+a],'~.'))
	bestModels[[a]]=do.call(funcs[a],c(list(form,clean.algae[,c(1:11,11+a)]),parSetts[[a]]))
}
clean.test=algae=knnImputation(test.algae,k=10,distData=algae[,1:11])
preds=matrix(ncol=7,nrow=140)
for( i in 1:nrow(clean.test.algae)){
	preds[i,]=sapply(1:7,
	                  function(x) predict(bestModels[[x]],clean.test.algae[i,]))
} 

avg.preds=apply(algae[,12:18],2,mean)
apply((algae.sol-preds)^2,2,mean)/apply((scale(algae.sols,avg.preds,F)^2),2,mean)
