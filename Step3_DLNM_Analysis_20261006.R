####DLNM Analysis####
#Madeline and Carina
#Updated 2023-03-20
#Note: This version includes precipitation waves.

####LIBRARIES AND FILE PATHS####
source(paste0('T:\\Analysis\\Madeline\\',"MS_Precip_Libraries.R"))
library(MASS)

####LOAD FILES####
#temp.RData is created in the step 2 program
load(paste0(f1,'temp.RData'))
#load additional files if previously created
#load(paste0(f1,'models.RData'))
#load('E://Madeline//cb1_ppt.RData')
#load('E://Madeline//cb1_tmax.RData')

####FUNCTIONS#####

#Function to create the precip and temperature crossbases piecemeal
#This loops through each set of 100,000 rows and can be used when the crossbasis function doesn't run on the entire exposure data set. This requires us to first specify internal and boundary knots, because the crossbasis won't otherwise have the full data set to calculate quantiles with. This also converts to 100ths of the value in integer format to save space, requiring the exposures for the final model interpretations to be multiplied by 100 to make them equivalent to these crossbasis data that the model was run with.
cbFxn<-function(expData=claims_master_agg,expNames=precipCols[1:21],cbArgVar=lsvar1,cbArgLag=lslag1) {
  mat1<-crossbasis(expData[1:100000,..expNames],argvar=cbArgVar,arglag = cbArgLag)
  pb<-txtProgressBar(min=100001,max=nrow(expData),style=3)
  for (i in seq(100001,nrow(expData),by=100000)) {
    setTxtProgressBar(pb,i)
    j<-min(nrow(expData),i+99999)
    mat2<-crossbasis(expData[i:j,..expNames],argvar=cbArgVar,arglag = cbArgLag)
    mat1<-rbind(mat1,mat2)
  }
  v4<-c('dimnames','df','range','lag','argvar','arglag','class')
  attributes(mat1)[v4]<-attributes(mat2)[v4]
  dt1<-as.integer(round(mat1*100,0)) #convert to 100ths of the value in integer format to save space
  attributes(dt1)<-attributes(mat1)
  dt1
}

##crosspred function for season-specific plots, taking into account the season-specific limits on the crossbases
cpSeasFxn<-function(cbName,seasCond=c1,modName,center=0,val=NULL,from=NULL,to=NULL) {
  v4<-c('dimnames','df','range','lag','argvar','arglag','class')
  cb1<-get(cbName)[seasCond,]
  attributes(cb1)[v4]<-attributes(get(cbName))[v4]
  v1<-grep(paste0(cbName,'\\['),names(models[[modName]]$coef),value=T)
  crosspred(basis=cb1,coef=models[[modName]]$coef[v1]*100,vcov=models[[modName]]$vcov[v1,v1]*100^2,model.link='log',cen=center,cumul=T,at=val,from=from,to=to)
}

##Function to create crossbasis for a given wave definition
#the data set, expData, must be ordered by super.zip and then CLM_FROM_DT within each super.zip, and there must be a CLM_FROM_DT for every possible date for every super.zip.
waveFxn<-function(expData=claims_master_agg[,.(super.zip,CLM_FROM_DT,PRISM_PPT)],waveCondition='(exp(PRISM_PPT/100)-1)/10/2.54==0',maxWave=9,cbArgVar=list('ns',knots=c(3,6)),cbArgLag=list('ns',knots=c(3,7))) {
  dat1<-expData[,waveName:=eval(str2lang(waveCondition))] #identify days meeting threshold, or condition
  dat1[,waveID:=c(NA,cumsum(abs(diff(waveName,lag=1))))] #create IDs for each set of waves
  dat1[,wave:=cumsum(waveName),by=waveID]
  #if day exceeds maxWave, set to maxWave value (interpretation is maxWave or higher)
  dat1[wave>maxWave,wave:=maxWave]
  #drop the Jan 2006 values up to 21 days (these may be from the prior super.zip)
  dat1[as.POSIXlt(CLM_FROM_DT)$year==106 & as.POSIXlt(CLM_FROM_DT)$yday<21,wave:=NA]
  cat('Percent of rows in each ')
  print(round(dat1[,table(wave,useNA='a')]/nrow(dat1)*100,7))
  mat1<-crossbasis(dat1[,wave],lag=20,argvar=cbArgVar,arglag = cbArgLag)
  v4<-c('dimnames','df','range','lag','argvar','arglag','class')
  dt1<-as.integer(round(mat1*100,0)) #convert to 100ths of the value in integer format to save space
  attributes(dt1)<-attributes(mat1)
  dt1
}

#Function for outputting model results into the models list
modFxn<-function(modChar='modelx') {
  list(call=get(modChar)$call,coef=get(modChar)$coefficients,vcov=summary(get(modChar))$cov.scaled,df.residual=get(modChar)$df.residual,df.null=get(modChar)$df.null,deviance=get(modChar)$deviance,null.deviance=get(modChar)$null.deviance,AIC=summary(get(modChar))$aic)
}

#Function for outputting model residuals
residModFxn<-function(modChar='modelx',cond=c1) {
  data.table(resid=resid(get(modChar)),super.zip=claims_master_agg[cond,super.zip],CLM_FROM_DT=claims_master_agg[cond,CLM_FROM_DT])
}

#Fucntion for outputting descriptives to evaluate for temporal correlation
residAcfFxn<-function(residMod=resid.modx,modChar='modelx') {
  dat1<-residMod[,mean(resid),by=CLM_FROM_DT]
  names(dat1)<-c('CLM_FROM_DT','resid.mean')
  dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
  par(mfrow=c(2,1))
  acf(dat2,main=modChar,ylim=c(-0.1,0.1))
  pacf(dat2,main=modChar)
  par(mfrow=c(1,1))
  cat('\n',modChar,'\nsum(acf):',
      round(sum(acf(dat2,plot=F)[[1]]),3),'     mean(acf):',
      round(mean(acf(dat2,plot=F)[[1]]),3),'\nsum(pacf):',
      round(sum(pacf(dat2,plot=F)[[1]]),3),'    mean(pacf):',
      round(mean(pacf(dat2,plot=F)[[1]]),3),'\n')
}



#####CREATE TMAX AND PRECIP AND PRECIP WAVE CROSSBASES####

##precip quantiles overall and max precip by county
quantile(claims_master_agg[,'PRISM_PPT'],c(0,0.5,0.75,0.9,0.99,1),na.rm=T)
#110 is 75th percentile, 350 is approx. 99th percentile and 250 is approx 90th percentile
min(claims_master_agg[,max(PRISM_PPT,na.rm=T),by=substr(super.zip,1,6)][,2])
#the lowest maximum in any county is 445, so all the counties have values above the overall 99th percentile

##precip crossbasis
cb1.ppt1<-cbFxn(expData=claims_master_agg,expNames=precipCols[1:21],cbArgVar=list('ns', knots=c(100,250,350)),cbArgLag=list('ns',knots=c(3,7)))

##temperature crossbasis:
cb1.tmax<-cbFxn(expData=claims_master_agg,expNames=tempCols[1:21],cbArgVar=list('ns', df = 4),cbArgLag=list('ns', knots=c(2,4,6)))

##delete the lag tmax and ppt variables in the master data set to save space
v1<-grep('\\_L',names(claims_master_agg),invert=T)
claims_master_agg<-claims_master_agg[,..v1]
rm(v1)

#save locally because it takes a few minutes to generate the above crossbases
#save(list='cb1.tmax',file='E://Madeline//cb1_tmax.RData')
#save(list='cb1.ppt1',file='E://Madeline//cb1_ppt.RData')

##precipitation waves crossbasis: number of consecutive days without rain

#dry spell, i.e., precip = 0
#2 and 5 are approximately equally spaced internal knots on the ln(#waves) scale of 0-12
cb1.ppt0in<-waveFxn(expData=claims_master_agg[,.(super.zip,CLM_FROM_DT,PRISM_PPT)],maxWave=12,waveCondition='(exp(PRISM_PPT/100)-1)/10/2.54==0',cbArgVar=list('ns',knots=c(2,5)),cbArgLag=list('ns',knots=c(3,7)))
#         0          1          2          3          4          5          6          7          8          9         10         11 
#35.7106510 18.1541795 13.1286391  9.3362322  6.4348833  4.5470500  3.2912776  2.4123491  1.8222844  1.3258906  0.9929952  0.6993960 
#         12       <NA> 
#  1.4262494  0.7179225 

#percentiles of precip in inches for guiding additional thresholds in commonly used thresholds on the local inch scale
claims_master_agg[,round(quantile((exp(PRISM_PPT/100)-1)/10/2.54,c(0.63,0.64,0.75,0.8,0.85,0.86,0.9,0.91,0.92,0.93,0.95,0.975,0.978,0.98)),2)]
#  63%   64%   75%   80%   85%   86%   90%   91%   92%   93%   95% 97.5% 97.8%   98% 
# 0.00  0.04  0.08  0.16  0.24  0.28  0.39  0.43  0.47  0.51  0.67  0.95  0.99  1.06
#0.25 inches is about the 85th percentile of ppt and 0.5 inches is about the 93rd percentile of ppt and 1 inch is approximately the 97.5th percentile

#0.25 inch threshold
#The lag period is again 21 days, and knots are the same as for continuous precip in the lag dimension.
cb1.ppt025in<-waveFxn(maxWave=3,expData=claims_master_agg[,.(super.zip,CLM_FROM_DT,PRISM_PPT)],waveCondition='(exp(PRISM_PPT/100)-1)/10/2.54>=0.25',cbArgVar=list('ns',knots=1),cbArgLag=list('ns',knots=c(3,7)))
#Percent of rows in each wave
#           0          1          2          3       <NA> 
#  85.2360446 10.7736557  2.6469627  0.6254146  0.7179225 

#0.5 inch threshold
cb1.ppt05in<-waveFxn(expData=claims_master_agg[,.(super.zip,CLM_FROM_DT,PRISM_PPT)],maxWave=2,waveCondition='(exp(PRISM_PPT/100)-1)/10/2.54>=0.5',cbArgVar=list('ns',knots=1),cbArgLag=list('ns',knots=c(3,7)))
#Percent of rows in eachwave
#           0          1          2       <NA> 
#  91.7888162  6.3930949  1.1001664  0.7179225 

#1 inch threshold
cb1.ppt1in<-waveFxn(maxWave=1,expData=claims_master_agg[,.(super.zip,CLM_FROM_DT,PRISM_PPT)],waveCondition='(exp(PRISM_PPT/100)-1)/10/2.54>=1',cbArgVar=list('strata',breaks=0.5),cbArgLag=list('ns',knots=c(3,7)))
#Percent of rows in each wave
#         0          1       <NA> 
#97.1166791  2.1653984  0.7179225

####CREATE THE TEMPORAL VARIABLES####
#create date spline with 6 df per year
ns.doy<-onebasis(as.POSIXlt(claims_master_agg$CLM_FROM_DT)$yday,'ns',df=6)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

##create a continuous year variable, assuming linear trend
claims_master_agg[,year:=as.integer(as.POSIXlt(CLM_FROM_DT)$year)]
#doing a spline just for date, rather than breaking up the effect into doy and year, creates a matrix that is too big for the model below. This approach has been used previously, e.g., Gasparrini & Hajat heat wave paper.

##create dow variable
claims_master_agg[,dow:=as.integer(as.POSIXlt(CLM_FROM_DT)$wday)]


####MERGE IN THE SPATIAL VARIABLES####
claims_master_agg<-superzip_to_1cnty[,.(super.zip,county)][claims_master_agg,on='super.zip']

####NAIVE GI DISEASE MODEL####
options(warn=1)

#models<-list()

#lag df=4, var df=4, lags=14
model1 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[,1:14]+cb1.tmax[,1:14]+ns.doy+year+as.factor(county), family=poisson(link='log'), data=claims_master_agg, model=F)
summary(model1)
#all the ppt crossbases have values when lag df=4 and var df=4. 3 of the 6 ppt crossbases are NA if var df=3 and lag df=3 or lag df=4. There are still significant effects at lag day 13, so need more lags.

models[['model1']]<-list(call=model1$call,coef=model1$coefficients,vcov=summary(model1)$cov.scaled,df.residual=model1$df.residual,df.null=model1$df.null,deviance=model1$deviance,null.deviance=model1$null.deviance,AIC=summary(model1)$aic)
rm(model1)

#lag df=4, var df=4, lags=21
model2 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1+cb1.tmax+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg, model=F)
summary(model2)
models[['model2']]<-list(call=model2$call,coef=model2$coefficients,vcov=summary(model2)$cov.scaled,df.residual=model2$df.residual,df.null=model2$df.null,deviance=model2$deviance,null.deviance=model2$null.deviance,AIC=summary(model2)$aic)
resid.mod2<-data.table(resid=resid(model2),super.zip=claims_master_agg[,super.zip],CLM_FROM_DT=claims_master_agg[,CLM_FROM_DT])
rm(model2)


####GI DISEASE SPATIAL DEPENDENCE####
#Examine semi-variogram to look for spatial correlation between ZIP codes
#This is done largely according to https://mgimond.github.io/Spatial/spatial-autocorrelation-in-r.html to analyze spatial dependence. This, in turn, is based on Roger S. Bivand, Virgilio Gomez-Rubio, Edzer Pebesma. 2013. Applied Spatial Data Analysis with R, Second Edition. New York: Springer.

dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][resid.mod2,on='super.zip']
dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
#there are still 824 unique X coordinates even if this is converted to an integer.
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]

#Only perform these diagnostics on a truly aggregated data set.
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2,xlim=c(0,3e4))
dat3<-vgm(psill=0.012, model="Exp", nugget=0.01, range=1.2e+04)
plot(dat2,model=dat3)
dat4<-fit.variogram(dat2,model=dat3)
plot(dat2,model=dat4)
plot(dat2,model=dat4,ylim=c(0,0.04),xlim=c(0,30000)) #Visually, there is autocorrelation up to about 10 km on the semivariogram.
#The range, or the distance up to which the values are autocorrelated, is 4146 m.
sum(dat2$np[dat2$dist<4146])/sum(dat2$np)*100
#This accounts for 0.17% of the data and likely is not substantially influencing the standard error.
#Look at Moran's I
#First, find the super zips within 30 km of each super zip. All super.zips have at least one zip within 30 km, and for nb2listw to work, there must be at least one neighbor.
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
length(nb1[!(nb1 %in% 0)]) #824 have neighbors within 30 km
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.19 (varies because the simulations draw random samples) when fully aggregated by super.zip only, so autocorrelation isn't significant and a spatial lag is not required. See DLNM_Analysis_211209.R for code for calculating a spatial lag.

####GI DISEASE TEMPORAL AUTOCORRELATION####
dat1<-resid.mod2[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #0.52 for model without dow, but there's still strong lag 1 correlation, and substantial lag 7, 14, etc. correlation. 3.1 with dow.
sum(pacf(dat2)[[1]]) #-0.047 for model without dow. 1.03 for model with dow.
#There is substantial positive residual correlation at lag days 1, 7, 14, 21, 28, etc. days without dow and negative correlation at lag days 2, 4, and 5. The rho for lag 1 was about 0.2. With dow, there is positive autocorrelation at lags 1, 3, 6 and 7, and the rho is about 0.17. Various efforts to increase df in the seasonal spline and change df in tmax and prcp exposure and lag dimensions failed. See the program Step3_DLNM_Analysis_210409.R for these details. 

#extract temporal lag
for (i in c(1:7)) {
  resid.mod2[,paste0('resid.lag',i):=stats::filter(resid,c(rep(0,i),1),sides=1)][CLM_FROM_DT %in% c(as.Date('2006-01-01'):as.Date(paste0('2006-01-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=resid]
}
resid.mod2[,resid.lag2_5:=rowMeans(.SD),.SDcols=paste0('resid.lag',2:5)]
resid.mod2[,resid.lag6_7:=rowMeans(.SD),.SDcols=paste0('resid.lag',6:7)]
resid.mod2[,resid.lag1.mean:=as.integer(mean(resid.lag1)*1000),by=CLM_FROM_DT]
resid.mod2[,resid.lag2_5.mean:=as.integer(mean(resid.lag2_5)*1000),by=CLM_FROM_DT]
resid.mod2[,resid.lag6_7.mean:=as.integer(mean(resid.lag6_7)*1000),by=CLM_FROM_DT]
for (i in c(1:10)) {resid.mod2[,paste0('resid.lag',i):=NULL]}
resid.mod2[,':='(resid.lag2_5=NULL,resid.lag6_7=NULL)]


model3<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1+cb1.tmax+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod2$resid.lag1.mean+resid.mod2$resid.lag2_5.mean+resid.mod2$resid.lag6_7.mean, family=poisson(link='log'),data=claims_master_agg,model=F)
models[['model3']]<-list(call=model3$call,coef=model3$coefficients,vcov=summary(model3)$cov.scaled,df.residual=model3$df.residual,df.null=model3$df.null,deviance=model3$deviance,null.deviance=model3$null.deviance,AIC=summary(model3)$aic)
resid.mod3<-data.table(resid=resid(model3),super.zip=claims_master_agg[,super.zip],CLM_FROM_DT=claims_master_agg[,CLM_FROM_DT])
rm(model3)
dat2<-resid.mod3[,mean(resid),by=CLM_FROM_DT]
names(dat2)<-c('CLM_FROM_DT','resid.mean')
dat3<-ts(data=dat2$resid.mean,start=min(dat2$CLM_FROM_DT),end=max(dat2$CLM_FROM_DT))
sum(acf(dat3,ylim=c(-0.05,0.05))[[1]]) #1.44
sum(pacf(dat3)[[1]]) #0.42
#No more significant residual temporal autocorrelation or obvious patterns.

v1<-grep('ppt1',names(models[['model2']]$coef),value=T)
cp.mod2<-crosspred(basis=cb1.ppt1,coef=models[['model2']]$coef[v1]*100,vcov=models[['model2']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=100,cumul=T)
v1<-grep('ppt1',names(models[['model3']]$coef),value=T)
cp.mod3<-crosspred(basis=cb1.ppt1,coef=models[['model3']]$coef[v1]*100,vcov=models[['model3']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=100,cumul=T)

save(models,file=paste0(f1,'models.RData'))
rm(resid.mod2,resid.mod3,cp.mod2,cp.mod3,dat2,dat3,v1,i)

####NAIVE RESPIRATORY DISEASE MODEL####
model4 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1+cb1.tmax+ns.doy+as.factor(dow)+log(year)+as.factor(county), family=quasipoisson(link='log'), data=claims_master_agg, model=F)
models[['model4']]<-list(call=model4$call,coef=model4$coefficients,vcov=summary(model4)$cov.scaled,df.residual=model4$df.residual,df.null=model4$df.null,deviance=model4$deviance,null.deviance=model4$null.deviance,AIC=summary(model4)$aic)
save(models,file=paste0(f1,'models.RData'))
resid.mod4<-data.table(resid=resid(model4),super.zip=claims_master_agg[,super.zip],CLM_FROM_DT=claims_master_agg[,CLM_FROM_DT])
rm(model4)
save(resid.mod4,file='resid_mod4.RData')

####RESPIRATORY DISEASE SPATIAL DEPENDENCE####
dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][resid.mod4,on='super.zip']

#aggregate residuals by super.zip and make the aggregated data set spatial
dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]

#plot the variograms to look for spatial autocorrelation
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2,xlim=c(0,3e4))
dat3<-vgm(psill=0.06, model="Exp", nugget=0.01, range=1.2e+04)
plot(dat2,model=dat3)
dat4<-fit.variogram(dat2,model=dat3)
plot(dat2,model=dat4)
plot(dat2,model=dat4,ylim=c(0,0.08),xlim=c(0,30000)) #There may be autocorrelation up to about 10 km on the semivariogram.
#The range, or the distance up to which the values are autocorrelated, is 4655.562
sum(dat2$np[dat2$dist<4652])/sum(dat2$np)*100
#This accounts for 0.26% of the data and likely is not substantially influencing the standard error.
#Look at Moran's I
#First, find the super zips within 30 km of each super zip. All super.zips have at least one zip within 30 km, and for nb2listw to work, there must be at least one neighbor.
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.97 (varies because the simulations draw random samples) when fully aggregated by super.zip only, so autocorrelation isn't significant. No spatial lag is necessary.
rm(dat1,dat2,dat3,dat4,ls2,nb1,nb3)

####RESPIRATORY TEMPORAL AUTOCORRELATION####
#Check the residuals for temporal autocorrelation.
dat1<-setDT(as.data.frame(do.call(rbind,lapply(claims_master_agg[,unique(super.zip)],FUN=function(x) {
  resid.mod4[claims_master_agg$super.zip==x,acf(resid,print=F)[[1]][,,1]]
}))))
names(dat1)<-paste0('lag',0:34)
dat1[,super.zip:=claims_master_agg[,unique(super.zip)]]
dat2<-melt(dat1[,-1],id.vars='super.zip',measure.vars=paste0('lag',1:34),variable='lag',value='r')
boxplot(r ~ lag,data=dat2)
lines(1:34,rep(0.036,34),col='blue',lty=2)
lines(1:34,rep(-0.036,34),col='blue',lty=2)
#Although none of the 75th percentiles for any lag exceed the critical value, on average, all the lags have medians > 0, so this suggests a small amount of autocorrelation.

dat3<-setDT(as.data.frame(do.call(rbind,lapply(claims_master_agg[,unique(super.zip)],FUN=function(x) {
  resid.mod4[claims_master_agg$super.zip==x,pacf(resid,plot=F)[[1]][,,1]]
}))))
names(dat3)<-paste0('lag',1:34)
dat3[,super.zip:=claims_master_agg[,unique(super.zip)]]
dat4<-melt(dat3,id.vars='super.zip',measure.vars=paste0('lag',1:34),variable='lag',value='r')
boxplot(r ~ lag,data=dat4)
lines(1:34,rep(0.036,34),col='blue',lty=2)
lines(1:34,rep(-0.036,34),col='blue',lty=2)
#Although none of the 75th percentiles for any lag exceed the critical value (except for lag 3), all the lags have medians > 0, so this suggests a small amount of autocorrelation. Therefore, we will work with the means of the residuals going forward, for which this autocorrelation is more detectable.
rm(dat3,dat4)

v2<-colMeans(dat1)
#the critical values for the individual acfs are -0.036 and +0.036
plot(0:34,v2,type='h',ylim=c(-0.1,0.1))
lines(0:34,rep(0,35))
lines(0:34,rep(0.036,35),col='blue',lty=2)
lines(0:34,rep(-0.036,35),col='blue',lty=2)

dat2<-resid.mod4[,mean(resid),by=CLM_FROM_DT]
names(dat2)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat2$resid.mean,start=min(dat2$CLM_FROM_DT),end=max(dat2$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #7.92
sum(pacf(dat2)[[1]]) #1.41
#There is substantial positive residual correlation and partial autocorrelation at lag days 1-7. Changing the ns.doy to week-of-year indicators resulted in even more temporal autocorrelation. Adding weekends, dropping ns.doy, changing ns.doy dfs, changing ppt1 dfs and tmax dfs, adding moving averages of NEW_COUNT or doy--nothing worked. Therefore, add residual lag.
v1<-pacf(dat2)[[1]][1:15,,1] #these are the partial autocorrelation coefficients for lags 1-15

#extract temporal lag by super.zip. Because the data set is ordered by super.zip and date, this code quickly accomplishes this.
#First, create individual residual lags, multiplied by their above corresponding pacf coefficient. Note that the values are all converted to integers to save space. They are first multiplied by 1e5 to preserve some of the precision.
for (i in c(1:15)) {
  resid.mod4[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1))*1e7][CLM_FROM_DT %in% c(as.Date('2006-01-01'):as.Date(paste0('2006-01-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}

#Create an inverse-day-weighted moving average of the residual lags.
resid.mod4[,resid.lag1_15:=as.integer(stats::filter(resid,1/exp(0:15/5),sides=1)*1e7)][is.na(resid.lag1_15),resid.lag1_15:=as.integer(resid*1e7)][,resid.lag1_15.mean:=as.integer(mean(resid.lag1_15)),by=CLM_FROM_DT]
#Drop unused residual lags to save space.
for (i in c('1_15')) {resid.mod4[,paste0('resid.lag',i):=NULL]}
rm(i)

#Redo the model with residual lags added.
#c1<-claims_master_agg[,county %in% c(26163,39035,42003)] #largest county in each state, useful for troubleshooting on a smaller set
c1<-1:nrow(claims_master_agg)
model5<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy[c1,]+log(year)+as.factor(county)+as.factor(dow)+resid.mod4[c1,resid.lag1_15.mean], family=poisson(link='log'),data=claims_master_agg[c1,],model=F)
models[['model5']]<-list(call=model5$call,coef=model5$coefficients,vcov=summary(model5)$cov.scaled,df.residual=model5$df.residual,df.null=model5$df.null,deviance=model5$deviance,null.deviance=model5$null.deviance,AIC=summary(model5)$aic)
save(models,file=paste0(f1,'models.RData'))
resid.mod5<-data.table(resid=resid(model5),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
rm(model5)

#Check for residual temporal autocorrelation.
dat2<-resid.mod5[,mean(resid),by=CLM_FROM_DT]
names(dat2)<-c('CLM_FROM_DT','resid.mean')
dat3<-ts(data=dat2$resid.mean,start=min(dat2$CLM_FROM_DT),end=max(dat2$CLM_FROM_DT))
sum(acf(dat3)[[1]]) #1.499. Still significant at lags 1 and 7.
sum(pacf(dat3)[[1]]) #0.387 Still significant at lags 1, 6 and 7.

#Re-run model with additional residual lags
for (i in c(1,6,7)) {
  resid.mod5[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1))*1e7][CLM_FROM_DT %in% c(as.Date('2006-01-01'):as.Date(paste0('2006-01-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',c(1,6,7))
resid.mod5[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod5[,(v1):=NULL]
rm(i,v1)

c1<-1:nrow(claims_master_agg)
model6<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy[c1,]+log(year)+as.factor(county)+as.factor(dow)+resid.mod4[c1,resid.lag1_15.mean]+resid.mod5[c1,resid.lag1]+resid.mod5[c1,resid.lag6]+resid.mod5[c1,resid.lag7], family=poisson(link='log'),data=claims_master_agg[c1,],model=F)
models[['model6']]<-list(call=model6$call,coef=model6$coefficients,vcov=summary(model6)$cov.scaled,df.residual=model6$df.residual,df.null=model6$df.null,deviance=model6$deviance,null.deviance=model6$null.deviance,AIC=summary(model6)$aic)
save(models,file=paste0(f1,'models.RData'))
resid.mod6<-data.table(resid=resid(model6),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
rm(model6)

#Check for residual temporal autocorrelation.
dat2<-resid.mod6[,mean(resid),by=CLM_FROM_DT]
names(dat2)<-c('CLM_FROM_DT','resid.mean')
dat3<-ts(data=dat2$resid.mean,start=min(dat2$CLM_FROM_DT),end=max(dat2$CLM_FROM_DT))
sum(acf(dat3)[[1]]) #1.32
mean(acf(dat3)[[1]]) #0.0368
sum(pacf(dat3)[[1]]) #0.299
mean(pacf(dat3)[[1]]) #0.00878
#These look really good. There is no remaining significant residual autocorrelation, except at lags 27-29, and the mean acf is low.
rm(dat2,dat3)


####GI WARM SEASON####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(5,6,7,8) #5=June, 6=July, 7=August, 8=September

#create new season spline with 3 df
ns.doy<-onebasis(sin(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday/365*2*3.14),'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model9 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod9<-data.table(resid=resid(model9),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model9']]<-list(call=model9$call,coef=model9$coefficients,vcov=summary(model9)$cov.scaled,df.residual=model9$df.residual,df.null=model9$df.null,deviance=model9$deviance,null.deviance=model9$null.deviance,AIC=summary(model9)$aic)
rm(model9)

#warm season spatial dependence
dat1<-data.table(resid=resid.mod9$resid,super.zip=claims_master_agg[c1,super.zip],year=claims_master_agg[c1,year],doy=as.POSIXlt(claims_master_agg$CLM_FROM_DT[c1])$yday,CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][dat1,on='super.zip']

dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2)
#Look at Moran's I
#First, find the super zips within 30 km of each super zip. All super.zips have at least one zip within 30 km, and for nb2listw to work, there must be at least one neighbor.
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.4, so spatial dependence not significant
rm(ls2,nb3,nb1,dat2,dat1)

#warm season temporal dependence
dat1<-resid.mod9[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #1.635
sum(pacf(dat2,main='model 10')[[1]]) #0.498
mean(pacf(dat2)[[1]]) #0.015
#Significant correlation at lag 2 (negative), 5, 6, and 13.
rm(dat1,dat2)

#re-run model with residual lags
for (i in c(2, 5, 6, 13)) {
  resid.mod9[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-06-01'):as.Date(paste0('2006-06-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',c(2,5,6,13))
resid.mod9[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod9[,(v1):=NULL]
rm(i,v1)

model10 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod9[,resid.lag2.mean]+resid.mod9[,resid.lag5.mean]+resid.mod9[,resid.lag6.mean]+resid.mod9[,resid.lag13.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod10<-data.table(resid=resid(model10),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model10']]<-list(call=model10$call,coef=model10$coefficients,vcov=summary(model10)$cov.scaled,df.residual=model10$df.residual,df.null=model10$df.null,deviance=model10$deviance,null.deviance=model10$null.deviance,AIC=summary(model10)$aic)
rm(model10)

dat1<-resid.mod10[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #1.41
sum(pacf(dat2,main='model 10')[[1]]) #0.36
mean(pacf(dat2)[[1]]) #0.011
#Looks pretty good.

#warm season effect estimates
v1<-grep('ppt1',names(models[['model10']]$coef),value=T)
cp.mod10<-crosspred(basis=cb1.ppt1,coef=models[['model10']]$coef[v1]*100,vcov=models[['model10']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=110,cumul=T)


####GI COLD SEASON####
#cold temperatures and ice or snow precipitation in November, December, January, and February
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(10,11,0,1)

#create new season spline with 3 df
#need the sin function in there because 10,11,0,1 is disjoint otherwise
ns.doy<-onebasis(sin(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday/354*3.14*2),'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model11<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod11<-data.table(resid=resid(model11),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model11']]<-list(call=model11$call,coef=model11$coefficients,vcov=summary(model11)$cov.scaled,df.residual=model11$df.residual,df.null=model11$df.null,deviance=model11$deviance,null.deviance=model11$null.deviance,AIC=summary(model11)$aic)
rm(model11)

#cold season spatial dependence
dat1<-data.table(resid=resid.mod11$resid,super.zip=claims_master_agg[c1,super.zip],year=claims_master_agg[c1,year],doy=as.POSIXlt(claims_master_agg$CLM_FROM_DT[c1])$yday,CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][dat1,on='super.zip']
dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2)
#Look at Moran's I
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.19, so spatial dependence not significant
rm(ls2,nb3,nb1,dat2,dat1)

#cold season temporal dependence
dat1<-resid.mod11[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #3.82
sum(pacf(dat2,main='model 10')[[1]]) #1.11
mean(pacf(dat2)[[1]]) #0.033
#Temporal dependence is present. There is significant autocorrelation at lag 1, 6, and 7.

#re-run model with residual lags
for (i in c(1,3,4,6,7,8,12)) {
  resid.mod11[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-01-01'):as.Date(paste0('2006-01-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',c(1,3,4,6,7,8,12))
resid.mod11[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod11[,resid.lag2_12.mean:=as.integer(stats::filter(resid,c(0,0,rep(1,11)),sides=1)*1e7)][is.na(resid.lag2_12.mean),resid.lag2_12.mean:=0L]
resid.mod11[,(v1):=NULL]
rm(i,v1)

model12<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+year+as.factor(county)+as.factor(dow)+resid.mod11[,resid.lag1.mean]+resid.mod11[,resid.lag2_12.mean]+resid.mod11[,resid.lag6.mean]+resid.mod11[,resid.lag7.mean]+resid.mod11[,resid.lag12.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod12<-data.table(resid=resid(model12),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model12']]<-list(call=model12$call,coef=model12$coefficients,vcov=summary(model12)$cov.scaled,df.residual=model12$df.residual,df.null=model12$df.null,deviance=model12$deviance,null.deviance=model12$null.deviance,AIC=summary(model12)$aic)
rm(model12)

#temporal dependence
dat1<-resid.mod12[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #1.28
sum(pacf(dat2,main='model 12')[[1]]) #0.257
mean(pacf(dat2)[[1]]) #0.007
rm(resid.mod12,dat1,dat2)
#Looks good.

#cold season effect estimates
v1<-grep('ppt1',names(models[['model12']]$coef),value=T)
cp.mod12<-crosspred(basis=cb1.ppt1,coef=models[['model12']]$coef[v1]*100,vcov=models[['model12']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=110,cumul=T)


####GI SPRING MELT SEASON####
#melt waters and spring rainy season
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(2,3,4) #March, April, May

#create new season spline with 3 df
ns.doy<-onebasis(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday,'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model13 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod13<-data.table(resid=resid(model13),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model13']]<-list(call=model13$call,coef=model13$coefficients,vcov=summary(model13)$cov.scaled,df.residual=model13$df.residual,df.null=model13$df.null,deviance=model13$deviance,null.deviance=model13$null.deviance,AIC=summary(model13)$aic)
rm(model13)
save(models,file=paste0(f1,'models.RData'))

#melt season spatial dependence
dat1<-data.table(resid=resid.mod13$resid,super.zip=claims_master_agg[c1,super.zip],year=claims_master_agg[c1,year],doy=as.POSIXlt(claims_master_agg$CLM_FROM_DT[c1])$yday,CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][dat1,on='super.zip']
dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2)
#Look at Moran's I
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.29, so spatial dependence not significant
rm(ls2,nb3,nb1,dat2,dat1)

#temporal dependence
dat1<-resid.mod13[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #1.66
sum(pacf(dat2,main='model 13')[[1]]) #0.52
mean(pacf(dat2)[[1]]) #0.0154
rm(dat1,dat2)
#significance at lags 3, 5, 6, 7, 8, and 15

#re-run with lags
for (i in c(15)) {
  resid.mod13[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-01-03'):as.Date(paste0('2006-03-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',c(15))
resid.mod13[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod13[,resid.lag1_8.mean:=as.integer(stats::filter(resid,c(0,rep(1,8)),sides=1)*1e7)][is.na(resid.lag1_8.mean),resid.lag1_8.mean:=0L]
resid.mod13[,(v1):=NULL]
rm(i,v1)

model14 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod13[,resid.lag1_8.mean]+resid.mod13[,resid.lag15.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod14<-data.table(resid=resid(model14),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model14']]<-list(call=model14$call,coef=model14$coefficients,vcov=summary(model14)$cov.scaled,df.residual=model14$df.residual,df.null=model14$df.null,deviance=model14$deviance,null.deviance=model14$null.deviance,AIC=summary(model14)$aic)
rm(model14)

#temporal dependence
dat1<-resid.mod14[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #1.145
sum(pacf(dat2,main='model 13')[[1]]) #0.16
mean(pacf(dat2)[[1]]) #0.0047
rm(dat1,dat2,resid.mod14)
#Looks good.

#melt season effect estimates
v1<-grep('ppt1',names(models[['model14']]$coef),value=T)
cp.mod14<-crosspred(basis=cb1.ppt1,coef=models[['model14']]$coef[v1]*100,vcov=models[['model14']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=110,cumul=T)


####RESP WARM SEASON####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(5,6,7,8) #5=June, 6=July, 7=August, 8=September

#create new season spline with 3 df
ns.doy<-onebasis(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday,'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model15 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod15<-data.table(resid=resid(model15),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model15']]<-list(call=model15$call,coef=model15$coefficients,vcov=summary(model15)$cov.scaled,df.residual=model15$df.residual,df.null=model15$df.null,deviance=model15$deviance,null.deviance=model15$null.deviance,AIC=summary(model15)$aic)
rm(model15)

#warm season spatial dependence
dat1<-data.table(resid=resid.mod15$resid,super.zip=claims_master_agg[c1,super.zip],year=claims_master_agg[c1,year],doy=as.POSIXlt(claims_master_agg$CLM_FROM_DT[c1])$yday,CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][dat1,on='super.zip']

dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2)
#Look at Moran's I
#First, find the super zips within 30 km of each super zip. All super.zips have at least one zip within 30 km, and for nb2listw to work, there must be at least one neighbor.
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.975, so spatial dependence not significant
rm(ls2,nb3,nb1,dat2,dat1)

#warm season temporal dependence
dat1<-resid.mod15[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #4.30
sum(pacf(dat2,main='model 15')[[1]]) #1.36
mean(pacf(dat2,ylim=c(-0.1,0.2))[[1]]) #0.04
#Significant correlation at lags 1-19.
rm(dat1,dat2)

#re-run model with residual lags
for (i in c(1)) {
  resid.mod15[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-06-01'):as.Date(paste0('2006-06-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',c(1))
resid.mod15[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod15[,resid.lag2_19:=as.integer(stats::filter(resid,c(0,0,2/2:19)*1e7))][is.na(resid.lag2_19),resid.lag2_19:=0L]
resid.mod15[,resid.lag2_19.mean:=mean(resid.lag2_19),by=CLM_FROM_DT]

model16 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod15[,resid.lag1.mean]+resid.mod15[,resid.lag2_19.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod16<-data.table(resid=resid(model16),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model16']]<-list(call=model16$call,coef=model16$coefficients,vcov=summary(model16)$cov.scaled,df.residual=model16$df.residual,df.null=model16$df.null,deviance=model16$deviance,null.deviance=model16$null.deviance,AIC=summary(model16)$aic)
rm(model16)

dat1<-resid.mod16[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #0.87
sum(pacf(dat2)[[1]]) #-0.12
mean(pacf(dat2)[[1]]) #-0.0034
#Negative correlation at lags 2, 7 and 8.

#re-run model with residual lags
v2<-c(2)
for (i in v2) {
  resid.mod16[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-06-01'):as.Date(paste0('2006-06-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',v2)
resid.mod16[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod16[,resid.lag3_8:=as.integer(stats::filter(resid,c(rep(0,2),-1/7:1/8),sides=1)*1e7)][is.na(resid.lag3_8),resid.lag3_8:=0L][,resid.lag3_8.mean:=mean(resid.lag3_8),by=CLM_FROM_DT]
resid.mod16[,(v1):=NULL]
rm(i,v1,v2)

model17 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod15[,resid.lag1.mean]+resid.mod15[,resid.lag2_19.mean]+resid.mod16[,resid.lag2.mean]+resid.mod16[,resid.lag3_8.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod17<-data.table(resid=resid(model17),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model17']]<-list(call=model17$call,coef=model17$coefficients,vcov=summary(model17)$cov.scaled,df.residual=model17$df.residual,df.null=model17$df.null,deviance=model17$deviance,null.deviance=model17$null.deviance,AIC=summary(model17)$aic)
rm(model17)

dat1<-resid.mod17[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #1.06
sum(pacf(dat2)[[1]]) #0.108
mean(pacf(dat2)[[1]]) #0.003
#Looks good.
rm(dat1,dat2)

save(models,file=paste0(f1,'models.RData'))

#warm season effect estimates
v1<-grep('ppt1',names(models[['model17']]$coef),value=T)
cp.mod17<-crosspred(basis=cb1.ppt1,coef=models[['model17']]$coef[v1]*100,vcov=models[['model17']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=110,cumul=T)

STOPPED HERE

####RESP COLD SEASON####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(10,11,0,1)

#create new season spline with 3 df
ns.doy<-onebasis(sin(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday/354*3.14*2),'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model18<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod18<-data.table(resid=resid(model18),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model18']]<-list(call=model18$call,coef=model18$coefficients,vcov=summary(model18)$cov.scaled,df.residual=model18$df.residual,df.null=model18$df.null,deviance=model18$deviance,null.deviance=model18$null.deviance,AIC=summary(model18)$aic)
rm(model18)

#cold season spatial dependence
dat1<-data.table(resid=resid.mod18$resid,super.zip=claims_master_agg[c1,super.zip],year=claims_master_agg[c1,year],doy=as.POSIXlt(claims_master_agg$CLM_FROM_DT[c1])$yday,CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][dat1,on='super.zip']
dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2)
#Look at Moran's I
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.97, so spatial dependence not significant
rm(ls2,nb3,nb1,dat2,dat1)

#cold season temporal dependence
dat1<-resid.mod18[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #8.4
sum(pacf(dat2,main='model 10')[[1]]) #1.44
mean(pacf(dat2)[[1]]) #0.042
#Temporal dependence is present. There is significant autocorrelation at lags 1-7.
rm(dat1,dat2)

resid.mod18[,resid.lag1_7:=as.integer(stats::filter(resid,c(0,1/1:7/2),sides=1)*1e7)][is.na(resid.lag1_7),resid.lag1_7:=0L][,resid.lag1_7.mean:=mean(resid.lag1_7),by=CLM_FROM_DT]

model19<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod18[, resid.lag1_7.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod19<-data.table(resid=resid(model19),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model19']]<-list(call=model19$call,coef=model19$coefficients,vcov=summary(model19)$cov.scaled,df.residual=model19$df.residual,df.null=model19$df.null,deviance=model19$deviance,null.deviance=model19$null.deviance,AIC=summary(model19)$aic)
rm(model19)

#redo with additional lags
v2<-c(1,5,7)
v1<-paste0('resid.lag',v2)
for (i in v2) {
  resid.mod19[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-01-01'):as.Date(paste0('2006-01-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
resid.mod19[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod19[,(v1):=NULL]
rm(i,v1,v2)

model20<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod18[, resid.lag1_7.mean]+resid.mod19[, resid.lag1.mean]+resid.mod19[, resid.lag5.mean]+resid.mod19[, resid.lag7.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod20<-data.table(resid=resid(model20),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model20']]<-list(call=model20$call,coef=model20$coefficients,vcov=summary(model20)$cov.scaled,df.residual=model20$df.residual,df.null=model20$df.null,deviance=model20$deviance,null.deviance=model20$null.deviance,AIC=summary(model20)$aic)
rm(model20)

save(models,file=paste0(f1,'models.RData'))

#effect estimates
v1<-grep('ppt1',names(models[['model20']]$coef),value=T)

cp.mod20<-crosspred(basis=cb1.ppt1,coef=models[['model20']]$coef[v1]*100,vcov=models[['model20']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=110,cumul=T)

####RESP MELT SEASON####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(2,3,4) #March, April, May

#create new season spline with 3 df
ns.doy<-onebasis(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday,'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model21<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod21<-data.table(resid=resid(model21),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model21']]<-list(call=model21$call,coef=model21$coefficients,vcov=summary(model21)$cov.scaled,df.residual=model21$df.residual,df.null=model21$df.null,deviance=model21$deviance,null.deviance=model21$null.deviance,AIC=summary(model21)$aic)
rm(model21)

#spatial dependence
dat1<-data.table(resid=resid.mod21$resid,super.zip=claims_master_agg[c1,super.zip],year=claims_master_agg[c1,year],doy=as.POSIXlt(claims_master_agg$CLM_FROM_DT[c1])$yday,CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
dat1<-superzip_to_1cnty[,c('super.zip','XCoord.mean','YCoord.mean')][dat1,on='super.zip']
dat1<-dat1[,mean(resid),by=.(super.zip,XCoord.mean,YCoord.mean)]
names(dat1)[names(dat1)=='V1']<-'resid.mean'
dat1[,c('XCoord.mean','YCoord.mean'):=lapply(.SD,FUN=as.integer),.SDcols=c('XCoord.mean','YCoord.mean')]
coordinates(dat1)<-dat1[,.(XCoord.mean,YCoord.mean)]
dat2<-variogram(resid.mean~1, data=dat1,width=100)
plot(dat2)
#Look at Moran's I
nb1<-dnearneigh(dat1,0,30000,row.names=dat1$super.zip)
nb3<-nbdists(nb1,coordinates(dat1)) #find the distances of each neighbor from that super.zip
ls2<-nb2listw(nb1,glist=nb3,style='B') #put in a list
moran.mc(dat1$resid.mean, ls2, nsim=599,zero.policy=T)
#p ~ 0.94, so spatial dependence not significant
rm(ls2,nb3,nb1,dat2,dat1)

#temporal dependence
dat1<-resid.mod21[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #5.1
sum(pacf(dat2,main='model 10')[[1]]) #1.27
mean(pacf(dat2)[[1]]) #0.037
#Temporal dependence is present. There is significant autocorrelation at lags 1-7.
rm(dat1,dat2)

#re-run model with residual lags
resid.mod21[,resid.lag1_7:=as.integer(stats::filter(resid,c(0,1/1:7/2),sides=1)*1e7)][is.na(resid.lag1_7),resid.lag1_7:=0L][,resid.lag1_7.mean:=mean(resid.lag1_7),by=CLM_FROM_DT]

model22<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod21[,resid.lag1_7.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod22<-data.table(resid=resid(model22),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model22']]<-list(call=model22$call,coef=model22$coefficients,vcov=summary(model22)$cov.scaled,df.residual=model22$df.residual,df.null=model22$df.null,deviance=model22$deviance,null.deviance=model22$null.deviance,AIC=summary(model22)$aic)
rm(model22)

#temporal dependence
dat1<-resid.mod22[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #1.99
sum(pacf(dat2,main='model 22')[[1]]) #0.77
mean(pacf(dat2)[[1]]) #0.022
rm(dat1,dat2)
#still residual lag 1, 3, and 7-8 autocorrelation.

#redo with additional lags
v2<-c(1,3,7,8)
v1<-paste0('resid.lag',v2)
for (i in v2) {
  resid.mod22[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-01-01'):as.Date(paste0('2006-01-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
resid.mod22[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod22[,(v1):=NULL]
rm(i,v1,v2)

model23<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+cb1.tmax[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod21[,resid.lag1_7.mean]+resid.mod22[,resid.lag1.mean]+resid.mod22[,resid.lag3.mean]+resid.mod22[,resid.lag7.mean]+resid.mod22[,resid.lag8.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod23<-data.table(resid=resid(model23),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model23']]<-list(call=model23$call,coef=model23$coefficients,vcov=summary(model23)$cov.scaled,df.residual=model23$df.residual,df.null=model23$df.null,deviance=model23$deviance,null.deviance=model23$null.deviance,AIC=summary(model23)$aic)
rm(model23)

#temporal dependence
dat1<-resid.mod23[,mean(resid),by=CLM_FROM_DT]
names(dat1)<-c('CLM_FROM_DT','resid.mean')
dat2<-ts(data=dat1$resid.mean,start=min(dat1$CLM_FROM_DT),end=max(dat1$CLM_FROM_DT))
sum(acf(dat2)[[1]]) #0.77
sum(pacf(dat2,main='model 12')[[1]]) #-0.203
mean(pacf(dat2)[[1]]) #-0.0060
rm(dat1,dat2)
#fine.

save(models,file=paste0(f1,'models.RData'))

#melt season effect estimates
v1<-grep('ppt1',names(models[['model23']]$coef),value=T)

cp.mod23<-crosspred(basis=cb1.ppt1,coef=models[['model23']]$coef[v1]*100,vcov=models[['model23']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=110,cumul=T)


####GI WARM SEASON WITH PRECIP WAVES####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(5,6,7,8) #5=June, 6=July, 7=August, 8=September

#create new season spline with 3 df
ns.doy<-onebasis(sin(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday/365*2*3.14),'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model31 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt0in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod31<-residModFxn('model31')
models[['model31']]<-modFxn('model31')
rm(model31)

#warm season temporal dependence
residAcfFxn(residMod=resid.mod31,'model31')
#sum(acf): 1.6      mean(acf): 0.046 
#sum(pacf): 0.501     mean(pacf): 0.015 
#Significant correlation at lag 2 (neg), 5, 6 and 13
#this does not change regardless of whether cb1.ppt1 is in the model or not or whether ns.doy is in the model or not or whether cb.ppt025in is in the model

#re-run model with residual lags
for (i in c(2,5,6,13)) {
  resid.mod31[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-06-01'):as.Date(paste0('2006-06-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',c(2,5,6,13))
resid.mod31[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod31[,(v1):=NULL]
rm(i,v1)

model45 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt0in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod31[,resid.lag2.mean]+resid.mod31[,resid.lag5.mean]+resid.mod31[,resid.lag6.mean]+resid.mod31[,resid.lag13.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod45<-residModFxn('model45')
models[['model45']]<-modFxn('model45')
rm(model45)

residAcfFxn(residMod=resid.mod45,'model45')
#sum(acf): 1.412      mean(acf): 0.04 
#sum(pacf): 0.366     mean(pacf): 0.011 
#Looks pretty good.


model46 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt025in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod46<-residModFxn('model46')
models[['model46']]<-modFxn('model46')
rm(model46)

residAcfFxn(resid.mod46,'model46')
#sum(acf): 1.522      mean(acf): 0.043 
#sum(pacf): 0.442     mean(pacf): 0.013 
#Fine considering that adding the residuals in the previous models didn't change the RRs or CIs.


model47 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt05in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod47<-residModFxn('model47')
models[['model47']]<-modFxn('model47')
rm(model47)

residAcfFxn(resid.mod47,'model47')
#sum(acf): 1.657      mean(acf): 0.047 
#sum(pacf): 0.541     mean(pacf): 0.016 
#Fine.
rm(resid.mod47)


model48 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt1in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod48<-residModFxn('model48')
models[['model48']]<-modFxn('model48')
rm(model48)

residAcfFxn(resid.mod48,'model48')
#sum(acf): 1.687      mean(acf): 0.048 
#sum(pacf): 0.56     mean(pacf): 0.016 
#Fine.

save(models,file=paste0(f1,'models.RData'))

wavemodels<-list(GI.Warm=c('model45','model46','model47','model48'))

####GI COLD SEASON WITH PRECIP WAVES####

#cold temperatures and ice or snow precipitation in November, December, January, and February
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(10,11,0,1)
c2<-!is.na(cb1.ppt0in[,1])

#create new season spline with 3 df
#need the sin function in there because 10,11,0,1 is disjoint otherwise
ns.doy<-onebasis(sin(as.POSIXlt(claims_master_agg[c1 & c2,CLM_FROM_DT])$yday/354*3.14*2),'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model32<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt0in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod32<-residModFxn('model32',cond=c1 & c2)
models[['model32']]<-modFxn('model32')
rm(model32)

#temporal dependence
residAcfFxn(resid.mod32,'model32')
#sum(acf): 3.699      mean(acf): 0.106 
#sum(pacf): 1.036     mean(pacf): 0.03
#Significant correlation at all lags and significant partial autocorrelation at lags 1-7

#re-run model with residual lags
for (i in c(1,7)) {
  resid.mod32[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-01-01'):as.Date(paste0('2006-01-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
v1<-paste0('resid.lag',c(1,7))
resid.mod32[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
rm(i,v1)
#create new condition with these residuals
test<-claims_master_agg[,.(super.zip,CLM_FROM_DT)]
test<-resid.mod32[test,on=c('super.zip','CLM_FROM_DT')]
c3<-test[,!is.na(resid.lag1) & !is.na(resid.lag7.mean)]
rm(test)

model33<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt0in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod32[,resid.lag1.mean]+resid.mod32[,resid.lag7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod33<-residModFxn('model33',cond=c1 & c2 & c3)
models[['model33']]<-modFxn('model33')
rm(model33)

#temporal dependence
residAcfFxn(resid.mod33,'model33')
#Looks better.
#sum(acf): 1.786      mean(acf): 0.051 
#sum(pacf): 0.593     mean(pacf): 0.017
rm(resid.mod33)


model49<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt025in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod32[,resid.lag1.mean]+resid.mod32[,resid.lag7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod49<-residModFxn('model49',cond=c1 & c2 & c3)
models[['model49']]<-modFxn('model49')
rm(model49)

#temporal dependence
residAcfFxn(resid.mod49,'model49')
#sum(acf): 1.786      mean(acf): 0.051 
#sum(pacf): 0.59     mean(pacf): 0.017
#Fine.
rm(resid.mod49)


model50<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt05in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod32[,resid.lag1.mean]+resid.mod32[,resid.lag7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod50<-residModFxn('model50',cond=c1 & c2 & c3)
models[['model50']]<-modFxn('model50')
rm(model50)

#temporal dependence
residAcfFxn(resid.mod50,'model50')
#sum(acf): 1.758      mean(acf): 0.05 
#sum(pacf): 0.577     mean(pacf): 0.017 
#Fine.
rm(resid.mod50)

model51<-glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt1in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod32[,resid.lag1.mean]+resid.mod32[,resid.lag7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod51<-residModFxn('model51',cond=c1 & c2 & c3)
models[['model51']]<-modFxn('model51')
rm(model51)

#temporal dependence
residAcfFxn(resid.mod51,'model51')
#sum(acf): 1.773      mean(acf): 0.051 
#sum(pacf): 0.582     mean(pacf): 0.017
#Fine.
rm(resid.mod51)

save(models,file=paste0(f1,'models.RData'))

wavemodels[['GI.Cold']]<-c('model33','model49','model50','model51')


####GI MELT SEASON WITH PRECIP WAVES####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(2,3,4) #March, April, May

#create new season spline with 3 df
ns.doy<-onebasis(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday,'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model34 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt0in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod34<-residModFxn('model34')
models[['model34']]<-modFxn('model34')
rm(model34)

residAcfFxn(resid.mod34,'model34')
#sum(acf): 1.615      mean(acf): 0.046 
#sum(pacf): 0.495     mean(pacf): 0.015 
#Fine.


model35 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt025in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod35<-residModFxn('model35')
models[['model35']]<-modFxn('model35')
rm(model35)

residAcfFxn(resid.mod35,'model35')
#sum(acf): 1.631      mean(acf): 0.047 
#sum(pacf): 0.501     mean(pacf): 0.015
#Fine.


model52 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt05in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod52<-residModFxn('model52')
models[['model52']]<-modFxn('model52')
rm(model52)

residAcfFxn(resid.mod52,'model52')
#sum(acf): 1.643      mean(acf): 0.047 
#sum(pacf): 0.515     mean(pacf): 0.015 
#Fine.


model53 <- glm(gi_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt1in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
#Warning message:
#In doTryCatch(return(expr), name, parentenv, handler) : restarting interrupted promise evaluation
resid.mod53<-residModFxn('model53')
models[['model53']]<-modFxn('model53')
rm(model53)

residAcfFxn(resid.mod53,'model53')
#sum(acf): 1.62      mean(acf): 0.046 
#sum(pacf): 0.493     mean(pacf): 0.014


wavemodels[['GI.Melt']]<-c('model34','model35','model52','model53')

save(models,file=paste0(f1,'models.RData'))

####RESP WARM SEASON WITH PRECIP WAVES####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(5,6,7,8) #5=June, 6=July, 7=August, 8=September

#create new season spline with 3 df
ns.doy<-onebasis(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday,'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model36 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt0in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod36<-residModFxn('model36')
models[['model36']]<-modFxn('model36')
rm(model36)

#temporal dependence
residAcfFxn(resid.mod36,'model36')
#sum(acf): 4.405      mean(acf): 0.126 
#sum(pacf): 1.38     mean(pacf): 0.041
#Significant correlation at all lags and pacf at lags 1-19.

#re-run model with residual lags
for (i in c(1,2)) {
  resid.mod36[,paste0('resid.lag',i):=as.integer(stats::filter(resid,c(rep(0,i),1),sides=1)*1e7)][CLM_FROM_DT %in% c(as.Date('2006-06-01'):as.Date(paste0('2006-06-',sprintf('%02.0f',i)))),paste0('resid.lag',i):=0L]
}
resid.mod36[,paste0(v1,'.mean'):=lapply(.SD,FUN=function(x) {as.integer(mean(x))}),by=CLM_FROM_DT,.SDcols=v1]
resid.mod36[,resid.lag1_19:=as.integer(stats::filter(resid,c(0,1/1:19)*1e7))][is.na(resid.lag1_19),resid.lag1_19:=0L]
resid.mod36[,resid.lag1_19.mean:=mean(resid.lag1_19),by=CLM_FROM_DT]
resid.mod36[,resid.lag1.mean:=mean(resid.lag1),by=CLM_FROM_DT]
resid.mod36[,resid.lag2.mean:=mean(resid.lag2),by=CLM_FROM_DT]

model37 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt0in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod36[,resid.lag1.mean]+resid.mod36[,resid.lag2.mean]+resid.mod36[,resid.lag1_19.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod37<-residModFxn('model37')
models[['model37']]<-modFxn('model37')
rm(model37)

residAcfFxn(resid.mod37,'model37')
#sum(acf): 1.3      mean(acf): 0.037 
#sum(pacf): 0.368     mean(pacf): 0.011
#Fine.


model38 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt025in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod36[,resid.lag1.mean]+resid.mod36[,resid.lag2.mean]+resid.mod36[,resid.lag1_19.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod38<-residModFxn('model38')
models[['model38']]<-modFxn('model38')

residAcfFxn(resid.mod38,'model38')
#sum(acf): 1.287      mean(acf): 0.037 
#sum(pacf): 0.341     mean(pacf): 0.01
#Fine.
rm(model38,resid.mod38)


model54 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt05in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod36[,resid.lag1.mean]+resid.mod36[,resid.lag2.mean]+resid.mod36[,resid.lag1_19.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod54<-residModFxn('model54')
models[['model54']]<-modFxn('model54')

residAcfFxn(resid.mod54,'model54')
#sum(acf): 1.643      mean(acf): 0.047 
#sum(pacf): 0.515     mean(pacf): 0.015 
#Fine.
rm(resid.mod54,model54)


model55 <- glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt1in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod36[,resid.lag1.mean]+resid.mod36[,resid.lag2.mean]+resid.mod36[,resid.lag1_19.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod55<-residModFxn('model55')
models[['model55']]<-modFxn('model55')

residAcfFxn(resid.mod55,'model55')
#sum(acf): 1.289      mean(acf): 0.037 
#sum(pacf): 0.337     mean(pacf): 0.01 
#Fine.
rm(resid.mod55,model55)


save(models,file=paste0(f1,'models.RData'))

wavemodels[['Resp.Warm']]<-c('model37','model38','model54','model55')


####RESP COLD SEASON WITH PRECIP WAVES####

c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(10,11,0,1)
c2<-!is.na(cb1.ppt0in[,1])

#create new season spline with 3 df
ns.doy<-onebasis(sin(as.POSIXlt(claims_master_agg[c1 & c2,CLM_FROM_DT])$yday/354*3.14*2),'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model39<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt0in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod39<-residModFxn('model39',cond=c1 & c2)
models[['model39']]<-modFxn('model39')

residAcfFxn(resid.mod39,'model39')
#sum(acf): 7.833      mean(acf): 0.224 
#sum(pacf): 1.345     mean(pacf): 0.04
#Temporal dependence is present. There is significant autocorrelation at lags 1-7.

resid.mod39[,resid.lag1:=as.integer(stats::filter(resid,c(0,1),sides=1)*1e7)][is.na(resid.lag1),resid.lag1:=0L][,resid.lag1.mean:=mean(resid.lag1),by=CLM_FROM_DT]
resid.mod39[,resid.lag1_7:=as.integer(stats::filter(resid,c(0,1/1:7/2),sides=1)*1e7)][is.na(resid.lag1_7),resid.lag1_7:=0L][,resid.lag1_7.mean:=mean(resid.lag1_7),by=CLM_FROM_DT]

test<-claims_master_agg[,.(super.zip,CLM_FROM_DT)]
test<-resid.mod39[test,on=c('super.zip','CLM_FROM_DT')]
c3<-test[,!is.na(resid.lag1_7.mean)]
rm(test)

model40<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt0in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod39[,resid.lag1.mean]+resid.mod39[, resid.lag1_7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod40<-residModFxn('model40',cond=c1 & c2 & c3)
models[['model40']]<-modFxn('model40')

residAcfFxn(resid.mod40,'model40')
#sum(acf): 1.4      mean(acf): 0.04 
#sum(pacf): 0.353     mean(pacf): 0.01
#Fine.
rm(resid.mod40,model40)


model41<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt025in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod39[,resid.lag1.mean]+resid.mod39[, resid.lag1_7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod41<-residModFxn('model41',cond=c1 & c2 & c3)
models[['model41']]<-modFxn('model41')

residAcfFxn(resid.mod41,'model41')
#sum(acf): 1.352      mean(acf): 0.039 
#sum(pacf): 0.312     mean(pacf): 0.009 
#Fine.
rm(resid.mod41,model41)


model56<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt05in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod39[,resid.lag1.mean]+resid.mod39[, resid.lag1_7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod56<-residModFxn('model56',cond=c1 & c2 & c3)
models[['model56']]<-modFxn('model56')

residAcfFxn(resid.mod56,'model56')
#sum(acf): 1.342      mean(acf): 0.038 
#sum(pacf): 0.303     mean(pacf): 0.009 
#Fine.
rm(resid.mod56,model56)


model57<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1 & c2,]+cb1.ppt1in[c1 & c2,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod39[,resid.lag1.mean]+resid.mod39[, resid.lag1_7.mean], family=poisson(link='log'), data=claims_master_agg[c1 & c2,], model=F)
resid.mod57<-residModFxn('model57',cond=c1 & c2 & c3)
models[['model57']]<-modFxn('model57')

residAcfFxn(resid.mod57,'model57')
#sum(acf): 1.316      mean(acf): 0.038 
#sum(pacf): 0.282     mean(pacf): 0.008  
#Fine.
rm(resid.mod57,model57)


save(models,file=paste0(f1,'models.RData'))

wavemodels[['Resp.Cold']]<-c('model40','model41','model56','model57')

####RESP MELT SEASON WITH PRECIP WAVES####
c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% c(2,3,4) #March, April, May

#create new season spline with 3 df
ns.doy<-onebasis(as.POSIXlt(claims_master_agg[c1,CLM_FROM_DT])$yday,'ns',df=3)
dt1<-as.integer(round(ns.doy*100,0))
attributes(dt1)<-attributes(ns.doy)
ns.doy<-dt1
rm(dt1)

model42<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt0in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod42<-residModFxn('model42')
models[['model42']]<-modFxn('model42')

residAcfFxn(resid.mod42,'model42')
#sum(acf): 4.861      mean(acf): 0.139 
#sum(pacf): 1.235     mean(pacf): 0.036
#Positive correlation, and positive correlation in pacf through lag 10
rm(model42)

#re-run model with residual lags
resid.mod42[,resid.lag1:=as.integer(stats::filter(resid,c(0,1),sides=1)*1e7)][is.na(resid.lag1),resid.lag1:=0L][,resid.lag1.mean:=mean(resid.lag1),by=CLM_FROM_DT]
resid.mod42[,resid.lag1_8:=as.integer(stats::filter(resid,c(0,1/1:8/2),sides=1)*1e7)][is.na(resid.lag1_8),resid.lag1_8:=0L][,resid.lag1_8.mean:=mean(resid.lag1_8),by=CLM_FROM_DT]

model43<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt0in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod42[,resid.lag1.mean]+resid.mod42[,resid.lag1_8.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod43<-residModFxn('model43')
models[['model43']]<-modFxn('model43')

residAcfFxn(resid.mod43,'model43')
#sum(acf): 1.676      mean(acf): 0.048 
#sum(pacf): 0.618     mean(pacf): 0.018 
#Fine.
rm(resid.mod43,model43)

model44<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt025in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod42[,resid.lag1.mean]+resid.mod42[,resid.lag1_8.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod44<-residModFxn('model44')
models[['model44']]<-modFxn('model44')

residAcfFxn(resid.mod44,'model44')
#sum(acf): 1.742      mean(acf): 0.05 
#sum(pacf): 0.641     mean(pacf): 0.019
#fine.
rm(resid.mod44,model44)


model58<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt05in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod42[,resid.lag1.mean]+resid.mod42[,resid.lag1_8.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod58<-residModFxn('model58')
models[['model58']]<-modFxn('model58')

residAcfFxn(resid.mod58,'model58')
#sum(acf): 1.755      mean(acf): 0.05 
#sum(pacf): 0.651     mean(pacf): 0.019
#Fine.


model59<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+cb1.ppt1in[c1,]+ns.doy+year+as.factor(county)+as.factor(dow)+resid.mod42[,resid.lag1.mean]+resid.mod42[,resid.lag1_8.mean], family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod59<-residModFxn('model59')
models[['model59']]<-modFxn('model59')

residAcfFxn(resid.mod59,'model59')
#sum(acf): 1.747      mean(acf): 0.05 
#sum(pacf): 0.64     mean(pacf): 0.019
#Fine.
rm(resid.mod59,model59)


save(models,file=paste0(f1,'models.RData'))

wavemodels[['Resp.Melt']]<-c('model43','model44','model58','model59')

####SENSITIVITY ANALYSES WITH TMAX AND SEASON####

#What if we remove tmax?
model21<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.ppt1[c1,]+ns.doy[c1,]+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
resid.mod21<-data.table(resid=resid(model21),super.zip=claims_master_agg[c1,super.zip],CLM_FROM_DT=claims_master_agg[c1,CLM_FROM_DT])
models[['model21']]<-list(call=model21$call,coef=model21$coefficients,vcov=summary(model21)$cov.scaled,df.residual=model21$df.residual,df.null=model21$df.null,deviance=model21$deviance,null.deviance=model21$null.deviance,AIC=summary(model21)$aic)
rm(model21)

v1<-grep('ppt1',names(models[['model21']]$coef),value=T)
cp.mod21<-crosspred(basis=cb1.ppt1,coef=models[['model21']]$coef[v1]*100,vcov=models[['model21']]$vcov[v1,v1]*100^2,from=0,to=400,by=10,model.link='log',cen=110,cumul=T)

par(mfrow=c(1,2),mai=c(0.1,0.5,0.1,0.1)); plot(cp.mod20,cum=T,zlim=c(0.95,1.25)); plot(cp.mod21,cum=T,zlim=c(0.95,1.25))
plot(cp.mod20,'slices',lag=20, main='Mod 20, Lag 20, Cum', cumul=T)
plot(cp.mod21,'slices',lag=20, main='Mod 20, Lag 20, Cum', cumul=T)
#In the model without temperature, effect estimates are high (1.2 vs. 1.05), and the effects at very high levels of precipitation are not diminshed but but instead flatten out.

#What if we remove precip and look at tmax?
model22<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+ns.doy[c1,]+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
models[['model22']]<-list(call=model22$call,coef=model22$coefficients,vcov=summary(model22)$cov.scaled,df.residual=model22$df.residual,df.null=model22$df.null,deviance=model22$deviance,null.deviance=model22$null.deviance,AIC=summary(model22)$aic)
rm(model22)

v1<-grep('tmax',names(models[['model20']]$coef),value=T)
cp.mod20<-crosspred(basis=cb1.tmax,coef=models[['model20']]$coef[v1]*100,vcov=models[['model20']]$vcov[v1,v1]*100^2,from=-17,to=15,by=1,model.link='log',cen=15,cumul=T)
v1<-grep('tmax',names(models[['model22']]$coef),value=T)
cp.mod22<-crosspred(basis=cb1.tmax,coef=models[['model22']]$coef[v1]*100,vcov=models[['model22']]$vcov[v1,v1]*100^2,from=-17,to=15,by=1,model.link='log',cen=15,cumul=T)
plot(cp.mod20,cum=T); plot(cp.mod22,cum=T)
#The effects are almost identical--with temperatures less than 15 being protective, i.e., a linear decrease in RR with decreasing temperature. Precip doesn't confound temperature, but temperature confounds precip.

#What if we just look at season?
model23<-glm(resp_disease ~ offset(log(NEW_COUNT))+ns.doy[c1,]+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
models[['model23']]<-list(call=model23$call,coef=model23$coefficients,vcov=summary(model23)$cov.scaled,df.residual=model23$df.residual,df.null=model23$df.null,deviance=model23$deviance,null.deviance=model23$null.deviance,AIC=summary(model23)$aic)
rm(model23)

v1<-grep('ns.doy',names(models[['model20']]$coef),value=T)
v2<-names(models[['model20']]$coef)[grepl('ns.doy',names(models[['model20']]$coef)) & !is.na(models[['model20']]$coef)]
v3<-c(1:length(v1))[v1 %in% v2]
v4<-models[['model20']]$coef[v1]; v4[is.na(v4)]<-0
dat1<-as.data.frame(models[['model20']]$vcov[v2,v2],row.names=v2)
dat1[,v1[!(v1 %in% v2)]]<-0
dat1[v1[!(v1 %in% v2)],]<-0
dat1<-as.matrix(dat1[order(row.names(dat1)),order(names(dat1))])
cp.mod20<-crosspred(basis=ns.doy,coef=v4*100,vcov=dat1*100^2,from=304,to=365,by=1,model.link='log',cen=304,cumul=T)

v1<-grep('ns.doy',names(models[['model22']]$coef),value=T)
v2<-names(models[['model22']]$coef)[grepl('ns.doy',names(models[['model22']]$coef)) & !is.na(models[['model22']]$coef)]
v3<-c(1:length(v1))[v1 %in% v2]
v4<-models[['model22']]$coef[v1]; v4[is.na(v4)]<-0
dat1<-as.data.frame(models[['model22']]$vcov[v2,v2],row.names=v2)
dat1[,v1[!(v1 %in% v2)]]<-0
dat1[v1[!(v1 %in% v2)],]<-0
dat1<-as.matrix(dat1[order(row.names(dat1)),order(names(dat1))])
cp.mod22<-crosspred(basis=ns.doy,coef=v4*100,vcov=dat1*100^2,from=304,to=365,by=1,model.link='log',cen=304,cumul=T)

v1<-grep('ns.doy',names(models[['model23']]$coef),value=T)
v2<-names(models[['model23']]$coef)[grepl('ns.doy',names(models[['model23']]$coef)) & !is.na(models[['model23']]$coef)]
v3<-c(1:length(v1))[v1 %in% v2]
v4<-models[['model23']]$coef[v1]; v4[is.na(v4)]<-0
dat1<-as.data.frame(models[['model23']]$vcov[v2,v2],row.names=v2)
dat1[,v1[!(v1 %in% v2)]]<-0
dat1[v1[!(v1 %in% v2)],]<-0
dat1<-as.matrix(dat1[order(row.names(dat1)),order(names(dat1))])
cp.mod23<-crosspred(basis=ns.doy,coef=v4*100,vcov=dat1*100^2,from=304,to=365,by=1,model.link='log',cen=304,cumul=T)

plot(cp.mod20,cumul=T); plot(cp.mod23,cumul=T)
#The pattern is the same, but the seasonal effects are much stronger in mod 20.

plot(cp.mod22,cumul=T); plot(cp.mod23,cumul=T)
#The pattern is the same, but the seasonal effects are stronger in mod22 (which has tmax but not precip), and mod22 is similar to mod20. Therefore, tmax is negatively confounding season in the cold season.

#What if we remove season from the models?
model24<-glm(resp_disease ~ offset(log(NEW_COUNT))+cb1.tmax[c1,]+year+as.factor(county)+as.factor(dow), family=poisson(link='log'), data=claims_master_agg[c1,], model=F)
models[['model24']]<-list(call=model24$call,coef=model24$coefficients,vcov=summary(model24)$cov.scaled,df.residual=model24$df.residual,df.null=model24$df.null,deviance=model24$deviance,null.deviance=model24$null.deviance,AIC=summary(model24)$aic)
rm(model24)

v1<-grep('tmax',names(models[['model22']]$coef),value=T)
cp.mod22<-crosspred(basis=cb1.tmax,coef=models[['model22']]$coef[v1]*100,vcov=models[['model22']]$vcov[v1,v1]*100^2,from=-17,to=15,by=1,model.link='log',cen=15,cumul=T)
v1<-grep('tmax',names(models[['model24']]$coef),value=T)
cp.mod24<-crosspred(basis=cb1.tmax,coef=models[['model24']]$coef[v1]*100,vcov=models[['model24']]$vcov[v1,v1]*100^2,from=-17,to=15,by=1,model.link='log',cen=15,cumul=T)
plot(cp.mod22,cum=T); plot(cp.mod24,cum=T)
plot(cp.mod22,'slices',lag=20, main='Mod 20, Lag 22, Cum', cumul=T)
plot(cp.mod24,'slices',lag=20, main='Mod 20, Lag 24, Cum', cumul=T)
#Temperature still behaves oddly. Without season, the temperature effect still increases from -17 to 0, contrary to expected, but then decreases from 0 to 15, as would be expected (one would expect warmer temperatures in this range to be protective).

####QUANTILES####
#Quantiles of precipitation. Note that the precipitation variable here is transformed to an uninterpretable quantity. 
quantile(claims_master_agg[,PRISM_PPT],p=c(0,0.25,0.5,0.75,0.9,0.95,.99,.995,.999))
mean(v1) #64.8
#   0%   25%   50%   75%   90%   95%   99% 99.5% 99.9% 
#   0     0     0   110   240   289   358   381   428

#Quantiles of precipitation converted back to mm of precip
quantile(exp(claims_master_agg[,PRISM_PPT]/100)-1,p=c(0,0.25,0.5,0.75,0.9,0.95,.99,.995,.999)) 
#       0%       25%       50%       75%       90%       95%       99%     99.5%     99.9% 
# 0.000000  0.000000  0.000000  2.004166 10.023176 16.993310 34.873541 44.150439 71.240440

#Quantiles of precipitation weighted by 100ths of population at risk
#all seasons
v1<-exp(claims_master_agg[,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]/100)-1
quantile(v1,p=c(0,0.25,0.5,0.75,0.9,0.95,.99,.995,.999,1))
#      0%        25%        50%        75%        90%        95%        99%      99.5%      99.9%       100% 
#0.000000   0.000000   0.000000   2.004166  10.023176  16.993310  34.873541  44.150439  71.240440 224.879123 

#Warm season
v1<-exp(claims_master_agg[as.POSIXlt(CLM_FROM_DT)$mon %in% 5:8,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]/100)-1
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#        0%        75%        90%        95%        99%      99.5%      99.9%       100% 
# 0.000000   2.004166  10.941264  19.085537  41.097990  53.054889  84.626944 197.343425 

#Cold season
v1<-exp(claims_master_agg[as.POSIXlt(CLM_FROM_DT)$mon %in% c(10,11,0,1),as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]/100)-1
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#        0%        75%        90%        95%        99%      99.5%      99.9%       100% 
#  0.000000   2.004166   8.025013  14.029276  28.964100  34.873541  54.146871 106.770073  

#Spring season
v1<-exp(claims_master_agg[as.POSIXlt(CLM_FROM_DT)$mon %in% 2:4,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]/100)-1
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#        0%        75%        90%        95%        99%      99.5%      99.9%       100% 
#  0.000000   2.004166  10.023176  16.993310  32.115452  39.044847  58.145470 134.639414 

#Quantiles of ln(precipitation+1) weighted by 100ths of population at risk
#All days
v1<-claims_master_agg[,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#      0%      75%      90%      95%      99%    99.5%    99.9%     100% 
#    -Inf 0.695228 2.304900 2.832820 3.551728 3.787603 4.266061 5.415563

#Warm season
v1<-claims_master_agg[as.POSIXlt(CLM_FROM_DT)$mon %in% 5:8,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#  0%   75%   90%   95%   99% 99.5% 99.9%  100% 
#  0   110   248   300   374   399   445   529 

#Cold season
v1<-claims_master_agg[as.POSIXlt(CLM_FROM_DT)$mon %in% c(10,11,0,1),as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#   0%   75%   90%   95%   99% 99.5% 99.9%  100% 
#   0   110   220   271   340   358   401   468 

#Spring season
v1<-claims_master_agg[as.POSIXlt(CLM_FROM_DT)$mon %in% 2:4,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#   0%   75%   90%   95%   99% 99.5% 99.9%  100% 
#   0   110   240   289   350   369   408   491 

#All seasons
v1<-claims_master_agg[,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]
quantile(v1,p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#    0%   75%   90%   95%   99% 99.5% 99.9%  100% 
#     0   110   240   289   358   381   428   542

#Quantiles of ln(precipitation) weighted by 100ths of population at risk
v1<-claims_master_agg[,as.integer(rep(PRISM_PPT,each=NEW_COUNT/100))]
quantile(log(exp(v1/100)-1+1),p=c(0,0.75,0.9,0.95,.99,.995,.999,1))
#      0%      75%      90%      95%      99%    99.5%    99.9%     100% 
#    -Inf 0.695228 2.304900 2.832820 3.551728 3.787603 4.266061 5.415563

#Quantiles of temperature weighted by 100ths of population at risk
v1<-claims_master_agg[,as.integer(rep(PRISM_TMAX,each=NEW_COUNT/100))]
quantile(v1,p=c(0,0.25,0.5,0.75,0.9,0.95,.99,.995,.999,1))
#   0%   75%   90%   95%   99% 99.5% 99.9% 100%
# -17    26    30    31    34    35    37   41

#associate model identifiers with plot names and labels
names.models<-data.table(model=c('model3','model6','model10','model12','model14','model17','model20','model23'),disease=c('GI','Resp','GI','GI','GI','Resp','Resp','Resp'),season=c('All seasons','All seasons','Warm season','Cold season','Melt season','Warm season','Cold season','Melt season'))
dat1<-data.table(season=c('Cold season','Melt season','Warm season','All seasons'),rbind(c(0,220,358),c(0,240,369),c(0,248,399),c(0,240,381)))
setnames(dat1,c('V1','V2','V3'),c('min','p90','p95'))
names.models<-dat1[names.models,on='season']

####FIGURES####

#Pick precip model to plot and plot it
for (i in 1:nrow(names.models)) {
  v2<-names.models[i,model]
  v1<-grep('ppt1',names(models[[v2]]$coef),value=T) #all-season precip crossbasis names
  #3-D and Lag 20 Slices Plot Titles
  v3<-paste('RR of',names.models[i,disease],'ED Visit\nCumulative\n',names.models[i,season])
  v4<-paste('RR of',names.models[i,disease],'ED Visit\nLag 2 Cumulative\n',names.models[i,season])
  v5<-paste('RR of',names.models[i,disease],'ED Visit\nLag 12 Cumulative\n',names.models[i,season])
  v6<-paste('RR of',names.models[i,disease],'ED Visit\nLag 20 Cumulative\n',names.models[i,season])
  #create crossprediction for min to 99.5th percentiles
  v7<-names.models[i,c('min','p90','p95')]
  cp1<-crosspred(basis=cb1.ppt1,coef=models[[v2]]$coef[v1]*100,vcov=models[[v2]]$vcov[v1,v1]*100^2,from=v7[1],to=v7[3],by=10,model.link='log',cen=100,cumul=T)
  cp1$predvar<-exp(cp1$predvar/100)-1 #change precip back to original scale in cm
  #actual 3-D plots and corresponding RRs
  png(paste0(f1,gsub(' ','_',paste(names.models[i,disease],names.models[i,season],'Lags_1_5_10_20')),'.png'),width=7.5*96,height=440)
  #par(mfrow=c(2,2),mar=c(0.5,0.1,0.2,0.1))
  #plot(cp1, cumul=T,xlab='Precipitation (mm)',zlab='RR',cex.lab=0.9,cex.axis=0.9,theta=315,phi=45)
  #title(main='All lags',cex.main=1,line=-1,adj=0.3)
  par(mfrow=c(2,2),mar=c(4,4,0.1,0.1))
  plot(cp1,'slices',lag=1, cumul=T,xlab='Precipitation (mm)',ylab='RR',log='x',xaxt='n')
  title(main='Lag 1',line=-2,cex.main=1,adj=0.5)
  axis(side=1,at=c(seq(0,exp(v7[3]/100)-1),by=1),log='x'); axis(side=1,at=c(0.1,0.5))
  plot(cp1,'slices',lag=5, cumul=T,xlab='Precipitation (mm)',ylab='RR',log='x',xaxt='n')
  title(main='Lag 5',line=-2,cex.main=1,adj=0.5)
  axis(side=1,at=c(seq(0,exp(v7[3]/100)-1,by=1)),log='x'); axis(side=1,at=c(0.1,0.5))
  plot(cp1,'slices',lag=10, cumul=T,xlab='Precipitation (mm)',ylab='RR',log='x',xaxt='n')
  title(main='Lag 10', line=-2,cex.main=1,adj=0.5)
  axis(side=1,at=c(seq(0,exp(v7[3]/100)-1,by=1)),log='x'); axis(side=1,at=c(0.1,0.5))
  plot(cp1,'slices',lag=20, cumul=T,xlab='Precipitation (mm)',ylab='RR',log='x',xaxt='n')
  title(main='Lag 20', line=-2, cex.main=1,adj=0.5)
  axis(side=1,at=c(seq(0,exp(v7[3]/100)-1,by=1)),log='x'); axis(side=1,at=c(0.1,0.5))
  dev.off()
}

#Figures with all seasons, and lags 1, 5, and 20, in a single figure for GI and a single figure for resp
for (i in c(1:nrow(names.models))[3:8]) {
  if (names.models[i,season]=='All seasons') next
  v2<-names.models[i,model]
  v1<-grep('ppt1',names(models[[v2]]$coef),value=T) #all-season precip crossbasis names
  #create crossprediction for min to 99.5th percentiles
  if (names.models[i,season]=='Cold season') {v7<-c(0,100,358); v10<-c('D) Cold,','E) Cold,','F) Cold,')
  } else if (names.models[i,season]=='Melt season') {v7<-c(0,100,369); v10<-c('G) Spring melt,','H) Spring melt,','I) Spring melt,')
  } else if (names.models[i,season]=='Warm season') {v7<-c(0,100,399); v10<-c('A) Warm,','B) Warm,','C) Warm,')
  } else v7<-c(0,100,381)
  cp1<-crosspred(basis=cb1.ppt1,coef=models[[v2]]$coef[v1]*100,vcov=models[[v2]]$vcov[v1,v1]*100^2,from=0,to=v7[3],by=10,model.link='log',cen=100,cumul=T)
  cp1$predvar<-exp(cp1$predvar/100)-1 #change precip back to original scale in cm
  
  #plots and corresponding RRs
  #png(paste0('C:\\Users\\gronlund\\Dropbox (University of Michigan)\\DocumentsC\\K99R00\\Madeline\\',gsub(' ','_',paste(names.models[i,disease],names.models[i,season],'Lags_1_5_20')),'.png'),width=7.5*96,height=220)
  par(mfrow=c(1,3),mar=c(4,4,1,0.1))
  plot(cp1,'slices',lag=1, cumul=T,xlab='Precipitation (mm)',ylab='RR',log='x',xaxt='n')
  title(main=paste(v10[1],'Day 1'),line=-2,cex.main=1,adj=0.5)
  axis(side=1,at=c(seq(0,exp(v7[3]/100)-1),by=1),log='x'); axis(side=1,at=c(0.1,0.5))
  plot(cp1,'slices',lag=5, cumul=T,xlab='Precipitation (mm)',ylab='RR',log='x',xaxt='n')
  title(main=paste(v10[2],'Day 5'),line=-2,cex.main=1,adj=0.5)
  axis(side=1,at=c(seq(0,exp(v7[3]/100)-1,by=1)),log='x'); axis(side=1,at=c(0.1,0.5))
  plot(cp1,'slices',lag=20, cumul=T,xlab='Precipitation (mm)',ylab='RR',log='x',xaxt='n')
  title(main=paste(v10[3],'Day 20'), line=-2, cex.main=1,adj=0.5)
  axis(side=1,at=c(seq(0,exp(v7[3]/100)-1,by=1)),log='x'); axis(side=1,at=c(0.1,0.5))
  #dev.off()
}

#Figure with waves (0 precip 7 days, 0 precip 11 days, 0 precip 14 days, 0.5 inch 2 days) at lag 20 for each disease and season
wavmodels<-list(model=c('model33','model35','model45','model41','model44','model38'),disease=c('GI','GI','GI','Resp','Resp','Resp'),season=c('Cold season','Melt season','Warm season','Cold season','Melt season','Warm season'),months=list(c(10,11,0,1),c(2,3,4),c(5,6,7,8),c(10,11,0,1),c(2,3,4),c(5,6,7,8)))
wavcbs<-list(cbName=c('cb1.ppt0in','cb1.ppt05in','cb1.ppt1in'),center=c(3,0,0),val=list(seq(0,13,by=0.1),1,1))


for (i in 1:6) {
  if (i %in% c(1,4)) {
    #png(paste0(f1,gsub(' ','_',paste('wav',wavmodels[['disease']][[i]],'Lags_1_5_20')),'.png'),width=7.5*96,height=220*3)
    par(mfrow=c(3,3),mar=c(4,4,1,0.1))
  }
  for (j in 1) {
    c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% wavmodels[['months']][[i]]
    cp1<-cpSeasFxn(cbName=wavcbs['cbName'][[1]][j],modName=wavmodels[['model']][[i]],seasCond=c1,center=wavcbs['center'][[1]][j],val=wavcbs[['val']][[j]])
    if (j==1) cp1$predvar<-cp1$predvar+1
    plot(cp1,'slices',lag=1,cumul=T,xlab='Consecutive Days',ylab='RR')
    title(main=paste(wavmodels[['season']][[i]],'Day 1'),line=-2,cex.main=1,adj=0.5)
    plot(cp1,'slices',lag=5,cumul=T,xlab='Consecutive Days',ylab='RR')
    title(main=paste(wavmodels[['season']][[i]],'Day 5'),line=-2,cex.main=1,adj=0.5)
    plot(cp1,'slices',lag=20,cumul=T,xlab='Consecutive Days',ylab='RR')
    title(main=paste(wavmodels[['season']][[i]],'Day 20'),line=-2,cex.main=1,adj=0.5)
  }
}


#Plot temperature effects
#Pick precip model to plot and plot it
for (i in 1:nrow(names.models)) {
  v2<-names.models[i,model]
  v1<-grep('tmax',names(models[[v2]]$coef),value=T) #all-season precip crossbasis names
  #3-D and Lag 20 Slices Plot Titles
  v3<-paste('RR of',names.models[i,disease],'ED Visit\nCumulative\n',names.models[i,season])
  v4<-paste('RR of',names.models[i,disease],'ED Visit\nLag 1 Cumulative\n',names.models[i,season])
  v5<-paste('RR of',names.models[i,disease],'ED Visit\nLag 4 Cumulative\n',names.models[i,season])
  v6<-paste('RR of',names.models[i,disease],'ED Visit\nLag 20 Cumulative\n',names.models[i,season])
  #create crossprediction
  if (names.models[i,season]=='Cold season') {v7<-c(-9,5,20)
  } else if (names.models[i,season]=='Melt season') {v7<-c(-2,22,31)
  } else if (names.models[i,season]=='Warm season') {v7<-c(16,26,36)
  } else v7<-c(-9,26,36)
  cp1<-crosspred(basis=cb1.tmax,coef=models[[v2]]$coef[v1]*100,vcov=models[[v2]]$vcov[v1,v1]*100^2,from=v7[1],to=v7[3],by=1,model.link='log',cen=v7[2],cumul=T)
  #actual 3-D plots
  png(paste0(f1,gsub(' ','_',paste(names.models[i,disease],names.models[i,season],'Temp_Lags_all_1_5_20')),'.png'),width=7.5*96,height=440)
  par(mfrow=c(2,2),mar=c(0.5,0.1,0.2,0.1))
  plot(cp1, cumul=T,xlab='Temperature',zlab='RR',cex.lab=0.9,cex.axis=0.9,theta=315,phi=45)
  title(main='All lags',cex.main=1,line=-1,adj=0.3)
  #plot(cp1, cumul=F,xlab='Temperature',zlab='RR',cex.lab=0.9,cex.axis=0.9,theta=315,phi=45)
  par(mar=c(4,4,0.1,0.1))
  plot(cp1,'slices',lag=1, cumul=T,xlab='Temperature',ylab='RR')
  title(main='Lag 1',line=-2,cex.main=1,adj=0.5)
  plot(cp1,'slices',lag=5, cumul=T,xlab='Temperature',ylab='RR')
  title(main='Lag 5', line=-2,cex.main=1,adj=0.5)
  plot(cp1,'slices',lag=20, cumul=T,xlab='Temperature',ylab='RR')
  title(main='Lag 20', line=-2, cex.main=1,adj=0.5)
  dev.off()
}


####BURDEN OF DISEASE CALCULATION####
#AF = 1 - sum(p.ci * (RR.i-1)/RR.i) where i indexes the exposure stratum and p.ci is the proportion of *cases* in the exposure stratum. here, precip > 1 in vs. precip <= 1 in are the exposure strata. Therefore, we first must average the ln(RR)s on the days where precip > 1 in and also average the lnRRs for days where precip <= 1 in.
#BOD = AF * cases
#BOD.annual.rate = BOD / annual-persons-at-risk

#1. extract the RR for each disease and average by stratum
#2. calculate the corresponding number of cases in each stratum, summing across days and zip codes
#3. calculate the corresponding AF for each stratum and calculate a case-weighted mean of the AF
#4. multiply AF and total # cases for present day
#5. For the precip > 1 in exposure stratum, increase the weight for this stratum by 2 additional days per year and decrease the weight of the precip <= 1 stratum by 2 fewer days per year and recalculate the AF.
#6. Multiply the new AF and total # cases to get BOD for future

#dat4<-data.table()
ls1<-list(All=0:11,Warm=5:8,Cold=c(10:11,0:1),Melt=2:4) #seasons
v1<-setNames(paste0('model',c(3,6,10,12,14,17,20,23)),c('GI.All','Resp.All','GI.Warm','GI.Cold','GI.Melt','Resp.Warm','Resp.Cold','Resp.Melt'))
mat5<-matrix(NA,8,9,dimnames=list(names(v1),c('AF.lt2mm','AF.gt2mmlt1in','AF.gt1in','AF.lt2mm.lower','AF.gt2mmlt1in.lower','AF.gt1in.lower','AF.lt2mm.upper','AF.gt2mmlt1in.upper','AF.gt1in.upper')))
mat6<-matrix(NA,8,6,dimnames=list(names(v1),paste0(c('lt2mm','gt2mmlt1in','gt1in'),rep(c('.annual.cases','.daily.rate'),each=3))))
set.seed(20100915)
for (m1 in names(v1)[2:8]) {
  cat('\n',m1,'\n')
  c2<-claims_master_agg[,as.POSIXlt(CLM_FROM_DT)$mon %in% ls1[[gsub('GI\\.|Resp\\.','',m1)]]] #restrict to the season
  c1<-claims_master_agg[c2,PRISM_PPT>326] #instances where precip > 1 in
  c3<-claims_master_agg[c2,PRISM_PPT<110] #instances where precip < 2 mm, or < 0.079 in
  m2<-v1[m1] #actual model name
  #first, calculate the log(RR) (allfit) for each superzip-day
  v2<-grep('ppt',names(models[[m2]]$coef),value=T)
  #mat3<-t(replicate(500,expr=rnorm(n=length(v2),mean=models[[m2]]$coef[v2],sd=sqrt(diag(models[[m2]]$vcov[v2,v2]))))) #randomly generate coefficients
  mat3<-mvrnorm(n=500,models[[m2]]$coef[v2],models[[m2]]$vcov[v2,v2])
  cb0<-crossbasis(rep(110,21),lag=20,argvar=attr(cb1.ppt1,'argvar'),arglag=attr(cb1.ppt1,'arglag'))[21,]*100 #variable values at ppt = 110
  #cp2<-rowSums((cb1.ppt1[1:10,]-matrix(cb0,10,16,byrow=T))*matrix(models[[m2]]$coef[v2],10,16,byrow=T)) #lnRR
  #mat2<-as.matrix(claims_master_agg[c2,.SD,.SDcols=precipCols])
  #cp1<-crosspred(basis=cb1.ppt1,coef=models[[m2]]$coef[v2]*100,vcov=models[[m2]]$vcov[v2,v2]*100^2,model.link='log',at=mat2[1:10,],cumul=T,cen=110)
  #all(abs((cp2-cp1$allfit)/cp1$allfit)<0.01) #all within 1% of each other, so verified that these two approaches produce about the same results
  mat4<-matrix(NA,500,3,dimnames=list(NULL,c('AF.lt2mm','AF.gt2mmlt1in','AF.gt1in')))
  #pull out cases and denominators for each super.zip-day and assign ppt categories
  #precip categories: -1 = < 2 mm, 1 = > 1 in, and 0 = all else, and 3 = all levels of precip
  dat1<-cbind(data.table(pptCat=as.integer(c1-c3)),claims_master_agg[c2,.SD,.SDcols=c(paste0(d1,'_disease'),'NEW_COUNT','super.zip')])
  setnames(dat1,paste0(d1,'_disease'),'cases')
  c4<-dat1[,pptCat==-1]; v4<-dat1[c4,sum(cases)]
  c5<-dat1[,pptCat==0]; v5<-dat1[c5,sum(cases)]
  c6<-dat1[,pptCat==1]; v6<-dat1[c6,sum(cases)]
  #calculate a case-weighted AF for each of the precip categories per Steenland and Armstrong 2006 (AF = sum(p.cases.i*(RR.i-1)/RR.i) = 1-sum(p.cases.i/RR.i)) where i indexes exposure level and p.cases.i is the proportion of cases in exposure level i. Monte Carlo over 500 iterations and take 2.5th and 97.5th percentiles to get CIs. Takes about 8 minutes for each model.
  pb<-txtProgressBar(min=1,max=100,style=3)
  for (i in 1:500) {
    setTxtProgressBar(pb,i/5)
    cp2<-rowSums((cb1.ppt1[c2,]-matrix(cb0,sum(c2),length(cb0),byrow=T))*matrix(mat3[i,],sum(c2),ncol(mat3),byrow=T)) #lnRR for each super.zip-day for replicate i
    mat4[i,'AF.lt2mm']<-1-dat1[c4,sum(cases/exp(cp2[c4]))]/v4
    mat4[i,'AF.gt2mmlt1in']<-1-dat1[c5,sum(cases/exp(cp2[c5]))]/v5
    mat4[i,'AF.gt1in']<-1-dat1[c6,sum(cases/exp(cp2[c6]))]/v6
  }
  v7<-setNames(apply(mat4,2,median),colnames(mat4))
  v8<-setNames(apply(mat4,2,quantile,0.025),paste0(colnames(mat4),'.lower'))
  v9<-setNames(apply(mat4,2,quantile,0.975),paste0(colnames(mat4),'.upper'))
  mat5[m1,names(c(v7,v8,v9))]<-c(v7,v8,v9)
  mat6[m1,paste0(c('lt2mm','gt2mmlt1in','gt1in'),'.annual.cases')]<-c(v4,v5,v6)/8
  #(note that each exposure stratum isn't over 365 days (except all) but is instead over the number of days in that stratum)
  #find the mean daily rate in each category
  mat6[m1,paste0(c('lt2mm','gt2mmlt1in','gt1in'),'.daily.rate')]<-
    c(dat1[pptCat==-1,mean(cases/NEW_COUNT)],dat1[pptCat==0,mean(cases/NEW_COUNT)],dat1[pptCat==1,mean(cases/NEW_COUNT)])
  write.csv(cbind(mat5,mat6),file='AFs_rates_20230320.csv')
}

#mat7<-read.csv('AFs_rates_20230320.csv')
#row.names(mat7)<-mat7[,'X']
#mat7[,'X']<-NULL
mat7<-cbind(mat5,mat6)
dat1<-setDT(as.data.table(mat7))[,mod:=row.names(mat7)]
dat2<-melt(dat1,id.vars='mod',measure.vars=list(c('AF.lt2mm','AF.gt2mmlt1in','AF.gt1in'),c('AF.lt2mm.lower','AF.gt2mmlt1in.lower','AF.gt1in.lower'),c('AF.lt2mm.upper','AF.gt2mmlt1in.upper','AF.gt1in.upper'),c('lt2mm.daily.rate','gt2mmlt1in.daily.rate','gt1in.daily.rate')),variable.name='pptCat',value.name=c('AF.med','AF.lower','AF.upper','daily.rate'))
dat2[,pptCat:=c('lt2mm','gt2mmlt1in','gt1in')[pptCat]]
dat2[,seas:=gsub('[[:alpha:]]{2,4}\\.','',mod)]
dat2<-data.table(seas=c('All','Warm','Cold','Melt'),nDaySeas=c(365.25,30+31+31+30,30+31+31+28.25,31+30+31))[dat2,on='seas']
dat2[,dz:=tolower(gsub('\\.All|\\.Warm|\\.Cold|\\.Melt','',mod))]

#Calculate the annual rates, based on the case-weighted number of days in each category
#First, Calculate number of days > 1 inch per year in 2006-2013, weighted by cases, and 2mm-1in and < 2mm
#length(seq.Date(from=as.Date('2006-01-01'),to=as.Date('2013-12-31'),by=1)) 2922 days in 2006-2013, or 365.25 days per year
claims_master_agg[,sum((PRISM_PPT<=log(2.54*10+1)*100 & PRISM_PPT>110)*(gi_disease+resp_disease)/sum(gi_disease+resp_disease))]*365.25 #78.33573
claims_master_agg[,sum((PRISM_PPT<=110)*(gi_disease+resp_disease)/sum(gi_disease+resp_disease))]*365.25 #279.2788
claims_master_agg[,sum((PRISM_PPT>log(2.54*10+1)*100)*(gi_disease+resp_disease)/sum(gi_disease+resp_disease))]*365.25 #7.64
#Therefore, if there are two more days per year in the future period, there will be 9.64 days per year on average of ppt > 1 in.

#estimate numerator as the sum of cases within pptCat, seas, and year
claims_master_agg[,sum(gi_disease,resp_disease)] #5091300 cases total
claims_master_agg[,pptCat:=cut(PRISM_PPT,c(0,110,326,542),include.lowest = T,labels=c('lt2mm','gt2mmlt1in','gt1in'))]
claims_master_agg[,seas:=c('Cold','Cold',rep('Melt',3),rep('Warm',4),'',rep('Cold',2))[as.POSIXlt(CLM_FROM_DT)$mon+1]]
claims_master_agg[,year:=as.POSIXlt(CLM_FROM_DT)$year+1900]
dat3<-claims_master_agg[,lapply(.SD,sum),.SDcols=c('gi_disease','resp_disease'),by=c('pptCat','seas','year')][seas!='',]
dat3<-melt(dat3,id.vars=c('pptCat','seas','year'),measure.vars=c('gi_disease','resp_disease'),variable.name='dz',value.name='cases')
dat3[,dz:=gsub('_disease','',dz)]

#estimate denominator for rate, first summing across super.zips and then taking the mean across all days within categories of pptCat,year, and seas
dat3<-claims_master_agg[,sum(NEW_COUNT),by=c('CLM_FROM_DT','pptCat','seas','year')][,mean(V1),by=c('pptCat','seas','year')][dat3,on=c('pptCat','seas','year')]
setnames(dat3,'V1','denom')

#estimate numerators and denominators for all seasons and bind to season-specific estimates
dat4<-claims_master_agg[,lapply(.SD,sum),.SDcols=c('gi_disease','resp_disease'),by=c('pptCat','year')][,seas:='All']
dat4<-melt(dat4,id.vars=c('pptCat','seas','year'),measure.vars=c('gi_disease','resp_disease'),variable.name='dz',value.name='cases')
dat4[,dz:=gsub('_disease','',dz)]
dat4<-claims_master_agg[,sum(NEW_COUNT),by=c('CLM_FROM_DT','pptCat','year')][,mean(V1),by=c('pptCat','year')][dat4,on=c('pptCat','year')]
setnames(dat4,'V1','denom')
dat4<-rbind(dat3,dat4)

#estimate sums of cases, means of denoms, and means of rates across all years
dat4[,rate.ann:=cases/denom]
dat5<-dat4[,lapply(.SD,mean),.SDcols=c('cases','denom','rate.ann'),by=c('pptCat','seas','dz')]
setnames(dat5,c('cases','denom','rate.ann'),c('cases.mean.ann','denom.mean.ann','rate.mean.ann'))

#merge into AFs
dat6<-dat5[dat2,on=c('pptCat','seas','dz')]

#calculate annual attributable rate (AAR) per 10,000
dat6[,paste0('AARe4.',c('med','lower','upper')):=lapply(.SD,'*',rate.mean.ann*10000),.SDcols=paste0('AF.',c('med','lower','upper'))]
dat6[,period:='present']

#calculate the increased proportion of effects for a change of 7.64 to 9.64 days per year
1+2/7.64
#increase of 1.26

#calculate future attributable rates
dat7<-dat6[seas=='All' & pptCat=='gt1in',]
dat7[,paste0('AARe4.',c('med','lower','upper')):=.SD*1.26,.SDcols=paste0('AARe4.',c('med','lower','upper'))][,period:='future']
dat7[,paste0('AF.',c('med','lower','upper')):=NULL]
dat8<-rbind(dat6,dat7,fill=T)

fwrite(dat8,file=paste0(f1,'bod_results_20230319.csv'))

####TABLES AND CALCULATIONS####

v9<-c('lag1','lag5','lag10','lag20')
mat1<-mat2<-mat3<-matrix(NA,8,4,dimnames=list(1:8,v9))
for (i in 1:nrow(names.models)) {
  v2<-names.models[i,model]
  v1<-grep('ppt1',names(models[[v2]]$coef),value=T) #all-season precip crossbasis names
  v7<-names.models[i,c('min','p90','p95')]
  cp2<-crosspred(basis=cb1.ppt1,coef=models[[v2]]$coef[v1]*100,vcov=models[[v2]]$vcov[v1,v1]*100^2,at=v7,model.link='log',cen=100,cumul=T) #crosspred at table values of 0, 90th, and 99.5th
  #tables of results
  ls1<-lapply(1:3,FUN=function(y) {mapply(cp2[c('cumRRfit','cumRRlow','cumRRhigh')],FUN=function(x) {sprintf('%.2f',round(x[y,v9],2))})})
  mat1[i,]<-paste0(ls1[[1]][,1],' (',ls1[[1]][,2],', ',ls1[[1]][,3],')')
  mat2[i,]<-paste0(ls1[[2]][,1],' (',ls1[[2]][,2],', ',ls1[[2]][,3],')')
  mat3[i,]<-paste0(ls1[[3]][,1],' (',ls1[[3]][,2],', ',ls1[[3]][,3],')')
}
vals.models<-data.table(names.models[,c('model','disease','season')],mat3,mat2,mat1) #mat1 is min, mat2 is 90th, and mat3 is 99.5th percentile
names(vals.models)[4:7]<-paste0(names(vals.models)[4:7],'.995th')
names(vals.models)[8:11]<-paste0(names(vals.models)[8:11],'.90th')
names(vals.models)[12:15]<-paste0(names(vals.models)[12:15],'.min')
fwrite(vals.models[order(disease,season),],paste0(f1,'model_values_995_90_0.csv'))

##Wave table
#wavemodels was constructed in the above wave code and is copied down here. For each unique combination of disease and season, there are 4 models corresponding to 0 in, 0.25 in, 0.5 in and 1 in. For the 0 in model, we will consider the effects of 1,4,8, and 12 consecutive days, for the 0.25 in model the effects of 1 and 3 consecutive days, for the 0.5 in model the effects of 1 and 2 consecutive days, and for the 1 in model the effects of 1 day.

wavemodels<-list(GI.Warm=c('model45','model46','model47','model48'))
wavemodels[['GI.Cold']]<-c('model33','model49','model50','model51')
wavemodels[['GI.Melt']]<-c('model34','model35','model52','model53')
wavemodels[['Resp.Warm']]<-c('model37','model38','model54','model55')
wavemodels[['Resp.Cold']]<-c('model40','model41','model56','model57')
wavemodels[['Resp.Melt']]<-c('model43','model44','model58','model59')


wavecbs<-data.table(cbName=c(rep('cb1.ppt0in',4),rep('cb1.ppt025in',2),rep('cb1.ppt05in',2),'cb1.ppt1in'),Amount=c(rep('0 in',4),rep('0.25 in',2),rep('0.5 in',2),'1 in'),Days=c(2,4,8,12,1,3,1,2,1),center=c(2,2,2,2,1,1,1,1,1),wavemod=c(1,1,1,1,2,2,3,3,4),refrow=c(1,1,1,1,5,5,7,7,9))

mat1<-mat2<-mat3<-mat4<-mat6<-mat7<-mat8<-matrix(NA,9,6)
for (i in 1:9) {
  k<-wavecbs[,'wavemod'][[1]][i] #the index in the dz-season-specific element of the wavemodels list corresponding to a crossbasis
  for (j in 1:6) {
    v1<-gsub('GI\\.|Resp\\.','',names(wavemodels))[j] #season
    v2<-list(Warm=c(5,6,7,8),Cold=c(10,11,0,1),Melt=c(2,3,4))[[v1]] #months in season
    c1<-as.POSIXlt(claims_master_agg$CLM_FROM_DT)$mon %in% v2 #select rows for that season
    cp1<-cpSeasFxn(cbName=wavecbs[,'cbName'][[1]][i],modName=wavemodels[[j]][[k]],seasCond=c1,center=0,val=wavecbs[,'Days'][[1]][i])
    mat1[i,j]<-round(cp1$allRRfit,2)
    mat2[i,j]<-round(cp1$allRRlow,2)
    mat3[i,j]<-round(cp1$allRRhigh,2)
    mat4[i,j]<-c(' ','*','**')[sum(abs(cp1$allfit/cp1$allse)>=1.96,abs(cp1$allfit/cp1$allse)>=2.57,1)]
    cp2<-cpSeasFxn(cbName=wavecbs[,'cbName'][[1]][i],modName=wavemodels[[j]][[k]],seasCond=c1,center=wavecbs[,'center'][[1]][i],val=wavecbs[,'Days'][[1]][i])
    mat6[i,j]<-cp2$allfit
    mat7[i,j]<-cp2$allse
    mat8[i,j]<-c('','\u2020','\u2021')[sum(abs(cp2$allfit/cp2$allse)>=1.96,abs(cp2$allfit/cp2$allse)>=2.57,1,na.rm=T)]
  }
}

exp(mat6)
mat5<-matrix(paste0(sprintf('%.2f',mat1),' (',sprintf('%.2f',mat2),', ',sprintf('%.2f',mat3),')',mat4,mat8),9,6)
dimnames(mat5)<-list(paste(wavecbs$Amount,'for',wavecbs$Days,'days'),paste(gsub('Resp','Respiratory',gsub('\\.',' ',names(wavemodels))),'Season'))
write.csv(mat5,paste0(f1,'wave_table.csv'))

#In the final manuscript, the 1-day durations are removed due to confusion over the differences between these and the independent effects, which have different reference categories (slight precip) and no constraints on length of the event.

##BOD Table
dat9<-fread(paste0(f1,'bod_results_20230319.csv'))
dat9[,AF.pretty:=paste0(sprintf('%.1f',AF.med*100),' (',sprintf('%.1f',AF.lower*100),', ',sprintf('%.1f',AF.upper*100),')')]
dat9[,AARe4.pretty:=paste0(sprintf('%.2g',AARe4.med),' (',sprintf('%.2g',AARe4.lower),', ',sprintf('%.2g',AARe4.upper),')')]
dat9[grepl('e\\+02',AARe4.pretty),AARe4.pretty:=gsub('e\\+02','0',gsub('\\.','',AARe4.pretty))]
dat9[AF.pretty=='NA (NA, NA)',AF.pretty:=gsub('\\(NA, NA\\)','',AF.pretty)]
dat9<-dcast(dat9,formula='period + seas + pptCat ~ dz',value.var=c('AF.pretty','AARe4.pretty'))
dat9[,period:=factor(period,levels=c('present','future'))][,pptCat:=factor(pptCat,levels=c('lt2mm','gt2mmlt1in','gt1in'))]
dat9<-dat9[order(period,pptCat,seas),.(period,seas,pptCat,AF.pretty_gi,AARe4.pretty_gi,AF.pretty_resp,AARe4.pretty_resp)]
fwrite(dat9,file=paste0(f1,'AF_BOD_table_20230320.csv'))

dat9[mod=='GI.All' & period=='present' & pptCat==2,'denom'] #1,381,023 persons

#from United States Census Bureau Quick Facts (https://www.census.gov/quickfacts/fact/table/PA,OH,MI,US/PST045222) accessed 3/10/2023:
12972008*.19 + 11756058*0.178 + 10034113*0.181 #6,373,434 persons (2,464,682 + 2,092,578 + 1,816,174)

#healthcare cost inflation https://www.usinflationcalculator.com/inflation/health-care-inflation-in-the-united-states/ based on Bureau of Labor Statistics Consumer Price Inflation Index Medical Care Index, inflation from 2015 dollars
2000*1.041*1.018*1.02*1.046*1.018 #$2302

(4.1+1.8+2+4.6+1.8)/5 #average of 2.86% per year

2000*(1+0.0286*5) #$2286, so similar using the 5-year average instead of each year

#2006-2013
6.37e+06*(37)/10000 #23,569 visits for respiratory annually for extreme precip
6.37e+06*(15)/10000 #9,555 visits for GI annually for extreme precip
6.37e+06*(37+15)/10000 #33,124 visits for resp + GI annually for extreme precip

6.37e+06*(37)/10000*2300 #$54,208,700 for resp annually for extreme precip
6.37e+06*(15)/10000*2300 #$21,976,500 for GI annually for extreme precip
6.37e+06*(37+15)/10000*2300 #$76,185,200 for GI+resp annually for extreme precip
2.47e+06*(37+15)/10000*2300 #$29,541,200 in PA
2.09e+06*(37+15)/10000*2300 #24,996,400 in OH
1.82e+06*(37+15)/10000*2300 #21,767,200 in MI
1.82e+06*(37+15)/10000*2000 #18,928,000 in MI in 2015 dollars

#2077-2099
6.37e+06*(24+59)/10000*2300 #$121,603,300 for GI+resp annually for extreme precip

####END####

