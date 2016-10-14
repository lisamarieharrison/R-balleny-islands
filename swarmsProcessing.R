exportDir='e:\\2015Tangaroa\\balleny\\20016dKrillExportAGG\\'
for(i in 1:length(EVFiles))
{
  message('Iteration ',i,' ; working with file ',raws[i])
  EVFile=EVOpenFile(EVAppObj,EVFiles[i])$EVFile
  Sys.sleep(0.1)
  EVSchoolsDetect(EVFile = EVFile,
                                acoVarName='120 7x7 convolution',
                                outputRegionClassName = 'agg',
                                deleteExistingRegions = TRUE,
                                distanceMode = "GPS distance",
                                maximumHorizontalLink = 30, #m
                                maximumVerticalLink = 10,#m
                                minimumCandidateHeight = 10, #m
                                minimumCandidateLength = 30, #m
                                minimumSchoolHeight = 10, #m
                                minimumSchoolLength = 30, #m
                                dataThreshold = -70)
  EVIntegrationByRegionsExport(EVFile=EVFile, acoVarName='120 7x7 convolution', regionClassName='agg', 
                               paste(exportDir,'120kHzAGG',exportNames[i],sep=''),
                                            dataThreshold = -80)
  Sys.sleep(0.1)
  EVIntegrationByRegionsExport(EVFile=EVFile, acoVarName='038 7x7 convolution', regionClassName='agg', 
                               paste(exportDir,'038kHzAGG',exportNames[i],sep=''),
                               dataThreshold = -80)
  EVSaveFile(EVFile=EVFile)
  Sys.sleep(0.25)
  EVCloseFile(EVFile)
  Sys.sleep(0.25)

  }

###Process schools
schoolsFiles=list.files(exportDir)
fnDF=cbind.data.frame(freq=sapply(strsplit(schoolsFiles,'kHzAGGtan1502-'),function(x) x[1]),
                 fn=sapply(strsplit(schoolsFiles,'kHzAGGtan1502-'),function(x) x[2]),fullname=schoolsFiles)
ufns=as.character(unique(fnDF$fn))
FirstRunFLAG=TRUE

for(i in ufns){
  message('Working on ',i)
  fnIND=which(fnDF$fn==i)
 if(length(fnIND)!=2){ warning('Incorrect number of unique files (2) for filename = ',i)
   next}
  freqs=as.numeric(as.character(fnDF$freq[fnIND]))
if(!all(freqs %in% c(38,120))) {warning('Incorrect frequencies detected (2) for filename = ', i)
  next}
  #hard wired:
  aggsL=vector(mode='list',length=2)
  aggsL[[1]]=read.csv(paste(exportDir,fnDF$fullname[fnIND[1]],sep=''))
  aggsL[[2]]=read.csv(paste(exportDir,fnDF$fullname[fnIND[2]],sep=''))
 if(nrow(aggsL[[1]])!=nrow(aggsL[[2]])){
   warning('Number of aggregations mismatch for files: ',fnDF$fullname[fnIND[1]],' ; ',fnDF$fullname[fnIND[2]])
   next
 }
  aggs=aggsL[[which(freqs==120)]]
  aggs=aggs[,c(1:34,52:71)]
  aggs$Sv_mean_038=aggsL[[which(freqs==38)]]$Sv_mean
  aggs$dBdiff=aggs$Sv_mean-aggs$Sv_mean_038
  aggs$filename=i
  write.table(aggs,paste(exportDir,'AllAggregations.csv',sep=''),row.names = F,sep=',',col.names=FirstRunFLAG,
            append=!FirstRunFLAG)
  FirstRunFLAG=FALSE
}  
###Work on processed aggregations:
aggs=read.csv(paste(exportDir,'AllAggregations.csv',sep=''))

#Check diner (2001) aquatic living resources beam corection:
Nb=function(Sv,L,P,Th,theta){
  #Sv=observed Sv
  #L=observed length
  #Th=minimum threshold
  #theta=beam width
  #P=mean depth of the echotrace.
  
  dST=Sv-Th
  
  B=0.44*theta * dST^0.45#detection angle
  
  B=pi*B/180
  Nb=L/(2*P*tan(B/2))
  return(Nb)}
aggs$Nbi=Nb(Sv=aggs$Sv_mean,L=aggs$Uncorrected_length,P=aggs$Depth_mean,Th=-70,theta=7)
#when the variable Nbi=NaN it means the acoustic energy in the aggregation is less than the minimum data threshold.
aggs$include=TRUE

aggs$include[is.nan(aggs$Nbi)]=FALSE
table(aggs$include)
aggs$include[aggs$Nbi<1.5]=FALSE #THis is a condition from the Diner 2001 paper
nrow(aggs)
aggs=aggs[aggs$include=TRUE,]
nrow(aggs)

#isolate krill based on 2 frequency dB difference:
mindB=0.37 ; maxdB=12
aggs=aggs[which(aggs$dBdiff>=mindB & aggs$dBdiff<=maxdB),]
nrow(aggs)
write.csv(aggs,paste(exportDir,'AllKrillSwarms.csv',sep=''),row.names=FALSE)

#scale 120 kHz Sv to volumetric density
TS=read.csv('file:///E:/2015Tangaroa/balleny/TS/TS.csv')

len=as.numeric(substr(names(TS)[-1],2,nchar(names(TS)[-1])))
#linear TS:
sigma_bs=10**(as.vector(unlist(TS[2,-1]))/10)
length(len)
length(sigma_bs)
sigmaBSDF=data.frame(len=len*1000,sigmaBS=sigma_bs)
#get trawl data:
LF=read.csv('file:///E:/2015Tangaroa/balleny/TS/tan1502_biologicals_to_43.csv')
summary(LF)

#krill only
nrow(LF)
LF=LF[LF$species=='EUS',]
nrow(LF)
KLF=LF$lgth*10
tKLF=table(KLF)
sigmaBSDF$density=0
sigmaBSDF$density[which(sigmaBSDF$len %in% as.numeric(names(tKLF)))]=as.vector(tKLF)/sum(tKLF)

TSbar=10*log10(weighted.mean(x=sigmaBSDF$sigmaBS,w=sigmaBSDF$density))
wetmassFunc=function(l)
  2.236e-6*l^3.314 #l in mm, wetmass g from jarvis et al 2010 dsrII p930
TSkg=10*log10(1000/weighted.mean(x=wetmassFunc(sigmaBSDF$len),w=sigmaBSDF$density))+TSbar

aggs$volDengPerm3=1000*(10**((aggs$Sv_mean-TSkg)/10))
write.csv(aggs,paste(exportDir,'AllKrillSwarms.csv',sep=''),row.names=FALSE)
