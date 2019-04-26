##1/4 degree daily precip recon from jan 1 1951 to dec 31 2015
#EOF data http://www.cpc.ncep.noaa.gov/products/janowiak/cmorph_description.html units mm/day
#download at ftp://ftp.cpc.ncep.noaa.gov/precip/CMORPH_V1.0/CRT/0.25deg-DLY_00Z/
#Obs data http://cdiac.ornl.gov/ftp/ushcn_daily/ units hundreths of inch/day



#########STATION DATA##############
library(ncdf4)
#dname='pobs'
obs=nc_open('ushcn_prcp.nc')
lon= ncvar_get(obs, 'LON')
#length 1218
lat=ncvar_get(obs,'LAT')
#length 1218
time=ncvar_get(obs,'TIME')
#length 44493371
prcp=ncvar_get(obs,'PRCP')
#length 44493371
ll=cbind(lat,lon)
#dim 1218 2
statind=ncvar_get(obs,'STATION_INDEX')
#length 44493371
library(lubridate) 
time=date_decimal(time)
sec=dseconds(1)
time=time + sec
time=format(time,"%Y%m%d")
tp=cbind(statind,time,prcp)
write.csv(tp,'tp.csv',row.names=F)

#must write out tp as csv, read back in, then use the read.zoo
#write.csv(tp,'tp.csv',row.names=F)
library(data.table)
tp=fread('tp.csv',header=T,sep=',',stringsAsFactors=F)

library(zoo)
prcpmat=read.zoo(tp, split=1, index=2, format= '%Y%m%d')
prcpmat=t(prcpmat)
#it is now station x time over the entire time, 1849-05-01 - 2014-12-31
#we want 1895-01-01 - 2014-12-31
grep('1895-01-01',colnames(prcpmat))
#column 13490 is 1895-01-01
prcpmat=prcpmat[,13490:57318]
prcpmat=cbind(ll,prcpmat) #attach lat/lon for each station
prcpmat=apply(prcpmat,2,FUN=as.numeric) #convert matrix to numeric


## stations active vs time
c=c()
for (i in 3:43831){
  y=complete.cases(prcpmat[,i])
  x=length(which(y))
  c=c(c,x)
}
t=seq(as.Date("1895/01/01"),as.Date("2014/12/31"),length.out=43829)
plot(t,c,type='l',xlab='Time',ylab='Active Stations',main='Active USHCN Stations 1895/01/01-2014/12/31')

## station monthly averages
ojan=grep('-01-',colnames(prcpmat))
ofeb=grep('-02-',colnames(prcpmat))
omar=grep('-03-',colnames(prcpmat))
oapr=grep('-04-',colnames(prcpmat))
omay=grep('-05-',colnames(prcpmat))
ojun=grep('-06-',colnames(prcpmat))
ojul=grep('-07-',colnames(prcpmat))
oaug=grep('-08-',colnames(prcpmat))
osep=grep('-09-',colnames(prcpmat))
ooct=grep('-10-',colnames(prcpmat))
onov=grep('-11-',colnames(prcpmat))
odec=grep('-12-',colnames(prcpmat))

janrm=rowMeans(prcpmat[,ojan],na.rm=T)
febrm=rowMeans(prcpmat[,ofeb],na.rm=T)
marrm=rowMeans(prcpmat[,omar],na.rm=T)
aprrm=rowMeans(prcpmat[,oapr],na.rm=T)
mayrm=rowMeans(prcpmat[,omay],na.rm=T)
junrm=rowMeans(prcpmat[,ojun],na.rm=T)
julrm=rowMeans(prcpmat[,ojul],na.rm=T)
augrm=rowMeans(prcpmat[,oaug],na.rm=T)
seprm=rowMeans(prcpmat[,osep],na.rm=T)
octrm=rowMeans(prcpmat[,ooct],na.rm=T)
novrm=rowMeans(prcpmat[,onov],na.rm=T)
decrm=rowMeans(prcpmat[,odec],na.rm=T)
statmean=cbind(janrm,febrm,marrm,aprrm,mayrm,junrm,julrm,augrm,seprm,octrm,novrm,decrm)

prcpmat[,ojan]=prcpmat[,ojan]-statmean[,1]
prcpmat[,ofeb]=prcpmat[,ofeb]-statmean[,2]
prcpmat[,omar]=prcpmat[,omar]-statmean[,3]
prcpmat[,oapr]=prcpmat[,oapr]-statmean[,4]
prcpmat[,omay]=prcpmat[,omay]-statmean[,5]
prcpmat[,ojun]=prcpmat[,ojun]-statmean[,6]
prcpmat[,ojul]=prcpmat[,ojul]-statmean[,7]
prcpmat[,oaug]=prcpmat[,oaug]-statmean[,8]
prcpmat[,osep]=prcpmat[,osep]-statmean[,9]
prcpmat[,ooct]=prcpmat[,ooct]-statmean[,10]
prcpmat[,onov]=prcpmat[,onov]-statmean[,11]
prcpmat[,odec]=prcpmat[,odec]-statmean[,12]

write.csv(prcpmat,'anomobs.csv',row.names=F)
anomobs=fread.csv('anomobs.csv',header=T,sep=',')


#### MAP STATION DATA ####


library(ggplot2)
library(ggmap)

myLocation <- c(-133, 24, -60, 60) #US bounds
#texas=c(-108, 25, -92, 37)
#lon-lat of lowerleft and lon-lat of upperright
#maptype = c("terrain", "toner", "watercolor")
#maptype = c("roadmap", "terrain", "satellite", "hybrid")
myMap = get_map(location = myLocation, source="google", maptype="roadmap", crop=TRUE)
#myMap = get_map(location = 'united states', zoom = 4, source = 'google', maptype="terrain")
#ggmap(myMap)

#cll=read.csv('cll3.csv',check.names=FALSE)
cll=data.frame(gridded[,2:3])
cll=data.frame(master[,1:2])
colnames(cll)=c('Lat','Lon')

mmDay=gridded$`2012-10-30`
mmDay=master[,31]

hi=max(mmDay,na.rm=T)
lo=min(mmDay,na.rm=T)
mid=(hi+lo)/2
range(mmDay,na.rm=T) #adjust individually
ggmap(myMap) + geom_point(data=cll, mapping=aes(x=cll$Lon, y=cll$Lat, colour=mmDay), size=1.5) +
  scale_colour_gradient2(limits=c(lo,hi),low="white",mid="blue", midpoint=mid, high = "red", space="rgb") #+ scale_x_continuous(limits=c(-125,-65))


############ CMORPH CREATION ################
#cmorph dim
#1440 480
lon=seq(0.125,360,0.25)
lat=seq(-59.875,59.875,0.25)

LAT=rep(lat,1440)
LON=rep(lon,each=480)
cll=cbind(LAT,LON)

##mask
#cmorph data is now oriented south to north
landmask=read.table('UMD60mask0.25.asc')
landmask=landmask[,3:5]
mask=which(landmask[,1]<60&landmask[,1]>-60)
landmask=landmask[mask,]
landmask=landmask[order(landmask$V4),]
lon2=seq(180.125,360,0.25)
lon2=rep(lon2,each=480)
landmask[1:345600,2]=lon2
landmask=landmask[order(landmask$V4),]
land=landmask[,3]
#mask now matches cmorph format
#add mask 
masksorted=cbind(land,LAT,LON)

top = 49.3457868 # north lat
left = 235.216 # west long
right = 293.049 # east long
bottom =  24.7433195 # south lat

#final mask is usland
usland=which(masksorted[,1]==1&masksorted[,2]>bottom&masksorted[,2]<top&masksorted[,3]>left&masksorted[,3]<right)

##master
#LAT=rep(lat,1440)
#LON=rep(lon,each=480)
FLAT=LAT[usland]
FLON=LON[usland]
master=cbind(FLAT,FLON)
master[,2]=master[,2]-360



#for each folder with cmorph data (can be any size you choose), change the "filenames" directory (make sure to do it chronologically)
filenames=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/1998')
#change the working directory to the same folder as above
setwd("C:/Users/Scott/Desktop/Sam Stuff/cmorph/1998") 
#do this until all cmorph data is read in (again, make sure to do all files in order)
#function for each bin
for (i in 1:length(filenames)){
  rawbin=readBin(filenames[i],numeric(),endian='little',n=691200,size=4)
  rawbin[rawbin==-999]=NA
  binvec=c(rawbin)
  binmat=matrix(binvec,nrow=1440,ncol=480)
  tbinmat=t(binmat)
  tbinvec=as.vector(tbinmat)
  landvec=tbinvec[usland]
  master=cbind(master,landvec)
}

#image(rotate(tbinmat),useRaster=T,axes=F,col=rainbow(20))

##missing date in 2014, 8/29

##remove mexico
mxco=which(master[,1]<31.2&master[,2]<(-106))
master=master[-mxco,]
mxco=which(master[,1]<28.5&master[,2]<(-100.5))
master=master[-mxco,]


#write.csv(master,file='C:/Users/Scott/Desktop/sum proj/cmorph19982016.csv')

#column names for master cmorph
#just use whatever folders contain all your cmorph data, in my case I had them by year
files98=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/1998')
files99=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/1999')
files00=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2000')
files01=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2001')
files02=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2002')
files03=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2003')
files04=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2004')
files05=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2005')
files06=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2006')
files07=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2007')
files08=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2008')
files09=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2009')
files10=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2010')
files11=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2011')
files12=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2012')
files13=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2013')
files14=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2014')
files15=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2015')
files16=list.files('C:/Users/Scott/Desktop/Sam Stuff/cmorph/2016')
namecol=c(files98,files99,files00,files01,files02,files03,files04,files05,files06,files07,files08,files09,files10,files11,files12,files13,files14,files15,files16)
namevec=c()
for (i in 1:length(namecol)){
  namevec=c(namevec,substr(namecol[i],33,40))
}
namevec=format(as.Date(namevec, format = "%Y%m%d"), "%Y-%m-%d")
col=c("Lat","Lon")
col=c(col,namevec)
colnames(master)=col
write.csv(master,file='C:/Users/Scott/Desktop/sum proj/cmorph19982014.csv')


#### CMORPH ANOM ####
#cmorph=fread.csv('C:/Users/Scott/Desktop/sum proj/cmorph19982014.csv',header=TRUE,sep=',')
#only use through 2014
cmorph=cmorph[,3:length(colnames(cmorph))] #leave out lat/lon
#get all months by column position
jan=grep('-01-',colnames(cmorph))
feb=grep('-02-',colnames(cmorph))
mar=grep('-03-',colnames(cmorph))
apr=grep('-04-',colnames(cmorph))
may=grep('-05-',colnames(cmorph))
jun=grep('-06-',colnames(cmorph))
jul=grep('-07-',colnames(cmorph))
aug=grep('-08-',colnames(cmorph))
sep=grep('-09-',colnames(cmorph))
oct=grep('-10-',colnames(cmorph))
nov=grep('-11-',colnames(cmorph))
dec=grep('-12-',colnames(cmorph))


janrm=rowMeans(cmorph[,jan],na.rm=T)
febrm=rowMeans(cmorph[,feb],na.rm=T)
marrm=rowMeans(cmorph[,mar],na.rm=T)
aprrm=rowMeans(cmorph[,apr],na.rm=T)
mayrm=rowMeans(cmorph[,may],na.rm=T)
junrm=rowMeans(cmorph[,jun],na.rm=T)
julrm=rowMeans(cmorph[,jul],na.rm=T)
augrm=rowMeans(cmorph[,aug],na.rm=T)
seprm=rowMeans(cmorph[,sep],na.rm=T)
octrm=rowMeans(cmorph[,oct],na.rm=T)
novrm=rowMeans(cmorph[,nov],na.rm=T)
decrm=rowMeans(cmorph[,dec],na.rm=T)
statmean=cbind(janrm,febrm,marrm,aprrm,mayrm,junrm,julrm,augrm,seprm,octrm,novrm,decrm)

cmorph[,jan]=cmorph[,jan]-statmean[,1]
cmorph[,feb]=cmorph[,feb]-statmean[,2]
cmorph[,mar]=cmorph[,mar]-statmean[,3]
cmorph[,apr]=cmorph[,apr]-statmean[,4]
cmorph[,may]=cmorph[,may]-statmean[,5]
cmorph[,jun]=cmorph[,jun]-statmean[,6]
cmorph[,jul]=cmorph[,jul]-statmean[,7]
cmorph[,aug]=cmorph[,aug]-statmean[,8]
cmorph[,sep]=cmorph[,sep]-statmean[,9]
cmorph[,oct]=cmorph[,oct]-statmean[,10]
cmorph[,nov]=cmorph[,nov]-statmean[,11]
cmorph[,dec]=cmorph[,dec]-statmean[,12]


write.csv(cmorph,'cmorphanom.csv',row.names=FALSE)

####monthly SVD####
cmorph=fread.csv('cmorphanom.csv',header=TRUE,sep=',')

cmorph[is.na(cmorph)]=0
#cll=read.csv('cll3.csv',check.names=FALSE)
#cmorph=cbind(cll[,1],cll[,2],cmorph)
jan=grep('-01-',colnames(cmorph))
feb=grep('-02-',colnames(cmorph))
mar=grep('-03-',colnames(cmorph))
apr=grep('-04-',colnames(cmorph))
may=grep('-05-',colnames(cmorph))
jun=grep('-06-',colnames(cmorph))
jul=grep('-07-',colnames(cmorph))
aug=grep('-08-',colnames(cmorph))
sep=grep('-09-',colnames(cmorph))
oct=grep('-10-',colnames(cmorph))
nov=grep('-11-',colnames(cmorph))
dec=grep('-12-',colnames(cmorph))


jansvd=svd(cmorph[,jan])
janeof=jansvd$u
janeig=jansvd$d #first 20 EOFs 34% var

febsvd=svd(cmorph[,feb])
febeof=febsvd$u
febeig=febsvd$d #first 20 EOFs 32% var

marsvd=svd(cmorph[,mar])
mareof=marsvd$u
mareig=marsvd$d #first 20 EOFs 28% var

aprsvd=svd(cmorph[,apr])
apreof=aprsvd$u
apreig=aprsvd$d #first 20 EOFs 24% var

maysvd=svd(cmorph[,may])
mayeof=maysvd$u
mayeig=maysvd$d #first 20 EOFs 19% var

junsvd=svd(cmorph[,jun])
juneof=junsvd$u
juneig=junsvd$d #first 20 EOFs 18% var

julsvd=svd(cmorph[,jul])
juleof=julsvd$u
juleig=julsvd$d #first 20 EOFs 15% var

augsvd=svd(cmorph[,aug])
augeof=augsvd$u
augeig=augsvd$d #first 20 EOFs 15% var

sepsvd=svd(cmorph[,sep])
sepeof=sepsvd$u
sepeig=sepsvd$d #first 20 EOFs 18% var

octsvd=svd(cmorph[,oct])
octeof=octsvd$u
octeig=octsvd$d #first 20 EOFs 22% var

novsvd=svd(cmorph[,nov])
noveof=novsvd$u
noveig=novsvd$d #first 20 EOFs 27% var

decsvd=svd(cmorph[,dec])
deceof=decsvd$u
deceig=decsvd$d #first 20 EOFs 33% var

#truncate at 480 b/c feb
eigmat=cbind(janeig[1:480],febeig,mareig[1:480],apreig[1:480],
             mayeig[1:480],juneig[1:480],juleig[1:480],augeig[1:480],
             sepeig[1:480],octeig[1:480],noveig[1:480],deceig[1:480])
colnames(eigmat)=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

eigcum=matrix(0,nrow=480,ncol=12)
for (i in 1:480){
  eigcum[i,12]=sum(eigmat[1:i,12])/sum(eigmat[,12])
}
eigcum=eigcum*100
colnames(eigcum)=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')



#### eig ####
eig=read.csv('cmorpheig.csv',check.names=FALSE)
#length 6938
eig=eig^(2/6938)
sum(eig[1:20])/sum(eig)
#for 95%:6600 90%:6300 75%:5200 50%:3500
#all eig are basically 1, seems pointless
#without ^2/6938, 95% is ~4000, 50% ~ 400


#### MAP EOF ####
eof=read.csv("new20.csv", check.names = FALSE) 
#eof1000=read.csv('first1000eof.csv',check.names=FALSE)
cll3=read.csv('cll3.csv',check.names=FALSE)
#cll[,2]=cll[,2]-360
rm(list=c("myMap"))
library(ggplot2)
library(ggmap)
myLocation <- c(-133, 24, -60, 60)
#lon-lat of lowerleft and lon-lat of upperright
#maptype = c("terrain", "toner", "watercolor")
#maptype = c("roadmap", "terrain", "satellite", "hybrid")
myMap = get_map(location = myLocation, source="google", maptype="terrain", crop=FALSE)
#ggmap(myMap)

eof=data.frame(eof)
eofcol=deceof[,1]
hi=max(eofcol)
lo=min(eofcol)
mid=(hi+lo)/2
range(eofcol)
ggmap(myMap) + geom_point(data=cll, mapping=aes(x=Lon, y=Lat, colour=eofcol), size=1) +
  scale_colour_gradient2(limits=c(lo,hi),low="red",mid="green", midpoint=mid, high = "blue", space="rgb") #+ scale_x_continuous(limits=c(-125,-65))
#5753

#### GRIDDING ####
#sorry about the confusing variable names in this section, just run it all and it works
#make sure you have the cmorph data read in correctly before this, as you will need the lat/lon from that
#set the cmorph lat/lon as variable "cll"
obs=fread('anomobs.csv',sep=',',header=T)
#cll=read.csv('cll3.csv',check.names=FALSE)


iround <- function(x, interval){
  ## Round numbers to desired interval
  ##
  ## Args:
  ##   x:        Numeric vector to be rounded
  ##   interval: The interval the values should be rounded towards.
  ## Retunrs:
  ##   a numeric vector with x rounded to the desired interval.
  ##
  interval[ifelse(x < min(interval), 1, findInterval(x, interval))]
}

x=sort(cll[,1])
y=sort(cll[,2])
gridlat=iround(obs$lat, x) #round obs lat to nearest 0.25x0.25 lat
gridlon=iround(obs$lon, y) #round obs lon to nearest 0.25x0.25 lon

gridll=cbind(gridlat,gridlon)
grid=c()
for (i in 1:14802){
  a=which((gridll[,1]==cll[i,1])&(gridll[,2]==cll[i,2]))
  grid=c(grid,a)
  }
gridb=c()
for (i in 1:1218){
     b=which(cll[,1]==gridll[i,1]&cll[,2]==gridll[i,2])
    gridb=c(gridb,b)
}

gridgrid=cbind(sort(grid),gridb)
gridgrid=as.data.frame(gridgrid)
colnames(gridgrid)=c('grid','gridb')

g=rep(0,1218)
gridded=cbind(g,obs)
gridded[gridgrid$grid,1]=gridgrid$gridb
colnames(gridded)[1]=c("grid")

#40 grids missing
missing=which(gridded[,1]==0)

#this loop will find the nearest grid point for the ungridded obs
#it takes a few minutes

nearest=c()
for (i in 1:length(missing)){ 
  m=c()
  for (j in 1:14802){
    dif=abs(cll[j,1]-gridll[missing[i],1]) + abs(cll[j,2]-gridll[missing[i],2])
    m=c(m,dif)
  }
  nearest=c(nearest,which.min(m))
}

#add in the newly gridded points
gridded[missing,1]=nearest

#now must deal with overlapping stations
#this loop finds duplicate grid points, averages those rows into a new row, removes the old rows
#and adds the new average row
notuni=which(duplicated(gridded$grid))
uninotuni=which(!duplicated(gridded[notuni,1]))
nu=notuni[uninotuni]
while (length(nu)!=0){
  
  same=which(gridded[,1]==as.numeric(gridded[nu[1],1]))
  rows=gridded[same,]
  newrow=as.matrix(colMeans(rows,na.rm=T))
  tnr=t(newrow)
  gridded=gridded[-same,]
  gridded=rbind(gridded,tnr)
  
  notuni=which(duplicated(gridded$grid))
  uninotuni=which(!duplicated(gridded[notuni,1]))
  nu=notuni[uninotuni]
}



#change units from .01s/inch/day to mm/day
gridded[,4:43832]=gridded[,4:43832]*(25.4/100) #1/100s inches to inches
range(gridded[,4:43832],na.rm=TRUE)
#-22.95022 2406.63876

write.csv(gridded,'anomobsgrid2.csv',row.names=FALSE)

#### RECON ####
#obs in 4 chunks, need gridll for last 3
library(data.table)
gridobs=fread('gridobs5080.csv', header=T, sep = ',')
dim(gridobs) #9520=1218 9136, 2050= 10961, 5080=10960, 8014= 12787
cll=read.csv('cll3.csv',check.names=FALSE)
#cll[,2]=cll[,2]-360
gridll=read.csv('gridll.csv',check.names=F)
gridobs=cbind(gridll,gridobs)

##monthly
ojan=grep('-01-',colnames(gridobs))
ofeb=grep('-02-',colnames(gridobs))
omar=grep('-03-',colnames(gridobs))
oapr=grep('-04-',colnames(gridobs))
omay=grep('-05-',colnames(gridobs))
ojun=grep('-06-',colnames(gridobs))
ojul=grep('-07-',colnames(gridobs))
oaug=grep('-08-',colnames(gridobs))
osep=grep('-09-',colnames(gridobs))
ooct=grep('-10-',colnames(gridobs))
onov=grep('-11-',colnames(gridobs))
odec=grep('-12-',colnames(gridobs))


janeof=fread('janeof.csv',header=T, sep = ',')
janeof=janeof[,1:220]
febeof=fread('febeof.csv',header=T, sep = ',')
febeof=febeof[,1:220]
mareof=fread('mareof.csv',header=T, sep = ',')
mareof=mareof[,1:220]
apreof=fread('apreof.csv',header=T, sep = ',')
apreof=apreof[,1:220]
mayeof=fread('mayeof.csv',header=T, sep = ',')
mayeof=mayeof[,1:220]
juneof=fread('juneof.csv',header=T, sep = ',')
juneof=juneof[,1:220]
juleof=fread('juleof.csv',header=T, sep = ',')
juleof=juleof[,1:220]
augeof=fread('augeof.csv',header=T, sep = ',')
augeof=augeof[,1:220]
sepeof=fread('sepeof.csv',header=T, sep = ',')
sepeof=sepeof[,1:220]
octeof=fread('octeof.csv',header=T, sep = ',')
octeof=octeof[,1:220]
noveof=fread('noveof.csv',header=T, sep = ',')
noveof=noveof[,1:220]
deceof=fread('deceof.csv',header=T, sep = ',')
deceof=deceof[,1:220]

colnames(noveof)=c('E1','E2','E3','E4','E5','E6','E7','E8','E9','E10',
                   'E11','E12','E13','E14','E15','E16','E17','E18','E19',
                   'E20','E21','E22','E23','E24','E25','E26','E27','E28',
                   'E29','E30','E31','E32','E33','E34','E35','E36','E37',
                   'E38','E39','E40','E41','E42','E43','E44','E45','E46',
                   'E47','E48','E49','E50','E51','E52','E53','E54','E55',
                   'E56','E57','E58','E59','E60','E61','E62','E63','E64',
                   'E65','E66','E67','E68','E69','E70','E71','E72','E73',
                   'E74','E75','E76','E77','E78','E79','E80','E81','E82',
                   'E83','E84','E85','E86','E87','E88','E89','E90','E91',
                   'E92','E93','E94','E95','E96','E97','E98','E99','E100',
                   'E101','E102','E103','E104','E105','E106','E107','E108','E109','E110',
                   'E111','E112','E113','E114','E115','E116','E117','E118','E119',
                   'E120','E121','E122','E123','E124','E125','E126','E127','E128',
                   'E129','E130','E131','E132','E133','E134','E135','E136','E137',
                   'E138','E139','E140','E141','E142','E143','E144','E145','E146',
                   'E147','E148','E149','E150','E151','E152','E153','E154','E155',
                   'E156','E157','E158','E159','E160','E161','E162','E163','E164',
                   'E165','E166','E167','E168','E169','E170','E171','E172','E173',
                   'E174','E175','E176','E177','E178','E179','E180','E181','E182',
                   'E183','E184','E185','E186','E187','E188','E189','E190','E191',
                   'E192','E193','E194','E195','E196','E197','E198','E199','E200',
                   'E201','E202','E203','E204','E205','E206','E207','E208','E209',
                   'E210','E211','E212','E213','E214','E215','E216','E217','E218',
                   'E219','E220')

recon=matrix(0,nrow=14802,ncol=10960)
recon[,2]=cll$Lat
recon[,3]=cll$Lon
for (i in onov) {
  y=complete.cases(gridobs[,i]) #obs with data
  v=which(y) #row of obs data
  u=gridobs[v,1] #grid numbers of obs data
  datr=gridobs[v,i] #obs data for column
  eofr=noveof[u,] #eofs at grids of obs
  df=data.frame(eofr,datr) #combine eofs with obs data
  #rs=rowSums(eofr)
  reg=lm(formula=datr~E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15
         +E16+E17+E18+E19+E20+E21+E22+E23+E24+E25+E26+E27+E28+E29+E30+E31
         +E32+E33+E34+E35+E36+E37+E38+E39+E40+E41+E42+E43+E44+E45+E46+E47
         +E48+E49+E50+E51+E52+E53+E54+E55+E56+E57+E58+E59+E60+E61+E62+E63
         +E64+E65+E66+E67+E68+E69+E70+E71+E72+E73+E74+E75+E76+E77+E78+E79
         +E80+E81+E82+E83+E84+E85+E86+E87+E88+E89+E90+E91+E92+E93+E94+E95
         +E96+E97+E98+E99+E100+E101+E102+E103+E104+E105+E106+E107+E108+E109
         +E110+E111+E112+E113+E114+E115+E116+E117+E118+E119+E120+E121+E122
         +E123+E124+E125+E126+E127+E128+E129+E130+E131+E132+E133+E134+E135
         +E136+E137+E138+E139+E140+E141+E142+E143+E144+E145+E146+E147+E148
         +E149+E150+E151+E152+E153+E154+E155+E156+E157+E158+E159+E160+E161
         +E162+E163+E164+E165+E166+E167+E168+E169+E170+E171+E172+E173+E174
         +E175+E176+E177+E178+E179+E180+E181+E182+E183+E184+E185+E186+E187
         +E188+E189+E190+E191+E192+E193+E194+E195+E196+E197+E198+E199+E200
         +E201+E202+E203+E204+E205+E206+E207+E208+E209+E210+E211+E212+E213
         +E214+E215+E216+E217+E218+E219+E220, data=df) #lm with obs data and eofs
  coe=reg$coefficients
  c1=rep(1,14802) #vector of 1's for each grid
  res=cbind(c1,noveof) #combine 1's with all eofs
  res=as.matrix(res)
  recon[,i]=res%*%coe # row dot prod
}
colnames(recon)=colnames(gridobs)


#### MAP RECON ####

cll=read.csv('cll3.csv',check.names=FALSE)
colnames(cll)=c('Latitude','Longitude')
#cll[,2]=cll[,2]-360
#rm(list=c("myMap"))
library(ggplot2)
library(ggmap)
myLocation <- c(-125, 24.5, -65, 49.5)
myLocation = c(lonmin-2,latmin-2,lonmax+2,latmax+2)
#texas=c(-108, 25, -92, 37)
#katrina=c(-95, 28, -80, 42)
#all_states=map_data()
#lon-lat of lowerleft and lon-lat of upperright
#maptype = c("terrain", "toner", "watercolor")
#maptype = c("roadmap", "terrain", "satellite", "hybrid")
myMap = get_map(location = myLocation, source="stamen", maptype='toner',color='bw')#, crop=FALSE)
#ggmap(myMap)
#plotname=paste('plot8014_', i, ".png", sep='')
#df=as.data.frame(recon)
#find day 2005-08-30(9377) 2012-10-30(11995) texas 6869
grep('1927-04-01',colnames(recon))
#mm=as.data.frame(recon[,2665])
mm=stormplot[,3]
days=colnames(recon)
day=days[2665]
#plothead=paste("United States Daily Precipitation Climatology [mm]:",day)
plothead=('May 24 - June 11 1903 Flood Total Precipitation [mm]')
ma=max(mm)
mi=min(mm)
range(mm)
#plot.new()
ggmap(myMap)+ geom_point(data=stormplot[,1:2], mapping=aes(x=stormplot[,2], y=stormplot[,1], colour=mm), size=2.5,alpha=.7) +
  scale_color_gradientn(limits=c(1,550),breaks=c(1,10,50,100,250,500),trans = 'log10',colours=c('grey','lightblue','blue',
                                                                                                'green','yellow','orange',
                                                                                                'red','pink','maroon'),guide=guide_colorbar(draw.ulim = T,draw.llim = T))+
  labs(x='Longitude',y='Latitude')+
  ggtitle(plothead)+
  theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"),legend.key.width=unit(.4,'inch'),legend.key.height=unit(1,'inch'))+
  theme(legend.text=element_text(size=12), legend.title=element_text(size=15))

ggsave('DecClimatology.png',width=16,height=9)
 #5753
#3/20/1899
##plot all save all
days=colnames(recon)
recon=as.data.frame(recon)
for (i in 7:9136){
  plotname=paste('plot8520_', i, ".png", sep='')
  mm=(recon[,i])
  day=days[i]
  plothead=paste("United States Daily Precipitation Anomalies [mm]:", day)
  ma=max(mm)
  mi=min(mm)
  range(mm)
  ggmap(myMap)+ geom_point(data=cll, mapping=aes(x=Longitude, y=Latitude, colour=mm), size=1.5,alpha=.7) +
    scale_color_gradientn(limits=c(1,500),breaks=c(1,10,50,100,250,500),trans = 'log10',colours=c('grey','lightblue','blue',
                                                    'green','yellow','orange',
                                                    'red','pink','maroon'),guide=guide_colorbar(draw.ulim = T,draw.llim = T))+
    labs(x='Longitude',y='Latitude')+
    ggtitle(plothead)+
    theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"),legend.key.width=unit(.4,'inch'),legend.key.height=unit(1,'inch'))+
    theme(legend.text=element_text(size=12), legend.title=element_text(size=15))
  ggsave(file=plotname,width=16,height=9)
}

##recon avg prcp time series
#first without cos weighting
cs=colMeans(recon[,4:23379])
ti=as.Date(colnames(gridobs[4:23379]))
plot(ti,cs,type='l',xlab='Time',ylab='Mean Prcp [mm]',main='Unweighted Means')

#cos weighted, have to change to matrix after using fread
crecon=matrix(0,nrow=14802,ncol=9133)
for (i in 1:14802){
  crecon[i,]=cos(recon[i,2]*pi/180)*recon[i,]}
ccs=colSums(crecon[,4:9133])/sum(cos(recon[,2]*pi/180))
ti=as.Date(colnames(recon[4:9133]))
plot(ti,ccs,type='l',xlab='Time',ylab='Mean Prcp [mm]',main='Cos-weighted Means')

#plot avg over time
c1=read.csv('cosweight9520.csv',check.names=F)
c2=read.csv('cosweight2050.csv',check.names=F)
c3=read.csv('cosweight5080.csv',check.names=F)
c4=read.csv('cosweight8014.csv',check.names=F)
cw=c(t(c1),t(c2),t(c3),t(c4))
t=seq(as.Date("1895/01/01"),as.Date("2014/12/31"),length.out=43829)
plot(t,cw,type='l',main='US Anomalous Precipitation Climatology 1895/01/01-2014/12/31',xlab='Time',ylab='Cos-weighted Mean Anomaly [mm]')
lm8514=lm(cw~t)
lm8520=lm(cw[1:9131]~t[1:9131])
lm2050=lm(cw[9131:20089]~t[9131:20089])
lm5080=lm(cw[20089:31046]~t[20089:31046])
lm8014=lm(cw[31046:43829]~t[31046:43829])
abline(lm8514,col='red')
clip(1,1000,-20,20)
abline(lm8520,col='blue')
abline(lm2050,col='green')
abline(lm5080,col='purple')
abline(lm8014,col='orange')


#### prcp volume & area####
recon=as.data.frame(recon)
stormday=grep('1913-01',colnames(recon))
stormday=colnames(recon[,3071:3089])
#recon[,2]=lat [,3]=lon
#ohio river valley coordinates jan 1937
latmin=36
latmax=41.5
lonmin=-90
lonmax=-80.5
#missouri/mississippi coords 5/24 - 6/11 1903 (use stormday=recon[,3071:3089])
latmin=38
latmax=43.5
lonmin=-103
lonmax=-89.5
stormbox=which(recon[,2]>latmin &recon[,2]<latmax & recon[,3]>lonmin &recon[,3]<lonmax)
storm=cbind(recon[stormbox,2],recon[stormbox,3],recon[stormbox,stormday])
totalVol=0
for (j in 1:length(stormday)){
  day=data.frame(storm[,j+2])
  great20=which(day[,1]>0)
  storm20=data.frame(cbind(storm[great20,1],storm[great20,2],day[great20,]))
  vol20=c()
  area20=c()
  for (i in 1:length(great20)){
    vol20[i]=cos(storm20[i,1]*pi/180)*(.25/180)^2*storm20[i,3]*.001*6371000^2
    area20[i]=cos(storm20[i,1]*pi/180)*(.25/180)^2*6371000^2
  }
  totalVol=totalVol + sum(vol20)
  sum(area20)
 
}

#sum each grid for plot
stormtime=storm[,3:(2+length(stormday))]
stormtotal=rowSums(stormtime)
stormplot=cbind(storm[,1:2],stormtotal)
#use 2.5 dot size


#### averages ####
##annual##
recon=as.data.frame(recon)
recon=recon+climm
avgs8014=c()
for (i in 1895:1919){
   w=grep(i,colnames(recon))
   avgs8014=c(avgs8014,mean(colMeans(recon[,w])))
}
avgs=c(avgs8014,avgs2050,avgs5080,avgs9520)
#daily
#recon=fread('recon5080EOF220.csv',sep=',',header=T)
recon=recon[,4:length(recon)]
#recon[recon[]<0]=0
daily5080=colMeans(recon)
daily=c(daily9520,daily2050,daily5080,daily8014)
x=1:length(daily)
slope=lm(daily~x)

#monthly
#recon=fread('recon8014EOF220.csv',sep=',',header=T)
#recon=recon[,7:length(recon)]
recon=as.data.frame(recon)
recon=recon+climm


w=grep('-01-',colnames(recon))
mon=recon[,w]
jan9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  jan9520=c(jan9520,mean(colMeans(mon[,w])))
}
#jan=c()
jan=c(jan,jan9520)
w=grep('-02-',colnames(recon))
mon=recon[,w]
feb9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  feb9520=c(feb9520,mean(colMeans(mon[,w])))
}
#feb=c()
feb=c(feb,feb9520)
w=grep('-03-',colnames(recon))
mon=recon[,w]
mar9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  mar9520=c(mar9520,mean(colMeans(mon[,w])))
}
#mar=c()
mar=c(mar,mar9520)
w=grep('-04-',colnames(recon))
mon=recon[,w]
apr9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  apr9520=c(apr9520,mean(colMeans(mon[,w])))
}
#apr=c()
apr=c(apr,apr9520)
w=grep('-05-',colnames(recon))
mon=recon[,w]
may9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  may9520=c(may9520,mean(colMeans(mon[,w])))
}
#may=c()
may=c(may,may9520)
w=grep('-06-',colnames(recon))
mon=recon[,w]
jun9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  jun9520=c(jun9520,mean(colMeans(mon[,w])))
}
#jun=c()
jun=c(jun,jun9520)
w=grep('-07-',colnames(recon))
mon=recon[,w]
jul9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  jul9520=c(jul9520,mean(colMeans(mon[,w])))
}
#jul=c()
jul=c(jul,jul9520)
w=grep('-08-',colnames(recon))
mon=recon[,w]
aug9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  aug9520=c(aug9520,mean(colMeans(mon[,w])))
}
#aug=c()
aug=c(aug,aug9520)
w=grep('-09-',colnames(recon))
mon=recon[,w]
sep9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  sep9520=c(sep9520,mean(colMeans(mon[,w])))
}
#sep=c()
sep=c(sep,sep9520)
w=grep('-10-',colnames(recon))
mon=recon[,w]
oct9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  oct9520=c(oct9520,mean(colMeans(mon[,w])))
}
#oct=c()
oct=c(oct,oct9520)
w=grep('-11-',colnames(recon))
mon=recon[,w]
nov9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  nov9520=c(nov9520,mean(colMeans(mon[,w])))
}
#nov=c()
nov=c(nov,nov9520)
w=grep('-12-',colnames(recon))
mon=recon[,w]
dec9520=c()
for (i in 1980:2014){
  w=grep(i,colnames(mon))
  dec9520=c(dec9520,mean(colMeans(mon[,w])))
}
#dec=c()
dec=c(dec,dec9520)

