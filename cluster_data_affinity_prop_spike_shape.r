require(e1071)
require(gplots)
require(apcluster)
require(cluster)
setwd("C:/Users/menonv/OneDrive - Howard Hughes Medical Institute/Misc/ct/20170727/")
dat=read.table("GLIF_param_plus_spike_features_7_27_17.csv",sep="\t",as.is=T,row.names=1,check.names=F,header=T)
metadata=dat[,1:2]
crecols=read.csv("../cre_colors.csv",as.is=T,header=F)
newcols=rgb(crecols[,2:4],maxColorValue = 255)
names(newcols)=crecols[,1]
colvec=newcols[match(metadata$cre,crecols[,1])]
fulldat=dat[,-c(1:2)]
for (ii in 1:ncol(fulldat)) {
  if (min(fulldat[,ii])*max(fulldat[,ii])>0) {
    if (min(fulldat[,ii])>0) {
      if (skewness(fulldat[,ii])>skewness(log10(fulldat[,ii]))) {
        fulldat[,ii]=log10(fulldat[,ii])
        print(ii)
      }
    } else {
      if (skewness(-fulldat[,ii])>skewness(log10(-fulldat[,ii]))) {
        fulldat[,ii]=log10(-fulldat[,ii])
        print(ii)
      }
    }
  }#
}

fulldat_all=fulldat


###features###
featdat=read.table("features_7_27_17.csv",as.is=T,row.names=1,check.names=F,sep=",",header=T)
featmetadata=featdat[,1:2]
featfulldat=featdat[,c("tau","ri","vrest","threshold_i_long_square","threshold_v_long_square","peak_v_long_square","fast_trough_v_long_square","trough_v_long_square","upstroke_downstroke_ratio_long_square","upstroke_downstroke_ratio_short_square","sag","f_i_curve_slope","latency","max_burstiness_across_sweeps")]
for (ii in 1:ncol(featfulldat)) {
  if (min(featfulldat[,ii])*max(featfulldat[,ii])>0) {
    if (min(featfulldat[,ii])>0) {
      if (skewness(featfulldat[,ii])>skewness(log10(featfulldat[,ii]))) {
        featfulldat[,ii]=log10(featfulldat[,ii])
        print(ii)
      }
    } else {
      if (skewness(-featfulldat[,ii])>skewness(log10(-featfulldat[,ii]))) {
        featfulldat[,ii]=log10(-featfulldat[,ii])
        print(ii)
      }
    }
  }
}
featfulldat_all=featfulldat

runaffprop=function(dat,k) {
  mmm=apclusterK(function (x){x=1-cor(t(x),method="pearson")},dat,K=k,seed=1)
  outvec=rep(0,nrow(dat))
  print(dim(ncol(dat)))
  for (ii in 1:length(mmm)) {
    outvec[mmm[[ii]]]=ii
  }
  outlist=list()
  outlist[['cluster']]=outvec
  return(outlist)
}


###Run affiinity propagation on model parameters###
allclustgap=list()
for (nameval in c("glif4_spike_shape","glif3_spike_shape","glif2_spike_shape","glif1_spike_shape")) {
  if (nameval=="glif1_spike_shape") {keepcols=c(1,3,4,5,8,13:16)}
  if (nameval=="glif2_spike_shape") {keepcols=c(1,3,4,5,8,9,10,13:16)}
  if (nameval=="glif3_spike_shape") {keepcols=c(2,3,4,5,6,7,8,13:16)}
  if (nameval=="glif4_spike_shape") {keepcols=c(2,3,4,5,6,7,8,9,10,13:16)}
  
  
  if (nameval %in% c("features","featuresnospike")) {
    startmat=featfulldat_all
  } else {
    startmat=fulldat_all
  }
  newstart=scale(startmat[,apply(startmat,2,var)>0])
  testnumclusts=1:25
  allclustgap[[nameval]]=clusGap(newstart[,keepcols],runaffprop,K.max=25)
  save(allclustgap,file="allclustgap_spike_shape.rda")
}
load("allclustgap_spike_shape.rda")

for (nameval in c("glif4_spike_shape","glif3_spike_shape","glif2_spike_shape","glif1_spike_shape")) {
  if (nameval=="glif1_spike_shape") {keepcols=c(1,3,4,5,8,13:16)}
  if (nameval=="glif2_spike_shape") {keepcols=c(1,3,4,5,8,9,10,13:16)}
  if (nameval=="glif3_spike_shape") {keepcols=c(2,3,4,5,6,7,8,13:16)}
  if (nameval=="glif4_spike_shape") {keepcols=c(2,3,4,5,6,7,8,9,10,13:16)}
  
  if (nameval %in% c("features","featuresnospike")) {
    startmat=featfulldat_all
  } else {
    startmat=fulldat_all
  }
  newstart=scale(startmat[,apply(startmat,2,var)>0])
  numclust=which.max(allclustgap[[nameval]][[1]][5:20,3])+4 
  print(c(nameval,numclust))
  if (nameval=="glif1_spike_shape") {numclust=10}
  if (nameval=="glif2_spike_shape") {numclust=10}
  if (nameval=="glif3_spike_shape") {numclust=16}
  if (nameval=="glif4_spike_shape") {numclust=15}

  clustout=apclusterK(function (x){x=1-cor(t(x),method="pearson")},newstart[,keepcols],K=numclust,seed=0,prc=0,verbose=T)
  #clustout=apcluster(function (x){x=1-cor(t(x),method="spearman")},newstart[,keepcols],p=-4.18)
  allclusts=rep(0,nrow(startmat))
  names(allclusts)=rownames(startmat)
  print(c(nameval,numclust,length(clustout)))
  for (ii in 1:length(clustout)) {
    allclusts[clustout[[ii]]]=ii
  }
  temptab=table(allclusts[intersect(names(allclusts),rownames(metadata))],metadata$cre[match(intersect(names(allclusts),rownames(metadata)),rownames(metadata))])
  colnames(temptab)=sapply(strsplit(colnames(temptab),"-"), `[`, 1)
  temptab=cbind(temptab,paste("Cluster ",rev(1:nrow(temptab)),sep=''))
  write.csv(temptab,file=paste0("affprop_cluster_composition_log_20170727_",nameval,".csv"))
  temptab=allclusts
  write.csv(temptab,file=paste0("affprop_cluster_ids_log_20170727_",nameval,".csv"))
  outtab2=read.csv(paste0("affprop_cluster_composition_log_20170727_",nameval,".csv"),as.is=T,row.names=1,check.names=F)
  outtab2=outtab2[,-ncol(outtab2)]
  colnames(outtab2)[grep("Scnn1a",colnames(outtab2))]=c("Scnn1a-Tg2","Scnn1a-Tg3")
  colnames(outtab2)[grep("Nkx2",colnames(outtab2))]="Nkx2-1"
  outtab2=outtab2[,c("Htr3a","Ndnf","Vip","Sst","Pvalb","Nkx2-1","Chat","Chrna2","Cux2","Nr5a1","Scnn1a-Tg2","Scnn1a-Tg3","Rorb","Rbp4","Ntsr1","Ctgf")]
  xvals=matrix(rep(1:ncol(outtab2),each=nrow(outtab2)),nrow=nrow(outtab2))
  yvals=matrix(rep(1:nrow(outtab2),ncol(outtab2)),nrow=nrow(outtab2))
  crecols=read.csv("../cre_colors.csv",as.is=T,header=F)
  newcols=rgb(crecols[,2:4],maxColorValue = 255)
  names(newcols)=crecols[,5]
  basecols=newcols[colnames(outtab2)]
  colvals=matrix(basecols[rep(1:ncol(outtab2),each=nrow(outtab2))],nrow=nrow(outtab2))
  outtab2=100*sweep(as.matrix(outtab2),2,colSums(outtab2),"/")
  dfx = data.frame(x=c(xvals), y=c(yvals), sizeval=sqrt(c(as.matrix(outtab2))),colsplot=c(colvals))
  dfx = dfx[dfx$sizeval>0,]
  pdf(paste0("sp_affprop_cluster_diagram_log_20170727_",nameval,".pdf"),useDingbats=F,width=12,height=10)
  par(fig=c(0.3,1,0,1), new=FALSE)
  plot(c(1,ncol(outtab2)),c(1,nrow(outtab2)),pch='',xlab='',ylab='',xaxt='n',yaxt='n')
  abline(h=1:nrow(outtab2),v=1:ncol(outtab2),col='grey')
  with(dfx, symbols(x=x, y=y, circles=sizeval, inches=1/4, ann=F, bg=as.character(colsplot), fg="black", xlab=colnames(outtab),add=T,xlim=c(1,ncol(outtab2)),ylim=c(1,nrow(outtab2)),xaxt='n',yaxt='n'))
  axis(1, at=1:ncol(outtab2),labels=colnames(outtab2),las=2,cex.axis=0.9)
  axis(2, at=1:nrow(outtab2),label=paste("Cluster ",rev(1:nrow(outtab2)),sep=''),las=2)
  dev.off()
}
