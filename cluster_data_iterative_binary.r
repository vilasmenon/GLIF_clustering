require(stats)
require(e1071)
require(gtools)
require(gplots)
require(ape)
setwd("C:/Users/menonv/OneDrive - Howard Hughes Medical Institute/Misc/ct/20170727/")
dat=read.table("GLIF_param_7_27_17.csv",sep="\t",as.is=T,row.names=1,check.names=F,header=T)
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
  }
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



###separate into two clusters and check for differences###
cluster_into_two=function(fulldat,startseed,meth='ward.D') {
  fulldat=scale(fulldat[,apply(fulldat,2,var)>0])
  hc=hclust(as.dist(1-cor(t(fulldat),method="pearson")),method=meth)
  clustids=cutree(hc,2)
  outlist=list()
  ###assess predictability using SVM prediction###
   fraction_incorrect=c()
   inds1=which(clustids==1)
   inds2=which(clustids==2)
   if (length(inds1)>5 & length(inds2)>5) {
     sampfrac1=round(0.5*length(inds1))
     sampfrac2=round(0.5*length(inds2))
     for (tt in 1:100) {
       set.seed(tt+startseed)
       sampvec=c(sample(inds1,sampfrac1),sample(inds2,sampfrac2))
       setcols=which(apply(fulldat[sampvec,],2,var)>0)
       svmpred=predict(svm(x=fulldat[sampvec,setcols],y=clustids[sampvec],type="C-classification"),fulldat[-sampvec,setcols])
       conf=table(svmpred,clustids[-sampvec])
       fraction_incorrect=c(fraction_incorrect,(conf[2,1]+conf[1,2])/sum(conf))
     }
   } else {
     fraction_incorrect=c(1,1)
     fraction_incorrect_rand=c(1,1)
   }
   outlist[['fraction_incorrect']]=fraction_incorrect
   outlist[['clustids']]=clustids
  return(outlist)
}



#####Run clustering with subsets of parameters###
methall='ward.D'
recursive_clustering=function(keepcols,fulldat_all,fraclim=0.2,splitlim=50,startseed,outlist) {
  clustmat=fulldat_all[,keepcols]
  tempout=cluster_into_two(clustmat,meth=methall,startseed)
  #splitval=length(which(tempout$fraction_incorrect<fraclim))
  #print(c(nrow(fulldat_all),splitval))
  #if (splitval>splitlim) {
  if (!is.na(tempout$fraction_incorrect)) {
    if (max(tempout$fraction_incorrect,na.rm=T)<=fraclim) {
      outlist$clustnames[names(tempout$clustids)]=paste(outlist$clustnames[names(tempout$clustids)],tempout$clustids,sep="_")
      outlist$fracmat=rbind(outlist$fracmat,tempout$fraction_incorrect)
      for (ii in 1:2) {
        if (length(which(tempout$clustids==ii))>=10) {
        outlist=recursive_clustering(keepcols,fulldat_all[names(tempout$clustids)[tempout$clustids==ii],],fraclim=fraclim,splitlim=splitlim,startseed+ii,outlist)
        }
      }
    }
  }
  return(outlist)
}

fraclimval=0.2
for (nameval in rev(c("features","featuresnospike","allparam","glif4","glif3","glif2","glif1","lif"))) {
#for (nameval in c("featuresnospike")) {
  if (nameval=="lif") {keepcols=c(1,3,4,5)}
  if (nameval=="glif1") {keepcols=c(1,3,4,5,8)}
  if (nameval=="glif2") {keepcols=c(1,3,4,5,8,9,10)}
  if (nameval=="glif3") {keepcols=c(2,3,4,5,6,7,8)}
  if (nameval=="glif4") {keepcols=c(2,3,4,5,6,7,8,9,10)}
  if (nameval=="allparam") {keepcols=c(1,2,3,4,6,7,8,9,10,11,12)}
  if (nameval=="glif1_lifopt_th") {keepcols=c(1,3,4,8,11)}
  if (nameval=="glif3_lifopt_th") {keepcols=c(2,3,4,6,7,8,12)}
  if (nameval=="features") {keepcols=1:ncol(featfulldat_all)}
  if (nameval=="featuresnospike") {keepcols=c(1,2,3,4,5,8,11,12,13,14)}

  if (nameval %in% c("features","featuresnospike")) {
    startmat=featfulldat_all
  } else {
    startmat=fulldat_all
  }
  print(colnames(startmat)[keepcols])
  startnames=rep("1",nrow(startmat))
  names(startnames)=rownames(startmat)
  outlist=list()
  outlist$clustnames=startnames
  outlist$fracmat=c()
  allclusts=recursive_clustering(keepcols,startmat,fraclim=fraclimval,splitlim=splitlimval,startseed=1,outlist=outlist)
  print(c(nameval,table(allclusts$clustnames)))
  temptab=table(allclusts$clustnames[intersect(names(allclusts$clustnames),rownames(metadata))],metadata$cre[match(intersect(names(allclusts$clustnames),rownames(metadata)),rownames(metadata))])
  colnames(temptab)=sapply(strsplit(colnames(temptab),"-"), `[`, 1)
  temptab=cbind(temptab,paste("Cluster ",rev(1:nrow(temptab)),sep=''))
  write.csv(temptab,file=paste0("cluster_composition_log_20170727_",nameval,".csv"))
  temptab=allclusts$clustnames
  write.csv(temptab,file=paste0("cluster_ids_log_20170727_",nameval,".csv"))
  
  outtab2=read.csv(paste0("cluster_composition_log_20170727_",nameval,".csv"),as.is=T,row.names=1,check.names=F)
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
  pdf(paste0("cluster_diagram_log_20170727_",nameval,".pdf"),useDingbats=F,width=12,height=10)
  par(fig=c(0.3,1,0,1), new=FALSE)
  plot(c(1,ncol(outtab2)),c(1,nrow(outtab2)),pch='',xlab='',ylab='',xaxt='n',yaxt='n')
  abline(h=1:nrow(outtab2),v=1:ncol(outtab2),col='grey')
  with(dfx, symbols(x=x, y=y, circles=sizeval, inches=1/4, ann=F, bg=as.character(colsplot), fg="black", xlab=colnames(outtab),add=T,xlim=c(1,ncol(outtab2)),ylim=c(1,nrow(outtab2)),xaxt='n',yaxt='n'))
  axis(1, at=1:ncol(outtab2),labels=colnames(outtab2),las=2,cex.axis=0.9)
  axis(2, at=1:nrow(outtab2),label=paste("Cluster ",rev(1:nrow(outtab2)),sep=''),las=2)
  if (nrow(outtab2)>1) {
    prevnam1=rownames(outtab2)[1]
    prevnam2=strsplit(rownames(outtab2)[1],"_")[[1]]
    textvec=paste0(rep("(",length(prevnam2)-1),collapse="")
    for (ii in 1:(nrow(outtab2)-1)) {
      currnam1=rownames(outtab2)[ii]
      currnam2=strsplit(currnam1,"_")[[1]]
      nextnam1=rownames(outtab2)[ii+1]
      nextnam2=strsplit(nextnam1,"_")[[1]]
      textvec=paste0(textvec,currnam1)
      if (currnam2[length(currnam2)]=="1") {
        textvec=paste0(textvec,",")
        if (length(currnam2)<length(nextnam2)) {
          textvec=paste0(textvec,paste0(rep("(",length(nextnam2)-length(currnam2)),collapse=''))
        }
      } else {
        closepar=length(currnam2)-which(currnam2!=nextnam2)[1]
        textvec=paste0(textvec,paste0(rep(")",closepar),collapse=''),",")
        if (nextnam2[length(nextnam2)]=="1") {
          openpar=length(nextnam2)-which(currnam2!=nextnam2)[1]
          textvec=paste0(textvec,paste0(rep("(",openpar),collapse=''))
        }
      }
    }
    textvec=paste0(textvec,nextnam1,paste0(rep(")",length(nextnam2)-1),collapse=""),";")
    dd2=read.tree(text=textvec)
    tipvec=dd2$tip.label
    dd2$tip.label=rep("",length(tipvec))
    par(fig=c(0,0.33,0,1), new=TRUE)
    plot(dd2)
    plotvals=apply(allclusts$fracmat,1,max)
    nodelabels(100-round(100*plotvals),cex=1.4)
  }
  dev.off()
}

