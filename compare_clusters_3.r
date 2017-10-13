####generate confusion matrices####
require(mclust)
require(gtools)
require(gplots)
setwd("C:/Users/menonv/OneDrive - Howard Hughes Medical Institute/Misc/ct/20170727/")

###read in Cre lines###
cellmeta=read.csv("cell_metadata.csv",as.is=T)
cellclass=read.csv("cluster_assignment_simple.csv",as.is=T)
allclust_rnaseq=unique(cellclass$primary[cellclass$coretype=="Core"])
allcre_rnaseq=unique(cellmeta$cre)

crefrac=t(table(cellclass$primary[cellclass$coretype=="Core"],cellmeta$cre[match(cellclass$cellid[cellclass$coretype=="Core"],cellmeta$long_name)]))
crefrac=crefrac[-grep("Neg",rownames(crefrac)),-grep("Astro|Endo|Olig|OPC|Micro",colnames(crefrac))]
crefrac=sweep(crefrac,1,rowSums(crefrac),"/")
credistmat=outer(1:nrow(crefrac),1:nrow(crefrac),FUN=Vectorize(function(x,y) sqrt(sum(crefrac[x,crefrac[x,]>0]-crefrac[y,crefrac[x,]>0])^2)))
credistmat=credistmat+t(credistmat)
rownames(credistmat)=rownames(crefrac)
colnames(credistmat)=rownames(credistmat)


#distmat=as.matrix(dist((crefrac)))
#plot(hclust(as.dist(distmat),method="average"))
#heatmap.2(credistmat,scale='none',trace='none',col=bluered)

###function to generate cluster distances based on hierarchical tree
dists_from_tree=function(xvec) {
  uniqclust=unique(xvec)
  uniqclust=uniqclust[order(uniqclust)]
  nbranch=max(nchar(uniqclust)+1)/2
  idmat=matrix(0,nrow=length(uniqclust),ncol=nbranch)
  distmat=matrix(0,nrow=nrow(idmat),ncol=nrow(idmat))
  for (ii in 1:length(uniqclust)) {
    addvec=strsplit(uniqclust[ii],"_")[[1]]
    idmat[ii,1:length(addvec)]=addvec
  }
  for (ii in 1:(nrow(distmat)-1)) {
    for (jj in (ii+1):nrow(distmat)) {
      mismatchind=which(idmat[ii,]!=idmat[jj,])[1]
      longestinds=grep(paste0("^",paste(idmat[ii,1:(mismatchind-1)],collapse="_")),uniqclust)
      longestval=max(which(idmat[longestinds,]!="0",arr.ind=T)[,2])
      distmat[ii,jj]=longestval-mismatchind+1
    }
  }
  distmat=distmat+t(distmat)
  rownames(distmat)=uniqclust
  colnames(distmat)=uniqclust
  distmat=2*distmat/max(distmat)
  return(distmat)
}

dist_from_name=function(x,y) {
  str1=strsplit(x,"_")
  str2=strsplit(y,"_")
  commonpoint=which(str1==str2)[1]
  dval=max(length(str1),length(str2))-commonpoint+1
  
}

###function to calculate Variation of Information or adjusted rand index
calc_cluster_diff=function(xvec,yvec,functype=1,credistmat=c(),clustdistmat=c()) {
  if (functype==1) {
    totaltab=table(xvec,yvec)
    rowmat=sweep(totaltab,1,rowSums(totaltab),"/")  
    colmat=sweep(totaltab,2,colSums(totaltab),"/")
    summat=(totaltab*(log(rowmat)+log(colmat)))
    sumval=sum(summat[totaltab>0])/length(xvec)
    return(-sumval) 
  } else {
    if (functype==2) {
      return(adjustedRandIndex(xvec,yvec))
    } else {
      return(calc_rar(xvec,yvec,credistmat,clustdistmat))
    }
  }
}

###function to calculate score based on 100 random permutations
rand_cluster_diff=function(xvec,yvec,functype=1,credistmat=c(),clustdistmat=c()) {
  allvals=rep(0,100)
  for (ii in 1:100) {
    set.seed(ii)
    allvals[ii]=calc_cluster_diff(xvec,sample(yvec),functype,credistmat,clustdistmat)
  }
  return(allvals)
}

# gmat1=credistmat[unique(featclust[,2]),unique(featclust[,2])]
# set.seed(0);scrambvec=sample(1:nrow(gmat1))
# gmat2=gmat1
# rownames(gmat2)=rownames(gmat1)[scrambvec]
# colnames(gmat2)=rownames(gmat2)[scrambvec]
# gmat2=gmat2[rownames(gmat1),rownames(gmat1)]
# c1=featclust[,2]
# c2=featclust[,2]
###function to calculate Rank-Adjusted Rand Index (Pinto et al. 2005)
calc_rar=function(c1,c2,gmat1=c(),gmat2=c()) {
    nn=length(c1)
    ng1=length(unique(c1))
    ng2=length(unique(c2))
    if (length(gmat1)==0) {
      gmat1=2-diag(ng1)
      c1=match(c1,unique(c1))
    } else {
      c1=match(c1,rownames(gmat1))
      for (ii in 1:ncol(gmat1)) {
        uniqcol=unique(gmat1[,ii])
        uniqcol=uniqcol[order(uniqcol)]
        gmat1[,ii]=match(gmat1[,ii],uniqcol)
      }
    }
    if (length(gmat2)==0) {
      gmat2=2-diag(ng2)
      c2=match(c2,unique(c2))
    } else {
      c2=match(c2,rownames(gmat2))
      for (ii in 1:ncol(gmat2)) {
        uniqcol=unique(gmat2[,ii])
        uniqcol=uniqcol[order(uniqcol)]
        gmat2[,ii]=match(gmat2[,ii],uniqcol)
      }
    }
    bothclass=matrix(0,nrow=ng1,ncol=ng2)
    for (ii in 1:length(c1)) {
      bothclass[c1[ii],c2[ii]]=bothclass[c1[ii],c2[ii]]+1
    }
    
    dim1=max(gmat1)
    dim2=max(gmat2)
    rmm=matrix(0,nrow=dim1,ncol=dim2)
    mismatchrank1=c(gmat1[c1,c1])
    mismatchrank2=c(gmat2[c2,c2])
    for (ii in 1:length(mismatchrank1)) {
      rmm[mismatchrank1[ii],mismatchrank2[ii]]=rmm[mismatchrank1[ii],mismatchrank2[ii]]+1
    }
    rmm[1,1]=rmm[1,1]-nn
    nprod1=rowSums(bothclass) %*% t(rowSums(bothclass))
    nprod2=colSums(bothclass) %*% t(colSums(bothclass)) 
    rmmi=matrix(0,nrow=dim1,ncol=dim2)
    for (ii in 1:dim1) {
      for (jj in 1:dim2) {
        d1=sum(nprod1[gmat1==ii])
        d2=sum(nprod2[gmat2==jj])
        rmmi[ii,jj]=d1*d2/(nn^2)
      }
    }
    rmmi[1,1]=rmmi[1,1]-nn
    nr=nrow(rmm)
    nc=ncol(rmm)
    devmat=matrix(rep((0:(nr-1))/(nr-1),nc),ncol=nc)-matrix(rep((0:(nc-1))/(nc-1),each=nr),nrow=nr)
    sumdev=sum(abs(devmat)*rmm)
    indsumdev=sum(abs(devmat)*rmmi)
    mddi=indsumdev/sum(rmmi)
    mdd=sumdev/sum(rmm)
    rar=(mddi-mdd)/mddi
    return(rar)
}


keepcells=intersect(grep("Cre",cellmeta$cre),which(cellclass$coretype=="Core"))
transcript_voi=-calc_cluster_diff(xvec = cellclass$primary[keepcells],yvec=cellmeta$cre[keepcells],1)+mean(rand_cluster_diff(xvec = cellclass$primary[keepcells],yvec=cellmeta$cre[keepcells],1))
transcript_ari=calc_cluster_diff(xvec = cellclass$primary[keepcells],yvec=cellmeta$cre[keepcells],2)-mean(rand_cluster_diff(xvec = cellclass$primary[keepcells],yvec=cellmeta$cre[keepcells],2))



###all clusters to features###
all_voi=c()
all_ari=c()
mean_voi=c()
sd_voi=c()
mean_ari=c()
sd_ari=c()
featclust=read.csv("cluster_ids_log_20170727_features.csv",as.is=T)
for (nameval in c("lif","glif1","glif2","glif3","glif4","featuresnospike","glif1_spike_shape","glif2_spike_shape","glif3_spike_shape","glif4_spike_shape")) {
  glifclust=read.csv(paste0("cluster_ids_log_20170727_",nameval,".csv"),as.is=T)
  glifclust=glifclust[match(featclust[,1],glifclust[,1]),]
  all_voi=c(all_voi,calc_cluster_diff(featclust[,2],glifclust[,2],1))
  all_ari=c(all_ari,calc_cluster_diff(featclust[,2],glifclust[,2],2))
  rand_voi=rand_cluster_diff(featclust[,2],glifclust[,2],1)
  mean_voi=c(mean_voi,mean(rand_voi))
  sd_voi=c(sd_voi,sd(rand_voi))
  rand_ari=rand_cluster_diff(featclust[,2],glifclust[,2],2)
  mean_ari=c(mean_ari,mean(rand_ari))
  sd_ari=c(sd_ari,sd(rand_ari))
  temptab=table(glifclust[,2],featclust[,2])
  glifcomp=read.csv(paste0("cluster_composition_log_20170727_",nameval,".csv"),as.is=T)
  rowvec=glifcomp[match(rownames(temptab),glifcomp[,1]),ncol(glifcomp)]
  featcomp=read.csv("cluster_composition_log_20170727_features.csv",as.is=T)
  colvec=featcomp[match(colnames(temptab),featcomp[,1]),ncol(featcomp)]
  temptab=cbind(temptab,rowvec)
  temptab=rbind(temptab,c(colvec,''))
  write.csv(temptab,file=paste0("confusion_matrix_log_20170727_features_",nameval,".csv"))
  
  ###make plot###
  plottab=temptab
  rownames(plottab)=temptab[,ncol(temptab)]
  colnames(plottab)=temptab[nrow(temptab),]
  plottab=plottab[-nrow(plottab),-ncol(plottab)]
  plottab=t(apply(plottab,1,as.numeric))
  if (nameval=="lif") {
    rownames(plottab)=paste0("LIF ",temptab[-nrow(temptab),ncol(temptab)])
  } else {
    rownames(plottab)=paste0("GLIF ",temptab[-nrow(temptab),ncol(temptab)])
  }
  rownames(plottab)=gsub("Cluster ","C",rownames(plottab))
  colnames(plottab)=paste0("Feature ",temptab[nrow(temptab),-ncol(temptab)])
  colnames(plottab)=gsub("Cluster ","C",colnames(plottab))
  plottab=plottab[mixedorder(rownames(plottab)),mixedorder(colnames(plottab))]
  texttab=plottab
  texttab[texttab==0]=''
  pdf(paste0("confusion_matrix_log_20170727_features_",nameval,".pdf"),useDingbats=F)
  heatmap.2(plottab,col=colorRampPalette(c("white","orange","red")),scale='none',trace='none',cellnote=texttab,notecol="black",Rowv=F,Colv=F,cexRow=0.9,cexCol=0.9,key=F)
  dev.off()
}


ap_voi=c()
ap_ari=c()
ap_mean_voi=c()
ap_sd_voi=c()
ap_mean_ari=c()
ap_sd_ari=c()
featclust=read.csv("affprop_cluster_ids_log_20170727_features.csv",as.is=T)
for (nameval in c("lif","glif1","glif2","glif3","glif4","featuresnospike","glif1_spike_shape","glif2_spike_shape","glif3_spike_shape","glif4_spike_shape")) {
  glifclust=read.csv(paste0("affprop_cluster_ids_log_20170727_",nameval,".csv"),as.is=T)
  glifclust=glifclust[match(featclust[,1],glifclust[,1]),]
  ap_voi=c(ap_voi,calc_cluster_diff(featclust[,2],glifclust[,2],1))
  ap_ari=c(ap_ari,calc_cluster_diff(featclust[,2],glifclust[,2],2))
  rand_voi=rand_cluster_diff(featclust[,2],glifclust[,2],1)
  ap_mean_voi=c(ap_mean_voi,mean(rand_voi))
  ap_sd_voi=c(ap_sd_voi,sd(rand_voi))
  rand_ari=rand_cluster_diff(featclust[,2],glifclust[,2],2)
  ap_mean_ari=c(ap_mean_ari,mean(rand_ari))
  ap_sd_ari=c(ap_sd_ari,sd(rand_ari))
  temptab=table(glifclust[,2],featclust[,2])
  glifcomp=read.csv(paste0("affprop_cluster_composition_log_20170727_",nameval,".csv"),as.is=T)
  rowvec=glifcomp[match(rownames(temptab),glifcomp[,1]),ncol(glifcomp)]
  featcomp=read.csv("affprop_cluster_composition_log_20170727_features.csv",as.is=T)
  colvec=featcomp[match(colnames(temptab),featcomp[,1]),ncol(featcomp)]
  temptab=cbind(temptab,rowvec)
  temptab=rbind(temptab,c(colvec,''))
  write.csv(temptab,file=paste0("affprop_confusion_matrix_log_20170727_features_",nameval,".csv"))
  
  ###make plot###
  plottab=temptab
  rownames(plottab)=temptab[,ncol(temptab)]
  colnames(plottab)=temptab[nrow(temptab),]
  plottab=plottab[-nrow(plottab),-ncol(plottab)]
  plottab=t(apply(plottab,1,as.numeric))
  if (nameval=="lif") {
    rownames(plottab)=paste0("LIF ",temptab[-nrow(temptab),ncol(temptab)])
  } else {
    rownames(plottab)=paste0("GLIF ",temptab[-nrow(temptab),ncol(temptab)])
  }
  rownames(plottab)=gsub("Cluster ","C",rownames(plottab))
  colnames(plottab)=paste0("Feature ",temptab[nrow(temptab),-ncol(temptab)])
  colnames(plottab)=gsub("Cluster ","C",colnames(plottab))
  plottab=plottab[mixedorder(rownames(plottab)),mixedorder(colnames(plottab))]
  texttab=plottab
  texttab[texttab==0]=''
  pdf(paste0("affprop_confusion_matrix_log_20170727_features_",nameval,".pdf"),useDingbats=F)
  heatmap.2(plottab,col=colorRampPalette(c("white","orange","red")),scale='none',trace='none',cellnote=texttab,notecol="black",Rowv=F,Colv=F,cexRow=0.9,cexCol=0.9,key=F)
  dev.off()
}


###all clusters to Cre lines###
cre_voi=c()
cre_ari=c()
cre_mean_voi=c()
cre_sd_voi=c()
cre_mean_ari=c()
cre_sd_ari=c()
cre_rar=c()
cre_mean_rar=c()
cre_sd_rar=c()
cre_rar_tree=c()
cre_mean_rar_tree=c()
cre_sd_rar_tree=c()
featclust=read.table("features_7_27_17.csv",as.is=T,sep=",",header=T)[,1:2]
for (nameval in c("lif","glif1","glif2","glif3","glif4","featuresnospike","features","glif1_spike_shape","glif2_spike_shape","glif3_spike_shape","glif4_spike_shape")) {
  glifclust=read.csv(paste0("cluster_ids_log_20170727_",nameval,".csv"),as.is=T)
  glifclust=glifclust[match(featclust[,1],glifclust[,1]),]
  cre_voi=c(cre_voi,calc_cluster_diff(featclust[,2],glifclust[,2],1))
  cre_ari=c(cre_ari,calc_cluster_diff(featclust[,2],glifclust[,2],2))
  rand_voi=rand_cluster_diff(featclust[,2],glifclust[,2],1)
  rand_ari=rand_cluster_diff (featclust[,2],glifclust[,2],2)
  cre_mean_voi=c(cre_mean_voi,mean(rand_voi))
  cre_mean_ari=c(cre_mean_ari,mean(rand_ari))
  cre_sd_voi=c(cre_sd_voi,sd(rand_voi))
  cre_sd_ari=c(cre_sd_ari,sd(rand_ari))
  cre_rar_tree=c(cre_rar_tree,calc_cluster_diff(featclust[,2],glifclust[,2],3,credistmat[unique(featclust[,2]),unique(featclust[,2])],dists_from_tree(glifclust[,2])))
  rand_rar=rand_cluster_diff(featclust[,2],glifclust[,2],3,credistmat[unique(featclust[,2]),unique(featclust[,2])],dists_from_tree(glifclust[,2]))
  cre_mean_rar_tree=c(cre_mean_rar_tree,mean(rand_rar))
  cre_sd_rar_tree=c(cre_sd_rar_tree,sd(rand_rar))
  cre_rar=c(cre_rar,calc_cluster_diff(featclust[,2],glifclust[,2],3,credistmat[unique(featclust[,2]),unique(featclust[,2])]))
  rand_rar=rand_cluster_diff(featclust[,2],glifclust[,2],3,credistmat[unique(featclust[,2]),unique(featclust[,2])])
  cre_mean_rar=c(cre_mean_rar,mean(rand_rar))
  cre_sd_rar=c(cre_sd_rar,sd(rand_rar))
}

cre_ap_voi=c()
cre_ap_ari=c()
cre_ap_mean_voi=c()
cre_ap_sd_voi=c()
cre_ap_mean_ari=c()
cre_ap_sd_ari=c()
cre_ap_rar=c()
cre_ap_mean_rar=c()
cre_ap_sd_rar=c()
featclust=read.table("features_7_27_17.csv",as.is=T,sep=",",header=T)[,1:2]
for (nameval in c("lif","glif1","glif2","glif3","glif4","featuresnospike","features","glif1_spike_shape","glif2_spike_shape","glif3_spike_shape","glif4_spike_shape")) {
  glifclust=read.csv(paste0("affprop_cluster_ids_log_20170727_",nameval,".csv"),as.is=T)
  glifclust=glifclust[match(featclust[,1],glifclust[,1]),]
  cre_ap_voi=c(cre_ap_voi,calc_cluster_diff(featclust[,2],glifclust[,2],1))
  cre_ap_ari=c(cre_ap_ari,calc_cluster_diff(featclust[,2],glifclust[,2],2))
  rand_voi=rand_cluster_diff(featclust[,2],glifclust[,2],1)
  rand_ari=rand_cluster_diff(featclust[,2],glifclust[,2],2)
  cre_ap_mean_voi=c(cre_ap_mean_voi,mean(rand_voi))
  cre_ap_sd_voi=c(cre_ap_sd_voi,sd(rand_voi))
  cre_ap_mean_ari=c(cre_ap_mean_ari,mean(rand_ari))
  cre_ap_sd_ari=c(cre_ap_sd_ari,sd(rand_ari))
  cre_ap_rar=c(cre_ap_rar,calc_cluster_diff(featclust[,2],glifclust[,2],3,credistmat[unique(featclust[,2]),unique(featclust[,2])]))
  rand_rar=rand_cluster_diff(featclust[,2],glifclust[,2],3,credistmat[unique(featclust[,2]),unique(featclust[,2])])
  cre_ap_mean_rar=c(cre_ap_mean_rar,mean(rand_rar))
  cre_ap_sd_rar=c(cre_ap_sd_rar,sd(rand_rar))
}

###random features###
cre_featrand_voi=c()
cre_featrand_ari=c()
cre_featrand_ap_voi=c()
cre_featrand_ap_ari=c()
cre_featrand_mean_voi=c()
cre_featrand_sd_voi=c()
cre_featrand_mean_ari=c()
cre_featrand_sd_ari=c()
cre_featrand_ap_mean_voi=c()
cre_featrand_ap_sd_voi=c()
cre_featrand_ap_mean_ari=c()
cre_featrand_ap_sd_ari=c()
featclust=read.table("features_7_27_17.csv",as.is=T,sep=",",header=T)[,1:2]
for (nameval in c("features_rand_g1","features_rand_g2","features_rand_g3","features_rand_g4")) {
  glifclust=read.csv(paste0("cluster_ids_log_20170727_",nameval,".csv"),as.is=T)
  glifclust=glifclust[match(featclust[,1],glifclust[,1]),]
  cre_featrand_voi=c(cre_featrand_voi,calc_cluster_diff(featclust[,2],glifclust[,2],1))
  cre_featrand_ari=c(cre_featrand_ari,calc_cluster_diff(featclust[,2],glifclust[,2],2))
  rand_voi=rand_cluster_diff(featclust[,2],glifclust[,2],1)
  rand_ari=rand_cluster_diff(featclust[,2],glifclust[,2],2)
  cre_featrand_mean_voi=c(cre_featrand_mean_voi,mean(rand_voi))
  cre_featrand_sd_voi=c(cre_featrand_sd_voi,sd(rand_voi))
  cre_featrand_mean_ari=c(cre_featrand_mean_ari,mean(rand_ari))
  cre_featrand_sd_ari=c(cre_featrand_sd_ari,sd(rand_ari))
  glifclust=read.csv(paste0("affprop_cluster_ids_log_20170727_",nameval,".csv"),as.is=T)
  glifclust=glifclust[match(featclust[,1],glifclust[,1]),]
  cre_featrand_ap_voi=c(cre_featrand_ap_voi,calc_cluster_diff(featclust[,2],glifclust[,2],1))
  cre_featrand_ap_ari=c(cre_featrand_ap_ari,calc_cluster_diff(featclust[,2],glifclust[,2],2))
  rand_voi=rand_cluster_diff(featclust[,2],glifclust[,2],1)
  rand_ari=rand_cluster_diff(featclust[,2],glifclust[,2],2)
  cre_featrand_ap_mean_voi=c(cre_featrand_ap_mean_voi,mean(rand_voi))
  cre_featrand_ap_sd_voi=c(cre_featrand_ap_sd_voi,sd(rand_voi))
  cre_featrand_ap_mean_ari=c(cre_featrand_ap_mean_ari,mean(rand_ari))
  cre_featrand_ap_sd_ari=c(cre_featrand_ap_sd_ari,sd(rand_ari))
}





pdf("cross_clustering_all_binary_splitting.pdf",useDingbats=F)
par(mar = c(5,5,2,5))
plot(1:length(cre_voi[2:7]),cre_mean_voi[2:7]-cre_voi[2:7],type="l",col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_mean_voi[2:7]-cre_voi[2:7])))
axis(side=1,at=1:length(cre_voi[2:7]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape", "Features"),las=2)
par(new = T)
plot(1:length(cre_ari[2:7]), cre_ari[2:7]-cre_mean_ari[2:7], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(cre_ari[2:7]-cre_mean_ari[2:7])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))

par(mar = c(5,5,2,5))
plot(1:length(cre_voi[2:7]),cre_mean_voi[2:7]-cre_voi[2:7],type="l",col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_mean_voi[2:7]-cre_voi[2:7])))
lines(1:length(cre_featrand_mean_voi),cre_featrand_mean_voi-cre_featrand_voi,type="l",col="grey")
axis(side=1,at=1:length(cre_voi[2:7]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape", "Features"),las=2)
par(new = T)
plot(1:length(cre_ari[2:7]), cre_ari[2:7]-cre_mean_ari[2:7], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(cre_ari[2:7]-cre_mean_ari[2:7])))
lines(1:length(cre_featrand_ari), cre_featrand_ari-cre_featrand_mean_ari, type="l", col="pink")
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index","Adjusted VOI,\nRandom Feature set","Adjusted Rand Index,\nRandom Feature set"),fill=c("black","red","grey","pink"))

par(mar = c(5,5,2,5))
plot(1:length(cre_voi[2:7]),cre_mean_voi[2:7]-cre_voi[2:7],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_mean_voi[2:7]-cre_voi[2:7])))
axis(side=1,at=1:length(cre_voi[2:7]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape", "Features"),las=2)
par(new = T)
rangeval=range(c(cre_ari[2:7]-cre_mean_ari[2:7],cre_rar_tree[2:7]-cre_mean_rar_tree[2:7],cre_rar[2:7]-cre_mean_rar[2:7]))
plot(1:length(cre_ari[2:7]), cre_ari[2:7]-cre_mean_ari[2:7], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(rangeval)))
lines(1:length(cre_rar_tree[2:7]), cre_rar_tree[2:7]-cre_mean_rar_tree[2:7], type="l",col="magenta")
#lines(1:length(cre_rar[2:7]), cre_rar[2:7]-cre_mean_rar[2:7], type="l",col="pink")
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Rand Index value',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index","Rank-Adjusted Rand Index"),fill=c("black","red","magenta"))


par(mar = c(5,5,2,5))
plot(1:length(all_voi[2:6]),mean_voi[2:6]-all_voi[2:6],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to clustering by features",ylim=c(0,max(mean_voi[2:6]-all_voi[2:6])))
axis(side=1,at=1:length(all_voi[2:6]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape"),las=2)
par(new = T)
plot(1:length(all_ari[2:6]), all_ari[2:6]-mean_ari[2:6], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(all_ari[2:6],mean_ari[2:6])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))
dev.off()


pdf("cross_clustering_all_ap.pdf",useDingbats=F)
par(mar = c(5,5,2,5))
plot(1:length(cre_ap_voi[2:7]),cre_ap_mean_voi[2:7]-cre_ap_voi[2:7],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_ap_mean_voi[2:7]-cre_ap_voi[2:7])))
axis(side=1,at=1:length(cre_ap_voi[2:7]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape", "Features"),las=2)
par(new = T)
plot(1:length(cre_ap_ari[2:7]), cre_ap_ari[2:7]-cre_ap_mean_ari[2:7], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(cre_ap_ari[2:7]-cre_ap_mean_ari[2:7])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))

par(mar = c(5,5,2,5))
plot(1:length(cre_ap_voi[2:7]),cre_ap_mean_voi[2:7]-cre_ap_voi[2:7],type="l",col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_ap_mean_voi[2:7]-cre_ap_voi[2:7])))
lines(1:length(cre_featrand_ap_mean_voi),cre_featrand_ap_mean_voi-cre_featrand_ap_voi,type="l",col="grey")
axis(side=1,at=1:length(cre_ap_voi[2:7]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape", "Features"),las=2)
par(new = T)
plot(1:length(cre_ap_ari[2:7]), cre_ap_ari[2:7]-cre_ap_mean_ari[2:7], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(cre_ap_ari[2:7]-cre_ap_mean_ari[2:7])))
lines(1:length(cre_featrand_ap_ari), cre_featrand_ap_ari-cre_featrand_ap_mean_ari, type="l", col="pink")
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("bottomright",c("Adjusted VOI","Adjusted Rand Index","Adjusted VOI,\nRandom Feature set","Adjusted Rand Index,\nRandom Feature set"),fill=c("black","red","grey","pink"))


par(mar = c(5,5,2,5))
plot(1:length(cre_ap_voi[2:7]),cre_ap_mean_voi[2:7]-cre_ap_voi[2:7],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_ap_mean_voi[2:7]-cre_ap_voi[2:7])))
axis(side=1,at=1:length(cre_ap_voi[2:7]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape", "Features"),las=2)
par(new = T)
rangeval=range(c(cre_ap_ari[2:7]-cre_ap_mean_ari[2:7],cre_ap_rar[2:7]-cre_ap_mean_rar[2:7]))
plot(1:length(cre_ap_ari[2:7]), cre_ap_ari[2:7]-cre_ap_mean_ari[2:7], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(rangeval)))
lines(1:length(cre_ap_rar[2:7]), cre_ap_rar[2:7]-cre_ap_mean_rar[2:7], type="l",col="magenta")
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Rand Index value',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index","Rank-Adjusted Rand Index"),fill=c("black","red","magenta"))

par(mar = c(5,5,2,5))
plot(1:length(ap_voi[2:6]),ap_mean_voi[2:6]-ap_voi[2:6],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to clustering by features",ylim=c(0,max(ap_mean_voi[2:6]-ap_voi[2:6])))
axis(side=1,at=1:length(ap_voi[2:6]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape"),las=2)
par(new = T)
plot(1:length(ap_ari[2:6]), ap_ari[2:6]-ap_mean_ari[2:6], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(ap_ari[2:6]-ap_mean_ari[2:6])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))
dev.off()




pdf("cross_clustering_all_spike_shape_binary_ap.pdf",useDingbats=F)
par(mar = c(5,5,2,5))
plot(1:length(cre_voi[2:11]),cre_mean_voi[2:11]-cre_voi[2:11],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_mean_voi[2:11]-cre_voi[2:11])))
axis(side=1,at=1:length(cre_voi[2:11]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape","Features","GLIF1+Spike Shape","GLIF2+Spike Shape","GLIF3+Spike Shape","GLIF4+Spike Shape"),las=2)
par(new = T)
plot(1:length(cre_ari[2:11]), cre_ari[2:11]-cre_mean_ari[2:11], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(cre_ari[2:11]-cre_mean_ari[2:11])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))

par(mar = c(5,5,2,5))
plot(1:length(cre_voi[2:11]),cre_mean_voi[2:11]-cre_voi[2:11],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_mean_voi[2:11]-cre_voi[2:11])))
axis(side=1,at=1:length(cre_voi[2:11]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape","Features","GLIF1+Spike Shape","GLIF2+Spike Shape","GLIF3+Spike Shape","GLIF4+Spike Shape"),las=2)
rangeval=range(c(cre_ari[2:11]-cre_mean_ari[2:11],cre_rar_tree[2:11]-cre_mean_rar_tree[2:11],cre_rar[2:11]-cre_mean_rar[2:11]))
par(new=T)
plot(1:length(cre_ari[2:11]), cre_ari[2:11]-cre_mean_ari[2:11], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(rangeval)))
lines(1:length(cre_rar_tree[2:11]), cre_rar_tree[2:11]-cre_mean_rar_tree[2:11], type="l",col="magenta")
axis(side = 4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Rand Index value',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index","Rank-Adjusted Rand Index"),fill=c("black","red","magenta"))

par(mar = c(5,5,2,5))
plot(1:length(all_voi[2:10]),mean_voi[2:10]-all_voi[2:10],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to clustering by features",ylim=c(0,max(mean_voi[2:10]-all_voi[2:10])))
axis(side=1,at=1:length(all_voi[2:10]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape","GLIF1+Spike Shape","GLIF2+Spike Shape","GLIF3+Spike Shape","GLIF4+Spike Shape"),las=2)
par(new = T)
plot(1:length(all_ari[2:10]), all_ari[2:10]-mean_ari[2:10], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(all_ari[2:10]-mean_ari[2:10])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))

par(mar = c(5,5,2,5))
plot(1:length(cre_ap_voi[2:11]),cre_ap_mean_voi[2:11]-cre_ap_voi[2:11],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_ap_mean_voi[2:11]-cre_ap_voi[2:11])))
axis(side=1,at=1:length(cre_ap_voi[2:11]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape","Features","GLIF1+Spike Shape","GLIF2+Spike Shape","GLIF3+Spike Shape","GLIF4+Spike Shape"),las=2)
par(new = T)
plot(1:length(cre_ap_ari[2:11]), cre_ap_ari[2:11]-cre_ap_mean_ari[2:11], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(cre_ap_ari[2:11]-cre_ap_mean_ari[2:11])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))

par(mar = c(5,5,2,5))
plot(1:length(cre_ap_voi[2:11]),cre_ap_mean_voi[2:11]-cre_ap_voi[2:11],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to Cre lines",ylim=c(0,max(cre_ap_mean_voi[2:11]-cre_ap_voi[2:11])))
axis(side=1,at=1:length(cre_ap_voi[2:11]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no\nspike-shape","Features","GLIF1+Spike Shape","GLIF2+Spike Shape","GLIF3+Spike Shape","GLIF4+Spike Shape"),las=2)
par(new = T)
rangeval=range(c(cre_ap_ari[2:11]-cre_ap_mean_ari[2:11],cre_ap_rar[2:11]-cre_ap_mean_rar[2:11]))
plot(1:length(cre_ap_ari[2:11]), cre_ap_ari[2:11]-cre_ap_mean_ari[2:11], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(rangeval)))
lines(1:length(cre_ap_rar[2:11]), cre_ap_rar[2:11]-cre_ap_mean_rar[2:11], type="l",col="magenta")
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index","Rank-Adjusted Rand Index"),fill=c("black","red","magenta"))

par(mar = c(5,5,2,5))
plot(1:length(ap_voi[2:10]),ap_mean_voi[2:10]-ap_voi[2:10],type="l", col="black",ylab="Adjusted VOI score",xaxt='n',xlab='',main="Comparison to clustering by features",ylim=c(0,max(ap_mean_voi[2:10]-ap_voi[2:10])))
axis(side=1,at=1:length(ap_voi[2:10]),labels=c("GLIF1","GLIF2","GLIF3","GLIF4","Features, no spike-shape","GLIF1+Spike Shape","GLIF2+Spike Shape","GLIF3+Spike Shape","GLIF4+Spike Shape"),las=2)
par(new = T)
plot(1:length(ap_ari[2:10]), ap_ari[2:10]-ap_mean_ari[2:10], type="l", col="red", axes=F, xlab=NA, ylab=NA,ylim=c(0,max(ap_ari[2:10]-ap_mean_ari[2:10])))
axis(side=4,labels=F)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "red", line = 1)
mtext(side = 4, line = 3, 'Adjusted Rand Index',col='red')
legend("topleft",c("Adjusted VOI","Adjusted Rand Index"),fill=c("black","red"))
dev.off()

