library(reshape2)
library(ggplot2)
library(readxl)
library(dplyr)
library(adegenet)
library(stringr)
library(data.table)

#lob_mac_bic<-read.delim(file="lob_mac_bic.perl.txt",na.strings = "NA",sep="\t")
##lob_mac_bic<-subset(lob_mac_bic, pos/10==round(pos/10,0))
##lob_mac_bic_copy<-lob_mac_bic
##for(i in 1:length(lob_mac_bic$p1)){
##  for(j in 1:526){
##    if(!is.na(lob_mac_bic[i,j+4])){
##      if(lob_mac_bic[i,j+4]==lob_mac_bic[i,3]){
##        lob_mac_bic_copy[i,j+4]<-1
##      }else{
##        if(lob_mac_bic[i,j+4]==lob_mac_bic[i,4]){
##         lob_mac_bic_copy[i,j+4]<-0
##        }else{
##          lob_mac_bic_copy[i,j+4]<-NA
##        }
##      }
##    }
##  }
##}
##lob_mac_bic_long<-melt(data = lob_mac_bic_copy,id.vars = c("scaffold","pos","p1","p2"))

#samp_meta<-read.delim("/hpc/group/manoslab/QMACROCARPA_GENOMES/26_MAC_PAINTING/SampSpUpdate.txt")
##update samp_meta to include for both copies
##update to loop through chromosomes to minimize memory.
#csomes<-c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12")
#j<-1
#for(i in csomes){
#	print(paste("beginning Chromosome",i))
#	lob_mac_ss<-subset(lob_mac_bic,scaffold==i)
#	lob_mac_bic_long<-reshape2::melt(data = lob_mac_ss,id.vars = c("scaffold","pos","p1","p2"))
#	lob_mac_bic_long$roundPos<-round(lob_mac_bic_long$pos/100000,0)*100000
#	lob_mac_bic_long<-merge(lob_mac_bic_long,samp_meta,by.x="variable",by.y="Bam",all.x = TRUE,all.y = FALSE)
#	print("merge complete")
#	lob_mac_bic_long$count<-1
#	lob_mac_bic_long$samp_chr_pos<-paste(lob_mac_bic_long$variable,lob_mac_bic_long$scaffold,lob_mac_bic_long$roundPos,sep="_")
#	lob_mac_bic_long_agg<-aggregate(as.numeric(lob_mac_bic_long$value),by=list(lob_mac_bic_long$samp_chr_pos,lob_mac_bic_long$collectionNumber,lob_mac_bic_long$roundPos,lob_mac_bic_long$latitude.orig,lob_mac_bic_long$Species_ept,lob_mac_bic_long$scaffold),function(x) mean(x, na.rm=TRUE))
#	lob_mac_bic_long_count<-aggregate(as.numeric(lob_mac_bic_long$count),by=list(lob_mac_bic_long$samp_chr_pos),FUN=sum)
#	lob_mac_bic_long_agg_count<-merge(lob_mac_bic_long_agg, lob_mac_bic_long_count,by="Group.1")
#	if(j==1){
#		write.table(lob_mac_bic_long_agg_count, "lob_mac_bic_long_agg.txt",sep="\t",quote=F,row.names=F,col.names=TRUE,append=TRUE)
#	}else{
#		write.table(lob_mac_bic_long_agg_count, "lob_mac_bic_long_agg.txt",sep="\t",quote=F,row.names=F,col.names=FALSE,append=TRUE)
#	}
#	j<-j+1
#	print(paste("finished Chromosome",i))
#}


lob_mac_bic_long_agg_count<-read.delim("lob_mac_bic_long_agg.txt",sep="\t")

windows<-aggregate(cbind(x.x,x.y)~Group.3 + Group.6 + Group.5,data=lob_mac_bic_long_agg_count,function(x) mean(x, na.rm=TRUE))
corLatIntro<-data.frame(Chr=character(),window=numeric(),species=character(),nsites=integer(),corgeneVsp=numeric(),r2=numeric(),mean=numeric(),median=numeric(),meanPVal=numeric(),pless.05=numeric(),sd=character())

for(j in 1:length(windows$Group.3)){
	tmpdata<-subset(lob_mac_bic_long_agg_count,Group.3==windows$Group.3[j]&Group.5=="macrocarpa"&Group.6==windows$Group.6[j])
	tmpdata<-tmpdata[rowSums(is.na(tmpdata[,7:8]))==0,]
	ncount=windows$x.y[j]
	meanval=mean(tmpdata$x.x,na.rm=TRUE)
	nSamps=length(tmpdata$x.x)
	if(!is.na(ncount)&ncount>20){
		tmpdata$pvals<-NA
		tmpdata$pvals[tmpdata$x.x<=meanval]<-pbinom(round(tmpdata$x.x[tmpdata$x.x<=meanval]*ncount,0),size=ncount, prob = meanval,lower.tail = T)
		tmpdata$pvals[tmpdata$x.x>meanval]<-pbinom(round(tmpdata$x.x[tmpdata$x.x>meanval]*ncount,0),size=ncount, prob = meanval,lower.tail = F)
		#tmpdata$pvals<-pbinom(tmpdata$x.x*ncount,size=ncount, prob = meanval,lower.tail = FALSE)
		pval<-mean(tmpdata$pvals)
		countpvals<-sum(tmpdata$pvals<.025/(378*7663))
		tmpdata$latCorrected<-(tmpdata$Group.4-29)/20
		model<-lm(x.x~latCorrected,tmpdata)
		modelsum<-summary(model)
		meanSp<-mean(tmpdata$x.x)
		medianSp<-median(tmpdata$x.x)
		sdSp<-sd(tmpdata$x.x)
		corLatIntro<-rbind(corLatIntro,c(Chr=windows$Group.6[j],window=windows$Group.3[j],species="lob0_bic1",nsites=ncount,corgeneVsp=model$coefficients[2],r2=modelsum$r.squared,mean=meanSp,median=medianSp,meanPVal=pval,pless.05=countpvals, sd=sdSp))
	}
}


write.table(corLatIntro, "lob_mac_bic_ItrgByLat_2.tsv",sep="\t",quote=F,row.names=F)


