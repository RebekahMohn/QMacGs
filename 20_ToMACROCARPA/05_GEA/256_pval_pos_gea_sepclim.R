library(ggplot2)
library(readxl)

ChrStr<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="Sheet3")
group.colors<-c(lyr="#000000",bic="#88CCEE",lob="#116644",mue="#AA5599",ste="#E69F00",alb="#332288",none="white",CenArray="gray20",InterArray="gray80", IntraSpInv="gray70",InterSpInv="gray35")
ChrStr$csome<-ChrStr$scaffold
pval_genome<-data.frame()

genomesize<-765592681



for(j in c(1:6)){
  pval_genome_env<-data.frame()
  for(i in c(1:12)){
    pval<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/C",i,"e",j,"_lfmmpval.txt",sep=""))
    position<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/MAC99_MACREF_HQSNP_C",i,"_maf01.map",sep=""),header=FALSE)
    pvalpos<-cbind(pval,position)
    write.table(pvalpos,file=paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/C",i,"e",j,"_lfmmpvalpos.txt",sep=""),quote=FALSE,row.names=FALSE)
    if(i<10){
      pvalpos$csome<-paste("Chr0",i,sep="")
    }else{
      pvalpos$csome<-paste("Chr",i,sep="")
      }
    pvalpos10<-pvalpos[seq(1,to=nrow(pvalpos),by=10),]
    pval_genome_env<-rbind(pval_genome_env,pvalpos10)
  }
 
   #pval_genome$pBH<-p.adjust(pval_genome$x,method=c("BH"),n=length(pval_genome$x))
  #pval_genome$pbonferroni<-p.adjust(pval_genome$x,method=c("bonferroni"),n=length(pval_genome$x))
  colnames(pval_genome_env)<-c("P","Count","uk","uk2","BP","CHR")
  pval_genome_env$SNP<-str_split_i(row.names(pval_genome_env)," ",2)
  #write.table(pval_genome_env,paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e",j,"_lfmmpval.txt",sep=""),quote = FALSE,row.names = FALSE) #these are output bellow
 # write.table(subset(pval_genome,pbonferroni<=0.05),paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e",j,"_lfmmpval_cor.txt",sep=""),sep="\t")
}



for(j in c(1:6)){
  pval_genome_env<-data.frame()
  for(i in c(1:12)){
    pvalpos<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/C",i,"e",j,"_lfmmpval_clump.clumped",sep=""),sep = " ")
    if(i<10){
      pvalpos$csome<-paste("Chr0",i,sep="")
    }else{
      pvalpos$csome<-paste("Chr",i,sep="")
    }
    pval_genome_env<-rbind(pval_genome_env,pvalpos)
  }
  ##pval_genome_env$pBH<-p.adjust(pval_genome_env$P,method=c("BH"),n=(6*genomesize/10000))
  ##pval_genome_env$pbonferroni<-p.adjust(pval_genome_env$P,method=c("bonferroni"),n=6*genomesize/10000)
  
  pval_genome_env$pBH<-p.adjust(pval_genome_env$P,method=c("BH"),n=(6*genomesize/10000)) #I believe this should be uncommented
  pval_genome_env$pbonferroni<-p.adjust(pval_genome_env$P,method=c("bonferroni"),n=6*genomesize/10000) #I believe this are uncommented
  write.table(pval_genome_env,paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e",j,"_lfmmpval.txt",sep=""),quote = FALSE,row.names = FALSE)
  write.table(subset(pval_genome_env,pbonferroni<0.05),paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e",j,"_lfmmpval_cor.txt",sep=""),quote = FALSE,row.names = FALSE)
}




gea1<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e1_lfmmpval.txt",sep=" ")
gea2<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e2_lfmmpval.txt",sep=" ")
gea3<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e3_lfmmpval.txt",sep=" ")
gea4<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e4_lfmmpval.txt",sep=" ")
gea5<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e5_lfmmpval.txt",sep=" ")
gea6<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_e6_lfmmpval.txt",sep=" ")

print(ggplot()+
        #geom_rect(data=ChrStr,mapping=aes(ymin=0,ymax=15,xmin=Start,xmax=Stop, fill=Type),alpha=0.5)+
        geom_point(data=gea6,mapping=aes(y=-log10(pbonferroni),x=BP),size=.1)+
        geom_hline(yintercept=-log10(0.001),color="red")+
        labs(x="position",y="-log10(p-value)")+
        scale_fill_manual(values=group.colors)+
        theme_light()+
        facet_grid(CHR~.))
print(ggplot()+
        #geom_rect(data=ChrStr,mapping=aes(ymin=0,ymax=15,xmin=Start,xmax=Stop, fill=Type),alpha=0.5)+
        geom_point(data=gea5,mapping=aes(y=-log10(pbonferroni),x=BP),size=.1)+
         geom_hline(yintercept=-log10(0.001),color="red")+
        labs(x="position",y="-log10(p-value)")+
        scale_fill_manual(values=group.colors)+
        theme_light()+
        facet_grid(CHR~.))
print(ggplot()+
        #geom_rect(data=ChrStr,mapping=aes(ymin=0,ymax=15,xmin=Start,xmax=Stop, fill=Type),alpha=0.5)+
        geom_point(data=gea4,mapping=aes(y=-log10(pbonferroni),x=BP),size=.1)+
        geom_hline(yintercept=-log10(0.001),color="red")+
        labs(x="position",y="-log10(p-value)")+
        scale_fill_manual(values=group.colors)+
        theme_light()+
        facet_grid(CHR~.))
print(ggplot()+
        #geom_rect(data=ChrStr,mapping=aes(ymin=0,ymax=15,xmin=Start,xmax=Stop, fill=Type),alpha=0.5)+
        geom_point(data=gea3,mapping=aes(y=-log10(pbonferroni),x=BP),size=.1)+
        geom_hline(yintercept=-log10(0.001),color="red")+
        labs(x="position",y="-log10(p-value)")+
        scale_fill_manual(values=group.colors)+
        theme_light()+
        facet_grid(CHR~.))
print(ggplot()+
        #geom_rect(data=ChrStr,mapping=aes(ymin=0,ymax=15,xmin=Start,xmax=Stop, fill=Type),alpha=0.5)+
        geom_point(data=gea2,mapping=aes(y=-log10(pbonferroni),x=BP),size=.1)+
        geom_hline(yintercept=-log10(0.001),color="red")+
        labs(x="position",y="-log10(p-value)")+
        scale_fill_manual(values=group.colors)+
        theme_light()+
        facet_grid(CHR~.))
print(ggplot()+
        #geom_rect(data=ChrStr,mapping=aes(ymin=0,ymax=15,xmin=Start,xmax=Stop, fill=Type),alpha=0.5)+
        geom_point(data=gea1,mapping=aes(y=-log10(pbonferroni),x=BP),size=.1)+
        geom_hline(yintercept=-log10(0.001),color="red")+
        labs(x="position",y="-log10(p-value)")+
        scale_fill_manual(values=group.colors)+
        theme_light()+
        facet_grid(CHR~.))

group.colors_e<-c(e1="skyblue", e2="salmon",e3="gold3",e4="skyblue4",e5="blue2",e6="red4",
                Centromere="#008b8b",Inversion="#5b0000",Chromosome="gray40")
invs<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2b.txt")
invs$csome<-invs$chromosome
cents<-subset(read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="CentRegion"))
cents$csome<-cents$scaffold
C_Coords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/mac_csome_dims.xlsx",sheet="Sheet1")
C_Coords$csome<-C_Coords$chr

subset(gea1,pbonferroni>=0.05)
pdf("C:/Users/rmohn/Desktop/20_WRITING_and_NOTES/Paper1/03_figures/GEA_CSOMES.pdf",width=7,height=5)
print(ggplot()+
        geom_rect(data=invs,mapping=aes(ymin=-2.5,ymax=-1.5,xmin=start,xmax=end, fill="Inversion"))+
        geom_rect(data=cents,mapping=aes(ymin=-1.25,ymax=-0.25,xmin=Start,xmax=Stop, fill="Centromere"))+
        geom_point(data=subset(gea1,pbonferroni<=0.05),mapping=aes(y=-log10(pbonferroni),x=BP,color="e1"),size=.5)+
        geom_point(data=subset(gea2,pbonferroni<=0.05),mapping=aes(y=-log10(pbonferroni),x=BP,color="e2"),size=.5)+
        geom_point(data=subset(gea3,pbonferroni<=0.05),mapping=aes(y=-log10(pbonferroni),x=BP,color="e3"),size=.5)+
        geom_point(data=subset(gea4,pbonferroni<=0.05),mapping=aes(y=-log10(pbonferroni),x=BP,color="e4"),size=.5)+
        geom_point(data=subset(gea5,pbonferroni<=0.05),mapping=aes(y=-log10(pbonferroni),x=BP,color="e5"),size=.5)+
        geom_point(data=subset(gea6,pbonferroni<=0.05),mapping=aes(y=-log10(pbonferroni),x=BP,color="e6"),size=.5)+
        geom_rect(data=C_Coords,mapping=aes(ymin=0,ymax=-log10(0.10),xmin=start,xmax=stop,fill="Chromosome"))+
        #geom_hline(yintercept=-log10(0.001),color="red")+
        labs(x="position",y="-log10(p-value)",colour=NULL,fill=NULL)+
        scale_color_manual(values=group.colors_e,labels=c("Mean diurnal range",
                          "Min temp. of coldest month","Mean temp. of wettest quarter",
                          "Precip. of wettest month","Precip. of driest month",
                          "Precip. of warmest quarter"))+
        scale_fill_manual(values=group.colors_e)+#,labels=c("Cent. array","Chromosome",
           #                                            "Inter−hap. inversion",
           #                                            "Cent. Inter−array region",
            #                                           "Interspecies inversion"))+
        theme_light()+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank(),text=element_text(size = 9),legend.position = "bottom",legend.key.size = unit(.25,'cm'))+
        
        facet_wrap(csome~.,nrow = 4,ncol = 3))
dev.off()
