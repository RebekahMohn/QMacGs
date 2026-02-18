#Testing relationship of inversions, gea, introgression, etc.

library(readxl)
library(dplyr)
library(ggplot2)
#genomic regions
ChrStr<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="Sheet3")
#assigning colors for use later
group.colors<-c(lyr="#000000",bic="#88CCEE",lob="#116644",mue="#AA5599",ste="#E69F00",alb="#332288",none="white",CenArray="gray20",InterArray="gray80", IntraSpInv="gray70",InterSpInv="gray35")

#renaming column for ease of use
ChrStr$csome<-ChrStr$scaffold

#importing chromosome dimensions
C_Coords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/mac_csome_dims.xlsx",sheet="Sheet1")


csomeSamp<-rep(C_Coords$chr,round(C_Coords$stop/100000,0))

gea_env<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")
##Read in gea significant

##read in all gea
pval_genome_env<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval.txt",sep=" ")

invs<-read.delim(file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/Invs2b.txt")

InvGEATest<-data.frame(invSize=integer(),samplenumber=integer(),invCsome=character(),invEnd=integer(),geaSite=integer())
for(j in invs$span_bp.x){
  randCsome<-sample(csomeSamp,200)
  for(i in 1:length(randCsome)){
    GEAs<-data.frame(e1=0,e2=0,e3=0,e4=0,e5=0,e6=0)
    randSamps<-round(runif(min=j,max=C_Coords[C_Coords[,1]==randCsome[i],3][[1]],1),0)
    geaSites<-subset(gea_env,csome==randCsome[i]&BP>randSamps-j&BP<randSamps)
    SNPsites<-length(subset(pval_genome_env,csome==randCsome[i]&V4>randSamps-j&V4<randSamps)$V4)
    if(length(geaSites$BP>0)){
      geaN<-length(geaSites$BP)
      GEAs<-as.data.frame(t(data.frame(table(geaSites$ENV))[,2]))
      colnames(GEAs)<-t(data.frame(table(geaSites$ENV))[,1])
      rownames(GEAs)<-NULL
    }else{
      geaN<-0
    }
    InvGEATest<-bind_rows(cbind(invSize=j,samplenumber=i,invCsome=randCsome[i],invEnd=randSamps,sites=SNPsites,geaSite=geaN,GEAs),InvGEATest) #sites=SNPsites,
  }
}

InvGEATest[is.na(InvGEATest[,6]),6]<-0
InvGEATest[is.na(InvGEATest[,7]),7]<-0
InvGEATest[is.na(InvGEATest[,8]),8]<-0
InvGEATest[is.na(InvGEATest[,9]),9]<-0
InvGEATest[is.na(InvGEATest[,10]),10]<-0
InvGEATest[is.na(InvGEATest[,11]),11]<-0
InvGEATest[is.na(InvGEATest[,12]),12]<-0
InvGEATest$pe1<-InvGEATest$e1/InvGEATest$sites
InvGEATest$pe2<-InvGEATest$e2/InvGEATest$sites
InvGEATest$pe3<-InvGEATest$e3/InvGEATest$sites
InvGEATest$pe4<-InvGEATest$e4/InvGEATest$sites
InvGEATest$pe5<-InvGEATest$e5/InvGEATest$sites
InvGEATest$pe6<-InvGEATest$e6/InvGEATest$sites


####


InvGEA<-data.frame(invSize=integer(),samplenumber=integer(),invCsome=character(),invEnd=integer(),geaSites=integer(),e1=integer(),e2=integer(),e3=integer(),e4=integer(),e5=integer(),e6=integer())
for(j in 1:length(invs$chr_start)){
  GEAs<-data.frame(e1=0,e2=0,e3=0,e4=0,e5=0,e6=0)
  geaSites<-subset(gea_env,csome==invs$chromosome[j]&BP>invs$start[j]&BP<invs$end[j])
  SNPsites<-length(subset(pval_genome_env,csome==invs$chromosome[j]&V4>invs$start[j]&V4<invs$end[j])$V4)
    if(length(geaSites$BP>0)){
      geaN<-length(geaSites$BP)
      GEAs$e1[1]<-length(subset(geaSites,ENV=="e1")$ENV)
      GEAs$e2[1]<-length(subset(geaSites,ENV=="e2")$ENV)
      GEAs$e3[1]<-length(subset(geaSites,ENV=="e3")$ENV)
      GEAs$e4[1]<-length(subset(geaSites,ENV=="e4")$ENV)
      GEAs$e5[1]<-length(subset(geaSites,ENV=="e5")$ENV)
      GEAs$e6[1]<-length(subset(geaSites,ENV=="e6")$ENV)
    }else{
      geaN<-0
    }
    InvGEA<-bind_rows(cbind(invSize=invs$span_bp.x[j],invnum=j,invCsome=invs$chromosome[j],invEnd=invs$end[j],sites=SNPsites,geaSite=geaN,GEAs,type=invs$INV[j]),InvGEA)
}





InvGEA$qe1<-NA
InvGEA$qe2<-NA
InvGEA$qe3<-NA
InvGEA$qe4<-NA
InvGEA$qe5<-NA
InvGEA$qe6<-NA

for(k in 1:length(InvGEA$invSize)){
  InvGEA$qe1[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&e1>=InvGEA$e1[k])$e1)/200
  InvGEA$qe2[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&e2>=InvGEA$e2[k])$e2)/200
  InvGEA$qe3[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&e3>=InvGEA$e3[k])$e3)/200
  InvGEA$qe4[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&e4>=InvGEA$e4[k])$e4)/200
  InvGEA$qe5[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&e5>=InvGEA$e5[k])$e5)/200
  InvGEA$qe6[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&e6>=InvGEA$e6[k])$e6)/200
}


#Correct for number of SNPs in an area.
InvGEA$Sitesqe1<-NA
InvGEA$Sitesqe2<-NA
InvGEA$Sitesqe3<-NA
InvGEA$Sitesqe4<-NA
InvGEA$Sitesqe5<-NA
InvGEA$Sitesqe6<-NA

for(k in 1:length(InvGEA$invSize)){
  InvGEA$Sitesqe1[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe1>=InvGEA$pe1[k])$e1)/200
  InvGEA$Sitesqe2[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe2>=InvGEA$pe2[k])$e2)/200
  InvGEA$Sitesqe3[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe3>=InvGEA$pe3[k])$e3)/200
  InvGEA$Sitesqe4[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe4>=InvGEA$pe4[k])$e4)/200
  InvGEA$Sitesqe5[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe5>=InvGEA$pe5[k])$e5)/200
  InvGEA$Sitesqe6[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe6>=InvGEA$pe6[k])$e6)/200
}


InvGEA$Siteslqe1<-NA
InvGEA$Siteslqe2<-NA
InvGEA$Siteslqe3<-NA
InvGEA$Siteslqe4<-NA
InvGEA$Siteslqe5<-NA
InvGEA$Siteslqe6<-NA

for(k in 1:length(InvGEA$invSize)){
  InvGEA$Siteslqe1[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe1<=InvGEA$pe1[k])$e1)/200
  InvGEA$Siteslqe2[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe2<=InvGEA$pe2[k])$e2)/200
  InvGEA$Siteslqe3[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe3<=InvGEA$pe3[k])$e3)/200
  InvGEA$Siteslqe4[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe4<=InvGEA$pe4[k])$e4)/200
  InvGEA$Siteslqe5[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe5<=InvGEA$pe5[k])$e5)/200
  InvGEA$Siteslqe6[k]<-length(subset(InvGEATest,invSize==InvGEA$invSize[k]&pe6<=InvGEA$pe6[k])$e6)/200
}


ggplot()+
  geom_boxplot(mapping=aes(x="e1",y=InvGEA$Sitesqe1))+
  geom_boxplot(mapping=aes(x="e2",y=InvGEA$Sitesqe2))+
  geom_boxplot(mapping=aes(x="e3",y=InvGEA$Sitesqe3))+
  geom_boxplot(mapping=aes(x="e4",y=InvGEA$Sitesqe4))+
  geom_boxplot(mapping=aes(x="e5",y=InvGEA$Sitesqe5))+
  geom_boxplot(mapping=aes(x="e6",y=InvGEA$Sitesqe6))+
  theme_light()

write.table(InvGEA,file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/InvGEA.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(InvGEATest,file="C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/InvGEATest.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

####Do for centromere interarray regions
library(readxl)
library(dplyr)
library(ggplot2)
#genomic regions
# ChrStr<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="Sheet3")
# #assigning colors for use later
# group.colors<-c(lyr="#000000",bic="#88CCEE",lob="#116644",mue="#AA5599",ste="#E69F00",alb="#332288",none="white",CenArray="gray20",InterArray="gray80", IntraSpInv="gray70",InterSpInv="gray35")
# 
# #renaming column for ease of use
# ChrStr$csome<-ChrStr$scaffold
# 
# #importing chromosome dimensions
# C_Coords<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/mac_csome_dims.xlsx",sheet="Sheet1")
# 
# 
# csomeSamp<-rep(C_Coords$chr,round(C_Coords$stop/100000,0))
# 
# gea_env<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval_cor.txt",sep=" ")
# ##Read in gea significant
# 
# ##read in all gea
# pval_genome_env<-read.delim("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/AllC_lfmmpval.txt",sep=" ")


# InvGEATest<-data.frame(invSize=integer(),samplenumber=integer(),invCsome=character(),invEnd=integer(),geaSite=integer())
# for(j in seq(from=50000,to=16000000,by=50000)){
#   randCsome<-sample(csomeSamp,100)
#   for(i in 1:length(randCsome)){
#     GEAs<-data.frame(e1=0,e2=0,e3=0,e4=0,e5=0,e6=0)
#     randSamps<-round(runif(min=j,max=C_Coords[C_Coords[,1]==randCsome[i],3][[1]],1),0)
#     geaSites<-subset(gea_env,csome==randCsome[i]&BP>randSamps-j&BP<randSamps)
#     SNPsites<-sum(pval_genome_env$csome==randCsome[i]&pval_genome_env$V4>randSamps-j&pval_genome_env$V4<randSamps)
#     if(length(geaSites$BP>0)){
#       geaN<-length(geaSites$BP)
#       GEAs<-as.data.frame(t(data.frame(table(geaSites$ENV))[,2]))
#       colnames(GEAs)<-t(data.frame(table(geaSites$ENV))[,1])
#       rownames(GEAs)<-NULL
#     }else{
#       geaN<-0
#     }
#     InvGEATest<-bind_rows(cbind(invSize=j,samplenumber=i,invCsome=randCsome[i],invEnd=randSamps,sites=SNPsites,geaSite=geaN,GEAs),InvGEATest) #sites=SNPsites,
#   }
# }
# 
# InvGEATest[is.na(InvGEATest[,6]),6]<-0
# InvGEATest[is.na(InvGEATest[,7]),7]<-0
# InvGEATest[is.na(InvGEATest[,8]),8]<-0
# InvGEATest[is.na(InvGEATest[,9]),9]<-0
# InvGEATest[is.na(InvGEATest[,10]),10]<-0
# InvGEATest[is.na(InvGEATest[,11]),11]<-0
# InvGEATest[is.na(InvGEATest[,12]),12]<-0
# InvGEATest$pe1<-InvGEATest$e1/InvGEATest$sites
# InvGEATest$pe2<-InvGEATest$e2/InvGEATest$sites
# InvGEATest$pe3<-InvGEATest$e3/InvGEATest$sites
# InvGEATest$pe4<-InvGEATest$e4/InvGEATest$sites
# InvGEATest$pe5<-InvGEATest$e5/InvGEATest$sites
# InvGEATest$pe6<-InvGEATest$e6/InvGEATest$sites


##############################
### Centromeres ##############
##############################

cents<-read_excel("C:/Users/rmohn/Desktop/10_Analysis/20_MacGenomeRegions/MacGenomeRegionsB.xlsx",sheet="InterArrays")

###test centromeres
centGEATest<-data.frame(centSize=integer(),samplenumber=integer(),centCsome=character(),centEnd=integer(),geaSite=integer())
for(j in unique(cents$Length)){
  randCsome<-sample(csomeSamp,200)
  for(i in 1:length(randCsome)){
    GEAs<-data.frame(e1=0,e2=0,e3=0,e4=0,e5=0,e6=0)
    randSamps<-round(runif(min=j,max=C_Coords[C_Coords[,1]==randCsome[i],3][[1]],1),0)
    geaSites<-subset(gea_env,csome==randCsome[i]&BP>randSamps-j&BP<randSamps)
    SNPsites<-length(subset(pval_genome_env,csome==randCsome[i]&V4>randSamps-j&V4<randSamps)$V4)
    if(length(geaSites$BP>0)){
      geaN<-length(geaSites$BP)
      GEAs<-as.data.frame(t(data.frame(table(geaSites$ENV))[,2]))
      colnames(GEAs)<-t(data.frame(table(geaSites$ENV))[,1])
      rownames(GEAs)<-NULL
    }else{
      geaN<-0
    }
    centGEATest<-bind_rows(cbind(centSize=j,samplenumber=i,centCsome=randCsome[i],centEnd=randSamps,sites=SNPsites,geaSite=geaN,GEAs),centGEATest) #sites=SNPsites,
  }
}

centGEATest[is.na(centGEATest[,6]),6]<-0
centGEATest[is.na(centGEATest[,7]),7]<-0
centGEATest[is.na(centGEATest[,8]),8]<-0
centGEATest[is.na(centGEATest[,9]),9]<-0
centGEATest[is.na(centGEATest[,10]),10]<-0
centGEATest[is.na(centGEATest[,11]),11]<-0
centGEATest[is.na(centGEATest[,12]),12]<-0
centGEATest$pe1<-centGEATest$e1/centGEATest$sites
centGEATest$pe2<-centGEATest$e2/centGEATest$sites
centGEATest$pe3<-centGEATest$e3/centGEATest$sites
centGEATest$pe4<-centGEATest$e4/centGEATest$sites
centGEATest$pe5<-centGEATest$e5/centGEATest$sites
centGEATest$pe6<-centGEATest$e6/centGEATest$sites


CentsGEA<-data.frame(CentsSize=integer(),samplenumber=integer(),CentsCsome=character(),CentsEnd=integer(),sites=integer(),geaSite=integer(),e1=integer(),e2=integer(),e3=integer(),e4=integer(),e5=integer(),e6=integer())
for(j in 1:length(cents$scaffold)){
  GEAs<-data.frame(e1=0,e2=0,e3=0,e4=0,e5=0,e6=0)
  geaSites<-subset(gea_env,csome==cents$scaffold[j]&BP>cents$Start[j]&BP<cents$Stop[j])
  SNPsites<-length(subset(pval_genome_env,csome==cents$scaffold[j]&V4>cents$Start[j]&V4<cents$Stop[j])$V4)
  if(length(geaSites$BP>0)){
    geaN<-length(geaSites$BP)
    GEAs$e1[1]<-length(subset(geaSites,ENV=="e1")$ENV)
    GEAs$e2[1]<-length(subset(geaSites,ENV=="e2")$ENV)
    GEAs$e3[1]<-length(subset(geaSites,ENV=="e3")$ENV)
    GEAs$e4[1]<-length(subset(geaSites,ENV=="e4")$ENV)
    GEAs$e5[1]<-length(subset(geaSites,ENV=="e5")$ENV)
    GEAs$e6[1]<-length(subset(geaSites,ENV=="e6")$ENV)
  }else{
    geaN<-0
  }
  CentsGEA<-bind_rows(cbind(CentsSize=cents$Length[j],Centsnum=j,CentsCsome=cents$scaffold[j],CentsEnd=cents$Stop[j],sites=SNPsites,geaSite=geaN,GEAs),CentsGEA)
}





CentsGEA$qe1<-NA
CentsGEA$qe2<-NA
CentsGEA$qe3<-NA
CentsGEA$qe4<-NA
CentsGEA$qe5<-NA
CentsGEA$qe6<-NA

for(k in 1:length(CentsGEA$CentsSize)){
  CentsGEA$qe1[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e1>=CentsGEA$e1[k])$e1)/200
  CentsGEA$qe2[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e2>=CentsGEA$e2[k])$e2)/200
  CentsGEA$qe3[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e3>=CentsGEA$e3[k])$e3)/200
  CentsGEA$qe4[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e4>=CentsGEA$e4[k])$e4)/200
  CentsGEA$qe5[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e5>=CentsGEA$e5[k])$e5)/200
  CentsGEA$qe6[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e6>=CentsGEA$e6[k])$e6)/200
}


CentsGEA$lqe1<-NA
CentsGEA$lqe2<-NA
CentsGEA$lqe3<-NA
CentsGEA$lqe4<-NA
CentsGEA$lqe5<-NA
CentsGEA$lqe6<-NA

for(k in 1:length(CentsGEA$CentsSize)){
  CentsGEA$lqe1[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e1<=CentsGEA$e1[k])$e1)/200
  CentsGEA$lqe2[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e2<=CentsGEA$e2[k])$e2)/200
  CentsGEA$lqe3[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e3<=CentsGEA$e3[k])$e3)/200
  CentsGEA$lqe4[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e4<=CentsGEA$e4[k])$e4)/200
  CentsGEA$lqe5[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e5<=CentsGEA$e5[k])$e5)/200
  CentsGEA$lqe6[k]<-length(subset(centGEATest,centSize==CentsGEA$CentsSize[k]&e6<=CentsGEA$e6[k])$e6)/200
}


ggplot()+
  geom_boxplot(mapping=aes(x="e1",y=CentsGEA$qe1))+
  geom_boxplot(mapping=aes(x="e2",y=CentsGEA$qe2))+
  geom_boxplot(mapping=aes(x="e3",y=CentsGEA$qe3))+
  geom_boxplot(mapping=aes(x="e4",y=CentsGEA$qe4))+
  geom_boxplot(mapping=aes(x="e5",y=CentsGEA$qe5))+
  geom_boxplot(mapping=aes(x="e6",y=CentsGEA$qe6))+
  theme_light()

