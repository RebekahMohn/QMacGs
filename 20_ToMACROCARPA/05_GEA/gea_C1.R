library(LEA)

###convert ped to lfmm
ped2lfmm(input.file="/hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C1_maf01.ped",force=FALSE)
lfmmInput<-read.lfmm("/hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C1_maf01.lfmm")
###import environmental daa
write.env(read.delim("BioClim6vars.txt",sep="\t")[,-1],"BC6.env")
lfmmEnv<-read.env(input.file="BC6.env")
###run lfmm
lfmmOut<-lfmm2(input = lfmmInput, env = lfmmEnv, K=2)
#write.lfmm(lfmmOut, "/hpc/group/manoslab/QMACROCARPA_GENOMES/24_VARIANTS_MACROCARPA/MAC99_MACREF_HQSNP_C1_maf01.lfmmOut.txt")


write.table(lfmmOut@U,"gea_out/C1_lfmmfactors.txt",sep="\t")
write.table(lfmmOut@V,"gea_out/C1_lfmmloadings.txt",sep="\t")
write.table(lfmmOut@K,"gea_out/C1_lfmmK.txt",sep="\t")
write.table(lfmmOut@lambda,"gea_out/C1_lfmmlambda.txt",sep="\t")
write.table(lfmmOut@B,"gea_out/C1_lfmmB.txt",sep="\t")
pdf("gea_out/C1_lfmmfactors.pdf",height=4,width=4)
plot(lfmmOut@U, col = "grey", pch = 19,xlab = "Factor 1",ylab = "Factor 2")
dev.off()
pv <- lfmm2.test(object = lfmmOut,input = lfmmInput,env = lfmmEnv,full = TRUE)

write.table(pv$pvalues,"gea_out/C1_lfmmpval.txt",sep="\t")
write.table(pv$zscores,"gea_out/C1_lfmmzsco.txt",sep="\t")
write.table(pv$fscores,"gea_out/C1_lfmmfsco.txt",sep="\t")
write.table(pv$adj.r.squared,"gea_out/C1_lfmmadjrsq.txt",sep="\t")
write.table(pv$gif,"gea_out/C1_lfmmgif.txt",sep="\t")
pdf("gea_out/C1_lfmmpval.pdf",height=4,width=4)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")
dev.off()
