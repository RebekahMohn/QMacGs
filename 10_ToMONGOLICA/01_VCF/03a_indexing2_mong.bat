#!/bin/bash
i=1
for file in /cwork/ram163/OAK_GENOMES/14_MAPPING/02_VARIANTS/*.gz;
do
        fq=`basename ${file} .gvcf.gz`
#	echo "#!/bin/bash" > /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}#
#	echo "#SBATCH --cpus-per-task=20" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
#	echo "#SBATCH --mem=65G"  >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
#	echo "#SBATCH --mail-type=ALL"  >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
#	echo "#SBATCH --mail-user=rmohn@mortonarb.org"  >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
	echo "cd /cwork/ram163/OAK_GENOMES/14_MAPPING/02_VARIANTS/" > /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
	echo "module load bcftools" >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
	echo "bcftools call ${file} -m -O z -g 5 -o /cwork/ram163/OAK_GENOMES/14_MAPPING/03_CALLED/${fq}.vcf.gz --threads 20"  >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
	echo "bcftools index /cwork/ram163/OAK_GENOMES/14_MAPPING/03_CALLED/${fq}.vcf.gz --threads 20 " >> /hpc/group/manoslab/QMACROCARPA_GENOMES/SCRIPTS/10_ToMONGOLICA/index2_S${i}
        ((i=i+1))
done
