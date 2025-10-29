#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=200G
#SBATCH --time=0-10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmohn@mortonarb.org
#SBATCH --job-name=albPaintPower


module load R

sed -i 's/\t\.\t\./\tNA\t\NA/g' alb_mac_bic.fixed2.txt
perl fix_num.pl ../alb_mac_bic.fixed2.txt > alb_mac_bic.perl.txt
Rscript splitR_postPerl_alb.R

#sed -i 's/\t\.\t\./\tNA\t\NA/g' lyr_mac_bic.fixed2.txt
#perl fix_num_nothin.pl ../lyr_mac_bic.fixed2.txt > lyr_mac_bic.perl.txt
#Rscript splitR_postPerl_lyr.R


sed -i 's/\t\.\t\./\tNA\t\NA/g' mue_mac_bic.fixed2.txt
perl fix_num.pl ../mue_mac_bic.fixed2.txt > mue_mac_bic.perl.txt
Rscript splitR_postPerl_2.R

sed -i 's/\t\.\t\./\tNA\t\NA/g' lob_mac_bic.fixed2.txt
perl fix_num.pl ../lob_mac_bic.fixed2.txt > lob_mac_bic.perl.txt
Rscript splitR_postPerl_lob.R

sed -i 's/\t\.\t\./\tNA\t\NA/g' ste_mac_bic.fixed2.txt
perl fix_num.pl ../ste_mac_bic.fixed2.txt > ste_mac_bic.perl.txt
Rscript splitR_postPerl_ste.R


Rscript AdIntWindows_mue.R
Rscript AdIntWindows_alb.R
Rscript AdIntWindows_ste.R
Rscript AdIntWindows_lob.R
