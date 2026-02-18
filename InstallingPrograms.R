#ClumpGEAvals
install.packages("genetics.binaRies", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
library(genetics.binaRies)
install.packages("devtools")
install.packages("remotes")
remotes::install_github("MRCIEU/gwasglue")
plink_bin <- genetics.binaRies::get_plink_binary()


