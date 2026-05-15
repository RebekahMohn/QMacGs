To run sweepfinder2, we first used VCFtools to extract each chromosome and each population of samples in the form of counts. (done in sweepfinder_mac_W_run).
Then we used an new Rscript (sweepformat.R) to reformat those counts into the input file for SweepFinder2 (done in sweepfinder_mac_W_run).

Then we calculated the population frequency using SweepFinder2 -f (done in sweepfinder_mac_*_run).

Lastly we calculated the sweep probability with SweepFinder2 -lg 2000 (done in sf_mac_*).

Since a few of the chromosomes timed out before SweepFinder2 finished, we wrote an R script to make a command file (sf2_catchup_scripts.R)
This command file did the same as above except that it trimmed the chromosomes to the relevant regions and used the previously calculated population frequency. (find in folder "catchup").