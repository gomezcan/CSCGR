# This script requared RDM, and peak summit file in order to calculate topics 
# run gem RMD to identify regions 
./gem RMD --g Athaliana_chr.bed --tf_peak_file List.TFs.txt --distance 50 --min_site 3 --out At_RDM_results 
# Run hdp to 
hdp --algorithm train --data 0_BS_clusters.At_RDM_results.d50.min3.HDP.txt --directory train_dir --eta 0.1 --random_seed 0 --init_topics 50 --max_iter 20 
