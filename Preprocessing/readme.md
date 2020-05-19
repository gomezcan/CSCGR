#### This folder does have descriptions of preprocessing steps required before any further analysis

1. Define local co-expression based on wPCC

- 1.1 Cluster samples to define co-expression: 

- 1.2 Using clusters, calculate weithed PCC for each cluster of samples

2. Define global co-expression based on Mutual information (MI)

3. Define local co-expression based on MI. Require cluster results from 1.1

4. Define targets genes based on peaks

5. Calculate partial correlacition of TFx-Targets conditioned by TFz. TFx and TFz should interect physically.

