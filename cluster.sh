#!/bin/bash

#k-medoids loop on number of clusters
#for k in {1..30}
#do
#   Rscript cluster.R kmedoids $k "/Users/aschils/Software/plagiat_detection/clustering/distances_Sector4_Level6.txt" "/Users/aschils/Software/plagiat_detection/code_hunt/levels/Sector4_Level6/" "results/Sector4_Level6/kmedoids/k$k"
#done

hclust loop on tree height
for h in {2..30}
do
   Rscript cluster.R hclust $h "/Users/aschils/Software/plagiat_detection/clustering/distances_Sector4_Level6.txt" "/Users/aschils/Software/plagiat_detection/code_hunt/levels/Sector4_Level6/" "results/Sector4_Level6/hclust/k$h"
done
