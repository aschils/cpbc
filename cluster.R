library(cluster)
library(MASS)

#PATH - parameters:
start.time <- Sys.time()

args = commandArgs(trailingOnly=TRUE)
algo_type = args[1]

hclust_type = "hclust"
kmedoids_type = "kmedoids"

dist_file_path = args[3]#"distances.txt"
files_dir = args[4]
out_dir = args[5]

########## parse the Jplag output ###########
data = readLines(dist_file_path)
#only keep file1-file2 and corresponding distance
data = data[grep("^Comparing", data)]
#Remove useless Comparing word at beginning of line
data = gsub("^.*? ", " ", data)

print("Removing spaces left at the beginning of each line")
for(i in 1:length(data)){
  data[i] = substring(data[i], 2)
}

data = strsplit(data, " ")

print("Retrieve files names in a vector")
files_names = c()
for(i in 1:length(data)){
  names_pair = unlist(strsplit(unlist(data[i]), "-"))
  for(j in 1:(length(names_pair)-1)){
    file_name = gsub(":", "", names_pair[j])
    if(!(file_name %in% files_names)){
        files_names = c(files_names, file_name)
    }
  }
}

#Build the distance matrix
#Position of the name of a file in the vector files_names is used as file index
n = length(files_names)

print("Create distances matrix")
dist_matrix = matrix(rep(0,n*n), nrow = n, ncol = n, byrow = TRUE)

get_file_idx <- function(file_name){
  which(files_names == file_name)
}

max_dist = 100 #With Jplag

for(i in 1:length(data)){
  files_pair_dist = unlist(strsplit(unlist(data[i]), "-")) #file_name1 file_name2 dist
  file_name_1 = files_pair_dist[1]
  file_name_2 = gsub(":", "", files_pair_dist[2])
  #print(c("f1 ", file_name_1, "f2 ", file_name_2))
  dist = max_dist-as.numeric(files_pair_dist[3])
  file_idx_1 = get_file_idx(file_name_1)
  file_idx_2 = get_file_idx(file_name_2)
  dist_matrix[file_idx_1, file_idx_2] = dist
  dist_matrix[file_idx_2, file_idx_1] = dist
  #print(c("idx1 ", file_idx_1, "idx2 ", file_idx_2))
  #print(c("dist ", dist, "in matrix ", toString(dist_matrix[file_idx_1, file_idx_2])))
}

# print("DIST MATRIX:")
# for(i in 1:n){
#   for(j in 1:n){
#     print(dist_matrix[i,j])
#   }
# }

write.matrix(dist_matrix, file = "dist_matrix.txt", sep = " ", blocksize=n/100)

#dist_matrix = replicate(n, rnorm(n))

#print(dist_matrix[get_file_idx("attempt002_20140920_112132.cs"), get_file_idx("attempt002_20140920_110800.cs")])

print("Retrieve low triangle of distances matrix")
low_tri_matrix_el_vec = dist_matrix[lower.tri(dist_matrix, diag = FALSE)]

#For the clustering algorithm to interpret input as a distance matrix,
#the input object mus be of type "dissimilatiry matrix"
dissimilarity_mat<-low_tri_matrix_el_vec
class(dissimilarity_mat)='dist'
#attr(dissimilarity_mat,"Size")<-length(low_tri_matrix_el_vec)
attr(dissimilarity_mat,"Size")<-n
dissimilarity_mat<-as.dist(dissimilarity_mat)

#print(dissimilarity_mat)

out_clusters_info <- function(out_dir, cluster_nbr, clustering_vec, medoid=FALSE){

  dir.create(out_dir, showWarnings = TRUE)

  #Display members in each cluster
  for(i in 1:cluster_nbr){

    cluster_name = toString(i)
    #print(c("Cluster ", cluster_name))
    #print("############")

    if(medoid){
      #print("Medoid:")
      medoid_idx = medoids[paste("medoids", cluster_name, sep="")]
      medoid_name = files_names[medoid_idx]
      #print(medoid_name)
      #print("-----------")
    }

    out_file = paste(out_dir, paste("/", paste(cluster_name, ".txt", sep=""), sep=""), sep="")
    cat(cluster_name, file=out_file, append=FALSE, sep = "\n")

    members_counter = 0
    for(j in 1:length(clustering_vec)){
      if(clustering_vec[j] == i){
        file_name = files_names[j]
        file_content = readLines(paste(files_dir, file_name, sep=""))
        cat(file_name, file=out_file, append=TRUE, sep = "\n")
        cat("################", file=out_file, append=TRUE, sep = "\n")
        cat(file_content, file=out_file, append=TRUE, sep = "\n")
        #print(files_names)
        members_counter = members_counter+1
      }
    }
    cat(paste("Number of samples in cluster: ", members_counter, sep=""), file=out_file, append=TRUE, sep = "\n")
  }
}

reconstruct_err <- function(cluster_info){
  cluster_info[,1] %*% cluster_info[,3]
}

#sum of the diameters
#(maximal dissimilarity between two observations of the cluster) of the clusters
sum_diameters <- function(cluster_info){
  sum(cluster_info[,4])
}

#sum of the separations (minimal dissimilarity between an observation of the
#cluster and an observation of another cluster).
sum_separations <- function(cluster_info){
  sum(cluster_info[,5])
}


#pam(x, k, diss = inherits(x, "dist"), metric = "euclidean",
#    medoids = NULL, stand = FALSE, cluster.only = FALSE,
#    do.swap = TRUE,
#    keep.diss = !diss && !cluster.only && n < 100,
#    keep.data = !diss && !cluster.only,
#    pamonce = FALSE, trace.lev = 0)
if(algo_type == kmedoids_type){
  print("Process clustering - k medoids")
  cluster_nbr = args[2]
  clusters = pam(dissimilarity_mat, cluster_nbr, diss=TRUE)
  clustering_vec = unlist(clusters["clustering"]) #association of each sample to a cluster
  medoids = unlist(clusters["medoids"]) #representative objects of the clusters "centers" that belong to the dataset
  out_clusters_info(out_dir, cluster_nbr, clustering_vec, medoid=TRUE)
  cluster_info = clusters$clusinfo
  reconstruct_err = reconstruct_err(cluster_info)
  #print(c("k", cluster_nbr, "reconstruct_err", reconstruct_err))
  write(reconstruct_err,file="/Users/aschils/Software/plagiat_detection/clustering/results/Sector4_Level6/kmedoids/a.txt",append=TRUE)
  sum_diameters = sum_diameters(cluster_info)
  #print(c("k", cluster_nbr, "sum diameters", sum_diameters))
  write(sum_diameters,file="/Users/aschils/Software/plagiat_detection/clustering/results/Sector4_Level6/kmedoids/b.txt",append=TRUE)
  sum_separations = sum_separations(cluster_info)
  #print(c("k", cluster_nbr, "sum separations", sum_separations(cluster_info)))
  write(sum_separations,file="/Users/aschils/Software/plagiat_detection/clustering/results/Sector4_Level6/kmedoids/c.txt",append=TRUE)

  #print(cluster_info)
}

if(algo_type == hclust_type){
  print("Process clustering - hclust")
  hclusters = hclust(dissimilarity_mat, method = "ward.D", members = NULL)
  print(hclusters)
  pdf("plots.pdf")
  plot(hclusters, labels=files_names)
  dev.off()
  tree_height = args[2]
  hclust_vec = cutree(hclusters, h=tree_height)

  k = max(hclust_vec)
  out_clusters_info(out_dir, k, hclust_vec, medoid=FALSE)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
