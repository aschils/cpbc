library(cluster)
library(MASS)
library(dendextend)

VERBOSE = FALSE
options(warn=-1)

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

if(VERBOSE)
  print("Removing spaces left at the beginning of each line")
for(i in 1:length(data)){
  data[i] = substring(data[i], 2)
}

data = strsplit(data, " ")

if(VERBOSE)
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

##################### Build the distance matrix #####################
#Position of the name of a file in the vector files_names is used as file index
n = length(files_names)

if(VERBOSE)
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

#write.matrix(dist_matrix, file = "dist_matrix.txt", sep = " ", blocksize=n/100)

#dist_matrix = replicate(n, rnorm(n))

#print(dist_matrix[get_file_idx("attempt002_20140920_112132.cs"), get_file_idx("attempt002_20140920_110800.cs")])

if(VERBOSE)
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

################### manually made clustering ######################

c1 = c("attempt001_20140920_022105.cs",
"attempt001_20140920_032927.cs",
"attempt001_20140920_102656.cs",
"attempt001_20140920_102739.cs",
"attempt001_20140920_104554.cs",
"attempt001_20140920_112533.cs",
"attempt001_20140920_164800.cs",
"attempt001_20140920_195521.cs",
"attempt001_20140920_203756.cs",
"attempt001_20140920_210358.cs",
"attempt002_20140920_112536.cs",
"attempt002_20140920_204122.cs",
"attempt003_20140920_174456.cs")

c2 = c("attempt002_20140920_022248.cs",
"attempt003_20140920_022256.cs",
"attempt004_20140920_022426.cs",
"attempt005_20140920_022441.cs")

c3 = c("attempt002_20140920_103730.cs",
"attempt003_20140920_103748.cs")

c4 = c("attempt002_20140920_104459.cs",
"attempt004_20140920_105646.cs",
"attempt004_20140920_175138.cs",
"attempt005_20140920_105658_winning3.cs",
"attempt005_20140920_175154_winning3.cs")

c5 = c("attempt012_20140920_213449_winning2.cs",
"attempt013_20140920_213548.cs",
"attempt014_20140920_213605.cs")

c6 = c("attempt002_20140920_212357.cs",
"attempt003_20140920_212539.cs",
"attempt004_20140920_212607.cs",
"attempt005_20140920_212629.cs",
"attempt006_20140920_212741.cs",
"attempt007_20140920_212855_winning2.cs",
"attempt011_20140920_213242_winning2.cs",
"attempt015_20140920_213730.cs")

c7 = c("attempt005_20140920_103956.cs",
"attempt006_20140920_104006.cs",
"attempt007_20140920_104048.cs",
"attempt008_20140920_104121.cs",
"attempt009_20140920_104258_winning3.cs",
"attempt015_20140920_023545_winning3.cs")

c8 = c("attempt006_20140920_022537.cs",
"attempt007_20140920_022553.cs")

c9 = c("attempt008_20140920_022724.cs",
"attempt009_20140920_022742.cs",
"attempt010_20140920_022752.cs")

c10 = c("attempt011_20140920_022855.cs",
"attempt012_20140920_022927.cs")

c11 = c("attempt012_20140920_085955.cs",
"attempt013_20140920_023107.cs",
"attempt013_20140920_090132.cs",
"attempt014_20140920_023256.cs")

c12 = c("attempt002_20140920_111009_winning3.cs")

ideal_clusters <- list()
ideal_clusters[[1]] <- c1
ideal_clusters[[2]] <- c2
ideal_clusters[[3]] <- c3
ideal_clusters[[4]] <- c4
ideal_clusters[[5]] <- c5
ideal_clusters[[6]] <- c6
ideal_clusters[[7]] <- c7
ideal_clusters[[8]] <- c8
ideal_clusters[[9]] <- c9
ideal_clusters[[10]] <- c10
ideal_clusters[[11]] <- c11
ideal_clusters[[12]] <- c12

max_score <- function(){

  max_score = 0
  for(cluster in ideal_clusters){
    max_score = max_score + length(cluster)*(length(cluster)-1)
  }
  max_score
}

get_user_assigned_class <- function(file_name) {
  class_id = 1
  for(cluster in ideal_clusters){
    if(file_name %in% cluster){
        return(class_id)
    }
    class_id = class_id+1
  }
  0
}

same_user_assigned_class <- function(file_name1, file_name2){

  for(cluster in ideal_clusters){
    if(file_name1 %in% cluster){
        return(file_name2 %in% cluster)
    }
  }
  FALSE
}

clustering_score <- function(clustering_vec){

  score = 0

  for(i in 1:length(clustering_vec)){
    file_name1 = files_names[i]
    cluster = clustering_vec[i]

    for(j in 1:length(clustering_vec)){

      if(i != j && clustering_vec[j] == cluster){

        file_name2 = files_names[j]
        if(same_user_assigned_class(file_name1, file_name2))
          score = score+1
        else
          score = score-1
      }
    }
  }

  score/max_score()
}

out_clusters_info <- function(out_dir, cluster_nbr, clustering_vec, medoid=FALSE){

  dir.create(out_dir, showWarnings = TRUE)

  if(medoid){
    score = clustering_score(clustering_vec)
    #print(c("Clustering score: ", score))
    print(score)
  }

  #Display members in each cluster
  for(i in 1:cluster_nbr){

    cluster_name = toString(i)
    #print(c("Cluster ", cluster_name))
    #print("############")
    out_file = paste(out_dir, paste("/", paste(cluster_name, ".txt", sep=""), sep=""), sep="")
    cat(cluster_name, file=out_file, append=FALSE, sep = "\n")

    if(medoid){
      #print("Medoid:")
      medoid_idx = medoids[paste("medoids", cluster_name, sep="")]
      medoid_name = files_names[medoid_idx]
      cat(paste("medoid name:", medoid_name, sep=" "), file=out_file, append=TRUE, sep = "\n")
      #print(medoid_name)
      #print("-----------")
    }


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
  if(VERBOSE)
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

cut_level_height <- function(hclusters, cut_level){
  tree_heights = unique(unlist(hclusters["height"]))

  if(cut_level <= length(tree_heights))
    tree_heights[cut_level]
  else
    tree_heights[length(tree_heights)]
}

if(algo_type == hclust_type){
  if(VERBOSE)
    print("Process clustering - hclust")
  hclusters = hclust(dissimilarity_mat, method = "ward.D", members = NULL)
  dendogram = as.dendrogram(hclusters)
  #print(unlist(hclusters["height"]))

  #Plot dendogram with at leaf user assigned classfication of samples
  #pdf("plots.pdf")
  # jpeg("plots.jpeg", width = 8, height = 6, units = 'in', res = 1000)
  # user_assigned_class = rep("",length(files_names))
  # for(i in 1:length(user_assigned_class)){
  #   user_assigned_class[i] = paste(" ", paste(toString(get_user_assigned_class(files_names[i])), " "))
  # }
  #
  # dendogram <- set(dendogram, "labels", user_assigned_class)
  # #dendogram <- set(dendogram, "labels_cex", 0.1)
  # par(cex=0.5)#mar=c(5, 8, 4, 1))
  # #leaflab="textlike",
  # plot(dendogram, xlab="hierarchical clustering", leaflab="perpendicular", yaxt='n', horiz=TRUE, type = "rectangle")# leaflab=leaflab)
  # dev.off()

  jpeg("plots.jpeg", width = 8, height = 6, units = 'in', res = 1000)
  user_assigned_class = rep("",length(files_names))
  for(i in 1:length(user_assigned_class)){
    user_assigned_class[i] = toString(get_user_assigned_class(files_names[i]))
  }

  par(cex=0.8)
  plot(hclusters, labels=user_assigned_class, xlab="hierarchical clustering", yaxt='n')
  dev.off()

  #tree_height_level = as.numeric(args[2])
  k = as.numeric(args[2])

  #hclust_vec = cutree(hclusters, h=tree_height)
  #cut_height = cut_level_height(hclusters, tree_height_level)
  #print(cut_height)
  #hclust_vec = cutree(hclusters, h=cut_height)

  hclust_vec = cutree(hclusters, k=k)
  print(clustering_score(hclust_vec))

  #k = max(hclust_vec)
  #print(c("k:", k))
  out_clusters_info(out_dir, k, hclust_vec, medoid=FALSE)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
if(VERBOSE)
  print(time.taken)

# medoid_idx = get_file_idx("attempt004_20140920_175138.cs")
# print("Dist to medoid:")
#
# print("attempt002_20140920_104459.cs")
# f_idx = get_file_idx("attempt002_20140920_104459.cs")
# print(dist_matrix[medoid_idx, f_idx])
#
# print("attempt004_20140920_105646.cs")
# f_idx = get_file_idx("attempt004_20140920_105646.cs")
# print(dist_matrix[medoid_idx, f_idx])
#
# print("attempt005_20140920_105658_winning3.cs")
# f_idx = get_file_idx("attempt005_20140920_105658_winning3.cs")
# print(dist_matrix[medoid_idx, f_idx])
#
# print("attempt005_20140920_175154_winning3.cs")
# f_idx = get_file_idx("attempt005_20140920_175154_winning3.cs")
# print(dist_matrix[medoid_idx, f_idx])
