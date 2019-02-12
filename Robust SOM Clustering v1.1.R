####################################################################################
#                                                                                  #
#         An R script to perform combined SOM/superSOM data clustering             #
#                   with robust initilizations  V1.1                               #
#                                                                                  #
# AuthorS: Joan Perez*                                                             #
#          E-mail: Joan.PEREZ@univ-cotedazur.fr                                    # 
#                                                                                  #
#          Giovanni Fusco*                                                         #
#          E-mail: Giovanni.FUSCO@univ-cotedazur.fr                                #
#                                                                                  #
#          *Univ. Cote d'Azur, CNRS, ESPACE, UMR7300                               #
#          98 Bd Herriot, BP 3209, 06204 Nice, France                              #
#                                                                                  #
####################################################################################


# R version 3.5.1 
# other attached packages:
# kohonen_2.0.19   dendextend_1.9.0 class_7.3-14     MASS_7.3-50      usethis_1.4.0  
# devtools_2.0.1 

####################################################################################
# 0.1 : Environment Preparation                                                    #
####################################################################################

# Packages loading
library(devtools)
install_version("kohonen", version = "2.0.19", repos = "http://cran.us.r-project.org")
library(MASS)
library(class)
library(dendextend)

####################################################################################
# 0.2. General Functions and Utilities                                             #
####################################################################################

## 0.2.1 Automated Number of Primes Selection / Initialization                       

select.primes <- function(nb)
{
  primes <- vector()
  t <- 0
  a <- 3
  while (length(primes)<nb) {for (i in a) { for (j in 2:(i-1)) {if (i %% j == 0){t <- (t+1)} } 
                                            if (t==0){primes<- c(primes,i)}
                                            t <- 0
  }
  a <- (a+1)
  }
  return(primes)
}

## 0.2.2 Calculation of the Best Similarity Score for each Seed

# function returns two arguments, first one is the best seed in chr
# second one is the average similarity score for each seed
FM.Similarity <- function(results.list, clusters)
{
FM.matrix <- matrix(nrow=length(results.list), ncol=length(results.list))
for(i in 1:length(results.list)) {for(j in 1:length(results.list)) 
{FM.matrix[i,j]<-FM_index_R(clusters[,i+1],clusters[,j+1],assume_sorted_vectors = TRUE) }}
rownames(FM.matrix) <- c(paste0("seed", primes))
colnames(FM.matrix) <- c(paste0("seed", primes))
average <- apply(FM.matrix, 1, mean)
FM.df <- as.data.frame(cbind(FM.matrix, average))
best.seed <- row.names (FM.df[which.max(FM.df$average),])
best.seed <- paste0("cluster.",best.seed)
output <- list(best.seed, average)
return(output)
}

####################################################################################
# 0.3 Prerequise on Df Origin                                                      #
####################################################################################

# Importation of the main database completed using bayesian inference, including NA
df.origin <- read.delim2("Data_India.txt")

# Base 10 log transformation of Air_Flows
# Replace -Inf by 0 within the dataframe
# Replace NA by 0 within MACRO_AREA_COMPACITY
df.origin$AIR_FLOWS <- log10(df.origin$AIR_FLOWS)
df.origin[ df.origin == "-Inf" ] = 0
df.origin$MACRO_AREA_COMPACITY <- replace(df.origin$MACRO_AREA_COMPACITY, 
                                          is.na(df.origin$MACRO_AREA_COMPACITY), as.numeric(0.0))

####################################################################################
# 0.4 : Data modification according to the SOM standards                           #
####################################################################################

# Removing 31 lines with NA in order to perform the SOM
df.origin.no.NA <- df.origin[rowSums(is.na(df.origin))==0, ]
row.names(df.origin.no.NA) <- df.origin.no.NA[,2]

# ID removing 
# the data-frame for analysis is now df.no.id
df.no.id <- df.origin.no.NA[,c(-1, -2)]

# default method to center and scale the columns of a dataframe
# the data-frame for analysis is now df.SOM.scale 
df.SOM.scale <- scale(df.no.id)

####################################################################################
# 0.5 : Data modification according to the superSOM standards                      #
####################################################################################

# ID removing 
# A dataframe ready for supersom
df.supersom <- df.origin
id <- df.supersom[2]
df.supersom.no.id <- df.supersom[,c(-1, -2)]

# default method to centers and scales the columns of a dataframe
df.supersom.scale <- as.data.frame(scale(df.supersom.no.id))
rm(df.supersom.no.id, df.supersom)

####################################################################################
# 0.6 Number of initializations selection                                          #
####################################################################################

# number of primes to be inserted in the vector primes
primes <- select.primes(20)

####################################################################################
#                                                                                  #
#                                   MAIN 1                                         #
#                               STANDARD SOM                                       #
#                                                                                  #
####################################################################################
####################################################################################
# 1.1 : SOM                                                                        #
####################################################################################

## 1.2.1 Self-organising maps of variables for mapping high-dimensional data to 2D grid using nb different seeds.
results.standard.SOM.list <- list()
for(i in primes)
{
  set.seed(i)
  som.seed<- som(data = df.SOM.scale, grid = somgrid(4, 4, "hexagonal"),
                 rlen = 10000, toroidal = FALSE, radius = 3)
  results.standard.SOM.list[[length(results.standard.SOM.list)+1]] <- som.seed
  print(i)
}
rm(som.seed)

names(results.standard.SOM.list) <- c(paste0("cluster.seed", primes))

## 1.2.2 Clustering results for each seed
# Combine the results for each seed in one data.frame

clusters.standard.SOM <- do.call(cbind, lapply(1:length(results.standard.SOM.list),
                                               function(i) cbind(results.standard.SOM.list[[i]][[8]])))
observation.names <- row.names(df.origin.no.NA)
clusters.standard.SOM <- as.data.frame(cbind(observation.names, clusters.standard.SOM))
names(clusters.standard.SOM) <- c("id_observations", (paste0("cluster.seed", primes)))
rm(observation.names)

####################################################################################
# 1.2 : Pick Up the Best Seed                                                      #
####################################################################################

best.seed.standard.SOM <- FM.Similarity(results.standard.SOM.list,clusters.standard.SOM)[[1]]

# Extract the best standardSOM
best.standard.SOM <- results.standard.SOM.list[[best.seed.standard.SOM]]

####################################################################################
#                                                                                  #
#                                   MAIN 2                                         #
#                         STANDARD SOM on factors + superSOM                       #
#                                                                                  #
####################################################################################
####################################################################################
# 2.1 : Data modification for Performing SOM on Factors                            #
####################################################################################

# transposing of df.scale.no.id in order to perform SOM on variables
# the data-frame for analysis is now df.scale.no.id.t
df.SOM.scale.t <- t(df.SOM.scale)

####################################################################################
# 2.2 : SOM on Variables                                                           #
####################################################################################

## 2.2.1 Self-organising maps of variables for mapping high-dimensional data to 2D grid using nb different seeds.
results.SOM.list <- list()
for(i in primes)
{
  set.seed(i)
  som.seed<- som(data = df.SOM.scale.t, grid = somgrid(4, 4, "hexagonal"),
                 rlen = 10000, toroidal = FALSE, radius = 3)
  results.SOM.list[[length(results.SOM.list)+1]] <- som.seed
  print(i)
}
rm(som.seed)
names(results.SOM.list) = c(paste0("som.seed", primes))

## 2.2.2 Clustering results for each seed
# Combine the results for each seed in one data.frame
clusters.SOM <- do.call(cbind, lapply(1:length(results.SOM.list), function(i) cbind(results.SOM.list[[i]][[8]])))
variable.names <- row.names(df.SOM.scale.t)
clusters.SOM <- as.data.frame(cbind(variable.names, clusters.SOM ))
names(clusters.SOM) <- c("id.factor", (paste0("cluster.seed", primes)))   

####################################################################################
# 2.3 : Pick Up the Best Seed                                                      #
####################################################################################

## Fowlkes-Mallows Similarity Index Calculation
best.seed.SOM <- FM.Similarity(results.SOM.list,clusters.SOM)[[1]]

####################################################################################
# 2.4 : Assign Factors                                                             #
####################################################################################

# cluster results for the best seed as a dataframe with first column as row names
rownames(clusters.SOM) <- clusters.SOM[,1]
clusters.SOM <- clusters.SOM[,-1]

# assign each row.names cluster result in a results.supersom.list of vectors
# results.SOM.list[[1]][[2]][[2]] and results.SOM.list[[1]][[2]][[3]] are the dimensions of the SOM grid which was calculated, and hence the nb of factors
Factor.lst.superSOM <- lapply(1:(results.SOM.list[[1]][[2]][[2]]*results.SOM.list[[1]][[2]][[3]]), 
                              function(i) row.names(subset(clusters.SOM , get(best.seed.SOM) ==i)))

####################################################################################
# 2.5 : SuperSOM Layers Preparation                                                #
####################################################################################

# Link the main database with the factors
# function Map('[', ...]) extracts data accordint fo Factot.lst from df.supersom.scale
# supersom.layers is a results.supersom.list of data.frames (one per factor) with the id of each spatial unit and values for the variables within the each factor
supersom.layers <- c(id, Map(`[`, list(df.supersom.scale), Factor.lst.superSOM))

#combine in a vector of character strings the names of the combined variables
Factor.titles <- sapply(Factor.lst.superSOM, paste, collapse=" ")

# Put a title for each factor
names(supersom.layers) <- c("ID",Factor.titles)
rm(Factor.titles)

# Define the weight of each layer according to the number of variables within it
supersom.w <- do.call(c,lapply(supersom.layers, NCOL))

# supersom.layers is transpormed in a results.supersom.list of matrices, as supersom function needs matrices
# lapply could also be used
supersom.layers <- lapply(supersom.layers, as.matrix)
#supersom.layers2 <- sapply(supersom.layers, as.matrix)

####################################################################################
# 2.6 : SuperSOM 3x3                                                               #
####################################################################################

results.supersom.list = list()
for(i in primes)
{
  set.seed(i)
  supersom.seed <- supersom(supersom.layers, somgrid(3, 3, "hexagonal"), 
                            rlen = 20000, alpha = c(0.05, 0.01), radius = 2,
           contin = TRUE, toroidal = FALSE, n.hood = "circular",
           weights = supersom.w/55 , maxNA.fraction = 10 , whatmap = 2:17 , keep.data = TRUE)
  results.supersom.list[[length(results.supersom.list)+1]] <- supersom.seed
  print(i)
}
names(results.supersom.list) = c(paste0("cluster.seed", primes))
rm(supersom.w)

# Combine the results for each seed in one data.frame
clusters.superSOM <- do.call(cbind, lapply(1:length(results.supersom.list),
                                           function(i) cbind(results.supersom.list[[i]][[4]])))
clusters.superSOM <- as.data.frame(cbind(supersom.layers[[1]], clusters.superSOM))
names(clusters.superSOM) <- c("id.factor", (paste0("seed", primes)))  

####################################################################################
# 2.7 : Pick Up the Best Seed                                                      #
####################################################################################

# 2.9.1 Fowlkes-Mallows Similarity Index Calculation
best.seed.superSOM <- FM.Similarity(results.supersom.list,clusters.superSOM)[[1]]

best.superSOM <- results.supersom.list[[best.seed.superSOM]]

####################################################################################
## 2.8 ANOVA By Variables                                                         #
####################################################################################

cluster.best.supersom <- best.superSOM[[4]]
results.best.superSOM.anova.ready <- cbind(df.supersom.scale, cluster.best.supersom)
anova.best.superSOM <- apply(results.best.superSOM.anova.ready, 2,
                             function(i) summary(aov(i ~ cluster.best.supersom,
                                                     data=results.best.superSOM.anova.ready)))
p.vector <- do.call(c, lapply(1:length(anova.best.superSOM),
                              function(i) cbind(anova.best.superSOM[[i]][[1]][[5]][1])))
p.vector <- p.vector [-length(p.vector)]
names(p.vector) <- variable.names
var.010significance.best.superSOM <- p.vector < 0.10
names(var.010significance.best.superSOM) <- variable.names
significant.variables <- names(subset(var.010significance.best.superSOM, var.010significance.best.superSOM == TRUE))
best.superSOM.anova.results <- cbind(p.vector,var.010significance.best.superSOM)
rm(var.010significance.best.superSOM,cluster.best.supersom,results.best.superSOM.anova.ready,
   anova.best.superSOM,p.vector)


####################################################################################
#                                                                                  #
#                                   MAIN 3                                         #
#                         STANDARD SOM on factors + superSOM + ANOVA               #
#                                                                                  #
####################################################################################
####################################################################################
## 3.1 Remove non Significant Variables in the Data                                #
####################################################################################

# remove non significant variables in SOM/Supersom data
df.SOM.scale.bis <- df.SOM.scale[,significant.variables]
df.SOM.scale.t.bis <- t(df.SOM.scale.bis)

####################################################################################
# 3.2 : Som on Factors                                                             #
####################################################################################

## 3.2.1 Self-organising maps of variables for mapping high-dimensional data to 2D grid using nb different seeds.
results.SOM.list.bis <- list()
for(i in primes)
{
  set.seed(i)
  som.seed<- som(data = df.SOM.scale.t.bis, grid = somgrid(4, 4, "hexagonal"),
                 rlen = 10000, toroidal = FALSE, radius = 3)
  results.SOM.list.bis[[length(results.SOM.list.bis)+1]] <- som.seed
  print(i)
}
rm(som.seed)

names(results.SOM.list.bis) = c(paste0("som.seed", primes))

## 3.2.2 Clustering results for each seed
# Combine the results for each seed in one data.frame
clusters.SOM.bis <- do.call(cbind, lapply(1:length(results.SOM.list.bis), 
                                          function(i) cbind(results.SOM.list.bis[[i]][[8]])))
variable.names.bis <- row.names(df.SOM.scale.t.bis)
clusters.SOM.bis <- as.data.frame(cbind(variable.names.bis, clusters.SOM.bis ))
names(clusters.SOM.bis) <- c("id.factor", (paste0("cluster.seed", primes)))    

####################################################################################
# 3.3 : Pick Up the Best Seed                                                     # 
####################################################################################

## Fowlkes-Mallows Similarity Index Calculation
best.seed.SOM.bis <- FM.Similarity(results.SOM.list.bis,clusters.SOM.bis)[[1]]

####################################################################################
# 3.4 : Assign Factors                                                            #
####################################################################################

# cluster results for the best seed as dataframe with first column as row names
rownames(clusters.SOM.bis) <- clusters.SOM.bis[,1]
clusters.SOM.bis <- clusters.SOM.bis[,-1]

# assign each row.names cluster result in a results.supersom.list of vectors
# results.SOM.list[[1]][[2]][[2]] and results.SOM.list[[1]][[2]][[3]] are the dimensions of the SOM grid which was calculated, and hence the nb of factors
Factor.lst.bis <- lapply(1:(results.SOM.list.bis[[1]][[2]][[2]]*results.SOM.list.bis[[1]][[2]][[3]]), 
                         function(i) row.names(subset(clusters.SOM.bis , get(best.seed.SOM.bis) ==i)))

####################################################################################
# 3.5 : superSOM Layer Preparation                                                 #
####################################################################################

df.supersom.scale.bis <- df.supersom.scale[,significant.variables]

# Link the main database with the factors
# function Map('[', ...]) extracts data accordint fo Factot.lst from df.supersom.scale
# supersom.layers is a results.supersom.list of data.frames (one per factor) with the id of each spatial unit and values for the variables within the each factor
supersom.layers.bis <- c(id, Map(`[`, list(df.supersom.scale.bis), Factor.lst.bis))

#combine in a vector of character strings the names of the combined variables
# Put a title for each factor
Factor.titles.bis <- sapply(Factor.lst.bis, paste, collapse=" ")
names(supersom.layers.bis) <- c("ID",Factor.titles.bis)

# Define the weight of each layer according to the number of variables within it
supersom.w.bis <- do.call(c,lapply(supersom.layers.bis, NCOL))

# supersom.layers is transpormed in a results.supersom.list of matrices, as supersom function needs matrices
supersom.layers.bis <- lapply(supersom.layers.bis, as.matrix)

####################################################################################
# 3.6 : SUPERSOM 3x3                                                              #
####################################################################################

## 3.6.1 superSOM
# The following primes can not be optimized and are thus replaced
primes[6] <- 79
primes[10]<- 89
primes[11]<- 101
primes[12]<- 103

results.supersom.list.bis = list()
for(i in primes)
{
  set.seed(i)
  supersom.seed <- supersom(supersom.layers.bis, somgrid(3, 3, "hexagonal"), rlen = 20000, 
                            alpha = c(0.05, 0.01), radius = 2,
                            contin = TRUE, toroidal = FALSE, n.hood = "circular",
                            weights = supersom.w.bis/47 , maxNA.fraction = 10 , whatmap = 2:17 , keep.data = TRUE)
  results.supersom.list.bis[[length(results.supersom.list.bis)+1]] <- supersom.seed
  print(i)
}
names(results.supersom.list.bis) = c(paste0("cluster.seed", primes))

## 3.6.2 combine the learning steps for each seed in one matriX
# Combine the results for each seed in one data.frame
clusters.superSOM.bis <- do.call(cbind, lapply(1:length(results.supersom.list.bis),
                                               function(i) cbind(results.supersom.list.bis[[i]][[4]])))
clusters.superSOM.bis <- as.data.frame(cbind(supersom.layers.bis[[1]], clusters.superSOM.bis))
names(clusters.superSOM.bis) <- c("id.factor", (paste0("seed", primes)))  

####################################################################################
# 3.7 : Pick Up the Best Seed                                                      #
####################################################################################

## Fowlkes-Mallows Similarity Index Calculation
best.seed.superSOM.bis <- FM.Similarity(results.supersom.list.bis,clusters.superSOM.bis)[[1]]

best.superSOM.bis <- results.supersom.list.bis[[best.seed.superSOM.bis]]

####################################################################################
## 3.8 Anova By variables                                                         #
####################################################################################

cluster.best.supersom.bis <- best.superSOM.bis[[4]]
results.best.superSOM.bis.anova.ready <- cbind(df.supersom.scale.bis, cluster.best.supersom.bis)
anova.best.superSOM.bis <- apply(results.best.superSOM.bis.anova.ready, 2,
                             function(i) summary(aov(i ~ cluster.best.supersom.bis,
                                                     data=results.best.superSOM.bis.anova.ready)))
p.vector <- do.call(c, lapply(1:length(anova.best.superSOM.bis),
                              function(i) cbind(anova.best.superSOM.bis[[i]][[1]][[5]][1])))
p.vector <- p.vector [-length(p.vector)]
names(p.vector) <- variable.names.bis
var.010significance.best.superSOM.bis <- p.vector < 0.10
names(var.010significance.best.superSOM.bis) <- variable.names.bis
significant.variables.bis <- names(subset(var.010significance.best.superSOM.bis,
                                          var.010significance.best.superSOM.bis == TRUE))
best.superSOM.anova.results.bis <- cbind(p.vector,var.010significance.best.superSOM.bis)

rm(var.010significance.best.superSOM.bis,cluster.best.supersom.bis,results.best.superSOM.bis.anova.ready,
   anova.best.superSOM.bis,p.vector)

####################################################################################
#                                                                                  #
#                                   APPENDIX                                       #
#                         Save the Data-Frame Results                              #
#                                                                                  #
####################################################################################

main.directory <- getwd()
dir.create("Dataframe_Results")
setwd(paste0(main.directory,'/Dataframe_Results'))

## Save the standardSOM
write.csv(cbind(best.standard.SOM[[1]], best.standard.SOM[[8]], best.standard.SOM[[9]]),
                "results.standard.SOM.csv")

## Save the SOM + superSOM 
write.csv(cbind(best.superSOM[[1]][[1]], df.supersom.scale, best.superSOM[[4]], best.superSOM[[5]]),
          "results.supersom.csv")

## Save the SOM + superSOM + Anova
write.csv(cbind(best.superSOM.bis[[1]][[1]], df.supersom.scale.bis, best.superSOM.bis[[4]], 
                best.superSOM.bis[[5]]), "results.supersom.bis.csv")

setwd(main.directory)


