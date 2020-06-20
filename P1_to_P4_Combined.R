library(tidyverse)
library(FossilSim)

my_foss <- read.table(".\\WrightLab\\shreya_foss.log", head = TRUE)


n.ages <- function(tree){
  
  depth = ape::node.depth.edgelength(tree)
  node.ages = max(depth) - depth
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))
  
  # adding possible offset if tree fully extinct
  if(!is.null(tree$root.time)) node.ages = node.ages + tree$root.time - max(node.ages)
  
  return(node.ages)
}



#sim.fbd.rateshift.taxa, but it outputs the fossil in a list with the SAtree object. It also adds an (h) column to the fossil dataframe so that the fossil count functions will work.

sfbdrt <- function (n, numbsim, lambda, mu, psi, times, complete = FALSE) 
{
  if (length(lambda) != length(times)) 
    stop("Length mismatch between rate shift times and birth rates")
  if (length(mu) != length(times)) 
    stop("Length mismatch between rate shift times and death rates")
  if (length(psi) != length(times)) 
    stop("Length mismatch between rate shift times and sampling rates")
  trees = TreeSim::sim.rateshift.taxa(n, numbsim, lambda, 
                                      mu, rep(1, length(times)), times, complete = TRUE)
  for (i in 1:length(trees)) {
    t = trees[[i]]
    origin = max(n.ages(t)) + t$root.edge
    horizons = c(times, origin)
    f <- sim.fossils.intervals(tree = t, interval.ages = horizons, 
                               rates = psi)
    tree = SAtree.from.fossils(t, f)
    node.ages = n.ages(tree)
    if (complete == FALSE) {
      fossil.tips = is.extinct(tree, tol = 1e-06)
      sa.tips = tree$tip.label[tree$edge[, 2][(tree$edge[, 
                                                         2] %in% 1:length(tree$tip.label)) & (tree$edge.length == 
                                                                                                0)]]
      unsampled.tips = fossil.tips[!(fossil.tips %in% 
                                       sa.tips)]
      tree = ape::drop.tip(tree, unsampled.tips)
      node.ages = n.ages(tree)
    }
    trees[[i]] = tree
    trees[[i]]$root.edge = origin - max(node.ages)
    trees[[i]] = SAtree(trees[[i]], complete)
  }
  
  f$h <- (f$hmin + f$hmax)/2
  
  return(list(f, trees))
}




#Run sfbdrt with parameter specifications


SFB <- function(n_extant, n_trees, df,  mu, lambda, psi, row, mers){
  
  num <- max(c(length(mu), length(lambda), length(psi)))
  
  if(length(mu) == 1){
    m1 <- rep(df[row, mu], times = num)} else{
      m1 <- c()
      for(i in mu){m1 <- append(m1, df[row, i])}}
  
  if(length(lambda) == 1){
    l1 <- rep(df[row, lambda], times = num)} else{
      l1 <- c()
      for(i in lambda){l1 <- append(l1, df[row, i])}}
  
  if(length(psi) == 1){
    p1 <- rep(df[row, psi], times = num)} else{
      p1 <- c()
      for(i in psi){p1 <- append(p1, df[row, i])}}
  
  modelledfoss <- sfbdrt(n_extant, n_trees, l1, m1, p1, times = mers, complete = TRUE)
  return(modelledfoss)
}


#Example list generation by SFB. 

SampList <- SFB(666, 1, my_foss, "mu", "lambda", c("psi.1.", "psi.2."), row = 2033, c(0, 61))
View(SampList)

#Summarizes the generated list. Interval ages must be specified in a vector. Output is also a vector with the layout: tree.origin(output), count.fossils(output), count.fossils.binned(output). I just used a vector to make manipulating the numbers simpler.

SummPhylo <- function(phylolist, interval.ages){
  tree_origin <- tree.max(phylolist[[2]][[1]])
  num_fossils <- count.fossils(phylolist[[1]])
  num_fossils_binned <- count.fossils.binned(phylolist[[1]], interval.ages)
  return(c(tree_origin, num_fossils, num_fossils_binned))
}

#Function to run multiple simulations

Summ_Samples <- function(data, extant, mus, lambdas, psis, intervals, by = FALSE, path = c()){
  
  num <- length(intervals) + 2 #length of time bins plus 2 (the mu and lambda)
  
  summlist <- list()
  for(i in c(1:num)){summlist[i + 1] <- c()}
  
  if(by == FALSE){
    for(i in c(1:nrow(data))){
      
      tempTree <- SFB(n_extant = extant, n_trees = 1, df = data, mu = mus, lambda = lambdas, psi = psis, row = i, mers = intervals)
      tree_stats <- SummPhylo(tempTree, intervals)
      
      for(i in c(1:length(tree_stats))){
        summlist[[i]] <- append(summlist[[i]], tree_stats[i])
      }
    }
  }
  
  else{
    for(i in which(c(1:nrow(data) %% by == 0))){
      
      tempTree <- SFB(n_extant = extant, n_trees = 1, df = data, mu = mus, lambda = lambdas, psi = psis, row = i, mers = intervals)
      tree_stats <- SummPhylo(tempTree, intervals)
      
      for(i in c(1:length(tree_stats))){
        summlist[[i]] <- append(summlist[[i]], tree_stats[i]) #add to summlist the tree_stats value
      }
    }
  }
  
  summlist <- as.data.frame(summlist)
  
  for(i in c(1:ncol(summlist))){
    names(summlist)[1] <- "OriginTime"
    names(summlist)[2] <- "Total_Foss_Count"
    
    if(i > 2){
      names(summlist)[i] <- paste("Int", i - 2, "Fossils" ,sep = "_")
    }
  }
  
  if(length(path) != 0){
    write.csv(summlist, file = path)
  } 
  
  return(summlist)
}


#Run Summ_Samples to simulate about 100 trees
##multiple_sim_df <- Summ_Samples(my_foss, extant = 666, "mu", "lambda", psis = c("psi.1.", "psi.2."), intervals = c(0,61), by = 20, path = c("./WrightLab/output/simulation1.csv"))



#Histogram of Origin Time
origin_plot <- ggplot(data = multiple_sim_df, aes(x = OriginTime)) + 
  geom_histogram(binwidth = 5, color = "darkblue", fill = "lightblue") + 
  geom_vline(xintercept = 131, color = "red") +
  scale_x_reverse() + #seq in the -x axis
  labs(title = "Origin Time Histogram", x = "Million of Years Ago", y= "Number of Samples") + #axis labels
  theme(plot.title = element_text(hjust = 0.5)) #centers title

#Histogram of Fossil Count
(total_fossils_plot <- ggplot(data = multiple_sim_df, aes(x = Total_Foss_Count)) + 
  geom_histogram(binwidth = 4, color = "darkblue", fill = "lightpink") + 
  labs(title = "Total Fossil Count Histogram", x = "Fossils Count", y= "Number of Sample") + #axis labels
  theme(plot.title = element_text(hjust = 0.5)) #centers title
)
