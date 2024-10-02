##########################################################################
##########################################################################
### Codes to run GAMPI
### Input:
# Y = ; X = ;
### Y: Primary variables of interest matrix: n * p
### X: Instrumental variable matrix: n * q


### Some other user options:
### method: method for tuning parameter selection; options include "EBIC" (extended Bayesian Information Criterion) and "CV" (cross-validation); the default method is EBIC
### family: a character string or a string vector of data types of primary/outcome variables Y, including "gaussian", "binomial", and "poisson"

### Main function    
# result = GAMPI_DRI(X,Y,method="EBIC",family="binomial")
### Output:
### result$causal_relation: estimatd DAG/causal relationships of primary variables
### result$interv_relation: estimated intervention relationships




##########################################################################
##########################################################################
### A motivating example with codes
source("GAMPI_functions.R")

### Simulate data
library(MASS)
library(igraph)


set.seed(1)
n = 500; p = 60; q = 60; chain_length = 4 # chain graphs, each chain has length 4
number_of_chains = p / chain_length

### Simulate Chain graph
U = matrix(rep(0,p*p),ncol = p)
for (i in 1:(p-1)){
  U[i,i+1] = 1;
}

for (k in 1:(number_of_chains-1)){
  idx = k*chain_length
  U[idx,idx+1] = 0
}

# Plot true DAG
g0 <- graph_from_adjacency_matrix(abs(U), mode="directed") # weights of adjacency matrix must be nonnegative
coords <- layout.circle(g0)
plot(g0,main="True DAG",layout = coords,vertex.size = 4,edge.arrow.size = 0.4,edge.color = 'black',edge.width = .5,vertex.label="",vertex.color="darkred")


### Intervention matrix
W = matrix(rep(0,q*p),ncol = p)
for (i in 1:p){
  W[i,i] = 1;
}

### Intervention variables
X = matrix(rnorm(n*q),nrow = q, ncol = n)
WtX = t(W) %*% X

### Primary variables
Y = matrix(rep(0,p*n),ncol = n)

### Parameters
beta_causal_parent_child = 2.5 ### Parent-child relationship
beta_interv_root = 5; beta_interv_nonroot = 3

### Confounders
Sigma <- matrix(1,ncol=p,nrow=p)
Sigma <- 0.95 * Sigma; diag(Sigma) = 1
confounder = mvrnorm(n = n, rep(0, p), Sigma, empirical = TRUE) # multivariate normal, features can be correlated
confounder[,seq(1,p,by = chain_length)] = beta_interv_root * confounder[,seq(1,p,by = chain_length)]


### Simulate primary variables Y
for (k in 1:number_of_chains){
  
  ### root node
  idx = (k-1)*chain_length + 1
  Y[idx,] = exp(beta_interv_root*WtX[idx,] + confounder[,idx])/(1 + exp(beta_interv_root*WtX[idx,] + confounder[,idx]))
  Y[idx,] = rbinom(n,1,Y[idx,])
  
  ### non-root node
  for (j in 2:chain_length){
    idx = (k-1)*chain_length + j
    beta1 = beta_causal_parent_child; beta2 = beta_interv_nonroot;
    
    Y[idx,] = exp(beta1*Y[idx-1,] + beta2*WtX[idx,] + confounder[,idx] - beta1/2)/ (1 + exp(beta1*Y[idx-1,] + beta2*WtX[idx,] + confounder[,idx] - beta1/2))
    Y[idx,] = rbinom(n,1,Y[idx,])
  }
}

X = t(X); Y = t(Y);
U = U * beta1


##########################################################################
### Run GAMPI (the proposed method is GAMPI_DRI)
# X = as.matrix(X); Y = as.matrix(Y) # for a dataframe
family_list = "binomial"
result = GAMPI_DRI(X,Y,method="EBIC",family=family_list)
U_est.ebic = result$causal_relation


##########################################################################
### Visualize DAG estimated by GAMPI-DRI
adjm = U_est.ebic
adjm[adjm!=0] = 1

if (!is.null(colnames(Y))){
  colnames(adjm) = colnames(Y); rownames(adjm) = colnames(Y)
}

g <- graph_from_adjacency_matrix(adjm, mode="directed") # weights of adjacency matrix must be nonnegative

### Plot code
coords <- layout.circle(g)
label.color <- c(rep("black",ncol(Y)))
plot(g,main="Estimated DAG by GAMPI-DRI",layout = coords,vertex.size = 4,edge.arrow.size = 0.4,edge.color = 'black',edge.width = .5,vertex.label="",vertex.color="darkred")

## Apply labels manually
# Specify x and y coordinates of labels, adjust outward as desired
x = coords[,1]*1.17
y = coords[,2]*1.17

# Create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(coords[,1]/coords[,2]))*(180/pi) < 0,  90 + atan(-(coords[,1]/coords[,2]))*(180/pi), 270 + atan(-coords[,1]/coords[,2])*(180/pi))

# Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=V(g)$name[i], adj=NULL, pos=NULL, cex=.7, col=label.color[i], srt=angle[i], xpd=T)
}


##########################################################################
### Comparison
# Method without adjusting for confounders
result_no_deconf = GAMPI_no_deconf(X,Y,method="EBIC",family=family_list)
U_no_deconf.ebic = result_no_deconf$causal_relation

# Evaluate and compare the performance in terms of F-score
paste("Fscore of GAMPI-DRI:",output_metrics(ifelse(U!=0,1,0),ifelse(U_est.ebic!=0,1,0))$F_score)
paste("Fscore of GAMPI without adjusting for confounders:",output_metrics(ifelse(U!=0,1,0),ifelse(U_no_deconf.ebic!=0,1,0))$F_score)
### GAMPI adjusting for confounders outperforms one that does not adjust for confounders


# To illustrate, the DAG estimated by the method without adjusting for confounders includes many FPs
adjm = U_no_deconf.ebic
adjm[adjm!=0] = 1
if (!is.null(colnames(Y))){
  colnames(adjm) = colnames(Y); rownames(adjm) = colnames(Y)
}
g2 <- graph_from_adjacency_matrix(adjm, mode="directed") # weights of adjacency matrix must be nonnegative
### Plot code
coords <- layout.circle(g)
label.color <- c(rep("black",ncol(Y)))
plot(g2,main="Estimated DAG without adjusting for confounders",layout = coords,vertex.size = 4,edge.arrow.size = 0.4,edge.color = 'black',edge.width = .5,vertex.label="",vertex.color="darkred")
## Apply labels manually
# Specify x and y coordinates of labels, adjust outward as desired
x = coords[,1]*1.17
y = coords[,2]*1.17
# Create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(coords[,1]/coords[,2]))*(180/pi) < 0,  90 + atan(-(coords[,1]/coords[,2]))*(180/pi), 270 + atan(-coords[,1]/coords[,2])*(180/pi))
# Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=V(g)$name[i], adj=NULL, pos=NULL, cex=.7, col=label.color[i], srt=angle[i], xpd=T)
}
### includes many false positives
