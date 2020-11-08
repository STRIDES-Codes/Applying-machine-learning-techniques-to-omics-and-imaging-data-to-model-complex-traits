###### Library Setrion#####################
library(rcc)
library(plsdepot)
library(RGCCA)

data(Russett)
head(Russett)
# Three input matrix
X_agric = as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("inst", "ecks", "death","demostab", "dictator")])
A = list(X_agric, X_ind, X_polit)
sapply(A, head)

A = lapply(A, function(x) scale2(x, bias = TRUE))
C = matrix(c(0, 0, 1,0, 0, 1,1, 1, 0), 3, 3)

rgcca_B_factorial = rgcca(A, C, tau = rep(0, 3), scheme = "factorial",scale = T, verbose = TRUE)
rgcca_B_factorial$a # weight vectors

Y = rgcca_B_factorial$Y
lapply(Y, head)
rgcca_B_factorial$AVE

# an ?optimal? shrinkage parameter
rgcca_optimal_factorial = rgcca(A, C, tau = "optimal", scheme = "factorial",
scale = FALSE, verbose = FALSE)

rgcca_optimal_factorial$tau


####################### IF based approach#######################

 library(RKUM)
 CCAIF1 <- ifcca(Datag1, Datag2, 0.00001, 2,2)
 




