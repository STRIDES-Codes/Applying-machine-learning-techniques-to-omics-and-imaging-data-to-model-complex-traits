
### ############################################################################
################################################################################
#                        Library Section                                       #
################################################################################
library(e1071)
library(class)
library(digest)
library(RKUM)
library(dplyr)
library(kernlab)



############################################################ Classification##############################


cv.folds <- function(n, nfolds)
## Randomly split the n samples into folds
## Returns a list of nfolds lists of indices, each corresponding to a fold
 {
 	return(split(sample(n),rep(1:nfolds,length=n)))
 }



   ###############################################################################
#       Function for Prediction error  for cross-validation                   #
###############################################################################
SVM2F<-function(X,Y)
 {
          folds=2
          tn<-nrow(X)
          em2<-rep(0,folds)
          f<- cv.folds(tn,folds)
          for(i in  seq(folds))
             {
            PE1<-ksvm(X[-f[[i]],],Y[-f[[i]]], kernel="vanilladot",type="C-svc",C=0.1)
             rr<-predict(PE1,X[f[[i]],])
             rc <-length(which (rr==Y[f[[i]]]))
            em2[i] <- (1- rc/length(rr))*100
             }
        #return(em2)
          list(Mean=mean(em2), stn = sd(em2))

 }


KNN2F<-function(X,Y)
 {
          folds=2
          tn<-nrow(X)
          em2<-rep(0,folds)
          f<- cv.folds(tn,folds)
          for(i in  seq(folds))
             {
             rr<-knn( CXCCA[-f[[i]],], CXCCA[f[[1]],], Y[-f[[i]]], k =2, prob=TRUE)
             rc <-length(which (rr==Y[f[[i]]]))
            em2[i] <- (1- rc/length(rr))*100
             }
            #return(em2)
            list(Mean=mean(em2), stn = sd(em2))
 }



 X <- CXCCA # set of selected features
 Y <-  factor(c(rep("1", 11), rep("2", 53)))
 

 c1 <- factor(c(rep(1, 22), rep(2, 40)))

 AllSVM<- SVM2F(CX, c1)
 AllKNN<- KNN2F(CX, c1)
 

