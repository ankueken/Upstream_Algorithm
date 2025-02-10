args=(commandArgs(TRUE))
library(R.matlab)
library(Matrix)
data=readMat(args)

S=data$model[,,1][[1]]
rev=rep(0,dim(S)[2])
dim(S)

source('deficiency_sparse')
library(igraph)

D=deficiency_sparse(S,rev)

print('deficiency is calculated')

AI = which(D$A!=0, arr.ind=T)
AI = cbind(AI, D$A[AI]) 

YI = which(D$Y!=0, arr.ind=T)
YI = cbind(YI, D$Y[YI]) 

complexes=rownames(D$A)

write.table(AI,file=paste(substring(args,1,nchar(args)-4),'_A.dat',sep=''),row.names=FALSE,col.names=FALSE)
write.table(YI,file=paste(substring(args,1,nchar(args)-4),'_Y.dat',sep=''),row.names=FALSE,col.names=FALSE)
writeMat(paste(substring(args,1,nchar(args)-4),'_complexes.mat',sep=''),complexes=complexes)
