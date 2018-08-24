#This file contains an Example workflow
#Just Run it
source("./Functions.R")

#Generate Case & Control expression Dataset
Case<-matrix(data=rnorm(n = 1200,mean = 100,sd=30),ncol=8)
colnames(Case)<-paste0("sample",1:8);rownames(Case)<-paste0("gene",1:150)

Control<-matrix(data=rnorm(n = 1050,mean = 110,sd=24),ncol=7)
colnames(Control)<-paste0("sample",1:7);rownames(Control)<-paste0("gene",1:150)

#Generate Example PPI data
PPI_test<-data.frame(Node1=sample(rownames(Case),size = 300,replace = T),Node2=sample(rownames(Case),300,T))

#Assume now we have get a DEG gene17

seed_Greedy("gene17",PPI_test,Case,Control,thr = 0.01)



