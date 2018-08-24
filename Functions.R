#######################################
#This File implements funtions used in Identifying module biomarker in type 
#2 diabetes mellitus by discriminative area of functional activity in BMC Bioinfomatics
#######################################



#Function act_vec
#get the avtivity vector
#Dataset is a matrix&data.frame whose rows denoted the genes
#module is a vector contains the module's genes
act_vec<-function(Dataset,module){
  k<-sqrt(length(module))
  apply(Dataset,2,sum)->S
  return(S/k)
}

#Calculate the Discriminative area of Two normal Distribution#
#u1 is the mean of Distribution1 and sig1 is the sd of Distribution1
#u1 < u2
Dis_Area<-function(u1,sig1,u2,sig2){
  c_down<-sig1^2 - sig2^2
  c_upper<-u2 * sig1^2 - sig2*(u1*sig2 + sig1*sqrt((u1-u2)^2 + 2*(sig1^2 -sig2^2)*log(sig1/sig2)))
  c<-(c_upper / c_down)
  Area<- 1 - pnorm(c,mean=u1,sd = sig1) + pnorm(c,mean=u2,sd=sig2)
  return(Area)
}


##Greedy Strategy###
#need to complete
###################

#

#

