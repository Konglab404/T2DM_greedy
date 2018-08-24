#######################################
#This File implements funtions used in Identifying module biomarker in type 
#2 diabetes mellitus by discriminative area of functional activity in BMC Bioinfomatics
#######################################
library(stringr)
library(magrittr)

#Function act_vec
#get the avtivity vector
#Dataset is a numeric matrix&data.frame whose rows denoted the genes
#rownames(Dataset) should be the gene names
#module is a vector contains the module's genes
act_vec<-function(Dataset,module){
  k<-sqrt(length(module))
  apply(Dataset[module,],2,sum)->S
  return(S/k)
}

#Function Dis_Area
#Calculate the Discriminative area of Two normal Distribution#
#u1 is the mean of Distribution1 and sig1 is the sd of Distribution1
#u1 < u2
#return value is always > 0
Dis_Area<-function(u1,sig1,u2,sig2){
  if(u1>u2){temp<-sig1;sig1<-sig2;sig2<-temp;temp<-u1;u1<-u2;u2<-temp}
  c_down<-sig1^2 - sig2^2
  c_upper<-u2 * sig1^2 - sig2*(u1*sig2 + sig1*sqrt((u1-u2)^2 + 2*(sig1^2 -sig2^2)*log(sig1/sig2)))
  c<-(c_upper / c_down)
  Area<- 1 - pnorm(c,mean=u1,sd = sig1) + pnorm(c,mean=u2,sd=sig2)
  return(abs(Area))
}


##Greedy Strategy###
#need to complete
###################
#Function Module_Add
#module is a vector contains the genes have been in module
#PPI is a 2-column dataframe/matrix denotes the interactions
#case_Data & contorl_Data have the same structure as Funtion act_vec's input Dataset 
Module_Add<-function(module,PPI,case_Data,control_Data){
  node_candi<-setdiff(PPI_test[PPI_test$Node1 %in% module,]$Node2,module) #list all the neighbor of module
  nowScore<-Dis_Area(act_vec(Case,module)%>% mean,act_vec(Case,module)%>%sd,
                     act_vec(Control,module)%>% mean, act_vec(Control,module) %>% sd)#calculate current score
    
  score<-c() #vector contains the node_candi's scores
  for(i in 1:length(node_candi)){
    newMod<-c(module,node_candi[i])
    score<-c(Dis_Area(act_vec(Case,newMod)%>% mean,act_vec(Case,newMod)%>%sd,
                    act_vec(Control,newMod)%>% mean, act_vec(Control,newMod) %>% sd),score)
  } #Calculate node_candi's scores
  names(score)<-node_candi
  if(min(score)<nowScore){return(c(module,score[order(score,decreasing = F)][1] %>% names))}
  else{return(module)} #decide wheater add the new node into current module
}

#

#

