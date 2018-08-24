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
#return a New module which have a smallest Dis_Area
Module_Add<-function(module,PPI,case_Data,control_Data){
  node_candi<-setdiff(PPI[PPI$Node1 %in% module,]$Node2 %>% as.character(),module) #list all the neighbor of module
  nowScore<-Dis_Area(act_vec(case_Data,module)%>% mean,act_vec(case_Data,module)%>%sd,
                     act_vec(control_Data,module)%>% mean, act_vec(control_Data,module) %>% sd)#calculate current score
  if(length(node_candi) == 0){return(list(module,nowScore))}  
  score<-c() #vector contains the node_candi's scores
  for(i in 1:length(node_candi)){
    newMod<-c(module,node_candi[i])
    score<-c(score,Dis_Area(act_vec(case_Data,newMod)%>% mean,act_vec(case_Data,newMod)%>%sd,
                    act_vec(control_Data,newMod)%>% mean, act_vec(control_Data,newMod) %>% sd))
  } #Calculate node_candi's scores
  names(score)<-node_candi
  if(min(score)<nowScore){return(list(c(module,score[order(score,decreasing = F)][1] %>% names),min(score)))}
  else{return(list(module,nowScore))} #decide wheater add the new node into current module
}

#Function seed_Greedy
#get a network from a seed by greedy algothrm
#seed is a character which is the seed to initiate search job
#thr is the threhold which determines whether terminate the search
seed_Greedy<-function(seed,PPI,case_Data,control_Data,thr){
  flag = T
  modl<-c(seed,PPI[PPI$Node1 == seed,"Node2"] %>% as.character())
  nowScore<-Dis_Area(act_vec(case_Data,modl)%>% mean,act_vec(case_Data,modl)%>%sd,
                     act_vec(control_Data,modl)%>% mean, act_vec(control_Data,modl) %>% sd)
  while(flag == T){
    res<-Module_Add(modl,PPI,case_Data,control_Data)
    modl<-res[[1]];score<-res[[2]]
    print(score)
    if((nowScore-score)<thr){break}
    nowScore<-score
  }
  return(modl)
}
########################################################

#

