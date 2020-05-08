#define the function of U and Q95 test
u20_q95_test=function(df,chr,windowsize,stepsize,popa,popb,popc,afa,afb,afc){
  Chr <- vector()
  Start <- vector()
  End <- vector()
  N_U20 <- vector()
  qts <- vector()
  N_sites<-vector()
  tardf=df[which( df$Chr == chr),] 
  n_window = ceiling(max(tardf$Pos)/windowsize) * (windowsize/stepsize)
  
  for(i in 1:n_window){
    window<-tardf[which(tardf$Pos >=1+(i-1)*stepsize & tardf$Pos <= windowsize+(i-1)*stepsize),]
    u20sites<-window[which(window[[popa]] <= afa & window[[popc]]>=afc & window[[popb]] >= afb),]
    
    Chr<-c(Chr,chr)
    Start<-c(Start,1+(i-1)*stepsize)
    End<-c(End,windowsize+(i-1)*stepsize)
    N_U20<-c(N_U20,length(u20sites$Pos))
    qts<-c(qts,quantile(u20sites[[popb]],probs=.95))
    N_sites<-c(N_sites,length(window$Pos))
  }
  odf <- cbind.data.frame(Chr, Start,End,N_U20,qts,N_sites)
  colnames(odf) <- c("Chr","Start","End","N_U20","Q95","N_sites")
  return(odf)
}

colnames<-c("Chr","Pos","N_alleles","N_chr","Ref","Alt_Major","Alt_Minor","Alt_Target") #VAF data should be formatted as such
VAFs<- read.table("path_to_VAF",header = F,skip = 1,col.names=colnames)

#run the test specify the Chr, window size, step size, VAF of the populations, and threshold of VAFs accordingly

output <- u20_q95_test(VAFs,"CHR_1",100000,20000,"ALT_Major","ALT_Target","ALT_Minor",0.01,0.2,1) 

t1<-quantile(output$N_U20,probs=0.99) #find the 99% quantile of U20 test
t2<-quantile(output$Q95,probs=0.99) #find the 99% quantile of Q95 test

output[which(output$N_U20 >=as.numeric(t1) & output$Q95 >=as.numeric(t2)),] #select those candidate regions

#one can also pull out introgressed variants that likely to be under selection
AI_sites<-VAFs[which(VAFs$Alt_Target >= as.numeric(t2) & VAFs$Alt_Major<=0.01 & VAFs$Alt_Minor ==1 ),]

#Then we can also utilize the VAFs to classify the SNPs of the target population based on the VAFs of Major and Minor parental populations
#based on the threshold defined
Shared_polymorphisms_of_three_groups<-VAFs[which(VAFs$Alt_Minor >0 & VAFs$Alt_Target > 0 & VAFs$Alt_Major>0),]

Shared_between_Target_and_Minor<-VAFs[which(VAFs$Alt_Minor >0 & AF_matrix$Alt_Target > 0 & AF_matrix$Alt_Major == 0),]

Shared_between_Target_and_Major<-VAFs[which(VAFs$Alt_Minor == 0 & AF_matrix$Alt_Target > 0 & AF_matrix$Alt_Major > 0),]

Shared_between_Minor_and_Major<-VAFs[which(VAFs$Alt_Minor > 0 & AF_matrix$Alt_Target == 0 & AF_matrix$Alt_Major > 0),]

Uniquely_in_Target<-VAFs[which(VAFs$Alt_Minor == 0 & AF_matrix$Alt_Target > 0 & AF_matrix$Alt_Major == 0),]
Uniquely_in_Major<-VAFs[which(VAFs$Alt_Minor == 0 & AF_matrix$Alt_Target == 0 & AF_matrix$Alt_Major > 0),]
Uniquely_in_Minor<-VAFs[which(VAFs$Alt_Minor > 0 & AF_matrix$Alt_Target == 0 & AF_matrix$Alt_Major == 0),]

#Then we figure out SNPs under selection based on the Q95 threshold and 99% quantile VAF of Major_Target shared variants

t3<-quantile(Shared_between_Target_and_Major$Alt_Target, probs=0.99)# get the threshold of selection for Major Target shared polymorphisms in the Target population

Selection_on_Target_Minor_shared<-VAFs[which(VAFs$Alt_Minor > 0 & AF_matrix$Alt_Target > as.numeric(t2) & AF_matrix$Alt_Major == 0),]

Selection_on_Target_Major_shared<-VAFs[which(VAFs$Alt_Minor == 0 & AF_matrix$Alt_Target > as.numeric(t3) & AF_matrix$Alt_Major > 0),]

#draw a Venn diagram
library(VennDiagram)

draw.triple.venn(area1=nrow(Uniquely_in_Minor), area3=nrow(Uniquely_in_Major), area2=nrow(Uniquely_in_Target), 
                 n13=nrow(Shared_between_Target_and_Minor) + nrow(Shared_polymorphisms_of_three_groups), 
                 n23=nrow(Shared_between_Target_and_Major) + nrow(Shared_polymorphisms_of_three_groups), 
                 n12=nrow(Shared_between_Minor_and_Major) + nrow(Shared_polymorphisms_of_three_groups), 
                 n123=nrow(Shared_polymorphisms_of_three_groups),
                 category=c("Minor","Major","Target"),print.mode=c("percent","raw"))
