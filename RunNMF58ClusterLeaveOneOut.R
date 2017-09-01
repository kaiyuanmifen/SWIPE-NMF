#!/usr/bin/Rscript

#this tries to parallel the works 
#this code is trying to exclude one dataset to get contribution of dataset to results
#This the code to run pNMF
source('NMFfunctions20.R')#load the functions 
library(Matrix)
library(methods)
#VecArg<-commandArgs(trailingOnly = T)
#CellLine<-as.character(VecArg[1])
#Chr<-as.character(VecArg[2])#Only focus on a chr for now 

load('LeaveOneOutParametersetPredefined.Rdata')
task.id <-as.integer(Sys.getenv("SGE_TASK_ID"))
#task.id <-1000
print(task.id)
Address='/broad/compbio/Dianbo/'
CellLine<-as.character(ParameterSet$CellType[task.id])
Chr<-as.character(ParameterSet$Chr[task.id])
NumberOfRepeat=as.integer(ParameterSet$Repeats[task.id])
#CellLine<-"E003"
#Chr<-'chr22'
#load(paste0("~/BoxSync/research/MIT/record/ProducedData/NMFInputMats/",CellLine,'_',Chr,'.Rdata'))
#NMFinputMats2 indicate only 4 data types 
load(paste0(Address,'datasets/ProducedData/NMFInputMats/',CellLine,'_',Chr,'.Rdata'))


#Kvec<-as.numeric(VecArg[3])
Kvec<-ParameterSet$K[task.id]
print(Kvec)


#select which information to include 
# InfortoInclude<-c(1,2,3,4)
# VecNames<-names(NMFInputMats$R)
# VecR<-unique(c(mapply(VecNames,FUN = function(x){as.numeric(strsplit(x,'')[[1]][c(2,3)])})))
# VecR<-VecR[!VecR%in%InfortoInclude]
# Rwhich<-which(!grepl(VecNames,pattern = VecR))


head(summary(NMFInputMats$R$R12))

InforExclude<-ParameterSet$LeaveOut[task.id]
print(paste('Excluding_',InforExclude))
if((length(InforExclude)>0)&InforExclude!=0){
NMFInputMats$DimensionOfEachSet<- NMFInputMats$DimensionOfEachSet[-InforExclude]       
NMFInputMats$Rvec<-NMFInputMats$Rvec[-InforExclude] 
NMFInputMats$Theta<-NMFInputMats$Theta[-InforExclude]
VecNames<-names(NMFInputMats$R)
Rwhich<-which(!grepl(VecNames,pattern =InforExclude))
NMFInputMats$R<-NMFInputMats$R[Rwhich]
}


#Slideing window 
#ParameterSet$WinStart[task.id]
Start=ParameterSet$WinStart[task.id]
End=ParameterSet$WinEnd[task.id]

NMFInputMats<-SlideWindow(NMFInput = NMFInputMats,Start = Start,End =End )
print('Sliding window made')

#
if (all(unlist(lapply(NMFInputMats$R,class))=="dgCMatrix")&(all(unlist(lapply(NMFInputMats$R,dim))>0))){#NMF will only be performed in region where all input data have some at least 1 link

for (Nrepeats in 1:NumberOfRepeat){

#initialization 


names(NMFInputMats)

G<-InitalizeG(k = Kvec,DimensionOfEachSet = NMFInputMats$DimensionOfEachSet,
              Rvec = NMFInputMats$Rvec,AcolN = 50 )

lapply(G,dim)
lapply(G,function(x){sum(x>0)})
print('G initialized')
#run the pnmf other parameters 

Output<-NMFDianbo(R =NMFInputMats$R,Theta = NMFInputMats$Theta,MaxTime =10*60,
                  G = G,TargetMatrixName = 'R12',MaxInteraction = 200,
                  StopThreshold = 1e-2,DimensionOfEachSet = NMFInputMats$DimensionOfEachSet)

#save(Output,file = paste0(CellLine,'_',Chr,'_Start=',Start,'_end=',End,'_PNMF_k=',signif(Kvec,digits = 3),'_','Remove=',paste(InforExclude,collapse = ''),'.Rdata'))


save(Output,file = paste0(Address,"datasets/ProducedData/LeaveOneOutNMFoutput/",CellLine,'_',Chr,
                          "_Start=",Start,"_end=",End,"_PNMF_k=",signif(Kvec,digits = 3),"_Remove=",paste(InforExclude,collapse = ''),'_Rep=',Nrepeats,"_.Rdata"))
}
} else { print ('data not available in the window')}

#load('E003_chr17_Start=2500000_end=7500000_PNMF_k=0.0098636170721809_Remove=0.Rdata')

