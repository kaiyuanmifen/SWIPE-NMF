#this version generate all cell lines 
#this version include DNase data 
#this version prepare for Jason's data instead of JW
#get for all Chr 
#this version take cell type and and chr
#this version makes thes the generation modulizable 
#this is the code to is to transform the eaw matrices foramts into formats ready to be taken by PNMF
#a lot of the works needs to be done manually 
library(Matrix)
library(methods)
#First of all prepare the matrices for PNMF

#Prepare R
#t3his is a implemenetaion of NMF
#source('NMFfunctions8.R')
#VecArg<-commandArgs(trailingOnly = T)
CellLine<-as.character(commandArgs(trailingOnly = T)[1])
#Chr<-as.character(VecArg[2])#Only focus on a chr for now 

Address='/broad/compbio/Dianbo/'

AllChr<-paste0('chr',c(1:22,'X'))
FileList<-list.files(path = paste0(Address,'datasets/ProducedData/RawMatrices/'))
#ALlCellLine=unique(mapply(FileList,FUN = function(x){strsplit(x,'[_|.]')[[1]][2]}))

#for (CellLine in ALlCellLine){



for (Chr in AllChr){
#Chr<-'chr15'#Only focus on a chr for now 


print(CellLine)
print(Chr)


FileList<-list.files(path = paste0(Address,'datasets/ProducedData/RawMatrices/'))
File<-FileList[grepl(FileList,pattern = CellLine)]
X<-load(paste0(Address,'datasets/ProducedData/RawMatrices/',File))
MatricesForNMF<-get(X)
#lapply(MatricesForNMF,dim)
#load('PromoterDisMat_E116.Rdata')
#load('PromoterDisMat_E116.Rdata')

#three data types 1=Enhancer, 2=promtoer ,3=HiC


MatricesForNMF<-lapply(MatricesForNMF, function(x){
        VecWhich=which(grepl(pattern = 'chr',strsplit(rownames(x)[1],'_')[[1]]))
        VecSub=strsplit(rownames(x),'_')
        VecChr=unlist(VecSub)[seq(from = VecWhich,to = length(unlist(VecSub)),by = length(strsplit(rownames(x)[1],'_')[[1]]))]
        VecRow=which(VecChr==Chr)
        
        VecWhich=which(grepl(pattern = 'chr',strsplit(colnames(x)[1],'_')[[1]]))
        VecSub=strsplit(colnames(x),'_')
        VecChr=unlist(VecSub)[seq(from = VecWhich,to = length(unlist(VecSub)),by = length(strsplit(colnames(x)[1],'_')[[1]]))]
        VecCol=which(VecChr==Chr)
        return(x[VecRow,VecCol])})
names(MatricesForNMF)
lapply(MatricesForNMF,dim)
#MatricesForNMF$JWEPremapped


R12=MatricesForNMF$JEEPremapped
R13=MatricesForNMF$Enhancer_HiC
R15=MatricesForNMF$Enhancer_TAD
R23=MatricesForNMF$Promoter_HiC
R25=MatricesForNMF$Promoter_TAD
R41=MatricesForNMF$eQTLSNP_RoadmapEnhancer
R42=MatricesForNMF$eQTLSNP_Gene
R61=MatricesForNMF$DHSThurman_Enhancer
R62=MatricesForNMF$DHSThurman_Gene

#R31=t(MatricesForNMF$Enhancer_HiC)


Nsources<-6

DimensionOfEachSet<-c(nrow(R12),nrow(R23),ncol(R13),nrow(R41),ncol(R25),nrow(R61))
#construct R

R<-list(R12=R12,R13=R13,R15=R15,R23=R23,R25=R25,R41=R41,R42=R42,R61=R61,R62=R62)

R<-lapply(R,function(x){x<-as.matrix(x)
x[x>0]<-1;return(x)})#make the R matrices binary 
R<-lapply(R,function(x){Matrix(x,sparse = T)})
#lapply(R,function(x){range(summary(x)$x)})
lapply(R,dim)

# head(summary(R$R23))

print('R done')



#Theta

#Prepare Theta
Theta33<-MatricesForNMF$HiC
Theta55<-MatricesForNMF$TAD_Dixon
dim(Theta55)

#install.packages('blockmatrix')

#build Matrix 



head(summary(Theta33))

Theta33<-as.matrix(Theta33)
Theta33[Theta33>0]<--0.2#make this binary 


Theta33<-Matrix(Theta33,sparse = T)

Theta55<-as.matrix(Theta55)
Theta55[Theta55>0]<--0.2#make this binary 

Theta55<-Matrix(Theta55,sparse=T)

Theta<-list()
for (i in 1:Nsources){
        
        Theta[[i]]<-list()
}




Theta[[1]][[1]]<-list()
Theta[[2]][[1]]<-list()
Theta[[3]][[1]]<-Theta33
Theta[[4]][[1]]<-list()
Theta[[5]][[1]]<-Theta55
Theta[[6]][[1]]<-list()
length(Theta)

print('Theta done')

#initlaize G

Rvec<-list(Rvec1=cbind(R$R12,R$R13,t(R$R41)),Rvec2=cbind(R$R23,t(R$R42)),
           Rvec3=cbind(t(R$R13),t(R$R23)),Rvec4=cbind(R$R41,R$R42),
           Rvec5=cbind(t(R$R15),t(R$R25)),
           Rvec6=cbind(R$R61,R$R62))
lapply(Rvec, dim)
lapply(R,dim)

#G<-InitalizeG(k = 0.01,DimensionOfEachSet = DimensionOfEachSet ,RVec = Rvec)

#Save the NMF input inforation 
NMFInputMats<-list(DimensionOfEachSet=DimensionOfEachSet,
                   Rvec=Rvec,R=R,Theta=Theta)




save(NMFInputMats,file = paste0(Address,'datasets/ProducedData/NMFInputMats/',CellLine,'_',Chr,'.Rdata'))

}
#}
