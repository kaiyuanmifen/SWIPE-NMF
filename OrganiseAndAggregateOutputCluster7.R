#this version add a part to get optimum ranl for each block 
#this version include 'Stick past togeter'
#this is to look at contribution of each set 
#this is to check the pnmf result 
library(methods)
library(Matrix)
X=load('AllCellTypes.Rdata')
AllCellType=get(X)
Address<-"/homes/dliu/datasets/ProducedData/NMFoutput"

FileList<-list.files(Address)
length(FileList)
#For each each line and each window
#mapply(FileList,FUN = function(x){strsplit(x,split = '[_|=]')})
CellLine=AllCellType[as.integer(as.character(commandArgs(trailingOnly = T)))]
task.id <-as.integer(Sys.getenv("SGE_TASK_ID"))
#for (VecChr in 1:23){


AllChromosomes<-paste0('chr',c(1:22,'X'))



#Aggreate repititions for each sliding window

Chr=paste0('_',AllChromosomes[task.id],'_')


Files<-FileList[grepl(FileList,pattern = paste0(CellLine,Chr))]

FileInfor<-data.frame(CellType=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][1]),
           Chromosome=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][2]),
           Start=as.integer(mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][4])),
           End=as.integer(mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][6])),
           k=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][9]),
           Remove=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][11]),
           Rep=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][13]))



#For each specific K
 IntermediateAddress<-paste0(paste0(strsplit(Address,'/')[[1]][-length(strsplit(Address,'/')[[1]])],collapse = '/'),'/OrganisedNMFOutput/')

for (Kvec in unique(FileInfor$k)){

        FileInforPerK<-FileInfor[FileInfor$k==Kvec,]
        for (StartVec in unique(FileInforPerK$Start)){

                FileInforSub<-FileInforPerK[FileInforPerK$Start==StartVec,]
                print(paste0('working on',CellLine,Chr,Kvec,"Start_",StartVec,'_End=',unique(FileInforSub$End)))
                ToLoad<-paste0(Address,rownames(FileInforSub))

                EPmat<-list()
                G1<-list()
                G2<-list()
                Smats<-list()
                ErrorMat<-NULL
                for (i in ToLoad){
                        X<-load(i)
                        Loaded=get(X)
                        EPmat[[length(EPmat)+1]]<-Loaded$Rc$R12
                        
                        rownames(Loaded$G$G1)<-rownames(Loaded$Rc$R12)
                        rownames(Loaded$G$G2)<-colnames(Loaded$Rc$R12)
                        G1[[length(G1)+1]]<-Loaded$G$G1
                        G2[[length(G2)+1]]<-Loaded$G$G2
                        ErrorMat<-rbind(ErrorMat,Loaded$Error)
                }
                EPmat<-Reduce(EPmat,f = '+')/length(EPmat)
                
                G1<-Reduce(G1,f = '+')/length(G1)
                G2<-Reduce(G2,f = '+')/length(G2)
                
                Output<-list(Mat=EPmat,Error=ErrorMat,G1=G1,G2=G2)

                save(Output,file =paste0(IntermediateAddress,CellLine,Chr,'k=',Kvec,'_Start_',StartVec,'_End_',
                                        unique(FileInforSub$End),'_total_rep=',max(as.integer(as.character(FileInforSub$Rep))),'.Rdata') )

        }





}





#select the best rank for each window 



FileList<-list.files(IntermediateAddress)
Files<-FileList[grepl(FileList,pattern = paste0(CellLine,Chr))]

#mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')})
FileInfor<-data.frame(CellType=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][1]),
                      Chromosome=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][2]),
                      Start=as.integer(mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][6])),
                      End=as.integer(mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][8])),
                      k=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][4]))
VecBestK<-NULL
for (StartVec in sort(unique(FileInfor$Start))){
        FileInforPerStart=FileInfor[FileInfor$Start==StartVec,]
        VecErrors<-NULL
        for (Kvec in sort(unique(FileInforPerStart$k))){
                X<-load(paste0(IntermediateAddress,rownames(FileInforPerStart)[FileInforPerStart$k==Kvec]))  
                Vec<-get(X)
                VecErrors<-rbind( VecErrors,apply(Vec$Error,MARGIN = 2,mean))
        }
        rownames(VecErrors)<-sort(unique(FileInforPerStart$k))
        VecPlot<-apply(VecErrors,MARGIN = 1,FUN = sum)
        Y=VecPlot
        X=as.numeric(names(VecPlot))
        
        VecDiff<-diff(Y)/diff(X)
        
        Vecwhich<-which.max(VecDiff)+1
        
        KKink=X[Vecwhich-1]
        VecBestK<-c(VecBestK,KKink)
        #plot(Y~X,type='l')
     
        
}
names(VecBestK)<-sort(unique(FileInfor$Start))

BestKAddress<-paste0(paste0(strsplit(Address,'/')[[1]][-length(strsplit(Address,'/')[[1]])],collapse = '/'),'/BestRanks/')
save(VecBestK,file = paste0(BestKAddress,CellLine,Chr,'.Rdata'))



# #stick sliding windows  on each chromosome togethr
FileList<-list.files(IntermediateAddress)
Files<-FileList[grepl(FileList,pattern = paste0(CellLine,Chr))]

#mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')})
FileInfor<-data.frame(CellType=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][1]),
                      Chromosome=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][2]),
                      Start=as.integer(mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][6])),
                      End=as.integer(mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][8])),
                      k=mapply(Files,FUN = function(x){strsplit(x,split = '[=|_]')}[[1]][4]))

##Contruct the matrix for each chrmosome
load("Chromosomes.Rdata")
ChrLength=ChromosomeInfor$Base.pairs[ChromosomeInfor$Chromosome==gsub(Chr,pattern = '_',replacement = '')]
Starts=50000*(c(1:(ChrLength%/%5e6+1))-1)
Ends=50000*c(1:(ChrLength%/%5e6+1))

InputAddress<-paste0(paste0(strsplit(Address,'/')[[1]][-length(strsplit(Address,'/')[[1]])],collapse = '/'),'/NMFInputMats/')
InPutFiles<-list.files(InputAddress)
InPutFiles<-InPutFiles[grepl(InPutFiles,pattern = CellLine)]
InPutFiles<-InPutFiles[mapply(InPutFiles,FUN = function(x){strsplit(x,split = '[.|_]')[[1]][2]})==gsub(Chr,pattern = '_',replacement = '')]
X=load(paste0(InputAddress,InPutFiles))
InputMats<-get(X)
JasonEP<-InputMats$R$R12

OutPutAddress<-paste0(paste0(strsplit(Address,'/')[[1]][-length(strsplit(Address,'/')[[1]])],collapse = '/'),'/FinalOrganisedMats/')

EPOutPut<-sparseMatrix(i = 1,j = 1,x = 0,dims = dim(JasonEP),dimnames = dimnames(JasonEP))
G1Output<-sparseMatrix(i = 1,j = 1,x = 0,dims = c(nrow(JasonEP),nrow(JasonEP)),
                       dimnames = list(rownames(JasonEP),rownames(JasonEP)))
G2Output<-sparseMatrix(i = 1,j = 1,x = 0,dims = c(ncol(JasonEP),ncol(JasonEP)),
                       dimnames = list(colnames(JasonEP),colnames(JasonEP)))
OutPutMat<-list()

for (StartVec in sort(unique(FileInfor$Start))){
        print(paste0('Working on',StartVec))
        
        VecFileBestK<-rownames(FileInfor)[FileInfor$Start==StartVec&FileInfor$k==VecBestK[as.numeric(names(VecBestK))==StartVec]]
        X<-load(paste0(IntermediateAddress,VecFileBestK))
        WindowMat<-get(X)
        VecEP<-WindowMat$Mat
        VecG1<-WindowMat$G1
        VecG2<-WindowMat$G2
        
        Target<-VecEP
        EPOutPut[rownames(Target),colnames(Target)]<-(EPOutPut[rownames(Target),colnames(Target)]+Target)/
                max(((EPOutPut[rownames(Target),colnames(Target)]>0)+(Target>0)),1)
        
        Target<-VecG1
        Target<-tcrossprod(Target,Target)
        G1Output[rownames(Target),colnames(Target)]<-(G1Output[rownames(Target),colnames(Target)]+Target)/
                max(((G1Output[rownames(Target),colnames(Target)]>0)+(Target>0)),1)
        
        Target<-VecG2
        Target<-tcrossprod(Target,Target)
        G2Output[rownames(Target),colnames(Target)]<-(G2Output[rownames(Target),colnames(Target)]+Target)/
                max(((G2Output[rownames(Target),colnames(Target)]>0)+(Target>0)),1)
        
}

OutPutMat<-list(EP=EPOutPut,G1=G1Output,G2=G2Output)
save(OutPutMat,file = paste0(OutPutAddress,CellLine,Chr,'.Rdata'))




# X<-load(paste0(OutPutAddress,CellLine,Chr,0.15,'.Rdata'))
# Miega<-get(X)
# dim(Miega)
# head(summary(Miega))
# dim(JasonEP)




#Agree gate results 

OutPutAddress<-"homes/dliu/datasets/ProducedData/FinalOrganisedMats/"

FileList<-list.files(OutPutAddress)

#Integrate different chromosome together 
CellTypes<-mapply(FileList,FUN = function(x){strsplit(x,split = '[_|=]')[[1]][1]})
Chromosome<-mapply(FileList,FUN = function(x){strsplit(x,split = '[_|=]')[[1]][2]})
FileInfor<-data.frame(CellTypes,Chromosome)
FileInfor<-FileInfor[FileInfor$CellTypes==CellLine,]


AllChromosomes<-as.character(FileInfor$Chromosome)[order(as.numeric(mapply(as.character(FileInfor$Chromosome),FUN = function(x){substr(x,4,nchar(x))})))]

AllMatrix<-list()
AllG1<-list()
AllG2<-list()
for( Chr in AllChromosomes){
        print(paste('working on', Chr))
        File<-rownames(FileInfor)[as.character(FileInfor$Chromosome)==Chr] 
        X<-load(paste0(OutPutAddress,File))
        OutPut<-get(X) 
        AllMatrix[[length(AllMatrix)+1]]<-OutPut$EP
        AllG1[[length(AllG1)+1]]<-OutPut$G1
        AllG2[[length(AllG2)+1]]<-OutPut$G2
        
}

names(AllMatrix)<-AllChromosomes
names(AllG1)<-AllChromosomes
names(AllG2)<-AllChromosomes

EPLinks<-bdiag(AllMatrix)
G1links<-bdiag(AllG1)
G2links<-bdiag(AllG2)

dimnames(EPLinks)<-list(unlist(lapply(AllMatrix,rownames)),
                        unlist(lapply(AllMatrix,colnames)))

dimnames(G1links)<-list(unlist(lapply(AllG1,rownames)),
                        unlist(lapply(AllG1,colnames)))
dimnames(G2links)<-list(unlist(lapply(AllG2,rownames)),
                        unlist(lapply(AllG2,colnames)))
head(summary(G2links))

Output<-list(EP=EPLinks,G1=G1links,G2=G2links)

AggregateAddress<-paste0(paste0(strsplit(OutPutAddress,'/')[[1]][-length(strsplit(OutPutAddress,'/')[[1]])],collapse = '/'),'/PredictedLinks/')

save(Output,file = paste0(AggregateAddress,CellLine,'.Rdata'))

