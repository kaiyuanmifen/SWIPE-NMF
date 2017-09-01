#This version let each task to run the number of repeations 
#tHis version let each 
#this version doesn't introduce noise to tanks 
#This version implement the sliding window 
#this is to give parameter sets to run in task array 

CellTypes<-c('E127')
NumberOfRepeat<-1
K<-c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)
#Chromosome<-paste0('chr',c(1:22,'X'))
Chromosome<-paste0('chr',c(1:18,20:22,"X"))
SizeOfSlidingWindow<-5e6
NumberOfRepeats<-20
#DataExclude<-0

#Get the sliding windows
load('Chromosomes.Rdata')
Windows<-NULL
for (ChrVec in Chromosome){
        VecLength<-ChromosomeInfor$Base.pairs[ChromosomeInfor$Chromosome==ChrVec]      
        
        Nb<-ceiling(VecLength/SizeOfSlidingWindow)
        Nb<-2*Nb-1
        Vecragnes<-c(1:Nb)
        WinStart<-(Vecragnes-1)*(SizeOfSlidingWindow/2)
        WinEnd<-WinStart+SizeOfSlidingWindow
        Windows<-rbind(Windows,data.frame(Chr=ChrVec,WinStart=WinStart,WinEnd=WinEnd))
        
}

#Generate the parameter list 
ParameterSet<-NULL



#Windows<-Windows[rep(seq_len(nrow(Windows)), each=NumberOfRepeats), ]
#nrow(Windows)
EachCellLine<-NULL
for (i in K){
        EachCellLine<-rbind(EachCellLine,data.frame(K=i,Windows))
}
##Add some noise to each k 
# EachCellLine$K<-rnorm(n = length(EachCellLine$K),mean =EachCellLine$K,sd =0.01*EachCellLine$K  )
# head(EachCellLine)
# range(EachCellLine$K)

for (j in CellTypes ){#Dulplicate it for each cell line 
        ParameterSet<-rbind(ParameterSet,data.frame(CellType=j,EachCellLine))
}

ParameterSet$Repeats=NumberOfRepeats

nrow(ParameterSet)
head(ParameterSet)

unique(ParameterSet$WinStart)



#save
save(ParameterSet,file = 'ParametersetNormal.Rdata')
