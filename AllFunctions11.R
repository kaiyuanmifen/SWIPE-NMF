#this is the file containing all functions 


#load libraries 

library("flux")


#Define functions to operates on Edgelist 

##reduce size of matrix 
Shrink<-function(EdgeList,SHNumber){
        EdgeList<-EdgeList[-sample(1:nrow(EdgeList),SHNumber),]
        EdgeList[EdgeList$RowID<EdgeList$ColumnID,c(1,2)]<-EdgeList[EdgeList$RowID<EdgeList$ColumnID,c(2,1)]
        return(EdgeList)
}

##rewire
RewireDis<-function(EdgeList,RewireNumber){
        
        AllNode<-unique(c(Target$RowID,Target$ColumnID))
        
        EdgeNumber<-sample(1:nrow(EdgeList),RewireNumber)
        VecRewire<-EdgeList[EdgeNumber,c(1,2)]
        
        #rewrie to opposite site 
        NewColID<-VecRewire$RowID-(VecRewire$ColumnID-VecRewire$RowID)
        #if the opposite is not available , use original -1 or +1
        NewColID[!NewColID%in%AllNode]<-VecRewire$ColumnID[!NewColID%in%AllNode]-1
        NewColID[!NewColID%in%AllNode]<-VecRewire$ColumnID[!NewColID%in%AllNode]+1
        sum(!NewColID%in%AllNode)
        
        EdgeList$ColumnID[EdgeNumber]<-NewColID
        EdgeList[EdgeList$RowID<EdgeList$ColumnID,c(1,2)]<-EdgeList[EdgeList$RowID<EdgeList$ColumnID,c(2,1)]
        return(EdgeList)
}





AddEdges<-function(EdgeList,AddNumber){
        
        RowID<-sample(NodeNumber,AddNumber,replace = T)
        AllDis<-EdgeList$ColumnID-EdgeList$RowID
        SampledID<-sample(1:nrow(EdgeList),AddNumber,replace = T)
        SampledDis<-AllDis[SampledID]
        
        ColID<-RowID+SampledDis
        ColID[!ColID%in%NodeNumber]<-(RowID-SampledDis)[!ColID%in%NodeNumber]
        sum(!ColID%in%NodeNumber)
        
        VecToAdd<-data.frame(RowID,ColID,EdgeList$Score[SampledID])
        
        names(VecToAdd)<-names(EdgeList)[c(1,2,8)]
        
        EdgeList<-rbind(EdgeList[,c(1,2,8)],VecToAdd)
        EdgeList[EdgeList$RowID<EdgeList$ColumnID,c(1,2)]<-EdgeList[EdgeList$RowID<EdgeList$ColumnID,c(2,1)]
        return(EdgeList)
}






#Add Gaussian noise 

AddGaussianNoise<-function(x){
        x$Score<-rnorm(n = nrow(x),mean = x$Score,sd =x$Score/10 )        
        x$Score[x$Score<0]<-0
        x$Score[x$Score>1]<-1
        #range(x$Score)
        return(x)
}




#rewire and add gaussian noise to sparse matrices 

SparseMatRewire<-function(VecMat,RewireNumber){
        library(Matrix)
        EdgeList<-summary(VecMat)
        head(EdgeList)
        AllNode<-unique(c(EdgeList$i,EdgeList$j))
        
        EdgeNumber<-sample(x = 1:nrow(EdgeList),size = RewireNumber,replace = F)
        VecRewire<-EdgeList[EdgeNumber,c(1,2)]
        
        #rewrie to opposite site 
        NewColID<-VecRewire$i-(VecRewire$j-VecRewire$i)
        #if the opposite is not available , use original -1 or +1
        NewColID[!NewColID%in%AllNode]<-VecRewire$j[!NewColID%in%AllNode]-1
        NewColID[!NewColID%in%AllNode]<-VecRewire$j[!NewColID%in%AllNode]+1
        
        
        
        
        if(sum(!NewColID%in%AllNode)>0){
                EdgeNumber<-EdgeNumber[NewColID%in%AllNode]
                NewColID<-NewColID[NewColID%in%AllNode]
        }
        EdgeList$j[EdgeNumber]<-NewColID
        EdgeList[EdgeList$i<EdgeList$j,c(1,2)]<-EdgeList[EdgeList$i<EdgeList$j,c(2,1)]
        
        sum(!EdgeList$j%in%AllNode)
        
        NewMat<-sparseMatrix(i = EdgeList$i,j = EdgeList$j,x = EdgeList$x,
                             dims = dim(VecMat),dimnames = dimnames(VecMat),
                             symmetric = T)
        
        return(NewMat)
}


AddGaussianNoiseSpraseMat<-function(VecMat){
        library(Matrix)
        EdgeList<-summary(VecMat)
        EdgeList$x<-rnorm(n = nrow(EdgeList),mean = EdgeList$x,sd =EdgeList$x/10 )  
        
        EdgeList$x[EdgeList$x<0]<-0
        EdgeList$x[EdgeList$x>1]<-1
        NewMat<-sparseMatrix(i = EdgeList$i,j = EdgeList$j,x = EdgeList$x,
                             dims = dim(VecMat),dimnames = dimnames(VecMat),
                             symmetric = T)
        head(summary(NewMat))
        return(NewMat)
}






##only looka t edges with data suport

FilterMat<-function(TargetMat,NewDataMat){
        #function to filter one mat based on another 
        OnlySup<-summary(TargetMat)
        OnlySup<-OnlySup[OnlySup$x>0,]
        nrow(OnlySup)
        
        NewVec<-OnlySup
        NewVec$x<-NewDataMat[as.matrix(OnlySup[,c(1,2)])]
        NewVec<-NewVec[NewVec$x>0,]
        nrow(NewVec)
        FilteredMat<-sparseMatrix(i = NewVec$i,j = NewVec$j,x = NewVec$x,
                                  dims = dim(TargetMat),
                                  dimnames = dimnames(TargetMat),
                                  symmetric = T)
        
        return(FilteredMat)
}






SortAcTGD<-function(VecMat){
        #sort the weight according to genomic distance (decreasing )
        MidPointPosition<-mapply(dimnames(VecMat)[[1]],FUN = function(x){mean(as.numeric(strsplit(x =x ,split = '_')[[1]][c(2,3)]))})
        
        DisMat<-mapply(MidPointPosition,FUN = function(x){abs(x-MidPointPosition)})
        class(DisMat)
        diag(DisMat)
        VecSum<-summary(VecMat)
        VecSum$Dis<-DisMat[as.matrix(VecSum[,c(1,2)])]
        
        head(VecSum)
        VecSum<-VecSum[order(VecSum$Dis,decreasing = F),]
        VecSum$x<-sort(VecSum$x,decreasing = T)
        VecSum[VecSum$j>VecSum$i,c(1,2)]<-VecSum[VecSum$j>VecSum$i,c(2,1)]
        
        NewMat<-sparseMatrix(i = VecSum$i,j = VecSum$j,x = VecSum$x,dims = dim(VecMat),dimnames = dimnames(VecMat),symmetric = T)
        
        return(NewMat)}










GetNearestMat<-function(VecMat){
        #Use VecMat a template to indentify which anchor to use 
        MidPointPosition<-mapply(dimnames(VecMat)[[1]],FUN = function(x){mean(as.numeric(strsplit(x =x ,split = '_')[[1]][c(2,3)]))})
        
        DisMat<-mapply(MidPointPosition,FUN = function(x){abs(x-MidPointPosition)})
        class(DisMat)
        diag(DisMat)
        
        VecSum<-summary(VecMat)
        VecSum<-VecSum[VecSum$x!=0,]
        AllInvolvedNodes<-sort(unique(c(VecSum$i,VecSum$j)))
        
        
        NewMat<-matrix(data = 0,nrow = nrow(VecMat),ncol = ncol(VecMat),dimnames = dimnames(VecMat))
        VecUpdate<-apply(DisMat[AllInvolvedNodes,AllInvolvedNodes],MARGIN = 1,FUN = function(x){which.min(x)})
        
        NewMat[cbind(AllInvolvedNodes,VecUpdate)]<-1
        NewMat[cbind(VecUpdate,AllInvolvedNodes)]<-1
        NewMat<-as(NewMat,'sparseMatrix')
        return(NewMat)
}




##function to get distances among nodes 
GetDisMat<-function(VecMat){
        #this method works for each chromosome 
       
    
        
        StartPosition<-as.numeric(unlist(lapply(Vecdata,FUN = function(x){x[2]})))
        EndPosition<-as.numeric(unlist(lapply(Vecdata,FUN = function(x){x[3]})))
        MidPoistion<-(EndPosition+StartPosition)/2
        #Chrmosome<-unlist(lapply(Vecdata,FUN = function(x){x[1]}))
        
        
        DisMat<-mapply(MidPoistion,FUN = function(x){abs(x-MidPoistion)})
        
        class(DisMat)
        diag(DisMat)
        return(DisMat)
}




GetAllDis<-function(VecMat){
        #this method works for each chromosome 
        
        Pos1=which(grepl(strsplit(rownames(VecMat)[1],'_')[[1]],pattern = 'chr'))
        Pos2=which(grepl(strsplit(colnames(VecMat)[1],'_')[[1]],pattern = 'chr'))
        
        Chr1=mapply(rownames(VecMat),FUN = function(x){strsplit(x,'_')[[1]][[Pos1]]})
        Start1=mapply(rownames(VecMat),FUN = function(x){as.numeric(strsplit(x,'_')[[1]][[Pos1+1]])})
        End1=mapply(rownames(VecMat),FUN = function(x){as.numeric(strsplit(x,'_')[[1]][[Pos1+2]])})
        Middle1=(Start1+End1)/2
        
        Chr2=mapply(colnames(VecMat),FUN = function(x){strsplit(x,'_')[[1]][[Pos2]]})
        Start2=mapply(colnames(VecMat),FUN = function(x){as.numeric(strsplit(x,'_')[[1]][[Pos2+1]])})
        End2=mapply(colnames(VecMat),FUN = function(x){as.numeric(strsplit(x,'_')[[1]][[Pos2+2]])})
        Middle2=(Start2+End2)/2
        
        VecChrMat=list()
        for (i in unique(Chr1)){
                VecChrMat[[length(VecChrMat)+1]]=mapply(Middle2[Chr2==i],FUN = function(x){abs(Middle1[Chr1==i]-x)})
                dimnames(VecChrMat[[length(VecChrMat)]])=dimnames(VecMat[Chr1==i,Chr2==i])
                        
        }
        DisMat=bdiag(VecChrMat)
        head(summary(DisMat))
        return(DisMat)
}






BipartieToAA<-function(VecMat){#function to convert a bipartie matrix to AA matrix 
        AllNodes<-unique(c(rownames(VecMat),colnames(VecMat)))
        VecSum<-summary(VecMat)
        tail(VecSum)
        VecSum$inames<-rownames(VecMat)[VecSum$i]
        VecSum$jnames<-colnames(VecMat)[VecSum$j]
        VecSum$Newi<-match(VecSum$inames, AllNodes)
        VecSum$Newj<-match(VecSum$jnames, AllNodes)
        VecSum[VecSum$Newi<VecSum$Newj,c(6,7)]<-VecSum[VecSum$Newi<VecSum$Newj,c(7,6)]
        OutPutMat<-sparseMatrix(i =VecSum$Newi,j = VecSum$Newj,x =VecSum$x, dims = c(length(AllNodes),length(AllNodes)),
                                dimnames = list(AllNodes,AllNodes),symmetric = T)
        head(summary(OutPutMat))
        return(OutPutMat)
        
}





RemoveLinkByDis<-function(VecMat){
        #this is the function to remove links <5k and >50K among EP 
        VecDisMat<-GetDisMat(EPMat)
        VecSum<-summary(VecMat)
        VecSum$dis<-VecDisMat[as.matrix(VecSum[,c(1,2)])]
        #VecSum<-VecSum[(VecSum$dis<=50000)&(VecSum$dis>5000),]
        VecSum<-VecSum[(VecSum$dis>5000),]
        NewMat<-sparseMatrix(i = VecSum$i,j = VecSum$j,x =VecSum$x,dims = dim(VecMat),
                             dimnames = dimnames(VecMat) ,symmetric = T)
        return(NewMat)
        
}







MakeTracks<-function(VecTrack){
        library(GenomicRanges)
        VecTrack<-strsplit(rownames(VecTrack),split = '_')
        VecTrack<-as.data.frame(Reduce(rbind,VecTrack))
        names(VecTrack)<-c('chr','start','end')
        head(VecTrack)
        VecTrack$start<-as.numeric(as.character(VecTrack$start))
        VecTrack$end<-as.numeric(as.character(VecTrack$end))
        return(makeGRangesFromDataFrame(df = VecTrack))
}



MapMatrieces<-function(MatFrom,MatTo){
        #This function needs to be debugged 
        #from is the data, to is the matrix I want to check whhich edge is supported 
        
        ToTrack<-MakeTracks(MatTo)
        head(summary(MatTo))
        head(ToTrack)
        length(ToTrack)
        FromTrack<-MakeTracks(MatFrom)
        length(FromTrack)
        
        ToTrack@ranges@start<-as.integer(floor(ToTrack@ranges@start+(ToTrack@ranges@width/2)))#consider overlap when center of EP overlap
        ToTrack@ranges@width<-rep(as.integer(1),length(ToTrack))
        
        length(ToTrack)
        VecOverLap<-findOverlaps(query =ToTrack ,subject = FromTrack)
        length(VecOverLap)
        tail(VecOverLap)
        length(unique(VecOverLap@from))
        length(unique(VecOverLap@to))
        (1:nrow(MatTo))[!((1:nrow(MatTo))%in%unique(VecOverLap@to))]
        
        
        SumVec<-summary(MatTo)
        tail(SumVec)
        nrow(SumVec)
        SumVec$x2<-0
        tail(SumVec)
        for (p in 1:nrow(SumVec)){
                
                Side1<-VecOverLap@to[which(VecOverLap@from==SumVec$i[p])]
                Side2<-VecOverLap@to[which(VecOverLap@from==SumVec$j[p])]
                if((length(Side1)>0)&((length(Side2)>0))){
                
                        SumVec$x2[p]<-  max(MatFrom[Side1,Side2])
                        }
              
    
        
        }
        
        sum(SumVec$x2>0)/length(SumVec$x2)
        
        head(SumVec)
        SumVec<-SumVec[SumVec$x2!=0,]
        MapedTo<-sparseMatrix(i = SumVec$i,j = SumVec$j,x = SumVec$x2,
                              dims = dim(MatTo),dimnames = dimnames(MatTo),use.last.ij = T,
                              symmetric = T)
        head(summary(MapedTo))
        head(summary(MatTo))
        return(MapedTo)
        
}


MapHighReToLowRe<-function(MatFrom,MatTo){
       #This is function to map high resolution networks to low one
        
        
        ToTrack<-MakeTracks(MatTo)
        head(summary(MatTo))
        head(ToTrack)
        length(ToTrack)
        FromTrack<-MakeTracks(MatFrom)
        length(FromTrack)
        
        FromTrack@ranges@start<-as.integer(floor(FromTrack@ranges@start+(FromTrack@ranges@width/2)))#consider overlap when center of EP overlap
        FromTrack@ranges@width<-rep(as.integer(1),length(FromTrack))
        
        length(ToTrack)
        VecOverLap<-findOverlaps(query =FromTrack ,subject = ToTrack)
        length(VecOverLap)
        tail(VecOverLap)
        length(unique(VecOverLap@from))
        length(unique(VecOverLap@to))
        (1:nrow(MatTo))[!((1:nrow(MatTo))%in%unique(VecOverLap@to))]
        
        
        SumVec<-summary(MatTo)
        tail(SumVec)
        nrow(SumVec)
        SumVec$x2<-0
        tail(SumVec)
        for (p in 1:nrow(SumVec)){
                
                Side1<-VecOverLap@from[which(VecOverLap@to==SumVec$i[p])]
                Side2<-VecOverLap@from[which(VecOverLap@to==SumVec$j[p])]
                if((length(Side1)>0)&((length(Side2)>0))){
                        
                        SumVec$x2[p]<-  max(MatFrom[Side1,Side2])
                }
                
                
                
        }
        
        sum(SumVec$x2>0)/length(SumVec$x2)
        
        head(SumVec)
        
        SumVec<-SumVec[SumVec$x2!=0,]
        nrow(SumVec)
        MapedTo<-sparseMatrix(i = SumVec$i,j = SumVec$j,x = SumVec$x2,
                              dims = dim(MatTo),dimnames = dimnames(MatTo),use.last.ij = T,
                              symmetric = T)
        head(summary(MapedTo))
        head(summary(MatTo))
        return(MapedTo)
        
}





RandomiszeEdge<-function(InPutMat){
        #Randomisze edges with genomic distance constraint but wiihin existing anchor 
        
        DisMat<-GetDisMat(InPutMat)
        VecSum<-summary(InPutMat)
        
        CurrentDis<-DisMat[as.matrix(VecSum[,c(1,2)])]
        
        GetNewj<-function(x){
                Newj<-names(sort(abs(DisMat[VecSum$i[x],-c(VecSum$i[x],VecSum$j[x])]-CurrentDis[x]))[1:3])
                Newj<-sample(Newj,1)
                Newj<-which(rownames(InPutMat)==Newj)
                return(Newj)
        }
        
        VecSum$j<-mapply(1:nrow(VecSum),FUN = GetNewj)
        VecSum[VecSum$i<VecSum$j,c(1,2)]<-VecSum[VecSum$i<VecSum$j,c(2,1)]
        
        #DisMat[VecSum$i[x],which(rownames(InPutMat)==Newj)]
        #DisMat[VecSum$i[x],268]
        OutPutMat<-sparseMatrix(i = VecSum$i,j = VecSum$j,x = VecSum$x,dims = dim(InPutMat),
                                dimnames = dimnames(InPutMat),symmetric = T,use.last.ij = T)
        return(OutPutMat)
}





GetDisMat<-function(VecMat){
        MidPointPosition<-mapply(dimnames(VecMat)[[1]],FUN = function(x){mean(as.numeric(strsplit(x =x ,split = '_')[[1]][c(2,3)]))})
        DisMat<-mapply(MidPointPosition,FUN = function(x){abs(x-MidPointPosition)})
        class(DisMat)
        diag(DisMat)
        return(DisMat)
}



RandomiszeEdge<-function(InPutMat){
        #Randomisze edges with genomic distance constraint but wiihin existing anchor 
        
        DisMat<-GetDisMat(InPutMat)
        VecSum<-summary(InPutMat)
        
        CurrentDis<-DisMat[as.matrix(VecSum[,c(1,2)])]
        
        GetNewj<-function(x){
                Newj<-names(sort(abs(DisMat[VecSum$i[x],-c(VecSum$i[x],VecSum$j[x])]-CurrentDis[x]))[1:3])
                Newj<-sample(Newj,1)
                Newj<-which(rownames(InPutMat)==Newj)
                return(Newj)
        }
        
        VecSum$j<-mapply(1:nrow(VecSum),FUN = GetNewj)
        VecSum[VecSum$i<VecSum$j,c(1,2)]<-VecSum[VecSum$i<VecSum$j,c(2,1)]
        
        #DisMat[VecSum$i[x],which(rownames(InPutMat)==Newj)]
        #DisMat[VecSum$i[x],268]
        OutPutMat<-sparseMatrix(i = VecSum$i,j = VecSum$j,x = VecSum$x,dims = dim(InPutMat),
                                dimnames = dimnames(InPutMat),symmetric = T,use.last.ij = T)
        return(OutPutMat)
}






UpdateEPWithHiC<-function(MatHiC,MatEP){
        #this is acutually using HiC data to filter EP links , those EP links in the same block are maintained and kept the score
        #the links across blocks but not detected  by HiC are removed (zero) and those supported by HiC are maintained and score using the HiC score
        # a new version to map HiC to EP  
        library(GenomicRanges)
        
        HiCTrack<-MakeTracks(MatHiC)
        
        EPTrack<-MakeTracks(MatEP)
        
        #Firstly, math the nodes 
        
        VecOverLap<-findOverlaps(query = EPTrack,subject = HiCTrack)
        length(VecOverLap)
       
        
    
        VecSumEP<-summary(MatEP)
        VecSumEP$HiCi<-VecOverLap@to[match(VecSumEP$i,VecOverLap@from)]        
        VecSumEP$HiCj<-VecOverLap@to[match(VecSumEP$j,VecOverLap@from)] 
        VecSumEP$HiCx<-MatHiC[as.matrix(VecSumEP[,c("HiCi","HiCj")])]
        #VecSumEP$x[VecSumEP$HiCi==VecSumEP$HiCj]<-VecSumEP$x[VecSumEP$HiCi==VecSumEP$HiCj]#Keep score within the same segment the same
        VecSumEP$x[VecSumEP$HiCi==VecSumEP$HiCj]<-1#Keep score within the same segment 1
        VecSumEP$x[(VecSumEP$HiCi!=VecSumEP$HiCj)&(VecSumEP$HiCx==0)]<-0
        VecSumEP$x[(VecSumEP$HiCi!=VecSumEP$HiCj)&(VecSumEP$HiCx!=0)]<-VecSumEP$HiCx[(VecSumEP$HiCi!=VecSumEP$HiCj)&(VecSumEP$HiCx!=0)]
        
        sum((VecSumEP$x!=0)&(VecSumEP$HiCx==0))
        sum(VecSumEP$HiCi==VecSumEP$HiCj)
        
    
        VecSumEP<-VecSumEP[VecSumEP$x!=0,]
        MapedTo<-sparseMatrix(i =VecSumEP$i,j =VecSumEP$j,x =VecSumEP$x,dims = dim(MatEP),
                              symmetric =T,dimnames = dimnames(MatEP)   )
        head(summary(MapedTo))
        return(MapedTo)
        
}




MapNetworkForCheck<-function(TestNet,GroundTrueth){
    
        library(GenomicRanges)
        
        TruethTrack<-MakeTracks(GroundTrueth)
        
        TestTrack<-MakeTracks(TestNet)
        
        
        #Firstly, math the nodes 
        
        VecOverLap<-findOverlaps(query = TestTrack,subject = TruethTrack)
        length(VecOverLap)
        
        
        
        VecSumEP<-summary(TestNet)
        VecSumEP$Truethi<-VecOverLap@to[match(VecSumEP$i,VecOverLap@from)]
        VecSumEP$Truethj<-VecOverLap@to[match(VecSumEP$j,VecOverLap@from)]
        head(VecSumEP)
        VecSumEP$Truex<-0
        VecSumEP$Truex[(is.na(VecSumEP$Truethi)|is.na(VecSumEP$Truethj))]<-0
        WhichVec<-which((!is.na(VecSumEP$Truethi))&(!is.na(VecSumEP$Truethj)))#Which position can be mapped 
        VecSumEP$Truex[WhichVec]<-GroundTrueth[as.matrix(VecSumEP[WhichVec,c("Truethi","Truethj")])]
        VecSumEP$Truex[VecSumEP$HiCi==VecSumEP$HiCj]<-VecSumEP$x[VecSumEP$HiCi==VecSumEP$HiCj]#Keep score within the same segment 1
        
        sum(VecSumEP$Truex>0)/length(sum(VecSumEP$Truex>0) )      
       tail(VecSumEP)
        VecSumEP<-VecSumEP[VecSumEP$Truex>0,]
        MapedTo<-sparseMatrix(i =VecSumEP$i,j =VecSumEP$j,x =VecSumEP$Truex,dims = dim(TestNet),
                              symmetric =T,dimnames = dimnames(TestNet)   )
        head(summary(MapedTo))
        return(MapedTo)
        
}




MapNetworkForEP<-function(EPLinks,TargetNet){
        #This funciton maps different links to EP
        library(GenomicRanges)
        
        PositionVec=which(grepl(strsplit(rownames(EPLinks)[1],'_')[[1]],pattern = 'chr'))
        EnhancerAnchors<-makeGRangesFromDataFrame(data.frame(chr=mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][PositionVec]}),
                                                             start=as.numeric(mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+1]})),
                                                             end=as.numeric(mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+2]}))))
        
        #EnhancerAnchors@ranges@start<-as.integer(round(EnhancerAnchors@ranges@start+EnhancerAnchors@ranges@width/2))
        #EnhancerAnchors@ranges@width<-as.integer(rep(1,length(EnhancerAnchors)))#map the center to my definition
        
        
        dim(EPLinks)
        
        colnames(EPLinks)
        PositionVec=which(grepl(strsplit(colnames(EPLinks)[1],'_')[[1]],pattern = 'chr'))
        
        PromoterAnchors<-makeGRangesFromDataFrame(data.frame(chr=mapply(colnames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][PositionVec]}),
                                                             start=as.numeric(mapply(colnames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+1]})),
                                                             end=as.numeric(mapply(colnames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+2]}))))
        
        #PromoterAnchors@ranges@start<-as.integer(round(PromoterAnchors@ranges@start+PromoterAnchors@ranges@width/2))
        #PromoterAnchors@ranges@width<-as.integer(rep(1,length(PromoterAnchors)))#map the center to my definition
        
        
        
        
        
        
        #ChiaPet to Look at EP
        library(GenomicRanges)
        
        PositionVec=which(grepl(strsplit(rownames(TargetNet)[1],'_')[[1]],pattern = 'chr'))
        
        TargetTrackRow<-makeGRangesFromDataFrame(data.frame(chr=mapply(rownames(TargetNet),FUN = function(x){strsplit(x,'_')[[1]][PositionVec]}),
                                                          start=as.numeric(mapply(rownames(TargetNet),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+1]})),
                                                          end=as.numeric(mapply(rownames(TargetNet),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+1]}))))
        
        PositionVec=which(grepl(strsplit(colnames(TargetNet)[1],'_')[[1]],pattern = 'chr'))
        TargetTrackCol<-makeGRangesFromDataFrame(data.frame(chr=mapply(colnames(TargetNet),FUN = function(x){strsplit(x,'_')[[1]][PositionVec]}),
                                                            start=as.numeric(mapply(colnames(TargetNet),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+1]})),
                                                            end=as.numeric(mapply(colnames(TargetNet),FUN = function(x){strsplit(x,'_')[[1]][PositionVec+1]}))))
        
        
        
        VecOverLapE<-findOverlaps(query = TargetTrackRow,subject =EnhancerAnchors)
        head(VecOverLapE)
        length(VecOverLapE)==length(unique(VecOverLapE@to))#check whether one to one match 
        VecOverLapE@from
        ECMat<-sparseMatrix(i = VecOverLapE@to,j = VecOverLapE@from,x = 1,
                            dims = c(length(EnhancerAnchors),length(TargetTrackRow)),
                            dimnames = list(names(EnhancerAnchors),names(TargetTrackRow)))
        
        
        VecOverLapP<-findOverlaps(query = TargetTrackCol,subject =PromoterAnchors)
        head(VecOverLapP)
        length(VecOverLapP)==length(unique(VecOverLapP@from))#check whether one to one match 
        VecOverLapP
        PCMat<-sparseMatrix(i = VecOverLapP@to,j = VecOverLapP@from,x = 1,
                            dims = c(length(PromoterAnchors),length(TargetTrackCol)),
                            dimnames = list(names(PromoterAnchors),names(TargetTrackCol)))
        
        #diag(TargetNet)<-0
        
        
        InterSupported<-tcrossprod(ECMat%*%TargetNet,PCMat)
        dim(InterSupported)
        return(InterSupported)
}








GetNearestEP<-function(VecMat){
        EMidPoint<-mapply(dimnames(VecMat)[[1]],FUN = function(x){mean(as.numeric(strsplit(x =x ,split = '_')[[1]][c(3,4)]))})
        PMidPoint<-mapply(dimnames(VecMat)[[2]],FUN = function(x){mean(as.numeric(strsplit(x =x ,split = '_')[[1]][c(4,5)]))})
        VecP<-mapply(EMidPoint,FUN = function(x){which.min(abs(x-PMidPoint))})
        VecSum<-data.frame(i=1:nrow(VecMat),j=VecP,x=1)
        VecSum$Dis<-mapply(EMidPoint,FUN = function(x){min(abs(x-PMidPoint))})
        VecSum<-VecSum[(VecSum$Dis>5000),]
        
        OutputMat<-sparseMatrix(i = VecSum$i,j = VecSum$j,x = VecSum$x,dims = dim(VecMat),dimnames = dimnames(VecMat),symmetric = F)
        return(OutputMat)
}



#Caculate TFrates and False Posittive rates 

GetTPR<-function(Prediction,GDstandard){
        #VecSum<-summary(Prediction)
        #VecSum<-VecSum[VecSum$x>0,]
        #VecReal<-GDstandard[as.matrix(VecSum[,c(1,2)])]
        VecReal=GDstandard[Prediction>0]
        return(sum(VecReal>0)/sum(Prediction>0))
}

GetFPR<-function(Prediction,GDstandard){#Specifically for Nearest 
        #VecSum<-summary(Prediction)
        #VecSum<-VecSum[VecSum$x>0,]
        #VecReal<-GDstandard[as.matrix(VecSum[,c(1,2)])]
        #VecAll<-GDstandard[unique(VecSum[,1]),unique(VecSum[,2])]
        VecReal=GDstandard[Prediction>0]
        return(sum(VecReal==0)/(sum(GDstandard==0)))
}



#Setting cutoff for enhancer promtoer networks 



