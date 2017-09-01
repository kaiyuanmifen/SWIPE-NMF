#this version get rid dependence of position on names 
#This version tries to get rid of inversion of matrices 
#tries to increase speed a little 
#this version trie to improve the G initialization 
#this version modifies G initialization 
#this version debugs what is needed 
#this version adpat to cluster 
#function related to NMF
library('Matrix')
library(corpcor)

#initialise G
InitalizeG<-function(k,DimensionOfEachSet,Rvec,AcolN=100){


Ranks<-round(DimensionOfEachSet*k)
Ranks[Ranks<5]<-5
#Initalizaing G using random Acol 
G<-list()
# for (i in 1:length(Rvec)){
#         G[[i]]<-matrix(data =0,ncol =Ranks[i],nrow = DimensionOfEachSet[i] )
#         
# }

#G<-lapply(G,FUN = function(x){Matrix(x,sparse = T)})

Rvec<-lapply(Rvec,as.matrix)

RandomAcol<-function(OriMat,ColN){
        
        #ColN<-round(ncol(OriMat)*AcolN)
        
        
        VecMat<-lapply(as.list(1:ColN),FUN = function(x){OriMat[,sample(ncol(OriMat))]})
        VecMat<-Reduce(VecMat,f = '+')/length(VecMat)   
                
                return(VecMat)
        }
        

G<-lapply(as.list(1:length(Rvec)),FUN =function(x){ RandomAcol(OriMat = Rvec[[x]],ColN =  AcolN)[,1:Ranks[x]]})

G<-lapply(G,FUN = function(x){Matrix(x,sparse = T)})

names(G)<-paste0('G',substr(names(Rvec),start = 5,stop = 5))
      
 #head(summary(G$G3))
 #lapply(G,FUN = function(x){sum(x>0)})
 #lapply(G, dim)
# dim(R$R23)
# dim(G$G3)

return(G)
}





#Slideing window 

SlideWindow<-function(NMFInput,Start,End){
        NMFInput$R<-lapply(NMFInput$R,FUN = function(x){GetWindowRange(x,Start,End)})
        #lapply(NMFInput$R,dim)
        NMFInput$Rvec<-lapply(NMFInput$Rvec,FUN = function(x){GetWindowRange(x,Start,End)})
        #lapply(NMFInput$Rvec,dim)
        
        for (i in 1:length(NMFInput$Theta)){
                NMFInput$Theta[[i]]<-lapply(NMFInput$Theta[[i]],FUN = function(x){GetWindowRange(x,Start,End)})
        }
        #lapply(NMFInput$Theta,length)
        
        NMFInput$DimensionOfEachSet<-unlist(lapply(NMFInput$Rvec,nrow))
        
        # NMFInput$Window<-list()
        # NMFInput$Window[[1]]<-Start
        # NMFInput$Window[[2]]<-End
        # names(NMFInput$Window)<-c('Start','End')
        return(NMFInput)
}

GetWindowRange<-function(VecMat,Start,End){
        if (length(VecMat)>0){
                Position<-mapply(rownames(VecMat),FUN = function(x){VecS<-strsplit(x,'_')[[1]]
                return(list(VecS[(length(VecS)-1):length(VecS)]))})
                Position<-lapply(Position,as.numeric)
                Position<-Reduce(Position,f = rbind)
                VecRow<-which(Position[,1]>=Start&Position[,2]<=End)
                
                Position<-mapply(colnames(VecMat),FUN = function(x){VecS<-strsplit(x,'_')[[1]]
                return(list(VecS[(length(VecS)-1):length(VecS)]))})
                Position<-lapply(Position,as.numeric)
                Position<-Reduce(Position,f = rbind)
                VecCol<-which(Position[,1]>=Start&Position[,2]<=End)
                
                return(VecMat[VecRow,VecCol])}
        else {return(VecMat)}
        
}




#run the pnmf
# R<-NMFInputMats$R
#  Theta<-NMFInputMats$Theta
#  G<-G
#  TargetMatrixName<-'R12'
#  StopThreshold<-1e-5
#  DimensionOfEachSet<-NMFInputMats$DimensionOfEachSet

NMFDianbo<-function(R,Theta,G,TargetMatrixName,StopThreshold,MaxInteraction,DimensionOfEachSet,MaxTime=NULL){
        library(Matrix)
        library(MASS)
        
        PTM<-sum(proc.time())
        TargetMatrix<-R[[which(names(R)==TargetMatrixName)]]
        
        Ranks<-sapply(G,ncol)
        #R ,Theta and G needs to be transform into full matrices 
        
        ##R
        aR<-sparseMatrix(i = 1,j = 1,x = 0,dims = c(sum(DimensionOfEachSet),sum(DimensionOfEachSet)))
        VecR<-cbind(as.integer(sapply(names(R),FUN = function(x){strsplit(x,split = '')[[1]][2]})),
                    as.integer(sapply(names(R),FUN = function(x){strsplit(x,split = '')[[1]][3]})))
      
        
        VecPosition<-matrix(mapply(VecR,FUN = function(x){match(x,sort(unique(c(VecR))))}),ncol = 2)##reassign position becasue I may remove some data
        
        for (i in 1:nrow(VecPosition)){
                
                p<-VecPosition[i,1]
                j<-VecPosition[i,2]
                aR[(cumsum(DimensionOfEachSet)[p]+1-DimensionOfEachSet[p]):cumsum(DimensionOfEachSet)[p],
                   (cumsum(DimensionOfEachSet)[j]+1-DimensionOfEachSet[j]):cumsum(DimensionOfEachSet)[j]]<-R[[i]] 
                
        }
        head(summary(aR))
        sum(unlist(lapply(R,function(x){sum(x>0)})))
        #R$R12
        lapply(R,function(x){sum(x>0)})
        
        
        ##G
        aG<-Matrix(bdiag(G),sparse = T)
#         range(aG)
#         dim(aG)
#         head(summary(G$G1))
#         head(summary(aG))
#         range(summary(aG)$x)
#         sum(aG>0)
#         diag(aG)
#         lapply(G,function(x){sum(x>0)})
#         lapply(G,dim)
#         length(G)
#         
        #Theta
        
        NumOfTheta<-max(unlist(lapply(Theta, length)))
        aTheta<-list()
        for (i in 1:NumOfTheta){
                
                VecDig<-list()
                for (m in 1:length(Theta)){
                        if (length(Theta[[m]][[i]])>0){
                                VecDig[[m]]<-Theta[[m]][[i]]
                        } else {VecDig[[m]]<-sparseMatrix(i= 1,j = 1,x = 0,
                                                          dims =c(DimensionOfEachSet[m],DimensionOfEachSet[m]))}
                        
                }
                aTheta[[i]]<-bdiag(VecDig)
                
                names(aTheta)<-paste0('t',i)
                
        }
#         sum(DimensionOfEachSet)
#         dim(Theta33)
#         length(aTheta)
#         dim(aTheta[[1]])
#         sum(DimensionOfEachSet)
#         range(summary(aTheta[[1]])$x)
#         
        
        
        
        
        
        
        #start of repeat 
        Step<-0
        
        repeat{
        Step<-Step+1
          print(Step)      
        #calculate S
        
        #S<-ginv(as.matrix(t(aG)%*%aG))%*%t(aG)%*%aR%*%aG%*%ginv(as.matrix(t(aG)%*%aG))
        #S<-ginv(as.matrix(t(aG)%*%aG))%*%t(aG)%*%aR%*%aG%*%ginv(as.matrix(t(aG)%*%aG))
        #S<-solve((t(aG)%*%aG))%*%t(aG)%*%aR%*%aG%*%solve(t(aG)%*%aG)
          
          VecDot<-crossprod(aG,aG) 
          VecDot[is.na(VecDot)]<-0
          
          #VecMat<-ginv(as.matrix(VecDot))
          VecMat<-pseudoinverse(as.matrix(VecDot))
          #VecMat<-pseudoinverse(as.matrix(crossprod(aG,aG)))
         
        VecMat<-Matrix(VecMat,sparse = T)
        head(summary(VecMat))
        
        #S<-VecMat%*%t(aG)%*%aR%*%aG%*%VecMat
        S<-tcrossprod(VecMat,aG)%*%aR%*%aG%*%VecMat
        S<-Matrix(S,sparse = T)
        
        #head(summary(S))
        #range(summary(S)$x)
        
        aS<-S
        #aS<-as.matrix(aS)
        Smats<-list()
        #names(R)
        for (i in 1:nrow(VecPosition)){
                cumsum(Ranks)
                p<-VecPosition[i,1]
                j<-VecPosition[i,2]
                Smats[[i]]<-S[(cumsum(Ranks)[p]+1-Ranks[p]):cumsum(Ranks)[p],
                              (cumsum(Ranks)[j]+1-Ranks[j]):cumsum(Ranks)[j]] 
                p<-VecR[i,1]
                j<-VecR[i,2]
                names(Smats)[i]<-paste0('S',p,j)
        }
        
        S<-Smats
        names(S)
        #S$S12
        #lapply(S,dim)
        #lapply(S,function(x){sum(x>0)})
        
        
        
        #head(summary(aS))
        #range(summary(aS)$x)
        ##check if each S mat has values 
        SmatValues<-unlist(lapply(S,function(x){sum(x>0)}))
        if(any(SmatValues==0)){stop('S matrix empty')}
        
        #set Ge GD
        SetToZero<-function(Mat){Mat[Mat>0]<-0
        return(Matrix(Mat,sparse = T))}
        
        Ge<-lapply(G,SetToZero)
        head(summary(Ge$G1))
        lapply(Ge,dim)
        
        Gd<-lapply(G,SetToZero)
        
        
        
        #Update Ge Gd using R 
        
        
        GetPus<-function(Mat){
                Mat[Mat<0]<-0
                return(Mat)
        }
        
        GetMinus<-function(Mat){
                Mat<-Mat
                Mat[Mat>0]<-0
                return(abs(Mat))
        }
        
        
        
        for (i in 1:nrow(VecR)){
                #update Ge and Gd
                p<-VecR[i,1]
                j<-VecR[i,2]
                
                Svec<-S[[which(names(S)==paste0('S',p,j))]]
                Rvec<-R[[which(names(R)==paste0('R',p,j))]]
                #head(summary(Rvec))
                #head(summary(Svec))
                
                #Vec1<-Rvec%*%G[[j]]%*%t(Svec)
                Vec1<-tcrossprod(Rvec%*%G[[which(grepl(j,x = names(G)))]],Svec)        
                #Vec2<-Svec%*%t(G[[j]])%*%G[[j]]%*%t(Svec)
                Vec2<-tcrossprod(tcrossprod(Svec,G[[which(grepl(j,x = names(G)))]])%*%G[[which(grepl(j,x = names(G)))]],Svec)
                
                #Vec3<-t(Rvec)%*%G[[p]]%*%Svec
                Vec3<-crossprod(Rvec,G[[which(grepl(p,x = names(G)))]])%*%Svec
                
                #Vec4<-t(Svec)%*%t(G[[p]])%*%G[[p]]%*%Svec      
                Vec4<-crossprod(Svec,t(G[[which(grepl(p,x = names(G)))]]))%*%G[[which(grepl(p,x = names(G)))]]%*%Svec 
                #Vec4<-crossprod(Svec,crossprod(G[[p]],G[[p]]))%*%Svec 
           
#                 head(summary(Vec4))
#                 
#                 Vec1Plus<-GetPus(Vec1)
#                 Vec1minus<-GetMinus(Vec1)
#                 
#                 Vec2Plus<-GetPus(Vec2)
#                 Vec2minus<-GetMinus(Vec2)
#                 
#                 Vec3Plus<-GetPus(Vec3)
#                 Vec3minus<-GetMinus(Vec3)
#                 
#                 Vec4Plus<-GetPus(Vec4)
#                 Vec4minus<-GetMinus(Vec4)
#                 
#                 
#                 dim(Vec2minus)
                
                
#                 Ge[[p]]<-Ge[[p]]+GetPus(Vec1)+G[[p]]%*%GetMinus(Vec2)
#                 Gd[[p]]<-Gd[[p]]+GetMinus(Vec1)+G[[p]]%*%GetPus(Vec2)
                
                Ge[[which(grepl(p,x = names(G)))]]<-Ge[[which(grepl(p,x = names(G)))]]+GetPus(Vec1)+G[[which(grepl(p,x = names(G)))]]%*%GetMinus(Vec2)
                Gd[[which(grepl(p,x = names(G)))]]<-Gd[[which(grepl(p,x = names(G)))]]+GetMinus(Vec1)+G[[which(grepl(p,x = names(G)))]]%*%GetPus(Vec2)
                
#                 Ge[[j]]<-Ge[[j]]+GetPus(Vec3)+G[[j]]%*%GetMinus(Vec4)
#                 Gd[[j]]<-Gd[[j]]+GetMinus(Vec3)+G[[j]]%*%GetPus(Vec4)
                
                Ge[[which(grepl(j,x = names(G)))]]<-Ge[[which(grepl(j,x = names(G)))]]+GetPus(Vec3)+G[[which(grepl(j,x = names(G)))]]%*%GetMinus(Vec4)
                Gd[[which(grepl(j,x = names(G)))]]<-Gd[[which(grepl(j,x = names(G)))]]+GetMinus(Vec3)+G[[which(grepl(j,x = names(G)))]]%*%GetPus(Vec4)
#                 
#                 lapply(Ge,function(x){sum(x==0)/length(x)})
#                 
#                 lapply(Gd,function(x){sum(x==0)/length(x)})
        }
        
        
        
        #apply constrain matrix theta 
        
        for (i in 1:length(Theta)){
                if(length(unlist(Theta[[i]]))>0){
                        for(m in 1:length(Theta[[i]])){
                                Ge[[i]]<- Ge[[i]]+GetMinus(Theta[[i]][[m]])%*%G[[i]]  
                                Gd[[i]]<- Gd[[i]]+GetPus(Theta[[i]][[m]])%*%G[[i]]
                        }
                            
                }
                
        }
#         lapply(Ge,function(x){range(summary(x)$x)})
#         lapply(Ge,function(x){sum(x>0)/length(x)})
#         
#         lapply(Gd,function(x){range(summary(x)$x)})
#         lapply(Gd,function(x){sum(x>0)/length(x)}) 
        
        
        #update G
        VecDiagMat<-list()
        for (i in 1:length(G)){
                Vec<-sqrt(Ge[[i]]/Gd[[i]])
                Vec[is.infinite(Vec)]<-0
                Vec[is.na(Vec)]<-0
                VecDiagMat[[i]]<-Vec
        }
        
        DiagVec<-bdiag(VecDiagMat)
        
#         head(summary(DiagVec))
#         range(summary(DiagVec)$x)
#         sum(is.na(summary(DiagVec)$x))
#         diag(DiagVec)
        
        aG<-aG*(DiagVec)
        
#         head(summary(aG))
#         range(summary(aG))
#         
        
        
        
        
        
        #Break G and S back to list and update  
       
        
        ##G
        VecG<-cbind(1:length(G),1:length(G))
                
        for (i in 1:nrow(VecG)){
                
                p<-VecG[i,1]
                j<-VecG[i,2]
                G[[i]]<-aG[(cumsum(DimensionOfEachSet)[p]+1-DimensionOfEachSet[p]):cumsum(DimensionOfEachSet)[p],
                              (cumsum(Ranks)[j]+1-Ranks[j]):cumsum(Ranks)[j]] 
                
        }
        
#         class(G)
#         G$G3
#         lapply(G, dim)
#         lapply(G, function(x){sum(x>0)})
#         head(summary(G$G3))
        
        
        ##S
#         
#         
#         for (i in 1:nrow(VecR)){
#                 cumsum(Ranks)
#                 p<-VecR[i,1]
#                 j<-VecR[i,2]
#                 S[[i]]<-aS[(cumsum(Ranks)[p]+1-Ranks[p]):cumsum(Ranks)[p],
#                               (cumsum(Ranks)[j]+1-Ranks[j]):cumsum(Ranks)[j]] 
#               
#         }
#         
#         lapply(S,dim)
#         lapply(S, function(x){sum(x>0)})
#         head(summary(S$S12))
#         S$S23
        
        if((Step==1)|(Step%%20)==0){
        p<-as.integer(strsplit(TargetMatrixName,split = '')[[1]][2])
        j<-as.integer(strsplit(TargetMatrixName,split = '')[[1]][3])
        VecS<-S[[which(names(S)==paste0('S',p,j))]]
        MatRecon<-G[[p]]%*%VecS %*%t(G[[j]])
        }
                
                
        if(Step==1){#initialize CostT
                CostT<-norm(TargetMatrix-MatRecon,type = 'F')}
        
        
        
        
        if((Step%%110)==0){
        #calculate the cost function 
        Cost<-0
        for (i in 1:nrow(VecPosition)){
                
                p<-VecPosition[i,1]
                j<-VecPosition[i,2]
                
                Cost<-Cost+(norm(R[[i]]-G[[p]]%*%S[[i]]%*%t(G[[j]]),type = 'F'))^2
#                 dim(R[[i]])
#                 dim(G[[p]]%*%S[[i]]%*%t(G[[j]]))
#                 names(S)
#                 lapply(G, dim)
        }
        
        for (i in 1:length(aTheta)){
                if(length(aTheta[[i]])>0){
                                Vec<-t(aG)%*%aTheta[[i]]%*%aG
                                Vec[which(is.infinite(Vec))]<-0
                                Cost<-Cost+sum(diag(Vec))
                        }
                }
                
        #cost of target matrix 
        
        #cost of target matrix 
        
        VecPre<-CostT
        CostT<-norm(TargetMatrix-MatRecon,type = 'F')
        VecDiff<-abs(VecPre-CostT)
        
        #track the progress
        print(paste('CostTarget=',CostT,'Cost=',Cost))
        if((VecDiff<StopThreshold)|(Step>=MaxInteraction)|((!is.null(MaxTime))&(sum(proc.time())-PTM>MaxTime))){
                ReconErrors<-NULL
                ReReconStructed<-list()
                for (i in 1:length(R)){#'total reconstruction error'
                        p<-as.integer(strsplit(names(R)[i],split = '')[[1]][2])
                        j<-as.integer(strsplit(names(R)[i],split = '')[[1]][3])
                        VecS<-S[[which(names(S)==paste0('S',p,j))]]
                        #ReReconStructed[[i]]<-G[[p]]%*%VecS%*%t(G[[j]])
                        ReReconStructed[[i]]<-G[[which(grepl(p,names(G)))]]%*%VecS%*%t(G[[which(grepl(j,names(G)))]])
                        VecError<-norm(R[[i]]-ReReconStructed[[i]],type = 'F')
                        ReconErrors<-c(ReconErrors,VecError)
                        names(ReconErrors)[i]<-paste0(p,j)
                        names(ReReconStructed)[i]<-paste0('R',p,j)
                        dimnames(ReReconStructed[[i]])<-dimnames(R[[i]])
                }
                
        
                
                
                break}
                }
        #break if < threshold 
        

        }#end of repeat loop
        
        return(list(Rc=ReReconStructed,G=G,S=S,Error=ReconErrors))
        
}#end of function 






