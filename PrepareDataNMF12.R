
#this  version check which version of HiC and eQTL to use 
#this include DHS and TFmotifs information 
#this version include Jason data 
#this version include TAD
#this version include eQTL data from GTEX
#this is preparing matrices for NMF

library(Matrix)
source("AllFunctions4.R")
CellLine<-'E096'



#load roadmap enhancer 

FileList<-list.files("~/BoxSync/research/MIT/datasets/ProducedData/EnhancerByCellTypeRoadMap/")
X<-load(paste0("~/BoxSync/research/MIT/datasets/ProducedData/EnhancerByCellTypeRoadMap/",FileList[grepl(FileList,pattern = CellLine)]))
#X<-load("RMState7_EnhancerLocation_E017.Rdata")
RoadmapEnhancer<-get(X) 
head(RoadmapEnhancer) 
unique(RoadmapEnhancer$Chr)
names(RoadmapEnhancer)<-paste0("Enhacer_",names(RoadmapEnhancer))
#RoadmapEnhancerVec<-RoadmapEnhancer#rearragnge the track for mapping purpose. center of the enhancer and promtoer anrused
#str(RoadmapEnhancerVec)
#RoadmapEnhancerVec@ranges@width<-as.integer(rep(1,length(RoadmapEnhancerVec)))
#head(RoadmapEnhancerVec)
#names(RoadmapEnhancerVec)<-paste0("Enhacer_",names(RoadmapEnhancerVec))


tail(RoadmapEnhancer)

#gencode 

sessionInfo()
#install.packages("refGenome")

myGTF <- "~/BoxSync/research/MIT/datasets/GencodeHg19/gencode.v19.annotation.gtf"

require(rtracklayer)
GenCodeHg19 <- readGFF(myGTF)

head(GenCodeHg19)


GenCodeHg19<-GenCodeHg19[GenCodeHg19$gene_type=="protein_coding" ,]##Only select protein coding 
GenCodeHg19<-GenCodeHg19[GenCodeHg19$type=="gene" ,]##Only select protein coding
nrow(GenCodeHg19)
tail(GenCodeHg19,10)
summary(GenCodeHg19$end-GenCodeHg19$start)


##look at gene level and only looks at chr1-22 x,y
unique(GenCodeHg19$seqid)
GenCodeHg19<-GenCodeHg19[GenCodeHg19$seqid!='chrM',]

##Only get the ealiest starting site for each gene 

GenCodeHg19$EnsembleID<-substring(GenCodeHg19$gene_id,first = 1,last = 15)
nrow(GenCodeHg19)==length(unique(GenCodeHg19$EnsembleID))#in this set, only focus on genes 

head(GenCodeHg19,20)
#Miega<-aggregate(start~EnsembleID,data=GenCodeHg19,min)


GenCodeHg19$TTS<-GenCodeHg19$start
GenCodeHg19$TTS[GenCodeHg19$strand=='-']<-GenCodeHg19$end[GenCodeHg19$strand=='-']

GenCodeHg19$PromoterStart[GenCodeHg19$strand=='+']<-GenCodeHg19$TTS[GenCodeHg19$strand=='+']-2000#TSS upstream 2k and downstream 500bp are considered promtoer
GenCodeHg19$PromoterStart[GenCodeHg19$PromoterStart<=0]<-1

GenCodeHg19$PromoterEnd[GenCodeHg19$strand=='+']<-GenCodeHg19$TTS[GenCodeHg19$strand=='+']+500



GenCodeHg19$PromoterStart[GenCodeHg19$strand=='-']<-GenCodeHg19$TTS[GenCodeHg19$strand=='-']-500#TSS upstream 2k and downstream 500bp are considered promtoer
GenCodeHg19$PromoterStart[GenCodeHg19$PromoterStart<=0]<-1

GenCodeHg19$PromoterEnd[GenCodeHg19$strand=='-']<-GenCodeHg19$TTS[GenCodeHg19$strand=='-']+2000


head(GenCodeHg19)

GencodeTrack<-makeGRangesFromDataFrame(data.frame(chr=GenCodeHg19$seqid,Start=GenCodeHg19$PromoterStart,End=GenCodeHg19$PromoterEnd))
GencodeTrack$GeneNames<-GenCodeHg19$EnsembleID
names(GencodeTrack)<-paste0('Promoter_',GenCodeHg19$EnsembleID,"_",GenCodeHg19$seqid,'_',
                                                           GenCodeHg19$PromoterStart,"_",GenCodeHg19$PromoterEnd)


head(GenCodeHg19)
# GencodeTrackVec<-GencodeTrack#Modify the track for mapping (center)
# GencodeTrackVec@ranges@width<-as.integer(rep(1,length(GencodeTrackVec)))
# names(GencodeTrackVec)<-paste0('Promoter_',names(GencodeTrackVec),"_",GenCodeHg19$seqid,'_',
#                                GenCodeHg19$PromoterStart,"_",GenCodeHg19$PromoterEnd)
# head(GenCodeHg19)
nrow(GenCodeHg19)==length(unique(GenCodeHg19$EnsembleID))
head(GenCodeHg19$EnsembleID)
GencodeTrack$TSS<-GenCodeHg19$TTS





#Get HiC data 


FileList<-list.files("~/BoxSync/research/MIT/datasets/ProducedData/HiCRao/")
FileVec<-FileList[grepl(FileList,pattern = CellLine)]

if(length(FileVec)==1){#this Rao data has this HiC
X<-load(file = paste0("~/BoxSync/research/MIT/datasets/ProducedData/HiCRao/",FileVec))
HiC<-get(X)
head(summary(HiC))
dimnames(HiC)


Chr<-mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][1]})
Start<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][2]}))
End<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][3]}))


library(GenomicRanges)

HiCTrack<-makeGRangesFromDataFrame(data.frame(Chr,Start,End))




#make a enhancer-Hic anchor , a promoter-HiC anchor mats 

Target<-RoadmapEnhancer
TargetCenter<-Target

TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
head(TargetCenter)

VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
head(VecOverLap)


Enhancer_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                dimnames = list(names(Target),names(HiCTrack)))
dim(Enhancer_HiC)



Target<-GencodeTrack
TargetCenter<-Target

TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
head(TargetCenter)

VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
head(VecOverLap)


Promoter_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                           dimnames = list(names(Target),names(HiCTrack)))

dim(Promoter_HiC)}


FileList<-list.files("~/BoxSync/research/MIT/datasets/ProducedData/HiC21/",pattern = 'Organised')
FileVec2<-FileList[grepl(FileList,pattern = CellLine)]

if(length(FileVec)==0&length(FileVec2)==1){#if Hi-C is available in 21HiC
        X<-load(file = paste0("~/BoxSync/research/MIT/datasets/ProducedData/HiC21/",FileVec2))
        HiC<-get(X)
        head(summary(HiC))
        HiC
        dimnames(HiC)
        
        
        Chr<-mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][1]})
        Start<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][2]}))
        End<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][3]}))
        
        
        library(GenomicRanges)
        
        HiCTrack<-makeGRangesFromDataFrame(data.frame(Chr,Start,End))
        
        
        
        
        #make a enhancer-Hic anchor , a promoter-HiC anchor mats 
        
        Target<-RoadmapEnhancer
        TargetCenter<-Target
        
        TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
        TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
        head(TargetCenter)
        
        VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
        head(VecOverLap)
        
        
        Enhancer_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                                   dimnames = list(names(Target),names(HiCTrack)))
        dim(Enhancer_HiC)
        
        
        
        Target<-GencodeTrack
        TargetCenter<-Target
        
        TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
        TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
        head(TargetCenter)
        
        VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
        head(VecOverLap)
        
        
        Promoter_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                                   dimnames = list(names(Target),names(HiCTrack)))  
        dim(Promoter_HiC)
}



if(length(FileVec)==0&length(FileVec2)==1){#if Hi-C is available in 21HiC
        X<-load(file = paste0("~/BoxSync/research/MIT/datasets/ProducedData/HiC21/",FileVec2))
        HiC<-get(X)
        head(summary(HiC))
        HiC
        dimnames(HiC)
        
        
        Chr<-mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][1]})
        Start<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][2]}))
        End<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][3]}))
        
        
        library(GenomicRanges)
        
        HiCTrack<-makeGRangesFromDataFrame(data.frame(Chr,Start,End))
        
        
        
        
        #make a enhancer-Hic anchor , a promoter-HiC anchor mats 
        
        Target<-RoadmapEnhancer
        TargetCenter<-Target
        
        TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
        TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
        head(TargetCenter)
        
        VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
        head(VecOverLap)
        
        
        Enhancer_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                                   dimnames = list(names(Target),names(HiCTrack)))
        dim(Enhancer_HiC)
        
        
        
        Target<-GencodeTrack
        TargetCenter<-Target
        
        TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
        TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
        head(TargetCenter)
        
        VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
        head(VecOverLap)
        
        
        Promoter_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                                   dimnames = list(names(Target),names(HiCTrack)))  
        dim(Promoter_HiC)
}


if(length(FileVec)==0&length(FileVec2)==0){#if Hi-C is unvailable for that cell type at all
        X<-load(file = paste0("~/BoxSync/research/MIT/datasets/ProducedData/HiC21/Organised_Union.Rdata"))
        HiC<-get(X)
        head(summary(HiC))
        HiC
        dimnames(HiC)
        
        
        Chr<-mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][1]})
        Start<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][2]}))
        End<-as.numeric(mapply(dimnames(HiC)[[1]],FUN = function(x){strsplit(x,split = '_')[[1]][3]}))
        
        
        library(GenomicRanges)
        
        HiCTrack<-makeGRangesFromDataFrame(data.frame(Chr,Start,End))
        
        
        
        
        #make a enhancer-Hic anchor , a promoter-HiC anchor mats 
        
        Target<-RoadmapEnhancer
        TargetCenter<-Target
        
        TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
        TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
        head(TargetCenter)
        
        VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
        head(VecOverLap)
        
        
        Enhancer_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                                   dimnames = list(names(Target),names(HiCTrack)))
        dim(Enhancer_HiC)
        
        
        
        Target<-GencodeTrack
        TargetCenter<-Target
        
        TargetCenter@ranges@start<-as.integer(round(Target@ranges@start+Target@ranges@width/2))
        TargetCenter@ranges@width<-as.integer(rep(1,length(Target)))##center of enhancer and promoter 
        head(TargetCenter)
        
        VecOverLap<-findOverlaps(query =TargetCenter,subject = HiCTrack )
        head(VecOverLap)
        
        
        Promoter_HiC<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims =c(length(Target),length(HiCTrack)),
                                   dimnames = list(names(Target),names(HiCTrack)))  
        dim(Promoter_HiC)
}


# #Map EP data to JW EP links 
# 
# RoadmapEnhancer
# head(RoadmapEnhancer)
# EntireEP<-sparseMatrix(i = 1,j = 1,x = 0,dims = c(length(RoadmapEnhancer),length(GencodeTrack)),
#                        dimnames =list(names(RoadmapEnhancer),names(GencodeTrack)) )
# 
# dim(EntireEP)
# 
# ##JW links 
# 
# load(file = paste0('~/BoxSync/research/MIT/datasets/ProducedData/JWlinks/',
#                    CellLine,'.Rdata'))
# 
# 
# EPLinks<-EPmats
# 
# #EPEdgeList<-JW
# #head(EPEdgeList)
# #length(EPLinks)
# dim(EPLinks)
# library(GenomicRanges)
# 
# 
# 
# 
# EnhancerAnchors<-makeGRangesFromDataFrame(data.frame(chr=mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][2]}),
#                                                      start=as.numeric(mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][3]})),
#                                                      end=as.numeric(mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][4]}))))
# 
# EnhancerAnchors@ranges@start<-as.integer(round(EnhancerAnchors@ranges@start+EnhancerAnchors@ranges@width/2))
# EnhancerAnchors@ranges@width<-as.integer(rep(1,length(EnhancerAnchors)))#map the center to my definition
# 
# 
# dim(EPLinks)
# 
# colnames(EPLinks)
# 
# 
# ##Map the anchors , map enhancer using overlap and genes using gene names , ie, which e-p links are supported by JW
# 
# VecOverLap<-findOverlaps(query = EnhancerAnchors,subject =RoadmapEnhancer )
# head(VecOverLap)
# length(VecOverLap)==length(unique(VecOverLap@from))#check whether one to one match 
# VecOverLap@from
# 
# 
# Pvec<-substr(colnames(EPLinks),start = 10,stop = 24)
# VecMatch<-match(Pvec,GencodeTrack$GeneNames)##there are some consistency in ensemnl IDs ignore those missing matches
# 
# dim(VecMatch)
# VecMatch<-data.frame(1:length(Pvec),VecMatch)
# VecMatch<-VecMatch[!is.na(VecMatch[,2]),]
# head(VecMatch)
# 
# 
# InterSupportedByEP<-EntireEP
# range(VecMatch[,1])
# dim(EPLinks)
# range(VecMatch[,2])
# range(VecOverLap@from)
# 
# InterSupportedByEP[VecOverLap@to,VecMatch[,2]]<-
#         EPmats[VecOverLap@from,VecMatch[,1]]
# head(summary(InterSupportedByEP))
# 
# 
# JWEPremapped<-InterSupportedByEP
# 
# head(summary(JWEPremapped))


#Map EP data to JE EP links 

RoadmapEnhancer
head(RoadmapEnhancer)
EntireEP<-sparseMatrix(i = 1,j = 1,x = 0,dims = c(length(RoadmapEnhancer),length(GencodeTrack)),
                       dimnames =list(names(RoadmapEnhancer),names(GencodeTrack)) )

##JE links 
 load(file = paste0('~/BoxSync/research/MIT/datasets/ProducedData/JElinks/',
                    CellLine,'.Rdata'))

EPLinks<-JElinks


dim(EPLinks)
colnames(EPLinks)
library(GenomicRanges)




EnhancerAnchors<-makeGRangesFromDataFrame(data.frame(chr=mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][1]}),
                                                     start=as.numeric(mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][2]})),
                                                     end=as.numeric(mapply(rownames(EPLinks),FUN = function(x){strsplit(x,'_')[[1]][3]}))))

EnhancerAnchors@ranges@start<-as.integer(round(EnhancerAnchors@ranges@start+EnhancerAnchors@ranges@width/2))
EnhancerAnchors@ranges@width<-as.integer(rep(1,length(EnhancerAnchors)))#map the center to my definition


dim(EPLinks)

colnames(EPLinks)


##Map the anchors , map enhancer using overlap and genes using gene names , ie, which e-p links are supported by JW

VecOverLap<-findOverlaps(query = EnhancerAnchors,subject =RoadmapEnhancer )
head(VecOverLap)
length(VecOverLap)==length(unique(VecOverLap@from))#check whether one to one match 
VecOverLap@from


Pvec<-substr(colnames(EPLinks),start = 1,stop = 15)
VecMatch<-match(Pvec,GencodeTrack$GeneNames)##there are some consistency in ensemnl IDs ignore those missing matches

dim(VecMatch)
VecMatch<-data.frame(1:length(Pvec),VecMatch)
VecMatch<-VecMatch[!is.na(VecMatch[,2]),]
head(VecMatch)


InterSupportedByEP<-EntireEP
range(VecMatch[,1])
dim(EPLinks)
range(VecMatch[,2])
range(VecOverLap@from)

InterSupportedByEP[VecOverLap@to,VecMatch[,2]]<-
        EPLinks[VecOverLap@from,VecMatch[,1]]
head(summary(InterSupportedByEP))


JEEPremapped<-InterSupportedByEP

head(summary(JEEPremapped))








#eQTL information 
##Map eQTL to promoters (by gene name)
FileList<-list.files(path = "~/BoxSync/research/MIT/datasets/ProducedData/GTEX_v6p_eQTL/")
VecFile<-FileList[grepl(FileList,pattern = CellLine)]
if(length(VecFile)==1){

File<-paste0("~/BoxSync/research/MIT/datasets/ProducedData/GTEX_v6p_eQTL/",VecFile)
X<-load(File)
eQTL_GTEX<-get(X)
head(summary(eQTL_GTEX))

VecSum<-summary(eQTL_GTEX)
VecSum<-VecSum[VecSum$x<=1e-5,]#only focus on those wiht p val <1e-5 because I am working on binary data 
eQTL_GTEX<-sparseMatrix(i = VecSum$i,j = VecSum$j,x = VecSum$x,dims = dim(eQTL_GTEX),dimnames = dimnames(eQTL_GTEX))

VecGene<-mapply(colnames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][2]})
head(GencodeTrack)
eQTLSNP_Gene<-sparseMatrix(i = 1,j = 1,x = 0,dims = c(nrow(eQTL_GTEX),length(GencodeTrack)),
                           dimnames = list(rownames(eQTL_GTEX),names(GencodeTrack)))
VecMatch<-match(GencodeTrack$GeneNames,VecGene)
eQTLSNP_Gene[,which(!is.na(VecMatch))]<-eQTL_GTEX[,VecMatch[!is.na(VecMatch)]]
colnames(eQTLSNP_Gene)
head(summary(eQTLSNP_Gene))


##map SNP to enhancer in the cell type 
library(GenomicRanges)
rownames(eQTL_GTEX)
chr<-mapply(rownames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][2]})
Start<-as.numeric(mapply(rownames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][3]}))
End<-as.numeric(mapply(rownames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][4]}))
        
eQTL_SNP_Track<-makeGRangesFromDataFrame(df = data.frame(chr,Start,End))       

VecOverLap<-findOverlaps(query = eQTL_SNP_Track,subject = RoadmapEnhancer)

eQTLSNP_RoadmapEnhancer<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,
                                      dim=c(length(eQTL_SNP_Track),length(RoadmapEnhancer)),
                                      dimnames = list(names(eQTL_SNP_Track),names(RoadmapEnhancer)))
head(summary(eQTLSNP_RoadmapEnhancer))
}

if(length(VecFile)!=1){#no match cell type in eQTL
        
        
        File<-paste0("~/BoxSync/research/MIT/datasets/ProducedData/GTEX_v6p_eQTL/UnionGTEX.RData")
        X<-load(File)
        eQTL_GTEX<-get(X)
        head(summary(eQTL_GTEX))
        
        VecSum<-summary(eQTL_GTEX)
        VecSum<-VecSum[VecSum$x<=1e-5,]#only focus on those wiht p val <1e-5 because I am working on binary data 
        eQTL_GTEX<-sparseMatrix(i = VecSum$i,j = VecSum$j,x = VecSum$x,dims = dim(eQTL_GTEX),dimnames = dimnames(eQTL_GTEX))
        
        VecGene<-mapply(colnames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][2]})
        head(GencodeTrack)
        eQTLSNP_Gene<-sparseMatrix(i = 1,j = 1,x = 0,dims = c(nrow(eQTL_GTEX),length(GencodeTrack)),
                                   dimnames = list(rownames(eQTL_GTEX),names(GencodeTrack)))
        VecMatch<-match(GencodeTrack$GeneNames,VecGene)
        eQTLSNP_Gene[,which(!is.na(VecMatch))]<-eQTL_GTEX[,VecMatch[!is.na(VecMatch)]]
        colnames(eQTLSNP_Gene)
        head(summary(eQTLSNP_Gene))
        
        
        ##map SNP to enhancer in the cell type 
        library(GenomicRanges)
        rownames(eQTL_GTEX)
        chr<-mapply(rownames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][2]})
        Start<-as.numeric(mapply(rownames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][3]}))
        End<-as.numeric(mapply(rownames(eQTL_GTEX),FUN = function(x){strsplit(x,split = '_')[[1]][4]}))
        
        eQTL_SNP_Track<-makeGRangesFromDataFrame(df = data.frame(chr,Start,End))       
        
        VecOverLap<-findOverlaps(query = eQTL_SNP_Track,subject = RoadmapEnhancer)
        
        eQTLSNP_RoadmapEnhancer<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,
                                              dim=c(length(eQTL_SNP_Track),length(RoadmapEnhancer)),
                                              dimnames = list(names(eQTL_SNP_Track),names(RoadmapEnhancer)))
        head(summary(eQTLSNP_RoadmapEnhancer))
        
}


#Dixon TAD data

File<- "~/BoxSync/research/MIT/datasets/ProducedData/TAD_Dixon/E017.RData"#use IMR90 for all as TAD are quite simiar 
X<-load(File)
TAD_Dixon<-get(X)
dim(TAD_Dixon)
head(TAD_Dixon)
VecTAD<-data.frame(Chr=mapply(rownames(TAD_Dixon),FUN = function(x){strsplit(x,split = '_')[[1]][2]}),
                   Start=as.numeric(mapply(rownames(TAD_Dixon),FUN = function(x){strsplit(x,split = '_')[[1]][3]})),
                   End=as.numeric(mapply(rownames(TAD_Dixon),FUN = function(x){strsplit(x,split = '_')[[1]][4]})))     
head(VecTAD)        
library(GenomicRanges)
TAD_Dixon_Track<-makeGRangesFromDataFrame(df = VecTAD)

##Map Enhancer to TAD
EnhancerCenter<-RoadmapEnhancer
EnhancerCenter@ranges@start<-as.integer(round(EnhancerCenter@ranges@start+EnhancerCenter@ranges@width/2))
EnhancerCenter@ranges@width<-as.integer(rep(1,length(EnhancerCenter)))

VecOverLap<-findOverlaps(query = EnhancerCenter,subject = TAD_Dixon_Track)

Enhancer_TAD<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims = c(length(EnhancerCenter),length(TAD_Dixon_Track)),
                           dimnames = list(names(EnhancerCenter),names(TAD_Dixon_Track)))
head(summary(Enhancer_TAD))


##Map promoter to TAD
Promoter_TSS<-GencodeTrack#use TSS for overlap
Promoter_TSS@ranges@start<-GencodeTrack$TSS
Promoter_TSS@ranges@width<-as.integer(rep(1,length(Promoter_TSS)))


VecOverLap<-findOverlaps(query = Promoter_TSS,subject = TAD_Dixon_Track)

Promoter_TAD<-sparseMatrix(i = VecOverLap@from,j = VecOverLap@to,x = 1,dims = c(length(Promoter_TSS),length(TAD_Dixon_Track)),
                           dimnames = list(names(Promoter_TSS),names(TAD_Dixon_Track)))
head(summary(Promoter_TAD))






#DNaseI hypersnesittivity 

#

##Thurman
X<-load('~/BoxSync/research/MIT/datasets/ProducedData/DHSThurman/DHSThurman.RData')
DHSThurman<-get(X)
range(summary(DHSThurman)$x)
nrow(DHSThurman)
##Make the matrix binary 
VecSum<-summary(DHSThurman)

DHSThurman<-sparseMatrix(i = VecSum$i,j = VecSum$j,x = 1,dims = dim(DHSThurman),
                                dimnames =dimnames(DHSThurman) )
head(summary(DHSThurman))
colnames(DHSThurman)<-mapply(colnames(DHSThurman),FUN = function(x){strsplit(x,'_')[[1]][1]})
rownames(DHSThurman)
##DHSThurman_Gene

VecGene<-colnames(DHSThurman)

head(summary(DHSThurman))

RowNodes<-rownames(DHSThurman)
ColNodes<-GencodeTrack$GeneNames

VecRow<-match(rownames(DHSThurman)[summary(DHSThurman)$i],RowNodes)
VecCol<-match(colnames(DHSThurman)[summary(DHSThurman)$j],ColNodes)
VecPos<-data.frame(VecRow,VecCol)
VecPos<-VecPos[(!is.na(VecPos$VecRow))&(!is.na(VecPos$VecCol)),]

DHSThurman_Gene<-sparseMatrix(i = VecPos$VecRow,j = VecPos$VecCol,x = 1,dims = c(nrow(DHSThurman),length(GencodeTrack)),
                           dimnames = list(rownames(DHSThurman),names(GencodeTrack)))
head(summary(DHSThurman_Gene))


##DHS_enhancer
TrackVec<-rownames(DHSThurman)
TrackVec<-strsplit(rownames(DHSThurman),split = '_')
TrackVec<-matrix(unlist(TrackVec),ncol = 4,byrow = T)
TrackVec<-as.data.frame(TrackVec)
names(TrackVec)<-c('Type','Chr','Start','End')
TrackVec$Start<-as.numeric(as.character(TrackVec$Start))
TrackVec$End<-as.numeric(as.character(TrackVec$End))

summary(TrackVec$End-TrackVec$Start)
DHSThurmanTrack<-makeGRangesFromDataFrame(df = TrackVec)
names(DHSThurmanTrack)<-rownames(DHSThurman)

DHSThurmanTrackCenter<-DHSThurmanTrack
DHSThurmanTrackCenter@ranges@start<-as.integer(round(DHSThurmanTrackCenter@ranges@start+DHSThurmanTrackCenter@ranges@width/2))
DHSThurmanTrackCenter@ranges@width<-as.integer(rep(1,length(DHSThurmanTrack)))


VecOverLap<-findOverlaps(query =  DHSThurmanTrackCenter,subject =RoadmapEnhancer)

VecPos<-data.frame(VecOverLap@from,VecOverLap@to)
VecPos<-VecPos[(!is.na(VecPos$VecOverLap.from))&(!is.na(VecPos$VecOverLap.to)),]
DHSThurman_Enhancer<-sparseMatrix(i =VecPos$VecOverLap.from,j = VecPos$VecOverLap.to,x = 1,
                                        dims = c(length(DHSThurmanTrack),length(RoadmapEnhancer)),
                                  dimnames = list(names(DHSThurmanTrack),names(RoadmapEnhancer)),use.last.ij = T)

head(summary(DHSThurman_Enhancer))



# 
# #Marbach TF motifs 
# 
# X<-load('~/BoxSync/research/MIT/datasets/ProducedData/MarbachTFmortifs/TFmotifs.Rdata')
# MarbachTF<-get(X)
# head(MarbachTF)
# 
# MarbachTF$name<-gsub(pattern='(.*)[_](.*)','\\1',MarbachTF$name)
# MarbachTF$name<-paste0(MarbachTF$name,'_',rep(MarbachTF@seqnames@values,MarbachTF@seqnames@lengths))
# head(MarbachTF)
# MarbachTFCenter<-MarbachTF
# MarbachTFCenter@ranges@start<-as.integer(round(MarbachTFCenter@ranges@start+MarbachTFCenter@ranges@width/2))
# MarbachTFCenter@ranges@width<-as.integer(rep(1,length(MarbachTF)))  
# 
# ##TF_Enhancer
# VecOverlap<-findOverlaps(query =MarbachTFCenter,subject = RoadmapEnhancer )
# RowNodes<-unique(MarbachTF$name)
# RowID<-match(MarbachTF$name[VecOverlap@from],RowNodes)
# 
# 
# MarbachTF_enhancer<-sparseMatrix(i =RowID,j = VecOverlap@to,x = 1,dims = c(length(RowNodes),length(RoadmapEnhancer)),
#                                  dimnames = list(RowNodes,names(RoadmapEnhancer)),use.last.ij = T)
# 
# head(summary(MarbachTF_enhancer))
# 
# 
# ##TF_Enhancer
# VecOverlap<-findOverlaps(query =MarbachTFCenter,subject = GencodeTrack )
# RowNodes<-unique(MarbachTF$name)
# RowID<-match(MarbachTF$name[VecOverlap@from],RowNodes)
# 
# 
# MarbachTF_Gene<-sparseMatrix(i =RowID,j = VecOverlap@to,x = 1,dims = c(length(RowNodes),length(GencodeTrack)),
#                                  dimnames = list(RowNodes,names(GencodeTrack)),use.last.ij = T)
# 
# head(summary(MarbachTF_Gene))



#Save all matrice

MatricesForNMF<-list(HiC,Promoter_HiC,Enhancer_HiC,JEEPremapped,eQTLSNP_Gene,eQTLSNP_RoadmapEnhancer,
                     TAD_Dixon,Enhancer_TAD,Promoter_TAD,DHSThurman_Gene,DHSThurman_Enhancer)

lapply(MatricesForNMF,dim)

names(MatricesForNMF)<-c('HiC','Promoter_HiC','Enhancer_HiC','JEEPremapped',
                         'eQTLSNP_Gene','eQTLSNP_RoadmapEnhancer','TAD_Dixon','Enhancer_TAD','Promoter_TAD',
                         'DHSThurman_Gene','DHSThurman_Enhancer')

save(MatricesForNMF,file = paste0('~/BoxSync/research/MIT/datasets/ProducedData/RawMatrices/MatricesForNMF_',CellLine,'.Rdata'))
