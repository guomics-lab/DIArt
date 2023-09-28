#file matrices####
diann_hye_pre<-read.csv("diann_hye_pre.csv",stringsAsFactors = F,row.names = 1)
diann_hyc_pre<-read.csv("diann_hyc_pre.csv",stringsAsFactors = F,row.names = 1)
diann_hye_pro<-read.csv("diann_hye_pro.csv",stringsAsFactors = F,row.names = 1)
diann_hyc_pro<-read.csv("diann_hyc_pro.csv",stringsAsFactors = F,row.names = 1)
diart_hye_pre<-read.csv("diart_hye_pre.csv",stringsAsFactors = F,row.names = 1)
diart_hyc_pre<-read.csv("diart_hyc_pre2.csv",stringsAsFactors = F,row.names = 1)
diart_hye_pro<-read.csv("diart_hye_pro.csv",stringsAsFactors = F,row.names = 1)
diart_hyc_pro<-read.csv("diart_hyc_pro.csv",stringsAsFactors = F,row.names = 1)

diart_hyc_pro_pe<-read.delim("prot2023-09-12.txt",stringsAsFactors = F,sep='\t',row.names = 1)
diart_hyc_pro<-2^diart_hyc_pro_pe
row.names(diart_hyc_pro)<-sapply(sapply(row.names(diart_hyc_pro_pe),strsplit,"\\|"),"[[",1)
diart_hyc_proIn3<-read.delim("orbitrap_protein_quant_fillna.tsv",stringsAsFactors = F,sep='\t',row.names = 1)
diart_hyc_pro<-diart_hyc_proIn3
row.names(diart_hyc_pro)<-sapply(sapply(row.names(diart_hyc_proIn3),strsplit,"\\|"),"[[",1)

CV_diann_hye_pre<-c(apply(diann_hye_pre[,c(1,3,5)],1,sd,na.rm=T)/apply(diann_hye_pre[,c(1,3,5)],1,mean,na.rm=T),
                    apply(diann_hye_pre[,c(2,4,6)],1,sd,na.rm=T)/apply(diann_hye_pre[,c(2,4,6)],1,mean,na.rm=T))
CV_diart_hye_pre<-c(apply(diart_hye_pre[,c(1,3,5)],1,sd,na.rm=T)/apply(diart_hye_pre[,c(1,3,5)],1,mean,na.rm=T),
                    apply(diart_hye_pre[,c(2,4,6)],1,sd,na.rm=T)/apply(diart_hye_pre[,c(2,4,6)],1,mean,na.rm=T))
CV_diann_hye_pro<-c(apply(diann_hye_pro[,c(1,3,5)],1,sd,na.rm=T)/apply(diann_hye_pro[,c(1,3,5)],1,mean,na.rm=T),
                    apply(diann_hye_pro[,c(2,4,6)],1,sd,na.rm=T)/apply(diann_hye_pro[,c(2,4,6)],1,mean,na.rm=T))
CV_diart_hye_pro<-c(apply(diart_hye_pro[,c(1,3,5)],1,sd,na.rm=T)/apply(diart_hye_pro[,c(1,3,5)],1,mean,na.rm=T),
                    apply(diart_hye_pro[,c(2,4,6)],1,sd,na.rm=T)/apply(diart_hye_pro[,c(2,4,6)],1,mean,na.rm=T))
CV_diann_hyc_pre<-c(apply(diann_hyc_pre[,c(1:3)],1,sd,na.rm=T)/apply(diann_hyc_pre[,c(1:3)],1,mean,na.rm=T),
                    apply(diann_hyc_pre[,c(4:6)],1,sd,na.rm=T)/apply(diann_hyc_pre[,c(4:6)],1,mean,na.rm=T))
CV_diart_hyc_pre<-c(apply(diart_hyc_pre[,c(1:3)],1,sd,na.rm=T)/apply(diart_hyc_pre[,c(1:3)],1,mean,na.rm=T),
                    apply(diart_hyc_pre[,c(4:6)],1,sd,na.rm=T)/apply(diart_hyc_pre[,c(4:6)],1,mean,na.rm=T))
CV_diann_hyc_pro<-c(apply(diann_hyc_pro[,c(1:3)],1,sd,na.rm=T)/apply(diann_hyc_pro[,c(1:3)],1,mean,na.rm=T),
                    apply(diann_hyc_pro[,c(4:6)],1,sd,na.rm=T)/apply(diann_hyc_pro[,c(4:6)],1,mean,na.rm=T))
CV_diart_hyc_pro<-c(apply(diart_hyc_pro[,c(1:3)],1,sd,na.rm=T)/apply(diart_hyc_pro[,c(1:3)],1,mean,na.rm=T),
                    apply(diart_hyc_pro[,c(4:6)],1,sd,na.rm=T)/apply(diart_hyc_pro[,c(4:6)],1,mean,na.rm=T))

library(vioplot)
CV_hye_pre<-list(CV_diann_hye_pre,CV_diart_hye_pre)
CV_hye_pro<-list(CV_diann_hye_pro,CV_diart_hye_pro)
CV_hyc_pre<-list(CV_diann_hyc_pre,CV_diart_hyc_pre)
CV_hyc_pro<-list(CV_diann_hyc_pro,CV_diart_hyc_pro)
names(CV_hye_pre)<-names(CV_hye_pro)<-names(CV_hyc_pre)<-names(CV_hyc_pro)<-c("","")

colToPlot<-c("#192e57","#f49d17")
pdf("Figure3A_CV_v10.pdf",width=12,height=3)
par(mfrow=c(1,4),mar=c(3,4,3,3))
vioplot(CV_hye_pre,col=colToPlot,las=2,main="TOF-HYE-Precursor",ylim=c(0,2))
text(c(1,2),1.85,round(sapply(CV_hye_pre,median,na.rm=T),3),col = 2)
vioplot(CV_hye_pro,col=colToPlot,las=2,main="TOF-HYE-Protein",ylim=c(0,2))
text(c(1,2),1.85,round(sapply(CV_hye_pro,median,na.rm=T),3),col = 2)
vioplot(CV_hyc_pre,col=colToPlot,las=2,main="Orbi-HYC-Precursor",ylim=c(0,2))
text(c(1,2),1.85,round(sapply(CV_hyc_pre,median,na.rm=T),3),col = 2)
vioplot(CV_hyc_pro,col=colToPlot,las=2,main="Orbi-HYC-Protein",ylim=c(0,2))
text(c(1,2),1.85,round(sapply(CV_hyc_pro,median,na.rm=T),3),col = 2)
dev.off()

#LFQ####
hye_spec<-read.table("ecolihumanyeast_concat_mayu_IRR_cons_openswath_64var_curated.tsv"
                     ,sep = '\t',header = T,row.names = 1)
hye_specSpecies<-unique(hye_spec[,c("PeptideGroupLabel","UniprotId")])
for(i in 1:nrow(hye_specSpecies)){
  sss<-unlist(strsplit(hye_specSpecies$UniprotId[i],"_"))
  ssss<-unlist(strsplit(hye_specSpecies$UniprotId[i],"\\|"))
  hye_specSpecies$ProteinName[i]<-ssss[1]
  hye_specSpecies$species[i]<-sss[length(sss)]
  print(i)
}
hye_trans_spec<-hye_specSpecies$species
names(hye_trans_spec)<-gsub("_","",hye_specSpecies$PeptideGroupLabel)
hye_trans_spec<-gsub("YEAS8","YEAST",hye_trans_spec)

hye_prot_spec<-hye_specSpecies$species
names(hye_prot_spec)<-hye_specSpecies$ProteinName
hye_prot_spec<-hye_prot_spec[!duplicated(names(hye_prot_spec))]
hye_prot_spec<-gsub("YEAS8","YEAST",hye_prot_spec)

orbit_spec<-read.table("../HYC_orbit_DIArt/orbitrap_speclib.tsv",stringsAsFactors = F,sep='\t',header=T)
orbit_specSpecies<-unique(orbit_spec[,c("transition_group_id","ProteinName","UniProtIds")])
for(i in 1:nrow(orbit_specSpecies)){
  sss<-unlist(strsplit(orbit_specSpecies$ProteinName[i],"_"))
  orbit_specSpecies$species[i]<-sss[length(sss)]
  print(i)
}
hyc_trans_spec<-orbit_specSpecies$species
names(hyc_trans_spec)<-paste0(sapply(sapply(orbit_specSpecies$transition_group_id,strsplit,"_"),"[[",2),
                              sapply(sapply(orbit_specSpecies$transition_group_id,strsplit,"_"),"[[",3))

hyc_prot_spec<-orbit_specSpecies$species
names(hyc_prot_spec)<-orbit_specSpecies$UniProtIds
hyc_prot_spec<-hyc_prot_spec[!duplicated(names(hyc_prot_spec))]

# diann_hye_pre[,c(2,4,6)]
colSpecies<-c(7,5,4,3)
names(colSpecies)<-c("CAEEL","ECOLI","YEAST","HUMAN")

library(ggplot2)
library(viridis)
library(ggpubr)

hyePreProList<-list(diann_hye_pre,diart_hye_pre[inter_hye_pre,],diann_hye_pro,diart_hye_pro[inter_hye_pro,])
g1<-g2<-g3<-list()
for(m in 1:4){
  readRat<-hyePreProList[[m]]
  if(m<3){
    readRat$specAll<-hye_trans_spec[row.names(readRat)]
    sizePoint=0.1
    binSize=300
    readRat$Amean<-apply(readRat[,c(1,3,5)],1,mean,na.rm=F)
    readRat$Bmean<-apply(readRat[,c(2,4,6)],1,mean,na.rm=F)
    readRat$ABratio<-readRat$Amean/readRat$Bmean
    
  }else{
    readRat$specAll<-hye_prot_spec[row.names(readRat)]
    sizePoint=1
    binSize=120
    readRat$Amean<-apply(readRat[,c(1,3,5)],1,mean,na.rm=F)
    readRat$Bmean<-apply(readRat[,c(2,4,6)],1,mean,na.rm=F)
    readRat$ABratio<-readRat$Amean/readRat$Bmean
  }
  
  g1[[m]]<-ggplot(NULL, aes(x, y)) +
    geom_hline(yintercept=1,color=c('#3CB875'))+
    geom_hex(data=readRat[readRat$specAll=="YEAST",]
             , mapping = aes(x=log10(Bmean), y=log2(ABratio))
             ,bins= binSize*0.6,size=sizePoint,linewidth=0.1,na.rm=T,show.legend=F) +
    scale_fill_viridis(option = "D",direction = -1,begin=0,end=0.8,alpha = 0.9)+
    scale_x_continuous(name="log10(B)", limits=c(2,8))+ 
    scale_y_continuous(name="log2(A/B)", limits=c(-4, 4))+ 
    theme_minimal()+theme(panel.grid.minor=element_blank())
  g2[[m]]<-ggplot(NULL, aes(x, y)) +
    geom_hline(yintercept=0,color=c('#DE4B68'))+
    geom_hex(data=readRat[readRat$specAll=="HUMAN",], mapping = aes(x=log10(Bmean), y=log2(ABratio))
             ,bins= binSize,size=sizePoint,linewidth=0.1,na.rm=T,show.legend=F) +
    scale_fill_viridis(option = "A",direction = -1,begin=0,end=0.85,alpha = 0.9)+
    scale_x_continuous(name="log10(B)", limits=c(2,8))+ 
    scale_y_continuous(name="log2(A/B)", limits=c(-4, 4))+ 
    theme_minimal()+theme(panel.grid.minor=element_blank())
  g3[[m]]<-ggplot(NULL, aes(x, y)) +
    geom_hline(yintercept=-2,color=c('#369DA9'))+
    geom_hex(data=readRat[readRat$specAll=="ECOLI",], mapping = aes(x=log10(Bmean), y=log2(ABratio))
             ,bins= binSize*0.6,size=sizePoint,linewidth=0.1,na.rm=T,show.legend=F) +#binwidth=c(0.02,0.04)
    scale_fill_viridis(option = "G",direction = -1,begin=0,end=0.85,alpha = 0.9)+
    scale_x_continuous(name="log10(B)", limits=c(2,8))+ 
    scale_y_continuous(name="log2(A/B)", limits=c(-4, 4))+ 
    theme_minimal()+theme(panel.grid.minor=element_blank())
}
pdf(paste0("hye_prepro-LFQ_v92.pdf"),width=15,height = 8)
ggarrange(plotlist=c(g1,g2,g3),ncol=4,nrow=3)
dev.off()


hycPreProList<-list(diann_hyc_pre,diart_hyc_pre[inter_hyc_pre,],diann_hyc_pro,diart_hyc_pro[inter_hyc_pro,])
g1<-g2<-g3<-list()
for(m in 1:4){
  readRat<-hycPreProList[[m]]
  if(m<3){
    readRat$specAll<-hyc_trans_spec[row.names(readRat)]
    sizePoint=0.1
    binSize=300
    readRat$Amean<-apply(readRat[,c(1:3)],1,mean,na.rm=F)
    readRat$Bmean<-apply(readRat[,c(4:6)],1,mean,na.rm=F)
    readRat$ABratio<-readRat$Amean/readRat$Bmean
    
  }else{
    readRat$specAll<-hyc_prot_spec[row.names(readRat)]
    sizePoint=1
    binSize=120
    readRat$Amean<-apply(readRat[,c(1:3)],1,mean,na.rm=F)
    readRat$Bmean<-apply(readRat[,c(4:6)],1,mean,na.rm=F)
    readRat$ABratio<-readRat$Amean/readRat$Bmean
    
  }

  g1[[m]]<-ggplot(NULL, aes(x, y)) +
    geom_hline(yintercept=1,color=c('#3CB875'))+
    geom_hex(data=readRat[readRat$specAll=="YEAST",]
             , mapping = aes(x=log10(Bmean), y=log2(ABratio))
             ,bins= binSize*0.6,size=sizePoint,linewidth=0.1,na.rm=T,show.legend=F) +
    scale_fill_viridis(option = "D",direction = -1,begin=0,end=0.8,alpha = 0.9)+
    scale_x_continuous(name="log10(B)", limits=c(3,10))+ 
    scale_y_continuous(name="log2(A/B)", limits=c(-4, 4))+ 
    theme_minimal()+theme(panel.grid.minor=element_blank())
  g2[[m]]<-ggplot(NULL, aes(x, y)) +
    geom_hline(yintercept=0,color=c('#DE4B68'))+
    geom_hex(data=readRat[readRat$specAll=="HUMAN",], mapping = aes(x=log10(Bmean), y=log2(ABratio))
             ,bins= binSize,size=sizePoint,linewidth=0.1,na.rm=T,show.legend=F) +
    scale_fill_viridis(option = "A",direction = -1,begin=0,end=0.85,alpha = 0.9)+
    scale_x_continuous(name="log10(B)", limits=c(3,10))+ 
    scale_y_continuous(name="log2(A/B)", limits=c(-4, 4))+ 
    theme_minimal()+theme(panel.grid.minor=element_blank())
  g3[[m]]<-ggplot(NULL, aes(x, y)) +
    geom_hline(yintercept=log2(1/1.3),color=c('#369DA9'))+
    geom_hex(data=readRat[readRat$specAll=="CAEEL",], mapping = aes(x=log10(Bmean), y=log2(ABratio))
             ,bins= binSize*0.6,size=sizePoint,linewidth=0.1,na.rm=T,show.legend=F) +#binwidth=c(0.02,0.04)
    scale_fill_viridis(option = "G",direction = -1,begin=0,end=0.85,alpha = 0.9)+
    scale_x_continuous(name="log10(B)", limits=c(3,10))+ 
    scale_y_continuous(name="log2(A/B)", limits=c(-4, 4))+ 
    theme_minimal()+theme(panel.grid.minor=element_blank())
}
pdf(paste0("hyc_prepro-LFQ_v9.pdf"),width=15,height = 8)
ggarrange(plotlist=c(g1,g2,g3),ncol=4,nrow=3)
dev.off()


#Cross-corr####
inter_hye_pre<-intersect(row.names(diann_hye_pre),row.names(diart_hye_pre))
inter_hyc_pre<-intersect(row.names(diann_hyc_pre),row.names(diart_hyc_pre))
inter_hye_pro<-intersect(row.names(diann_hye_pro),row.names(diart_hye_pro))
inter_hyc_pro<-intersect(row.names(diann_hyc_pro),row.names(diart_hyc_pro))


cor_hye_pre<-cor((unlist(diann_hye_pre[inter_hye_pre,])),
(unlist(diart_hye_pre[inter_hye_pre,])),use='na.or.complete',method='spearman')
cor_hyc_pre<-cor((unlist(diann_hyc_pre[inter_hyc_pre,])),
(unlist(diart_hyc_pre[inter_hyc_pre,])),use='na.or.complete',method='spearman')
cor_hye_pro<-cor((unlist(diann_hye_pro[inter_hye_pro,])),
(unlist(diart_hye_pro[inter_hye_pro,])),use='na.or.complete',method='spearman')
cor_hyc_pro<-cor((unlist(diann_hyc_pro[inter_hyc_pro,])),
(unlist(diart_hyc_pro[inter_hyc_pro,])),use='na.or.complete',method='spearman')

library(RColorBrewer)
colGrad<-c("white","#ebebeb", brewer.pal(3,"Blues"), rev(brewer.pal(3,"Oranges")))

png("Figure3C_XCorr_v10.png",width = 210,height = 55,units = "mm",res = 300)
par(mfrow=c(1,4),mar=c(5,4,3,3))
smoothScatter(log10(unlist(diann_hye_pre[inter_hye_pre,])),
              log10(unlist(diart_hye_pre[inter_hye_pre,]))
              ,colramp = colorRampPalette(colGrad,space = "Lab"),nrpoints = 50
              ,main=paste("r =",round(cor_hye_pre,3)),xlab = "DIA-NN",ylab="DIArt")
smoothScatter(log10(unlist(diann_hye_pro[inter_hye_pro,])),
              log10(unlist(diart_hye_pro[inter_hye_pro,]))
              ,colramp = colorRampPalette(colGrad,space = "Lab"),nrpoints = 50
              ,main=paste("r =",round(cor_hye_pro,3)),xlab = "DIA-NN",ylab="DIArt")
smoothScatter(log10(unlist(diann_hyc_pre[inter_hyc_pre,])),
              log10(unlist(diart_hyc_pre[inter_hyc_pre,]))
              ,colramp = colorRampPalette(colGrad,space = "Lab"),nrpoints = 50
              ,main=paste("r =",round(cor_hyc_pre,3)),xlab = "DIA-NN",ylab="DIArt")
smoothScatter(log10(unlist(diann_hyc_pro[inter_hyc_pro,])),
              log10(unlist(diart_hyc_pro[inter_hyc_pro,]))
              ,colramp = colorRampPalette(colGrad,space = "Lab"),nrpoints = 50
              ,main=paste("r =",round(cor_hyc_pro,3)),xlab = "DIA-NN",ylab="DIArt")
dev.off()



