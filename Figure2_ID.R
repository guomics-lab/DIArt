diann_hye_preIn<-read.delim("diann_report_230829/hye/report.pr_matrix.tsv",header = T,sep='\t')
diann_hye_pre<-diann_hye_preIn[,11:16]
row.names(diann_hye_pre)<-diann_hye_preIn$Precursor.Id
write.csv(diann_hye_pre,"diann_hye_pre.csv")

diann_hye_proIn<-read.delim("diann_report_230829/hye/report.pg_matrix.tsv",header = T,sep='\t')
diann_hye_pro<-diann_hye_proIn[,6:11]
row.names(diann_hye_pro)<-diann_hye_proIn$Protein.Ids
diann_hye_pro<-diann_hye_pro[!grepl(";",row.names(diann_hye_pro)),]
diann_hye_pro<-diann_hye_pro[!grepl("-",row.names(diann_hye_pro)),]
row.names(diann_hye_pro)<-sapply(strsplit(row.names(diann_hye_pro),"\\|"),"[[",2)
write.csv(diann_hye_pro,"diann_hye_pro.csv")

diann_hyc_preIn<-read.delim("diann_report_230829/orbitrap/report.pr_matrix.tsv",header = T,sep='\t')
diann_hyc_pre<-diann_hyc_preIn[,11:16]
row.names(diann_hyc_pre)<-diann_hyc_preIn$Precursor.Id
write.csv(diann_hyc_pre,"diann_hyc_pre.csv")

diann_hyc_proIn<-read.delim("diann_report_230829/orbitrap/report.pg_matrix.tsv",header = T,sep='\t')
diann_hyc_pro<-diann_hyc_proIn[,6:11]
row.names(diann_hyc_pro)<-diann_hyc_proIn$Protein.Ids
diann_hyc_pro<-diann_hyc_pro[!grepl(";",row.names(diann_hyc_pro)),]
diann_hyc_pro<-diann_hyc_pro[!grepl("-",row.names(diann_hyc_pro)),]
write.csv(diann_hyc_pro,"diann_hyc_pro.csv")

#diart####
diart_hye_preIn<-read.delim("diart_report_20230915/hye_precursor_quant.txt",sep='\t',row.names = 1)
diart_hye_pre<-diart_hye_preIn[,1:6]
diart_hye_pre<-diart_hye_pre[apply(is.na(diart_hye_pre),1,sum)!=6,]
row.names(diart_hye_pre)<-paste0(sapply(sapply(row.names(diart_hye_pre),strsplit,"_"),"[[",2),
                                 sapply(sapply(row.names(diart_hye_pre),strsplit,"_"),"[[",3))
write.csv(diart_hye_pre,"diart_hye_pre.csv")

diart_hye_proIn<-read.delim("diart_report_20230915/hye_protein_quant.txt",sep='\t'
                            ,row.names = 1)
# diart_hye_pro<-diart_hye_proIn[,1:6]
# diart_hye_pro<-diart_hye_pro[apply(is.na(diart_hye_pro),1,sum)!=6,]
hye_protID<-unique(hye_spec$ProteinId)
hye_uniID<-c()
for(i in 1:length(hye_protID)){
  hye_uniID[i]<-unlist(strsplit(hye_protID[i],"\\|"))[1]
  names(hye_uniID)[i]<-unlist(strsplit(hye_protID[i],"\\|"))[2]
}
hye_uniIDmap<-hye_uniID[row.names(diart_hye_proIn)]
diart_hye_pro<-diart_hye_proIn[which(!duplicated(hye_uniIDmap)&!is.na(hye_uniIDmap)),1:6]
row.names(diart_hye_pro)<-hye_uniIDmap[!duplicated(hye_uniIDmap)&!is.na(hye_uniIDmap)]
# sapply(strsplit(row.names(diart_hye_pro),"\\|"),"[[",1)
write.csv(diart_hye_pro,"diart_hye_pro.csv")

diart_hyc_preIn<-read.delim("diart_report_20230915/orbitrap_precursor_quant.txt",sep='\t',header = T
                          ,row.names = 1)
diart_hyc_pre<-diart_hyc_preIn[,1:6]
row.names(diart_hyc_pre)<-paste0(sapply(sapply(row.names(diart_hyc_preIn),strsplit,"_"),"[[",2),
       sapply(sapply(row.names(diart_hyc_preIn),strsplit,"_"),"[[",3))
write.csv(diart_hyc_pre,"diart_hyc_pre.csv")

diart_hyc_proIn<-read.table("diart_report_20230915/orbitrap_protein_quant.txt",sep='\t',header = T
                            ,row.names = 1)
diart_hyc_proIn2<-diart_hyc_proIn[!grepl(";",row.names(diart_hyc_proIn)),]
diart_hyc_proIn3<-diart_hyc_proIn2[!grepl(";",diart_hyc_proIn2$ProteinGroups),]
diart_hyc_proIn4<-diart_hyc_proIn3[diart_hyc_proIn3$ProteinGroups!="",]
diart_hyc_pro<-diart_hyc_proIn4[,2:7]
row.names(diart_hyc_pro)<-diart_hyc_proIn4$ProteinGroups
write.csv(diart_hyc_pro,"diart_hyc_pro.csv")

diart_hyc_pre_0.01_1In<-read.delim("diart_report_230829/orbitrap/1%_fdr/MP-LFC-MS1var-OT-S1-120kMS1_MHRM_R01_precursor.tsv",stringsAsFactors = F,sep='\t',header = T)
diart_hyc_pre_0.01_1<-diart_hyc_pre_0.01_1In[,c("transition_group_id","q_value")]
for(i in 1:nrow(diart_hyc_pre_0.01_1)){
  diart_hyc_pre_0.01_1$transition_group_id[i]<-paste0(unlist(strsplit(diart_hyc_pre_0.01_1$transition_group_id[i],"_"))[2:3],collapse = "")
}
write.csv(diart_hyc_pre_0.01_1,"diart_hyc_pre_0.01_1.csv")

diann_hye_pre<-read.csv("diann_hye_pre.csv",stringsAsFactors = F,row.names = 1)
diann_hye_pro<-read.csv("diann_hye_pro.csv",stringsAsFactors = F,row.names = 1)
diann_hyc_pre<-read.csv("diann_hyc_pre.csv",stringsAsFactors = F,row.names = 1)
diann_hyc_pro<-read.csv("diann_hyc_pro.csv",stringsAsFactors = F,row.names = 1)


#Barplots####
colToPlot<-c("#192e57","#f49d17")
pdf("Figure2A_barplots_v6.pdf",width=9,height=4.5)
par(mfrow=c(1,4),mar=c(3,4,3,3))
barplot(c(nrow(diann_hye_pre),nrow(diart_hye_pre)),col=colToPlot,main="TOF-HYE-Precursor",las=2,space=0,ylim=c(0,50000))
text(c(0.5,1.5),c(nrow(diann_hye_pre),nrow(diart_hye_pre))+2000,c(nrow(diann_hye_pre),nrow(diart_hye_pre)))
barplot(c(nrow(diann_hye_pro),nrow(diart_hye_pro)),col=colToPlot,main="TOF-HYE-Protein",las=2,space=0,ylim=c(0,7000))
text(c(0.5,1.5),c(nrow(diann_hye_pro),nrow(diart_hye_pro))+200,c(nrow((diann_hye_pro)),nrow(diart_hye_pro)))
barplot(c(nrow(diann_hyc_pre),nrow(diart_hyc_pre)),col=colToPlot,main="Orbi-HYC-Precursor",las=2,space=0,ylim=c(0,250000))
text(c(0.5,1.5),c(nrow(diann_hyc_pre),nrow(diart_hyc_pre))+10000,c(nrow(diann_hyc_pre),nrow(diart_hyc_pre)))
barplot(c(nrow(diann_hyc_pro),nrow(diart_hyc_pro)),col=colToPlot,main="Orbi-HYC-Protein",las=2,space=0,ylim=c(0,12000))
text(c(0.5,1.5),c(nrow(diann_hyc_pro),nrow(diart_hyc_pro))+350,c(nrow(diann_hyc_pro),nrow(diart_hyc_pro)))
dev.off()

#Vennplots####
library(VennDiagram)
pdf("Figure2B_vennplots_v6.pdf",width=3,height = 3)
grid.newpage()
draw.pairwise.venn(nrow(diann_hye_pre),nrow(diart_hye_pre)
                                ,length(intersect(row.names(diann_hye_pre),row.names(diart_hye_pre)))
                       ,col = colToPlot,c("DIA-NN", "DIArt"),fill=c(0,0))
grid.newpage()
draw.pairwise.venn(nrow(diann_hye_pro),nrow(diart_hye_pro)
                   ,length(intersect(row.names(diann_hye_pro),row.names(diart_hye_pro)))
                   ,col = colToPlot,c("DIA-NN", "DIArt"),fill=c(0,0))
grid.newpage()
draw.pairwise.venn(nrow(diann_hyc_pre),nrow(diart_hyc_pre)
                   ,length(intersect(row.names(diann_hyc_pre),row.names(diart_hyc_pre)))
                   ,col = colToPlot,c("DIA-NN", "DIArt"),fill=c(0,0))
grid.newpage()
draw.pairwise.venn(nrow(diann_hyc_pro),nrow(diart_hyc_pro)
                   ,length(intersect(row.names(diann_hyc_pro),row.names(diart_hyc_pro)))
                   ,col = colToPlot,c("DIA-NN", "DIArt"),fill=c(0,0))
dev.off()


#Lineplots####
pdf("Figure2C_20230830.pdf",width=10,height = 5)
par(mfrow=c(1,2))
plot(names(table(diart_hye_pre_0.01_1$q_value)),cumsum(table(diart_hye_pre_0.01_1$q_value)),ylim=c(0,45000),las=2,type = 'l'
     ,col=colToPlot[2],lwd=2,ylab="",xlab="q value",main="TOF-HYE-Precursor")
lines(names(table(diann_hye_pre_0.01_1$Q.Value)),cumsum(table(diann_hye_pre_0.01_1$Q.Value)),col=colToPlot[1],lwd=2)

plot(names(table(diart_hyc_pre_0.01_1$q_value)),cumsum(table(diart_hyc_pre_0.01_1$q_value)),ylim=c(0,150000),las=2,type = 'l'
     ,col=colToPlot[2],lwd=2,ylab="",xlab="q value",main="Orbi-HYC-Precursor")
lines(names(table(diann_hyc_pre_0.01_1$Q.Value)),cumsum(table(diann_hyc_pre_0.01_1$Q.Value)),col=colToPlot[1],lwd=2)
dev.off()

