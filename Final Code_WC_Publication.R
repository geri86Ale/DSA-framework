### Final Code
### Finalized June 2, 2014

setwd("~/Desktop/Mirjana's Lab/DSA Project/Code Final/")
load("~/Desktop/Mirjana's Lab/DSA Project/Code Final/Raw Data/Server Workspace.RData")

##  Cleaning workspace and loading functions, 
##  packages, DSA functions, and microarray data

rm(list=ls())
load("step1-signal-extraction.rdata")
library(gplots)
library(matrixStats)
source('functions.R')

anno=output[,"geneassignment"]
names(anno)=rownames(rmaData)


## Removing annotations that have NA and/or dvec>=1.2
anno1 = anno
anno1[!is.na(anno1)]=1
anno1[is.na(anno)]=0
dvec=apply(rmaData,1,function(r){max(r)-min(r)})
#crit = anno1 > 0 & 2^dvec >= 1.2
crit = anno1 > 0 
anno = anno[crit]  ## anno is the annotation of 170

data = as.matrix(rmaData)
rmaData = data[crit,]
dvec=apply(rmaData,1,function(r){max(r)-min(r)})

## setting up training set
training=c()
training$NSC = c("Nes", "Sox2", "Gfap", "Thrsp", "Vim",
                 "Fabp7", "NM_009233", "Ascl1", "Hes5", 
                 "Pax6", "NM_010634","Pdgfra")  
                #"Fabp3", "Lfng", "Lpar1", "Nr2e1" fabp5 = NM_010634, sox1 = NM_009233
training$NB = c("Dcx","Pomc","Neurod1","Dpysl3","Eomes",
                "Sox4","Sox11","Prox1", "Neurog2", 
                "NM_009237", "Foxg1","Reln","Tubb4", "Disc1") 
                #Dlx2, Mcm2, Foxg1, sox3 = NM_009237
training$OMEG = c("Calb1","Csf1r","Cd68","Itgam","Pecam1",
                  "Cspg4","S100b","Cd3g","NM_010814", "NM_019467",
                  "Flt1", "Mecp2")#"Vegfa","Vegfb","Vegfc", NM_010814 = Mog, NM_019467 = Aif1

trainingMK=lapply(training,function(r){
  index = unlist(sapply(r, function(i){grep(i,anno)}))
  rownames(rmaData)[index]
})

# remove wrong probe 10422272 (Sox22)
trainingMK$NSC = trainingMK$NSC[-2]
#probeMK$OMEG = probeMK$OMEG[-2]
#probeMK$OMEG = probeMK$OMEG[-2]
#probeMK$ANP = probeMK$ANP[-3]

## Remove wrong probes 10375396 (Sox30)
#trainingMK$NB = trainingMK$NB[-10]

# remove wrong probes 10393926 (Dcxr)
trainingMK$NB = trainingMK$NB[-1]
#probeMK$pomc = probeMK$pomc[-1]

short_anno = match(unlist(trainingMK),rownames(rmaData))
training_rmaData = rmaData[short_anno,]
training_anno = as.data.frame(anno)[short_anno,]

### Marker set 1
markers=c()

markers$NSC = c("Nes","Pax6","Ascl1")
markers$NB = c("Sox3","Pomc","Disc1")
markers$OMEG =c("Calb1","Aif1","Pecam1","Cspg4","S100b")

probeMK=lapply(markers,function(r){
  index = unlist(sapply(r, function(i){grep(i,training_anno)}))
  rownames(training_rmaData)[index]
})

### Probe ID Check
### probeMK$NSC:  Nes (10493114), Pax6 (10474312), Ascl1 (10371578)
### probeMK$NB:  Sox3 (10604837), Pomc (10394240), Disc1 (10576538)
### probeMK$OMEG:  Calb1 (10503416), Aif1 (10450484), Pecam1 (10392221), Cspg4 (10585588), S100b (10364149)

probeMK

indexMK=lapply(probeMK,function(r){
  tmp=rownames(training_rmaData)
  tmp= tmp %in% r
})

###  W ~ 1, weight estimation using marker probes from training
model_w=weight_estimate_1(2^training_rmaData,indexMK, method="Lm")

###  estData estimation -- Training set
training_estData =log2(deconvoltion_general(2^training_rmaData,t(model_w),method="QP_LM"))
rownames(training_estData) = rownames(training_rmaData)
colnames(training_estData) = c("NPC", "IN", "OMEGA")

###  estData estimation -- Test set
test_estData = log2(deconvoltion_general(2^rmaData,t(model_w),method="QP_LM"))
rownames(test_estData) = rownames(rmaData)
colnames(test_estData) = c("NPC", "IN", "OMEGA")

print(sapply(markers$NSC,function(r){test_estData[grep(r,anno),]}))
print(sapply(markers$NB,function(r){test_estData[grep(r,anno),]}))
print(sapply(markers$OMEG,function(r){test_estData[grep(r,anno),]}))
print(test_estData[grep("Lpar1",anno),])

#write.csv(test_estData, "test_estData.csv")

##### stability
estData3=c()

library(multicore)

# for(i in 1:50){
#   print(i)
#   mix_ob = 2^rmaData[,sample(1:ncol(rmaData),5,replace=F)] #tried 6
#   # gene_list = indexMK
#   
#   w2=weight_estimate_1(mix_ob,indexMK, method="Lm")
#   temp = log2(deconvoltion_general(mix_ob,t(w2),method="QP_LM"))
#   estData3=cbind(estData3,temp) 
# }

TN=50
estData3 = mclapply(1:TN,function(r){
  mix_ob = 2^rmaData[,sample(1:ncol(rmaData),8,replace=F)] #tried 6
  # gene_list = indexMK
  train_ob = 2^training_rmaData[,colnames(mix_ob)]
  
  w2=weight_estimate_1(train_ob,indexMK, method="Lm")
  temp = log2(deconvoltion_general(mix_ob,t(w2),method="QP_LM"))
  
},mc.cores=20)

tt=c()
for(i in 1:length(estData3)){
  tt=cbind(tt,estData3[[i]])
}
estData3=tt

print("Multicore DONE!!!!!!!!!!!")

#write.csv(estData3, "estData3.csv")

estData4 = estData3[,c(seq(1,ncol(estData3),by=3),seq(2,ncol(estData3),by=3),seq(3,ncol(estData3),by=3))]

#write.csv(estData4, "estData4.csv")

#estData4 = read.csv("estData4.csv", header = TRUE)
#estData4 = estData4[,-1]
#estData5 = estData4[,2:151]

estData5 = estData4

method = "bonferroni"
pvec=unlist(apply(estData5,1,function(r){t.test(r[1:TN],r[(TN+1):(2*TN)],paired=T)$p.value}))
qvec = p.adjust(pvec,method=method)
dvec=apply(rmaData,1,function(r){max(r)-min(r)})


pvec2=unlist(apply(estData5,1,function(r){t.test(r[1:TN],r[(2*TN+1):(3*TN)],paired=T)$p.value}))
qvec2 = p.adjust(pvec2,method=method)
dvec=apply(rmaData,1,function(r){max(r)-min(r)})


pvec3=unlist(apply(estData5,1,function(r){t.test(r[(TN+1):(2*TN)],r[(2*TN+1):(3*TN)],paired=T)$p.value}))
qvec3 = p.adjust(pvec3,method=method)
dvec=apply(rmaData,1,function(r){max(r)-min(r)})


estData = test_estData

estData = cbind(apply(estData5[,1:TN],1,mean),
                apply(estData5[,(TN+1):(2*TN)],1,mean),
                apply(estData5[,(2*TN+1):(3*TN)],1,mean))

datafile2 = cbind(estData[,1:3], anno, 2^dvec, pvec, qvec, pvec2, qvec2, pvec3, qvec3, estData[,1]-estData[,2], estData[,1]-estData[,3], estData[,2]-estData[,3])
colnames(datafile2) = c("NPC","NB-IN","OMEGA","anno","2^dvec","pvec","qvec","pvec2", "qvec2", "pvec3", "qvec3",
                        "NPC-NB","NPC-OMEGA", "NB-OMEGA")

#write.csv(cbind(test_estData,datafile2), "datafile2.csv")
#write.csv(cbind(test_estData, datafile2), "datafile2-bonferroni.csv")

### Top gene NPC list
anno1 = anno
anno1[!is.na(anno1)]=1
anno1[is.na(anno)]=0
#high = (apply(estData4[,1:50],1,mean))
#low = (apply(estData4[,61:100],1,mean))
crit.npc = 2^dvec>1.4 & qvec < 0.001 & anno1 > 0 & qvec2 < 0.001 & (estData[,1]-estData[,2]) > 6 & (estData[,1]-estData[,3]) > 6  
crit.in = 2^dvec>1.4 & qvec < 0.001 & anno1 > 0 & qvec3 < 0.001 & (estData[,2]-estData[,1]) > 6 & (estData[,2]-estData[,3]) > 6  
crit.omeg = 2^dvec>1.4 & qvec2 < 0.001 & anno1 > 0 & qvec3 < 0.001 & (estData[,3]-estData[,1]) > 6 & (estData[,3]-estData[,2]) > 6  
#crit = 2^dvec>1.2 & qvec < 0.001 & anno1 > 0 & qvec3 < 0.001 & (estData[,2]-estData[,1]) > 6 & (estData[,2]-estData[,3]) > 6  
#crit = 2^dvec>1.2 & qvec < 0.01 & anno1 > 0  & abs(estData[,1]-estData[,2]) > 7 & high-low > 0 #& high > 10 | low > 10 #& (estData[,1]-estData[,3]) > 7



topestData.npc = as.data.frame(estData5[crit.npc,1:150])
topanno.npc = anno[crit.npc]
topgene.npc = cbind(topanno.npc, topestData.npc)

topestData.in = as.data.frame(estData5[crit.in,1:150])
topanno.in = anno[crit.in]
topgene.in = cbind(topanno.in, topestData.in)

topestData.omeg = as.data.frame(estData5[crit.omeg,1:150])
topanno.omeg = anno[crit.omeg]
topgene.omeg = cbind(topanno.omeg, topestData.omeg)


topgene = topgene.npc
dim(topgene)


dim(topgene.npc)




#write.csv(datafile2[crit,], "Topgenes.csv")
#write.csv(datafile2[crit.npc,], "Topgenes 145 NPC.csv")






library("gplots")
library("RColorBrewer")

#pdf('2dPlot-2.pdf')
par(mfrow=c(2,2))
q = test_estData
color = rgb(abs(q[,3]),abs(q[,1]),abs(q[,2]),maxColorValue=18)
pointsval = 0.5
plot(test_estData[,1],test_estData[,2], col = color, cex= pointsval, pch = 16, xlab = "NPC DSA Expression", ylab = "IN DSA Expression")
abline(0,1)
#points(test_estData[crit.npc,1], test_estData[crit.npc,2], col = "red", pch = 20)
plot(test_estData[,1],test_estData[,3], col = color, cex= pointsval, pch = 16,
     xlab = "NPC DSA Expression", ylab = "OMEGA DSA Expression")
abline(0,1)
#points(test_estData[crit.npc,1], test_estData[crit.npc,3], col = "red", pch = 20)
plot(test_estData[,2],test_estData[,3], col = color, cex= pointsval, pch = 16,
     xlab = "IN DSA Expression", ylab = "OMEGA DSA Expression")
abline(0,1)
#points(test_estData[crit.npc,2], test_estData[crit.npc,3], col = "red", pch = 20)
#abline(0,1, col = "red", pch = 10)
#dev.off()






###  NPC Genes Heatmap

#pdf('NPCheatmap.pdf')

x = topgene[,2:151]
color = brewer.pal(11, "Spectral")
#color = redgreen(20)
heatmap.2(as.matrix(x), Colv = F, col=color, scale="row", labRow = T,
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5,)
#dev.off()

#pdf('NPC group heatmap.pdf')

x = test_estData[crit.npc,]
color = brewer.pal(11, "Spectral")
#color = redgreen(20)
heatmap.2(as.matrix(x), Colv = F, col=color, scale="row", labRow = T,
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5,)
#dev.off()



par(mfrow=c(2,2))
q = test_estData
color = rgb(abs(q[,3]),abs(q[,1]),abs(q[,2]),maxColorValue=18)
pointsval = 0.5
plot(test_estData[,1],test_estData[,2], col = color, cex= pointsval, pch = 16, xlab = "NPC DSA Expression", ylab = "IN DSA Expression")
abline(0,1)
points(test_estData[crit.npc,1], test_estData[crit.npc,2], col = "red", pch = 20)
plot(test_estData[,1],test_estData[,3], col = color, cex= pointsval, pch = 16,
     xlab = "NPC DSA Expression", ylab = "OMEGA DSA Expression")
abline(0,1)
points(test_estData[crit.npc,1], test_estData[crit.npc,3], col = "red", pch = 20)
plot(test_estData[,2],test_estData[,3], col = color, cex= pointsval, pch = 16,
     xlab = "IN DSA Expression", ylab = "OMEGA DSA Expression")
abline(0,1)
points(test_estData[crit.npc,2], test_estData[crit.npc,3], col = "red", pch = 20)


x = topgene.npc[,2:151]
color = brewer.pal(11, "Spectral")
#color = redgreen(10)
heatmap.2(as.matrix(x), Colv = FALSE, col=color, scale="row", labRow = F, labCol = F,
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5)

x = topgene.in[,2:151]
color = brewer.pal(11, "Spectral")
#color = redgreen(20)
#color = rainbow(11)
heatmap.2(as.matrix(x), Colv = FALSE, col=color, scale="row", labRow = T,
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5)

x = topgene.omeg[,2:151]
color = brewer.pal(11, "Spectral")
#color = redgreen(20)
heatmap.2(as.matrix(x), Colv = FALSE, col=color, scale="row", labRow = T,
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5)

together = rbind(topgene.npc[,2:151], topgene.in[,2:151], topgene.omeg[,2:151])
x = together
color = brewer.pal(5, "Spectral")
#color = redgreen(20)
heatmap.2(as.matrix(x), Colv = FALSE, col=color, scale="row", labRow = T,
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5)


### heatmap of weights
#pdf("heatmap-weights.pdf")
weights = model_w
rownames(weights)  = c("NPC", "IN", "OMEGA")
color = brewer.pal(11, "Spectral")
#color = redgreen(20)
heatmap.2(t(weights), Colv = T, Rowv = F, col=color, scale="row", labRow = T, dendrogram=c("column"),
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5, cexCol = 0.5)
#dev.off()






x = rbind(topgene.npc, topgene.in, topgene.omeg)
color = brewer.pal(11, "Spectral")
#color = redgreen(20)
heatmap.2(as.matrix(x), Colv = FALSE, col=color, scale="row", labRow = F,
          key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5)




q = test_estData
color = rgb(abs(q[,3]),abs(q[,1]),abs(q[,2]),maxColorValue=18)
plot(test_estData[,1], test_estData[,2], col = color)

color = rgb(abs(q[,3]),abs(q[,1]),abs(q[,2]),maxColorValue=18)
plot(test_estData[,1], test_estData[,3], col = color)

color = rgb(abs(q[,3]),abs(q[,1]),abs(q[,2]),maxColorValue=18)
plot(test_estData[,2], test_estData[,3], col = color)


library(rgl)

#q = read.csv("07022013 estData4 Trial 59.csv", header = TRUE, sep = ",", row.names = 1)

mycol = cbind(abs(q[,1]-q[,2]),abs(q[,2]-q[,3]),abs(q[,1]-q[,3]))
mycol2=  apply(mycol,1,function(r){(r-min(r)/(max(r)-min(r)))})

mycol3=  t(apply(t(mycol2),1,function(r){r[r!=max(r)]=0;r}))

q = estData
plot3d(q, col = rgb(abs(q[,3]-q[,2]),abs(q[,1]-q[,2]),abs(q[,1]-q[,1]),maxColorValue=16),xlab="",ylab="",zlab="OMEGA")

pdf("3dplot.pdf")
q = test_estData
plot3d(q, col = rgb(abs(q[,3]),abs(q[,1]),abs(q[,2]),maxColorValue=18),xlab="",ylab="",zlab="", box = T, axes = T)
points3d(test_estData[crit.npc,1], test_estData[crit.npc,2], test_estData[crit.npc,3], col = "red", pch = 20)
dev.off()

library(scatterplot3d)
col = rgb(abs(q[,3]),abs(q[,1]),abs(q[,2]),maxColorValue=18)
scatterplot3d(q, color = col)





### Gensat
geneSAT=c()
geneSAT$NSC = scan("GENESAT-NSC",what="character")
geneSAT$NB = scan("GENESAT-NB",what="character")
geneSAT$IN = scan("GENESAT-IN.txt", what = "character")
g=c()
g$NSC=sapply(geneSAT$NSC,function(r){grep(r,anno)})
g$NB=sapply(geneSAT$NB,function(r){grep(r,anno)})
g$IN=sapply(geneSAT$IN, function(r){grep(r,anno)})


###### intersect with the list
m1=sum(test_estData[unlist(g$NSC),1]-test_estData[unlist(g$NSC),2] >5) / length(unlist(g$NSC))
m2=sum(test_estData[unlist(g$NB),2]-test_estData[unlist(g$NB),1] >5) / length(unlist(g$NB))

which(test_estData[unlist(g$NSC),1]-test_estData[unlist(g$NSC),2] >5)

intersect_npc = match(rownames(test_estData[crit.npc,]), rownames(test_estData[unlist(g$NSC),]))
topgene.npc[!is.na(intersect_npc),1]

intersect_nb = match(rownames(test_estData[crit.npc,]), rownames(test_estData[unlist(g$NB),]))
topgene.npc[!is.na(intersect_nb),1]

intersect_in = match(rownames(test_estData[crit.npc,]), rownames(test_estData[unlist(g$IN),]))
topgene.npc[!is.na(intersect_in),1]


intersect_npc = match(rownames(test_estData[crit.in,]), rownames(test_estData[unlist(g$NSC),]))
topgene.npc[!is.na(intersect_npc),1]

intersect_nb = match(rownames(test_estData[crit.in,]), rownames(test_estData[unlist(g$NB),]))
topgene.npc[!is.na(intersect_nb),1]

intersect_in = match(rownames(test_estData[crit.in,]), rownames(test_estData[unlist(g$IN),]))
topgene.npc[!is.na(intersect_in),1]





### Group corelation

nes = estData4[grep("Nes",anno),]
groupcor = function(topgene, crit_index, group = c(1,2,3)){
  result=c()
  topgene_index = topgene
  crit_index = crit_index
  group = group
  for(i in 1:dim(topgene_index)[1]){
    tmp = cor(model_w[group,], rmaData[crit_index,][i,])
    #tmp = cor(rmaData[grep("Nes",anno),], rmaData[crit.npc,][i,])
    #tmp = cor(nes,estData4[crit.npc,][i,] )
    result = c(result, tmp)
  }
  return(result)
}

groupcor.nsc = groupcor(topgene.npc, crit.npc, 1)
groupcor.nb = groupcor(topgene.in, crit.in, 1)
groupcor.omeg = groupcor(topgene.omeg, crit.omeg, 1)

group_cor = cbind(topanno.npc, as.numeric(result))
#pdf("histogramnpc.pdf")
color = brewer.pal(11, "Spectral")
hist(as.numeric(group_cor[,2]), col = color, breaks = 12, 
     xlab = "Correlation to NPC Group", ylab = "Number of Probes", main = "")
#dev.off()

# write.csv(group_cor, "group_cor_nes1.csv")


cor_matrix = cor(t(estData4[crit.npc,]), t(estData4[crit.npc,]))
rownames(cor_matrix) = rownames(estData4[crit.npc,])
colnames(cor_matrix) = rownames(estData4[crit.npc,])

#write.csv(cor_matrix, "cor_matrix.csv")

# for(i in 1:50){
#   print(i)
#   mix_ob = 2^rmaData[,sample(1:ncol(rmaData),5,replace=F)] #tried 6
#   # gene_list = indexMK
#   
#   w2=weight_estimate_1(mix_ob,indexMK, method="Lm")
#   temp = log2(deconvoltion_general(mix_ob,t(w2),method="QP_LM"))
#   estData3=cbind(estData3,temp) 
# }
