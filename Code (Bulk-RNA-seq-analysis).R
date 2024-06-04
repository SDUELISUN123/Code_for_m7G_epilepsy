
#>Tittle:
#>Integrative Role of RNA N7-methylguanosine in Epilepsy: 
#>Regulation of Neuronal Oxidative Phosphorylation, 
#>Programmed Death and Immune Microenvironment
#>
#>
#>Author:
#>Qingyuan Sun (孙清源)
#>
#>
#>Institute:
#>Shandong University
#>
#>
#>Date:
#>June 4th, 2024





###### 1 Normalize bulk RNA-seq data #####
library(limma)               #引用包
expFile="geneMatrix.txt"     #表达数据文件
conFile="s1.txt"             #对照组样品文件
treatFile="s2.txt"           #实验组样品文件
setwd("")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
#data=log2(data+1)     #如果表达数值较大，可以把这一行前面的#号去掉

#读取对照组样品文件,提取对照组的表达数据
s1=read.table(conFile, header=F, sep="\t", check.names=F)
sampleName1=as.vector(s1[,1])
conData=data[,sampleName1]

#读取实验组样品文件,提取实验组的表达数据
s2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName2=as.vector(s2[,1])
treatData=data[,sampleName2]

#数据矫正
data=cbind(conData, treatData)
data=normalizeBetweenArrays(data)

#对数据进行矫正，输出矫正后的表达量
conNum=ncol(conData)
treatNum=ncol(treatData)
Type=c(rep("con",conNum),rep("treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)



### 2 ROC #####
#输入文件行名基因，列名样本
bc<-read.table("diffGeneExp.txt",sep='\t',header=TRUE,row.names = 1)
bc<-t(bc)
Type<-gsub("(.*)\\_(.*)", "\\2", rownames(bc))

bc<-cbind("Type"=Type,bc)
bc<-as.data.frame(bc)

pdf("ROC.pdf",width = 6,height = 6)
roc1<- roc(bc$Type, as.numeric(bc[,2]))
plot(roc1,#可以用smooth函数变得光滑
     col=cors[1],
     #legacy.axes=T,#更改y轴格式
     #print.auc=TRUE,#显示auc面积
     #print.thres=TRUE,#添加节点和95ci
     
)
legend<-vector()
for (i in 3:ncol(bc)) {
  roc<- roc(bc$Type, as.numeric(bc[,i]))
  plot(roc,#可以用smooth函数变得光滑
       col=cors[i-1],
       #legacy.axes=T,#更改y轴格式
       #print.auc=TRUE,#显示auc面积
       #print.thres=TRUE,#添加节点和95ci
       add=TRUE
  )
  legend[i-1]<-paste0(colnames(bc)[i],"_AUC=",round(auc(roc),3))
}
legend[1]<-paste0(colnames(bc)[2],"_AUC=",round(auc(roc1),3))
legend("bottomright", legend=legend,
       col=cors,lty=1)
dev.off()



###### 3 Consensus clustering ######
library(ConsensusClusterPlus)      #引用包
expFile="diffGeneExp.txt"          #表达数据文件
workDir=""     #工作目录
setwd(workDir)       #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#删掉对照组样品，只保留实验组样品
group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="treat"]

#聚类
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


#输出分型结果
clusterNum=2        #分几类，根据前面的图形判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("m7Gcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$m7Gcluster))
cluster$m7Gcluster=letter[match(cluster$m7Gcluster, uniqClu)]
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="m7Gcluster.txt", sep="\t", quote=F, col.names=F)



###### 4 PCA #####
library(ggplot2)          
inputFile="input.txt"     
outFile="PCA.pdf"         
setwd("")   

rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
data=rt[,c(2:ncol(rt))]
Type=rt[,1]
var=colnames(rt)[1]

data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
PCA = data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Type)

col=c("blue","red")
if(length(levels(factor(Type)))>2){
  col=rainbow(length(levels(factor(Type))))}

pdf(file=outFile, height=2.5, width=3.5) 
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
  scale_colour_manual(name=var,  values =col)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()



###### 5 Calculation of m7G Score #####
expFile="diffGeneExp.txt"      #表达数据文件
setwd("")     #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#PCA分析
pca=prcomp(data, scale=TRUE)
value=predict(pca)
m7Gscore=value[,1]+value[,2]
m7Gscore=as.data.frame(m7Gscore)
scoreOut=rbind(id=colnames(m7Gscore), m7Gscore)
write.table(scoreOut, file="m7Gscore.txt", sep="\t", quote=F, col.names=F)



###### 6 GSVA ######
#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="normalize.txt"              #表达数据文件
clusterFile="m7Gcluster.txt"            #分型的结果文件
gmtFile="c2.cp.kegg.symbols.gmt"     #基因集文件
setwd("")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="treat",drop=F]

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GSVA分析
ssgseaScore=gsva(data, geneSets, method='gsva')
#对GSVA的打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=ssgseaScore[,nameC1,drop=F]
dataC2=ssgseaScore[,nameC2,drop=F]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
data=cbind(dataC1, dataC2)
Type=c(rep("C1",conNum), rep("C2",treatNum))

#通路差异分析
outTab=data.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=test$p.value
  t=test$statistic
  if(pvalue<0.05){
    Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
    outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
  }
}

#绘制柱状图
termNum=10      #展示通路的数目
outTab=outTab[order(outTab$t),]
outTab=outTab[c(1:termNum,(nrow(outTab)-termNum):nrow(outTab)),]
pdf(file="barplot.pdf", width=10, height=5)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("blue3", "red3"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title="",
              xlab="Term", ylab="t value of GSVA score, m7Gcluster B vs A",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()


###### 7 GSEA ######
#> Conducted by GSEA software (version 4.2.3)
#> Plotted via ggplot2 (version 3.4.4)
#> 
#> Here is the code for visualization


#引用包
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("")      #设置工作目录
files=grep(".tsv", dir(), value=T)      #获取目录下所有tsv结尾的文件
data=lapply(files, read.delim)          #读取每个tsv文件
names(data)=files

dataSet=ldply(data, data.frame)
dataSet$pathway = gsub(".tsv", "", dataSet$.id)     #将文件后缀.tsv删掉

gseaCol=c("#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
  geom_line(size = 1.5) + scale_color_manual(values = gseaCol[1:nrow(dataSet)]) +   
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits = c(min(dataSet$RUNNING.ES - 0.02), max(dataSet$RUNNING.ES + 0.02))) +   
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) + 
  guides(colour = guide_legend(title = NULL)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank())+theme(legend.key.size=unit(0.5,'cm'))
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "m7G cluster A<---------------------->m7G cluster B", y = "", title = "Pathways enriched in m7G cluster A (P<0.05)") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

#将图形可视化，保存在"multipleGSEA.pdf"
pdf(file="multipleGSEA.pdf",    #输出图片的文件
    width=9,                    #设置输出图片高度
    height=4.5
)                 #设置输出图片高度
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea, gGene, nrow=2, heights=c(.8,.25)))
dev.off()



###### 8 Microenvironemnt analysis ######
#>Immune Cells' abundance were calculated by TIMER2, an online tool.
#>
#>Here are code for "ESTIMATE" algorithm

library(estimate)
setwd("")
filterCommonGenes(input.f="GSE116174.txt", output.f="GSE116174.gct", id="EntrezID")
estimateScore(input.ds="GSE116174.gct", output.ds="GSE116174_estimate_score.gct", platform="affymetrix")
estimate_score <- read.table("GSE116174_estimate_score.gct", skip = 2, header = TRUE)
write.csv(estimate_score,"GSE116174_est.csv",row.names = FALSE)



###### 9 Identification of DEGs via limma ######
#引用包
library(limma) 
library(VennDiagram)

logFCfilter=0.5               #logFC的过滤条件
adj.P.Val.Filter=0.05         #矫正后p值的过滤条件
expFile="normalize.txt"       #表达数据文件
cluFile="m7Gcluster.txt"      #分型的结果文件
setwd("F:\\科研2\\11.11 M7G癫痫\\21 clusterdiff")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#读取分型的结果文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#数据整理
sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample, "m7Gcluster"]

#设置比较组
geneList=list()
Type=as.vector(cluster)
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)
allDiffGenes=c()

#对比较组进行差异分析
for(i in 1:ncol(comp)){
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  #print(contrast)
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #输出所有基因的差异情况
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #输出差异结果
  diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  geneList[[contrast]]=row.names(diffSig)
}

#绘制venn图
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#输出交集基因的名称
interGenes=Reduce(intersect,geneList)
write.table(file="interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

#输出交集基因的表达量
interGeneExp=data[interGenes,]
interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
write.table(interGeneExp, file="interGeneExp.txt", sep="\t", quote=F, col.names=F)



###### 10 Machine learning: GLM, RF, and SVM #####
#引用包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)      #设置种子
inputFile="diffGeneExp.txt"      #表达数据文件
geneFile="gene.txt"      #基因列表文件
setwd("")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#读取基因列表文件,提取交集核心基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

#获取样品分组信息
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

#对数据进行分组
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#RF随机森林树模型
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)

#SVM机器学习模型
mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)

#GLM模型
mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)


#定义预测函数
p_fun=function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="Control", 0, 1)

#RF随机森林树模型预测结果
explainer_rf=explain(mod_rf, label = "RF",
                     data = test, y = yTest,
                     predict_function = p_fun,
                     verbose = FALSE)
mp_rf=model_performance(explainer_rf)
#SVM机器学习模型预测结果
explainer_svm=explain(mod_svm, label = "SVM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_svm=model_performance(explainer_svm)
#GLM模型预测结果
explainer_glm=explain(mod_glm, label = "GLM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_glm=model_performance(explainer_glm)

#绘制四种方法的残差反向累计分布图
pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_glm)
print(p1)
dev.off()

#绘制四种方法的残差箱线图
pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_glm, geom = "boxplot")
print(p2)
dev.off()


#绘制ROC曲线
pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_glm, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="orange", add=T)
legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc)),
         paste0('GLM: ',sprintf("%.03f",roc3$auc))),
       col=c("red","blue","orange"), lwd=2, bty = 'n')
dev.off()

#对四种方法进行基因的重要性分析,得到四种方法基因重要性评分
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)

#绘制基因重要性图形
pdf(file="importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()
#输出重要性评分最高的基因
geneNum=9     #设置基因的数目
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)



###### 11 Nomogram model ######
#引用包
library(rms)
library(rmda)

inputFile="diffGeneExp.txt"       #输入文件
setwd("")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

#数据打包
ddist=datadist(rt)
options(datadist="ddist")

#构建模型，绘制列线图
lrmModel=lrm(Type~ METTL1+LSM1+SNUPN+NUDT3+EIF4E3, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.001,0.01,0.1,0.3,0.5,0.7,0.9,0.99,0.999,0.9999),
              lp=F, funlabel="Risk of Epilepsy")
#输出列线图
pdf("Nom.pdf", width=6, height=4)
plot(nomo)
dev.off()

#绘制校准曲线
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5, height=5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=T)
dev.off()

#绘制决策曲线
rt$Type=ifelse(rt$Type=="con", 0, 1)
dc=decision_curve(Type ~ METTL1+LSM1+SNUPN+NUDT3+EIF4E3, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
#输出DCA图形
pdf(file="DCA.pdf", width=5, height=5)
plot_decision_curve(dc,
                    curve.names="m7G genes",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()

#绘制临床影响曲线
pdf(file="clinical_impact.pdf", width=5, height=5)
plot_clinical_impact(dc,
                     confidence.intervals=T,
                     col = c("red", "blue"),
                     legend="bottomleft")
dev.off()



####### Codes for single cell seq analysis were in another file


##>Last version: June 4th, 2024
##>By: Qingyuan Sun
##>Shandong University

