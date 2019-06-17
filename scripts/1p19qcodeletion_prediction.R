tcga = read.delim("~/Documents/glioma_genefusion/projects/code/codel_from_RNA.txt", stringsAsFactors = F)
cyto = read.delim("~/Documents/codel/CNAPE/data/Hs.egSYMBOL.gene.cytoband.txt", stringsAsFactors = F)


###
library(caret)
library(glmnet)
library(pROC)
dtx = tcga[,1:20501]
dty = ifelse(tcga$cdl=="codel",1,0)
set.seed(111)
idx = createDataPartition(dty,p = 0.8,list = F)[,1]

xtrn = as.matrix(dtx[idx,])
ytrn = dty[idx]

xtst = as.matrix(dtx[-idx,])
ytst = dty[-idx]

md = cv.glmnet(x= xtrn ,y=ytrn,family="binomial",nfolds=3,alpha = 0.1)
confusionMatrix(as.factor(predict(md,newx = xtst, type = "class")[,1]),as.factor( ytst),positive = '1')

#apply to all samples
confusionMatrix(as.factor(predict(md,newx = as.matrix(dtx), type = "class")[,1]),as.factor( dty),positive = '1')
auroc = roc(as.factor(dty),predict(md,newx = as.matrix(dtx))[,1])
plot(auroc)