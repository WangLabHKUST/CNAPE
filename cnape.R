#check and install necessary packages

if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('pheatmap')) install.packages('pheatmap'); library('pheatmap')


load("model/tcga_aneuploidy_glmnetmodels.RData")
xx <- coef(glmnetmodels$glmnet_chr1p, s = "lambda.min")
gns = xx$`0`@Dimnames[[1]][-1] #extract the gene list


#check the input matrix
inputCheck <-function(df){
  inputGenes = rownames(df)
  commonGenes = intersect(inputGenes, cyto$symbol)
  try(if(length(commonGenes) < 100) stop("ERROR: failed to detect necessary genes. Please check your input."))
  allNumric = sum(apply(df,2,function(x) !is.numeric(x)))
  try(if(allNumric !=0) stop("ERROR: Non-numeric values for gene expression level detected. Please check your input."))
  m = max(inputGenes, na.rm = T)
  n = min(inputGenes, na.rm = T)
  try(if(n<0) stop("ERROR: Negative values detected. Please check your input file."))
  try(if(m>1000000) stop("ERROR: Come on, why is your RPKM greater than 1M???"))
  if (m < 50  ){
    print ("WARNING: Max value in expression data < 50. We will continue treating the input as log2 transformed.")
  }
  return(0)
}

#data preprocessing
renameNewSample <-function(df){ #change the names of the df, and record the name map in a dataframe
  inputNames = names(df)
  newNames = paste("tmp", 1:length(inputNames))
  names(df) <-newNames
  return(data.frame(inputNames = inputNames, newNames = newNames, stringsAsFactors = F))
}

dataPreprocess <-function(df, modelData){ #process the input data to get z scores
  df.log2 = log2(df + 0.1)
  commonGenes= intersect(rownames(df), rownames(modelData))
  df = df[rownames(df) %in% commonGenes,]
  df = df[order(rownames(df)),]
  modelData = modelData[rownames(modelData) %in% commonGenes,]
  modelData = modelData[order(rownames(modelData)),]

  mergeData = cbind(df, modelData)
  rmean = apply(mergeData, 1, mean)
  mergeData = mergeData[rmean >0,]

  mergeData.qt = normalize.quantiles(mergeData)
  rownames(mergeData.qt) = rownames(mergeData)
  colnames(mergeData.qt) = colnames(mergeData)

  modelData.qt.z = t(scale(t(mergeData[,!startsWith(colnames(mergeData.qt), "tmp")])))
  df.qt.z = t(scale(t(mergeData.qt[,startsWith(colnames(mergeData.qt), "tmp")])))

  mergeData.qt.z = cbind(modelData.qt.z, df.qt.z)

  return(mergeData.qt.z)
}

dataPreparation(df, target){
  idx = cyto[,target]
  genes = cyto$symbol[idx]
  df = df[rownames(df) %in% genes,]
  knowns = t(df[,!startsWith(colnames(df), "tmp")])
  unknowns = t(df[,startsWith(colnames(df), "tmp")])
  labs = rep("", nrow(knowna))
  for (i in 1:length(labs)){
    labs[i] = codel.status[which(codel.status$ID == rownames(knowns)[i]),target]
  }
  return(list(modelData = knowns, predictData = unknowns, modelLabel = labs))
}

##start prediction

for (i in 1:length(glmnetmodels)){
  md = glmnetmodels[[i]]
  m = predict(md, newx = as.matrix(dtexp))[,,1]
  m = exp(m)/rowSums(exp(m))
  lb = as.integer(dimnames(m)[[2]][unlist(apply(m, 1, which.max))])
  
  if (i ==1){
    mout = data.frame(Sample = unlist(lapply(dimnames(m)[[1]], function(x) substr(x,1,15))), cna = lb, stringsAsFactors = F)
    names(mout)[2] = gsub(pattern = "glmnet_","",names(glmnetmodels)[i])
  } else{
    tmp = data.frame(Sample = unlist(lapply(dimnames(m)[[1]], function(x) substr(x,1,15))), cna = lb, stringsAsFactors = F)
    names(tmp)[2] = gsub(pattern = "glmnet_","",names(glmnetmodels)[i])
    mout = merge(mout, tmp)
  }
}

rownames(mout) = mout$Sample
#mout = mout[,-1]

write.table(mout[1:40], file = "example/example.arm_level.cna.txt", quote =F, row.names = F,sep = "\t" )
moutc = mout[,c(1,41:52,26:28,53:57,39:40)]
names(moutc)[2:23] = gsub(pattern = "q",replacement = "",x =names(moutc)[2:23] )
write.table(mout[1:40], file = "example/example.chromosome_level.cna.txt", quote =F, row.names = F, sep = "\t" )

