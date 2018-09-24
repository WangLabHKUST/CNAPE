require(argparse)
require(preprocessCore)
require(caret)
require(ROSE)

#load the workspace
#in the image: tcga expression data (may also include other data), cyto-gene, codel.status
load("cnape.model.RData")

#check the input matrix
inputCheck <-function(df){
  inputGenes = rownames(df)
  commonGenes = intersect(inputGenes, cyto$symbol)
  try(if(length(commonGenes) < 100) stop("ERROR: failed to detect necessary genes. Please check the documentation and the example input file."))
  allNumric = sum(apply(df,2,function(x) !is.numeric(x)))
  try(if(allNumric !=0) stop("ERROR: Non-numeric values for gene expression level detected. Please check your input."))
  m = max(inputGenes, na.rm = T)
  n = min(inputGenes, na.rm = T)
  try(if(n<0) stop("ERROR: Negative values detected. Please check your input file."))
  try(if(m>1000000) stop("ERROR: Come on, why is your RPKM greater than 1M???"))
  if (m < 50  ){
    print ("WARNING: Max value in expression data < 50. We will continue but has it been log2 transformed?")
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

plsModeling <-function(df){ #df should be made by cbind(modelData,labs)
  idx = createDataPartition(y= df$labs, p = 0.75, list =F)
  train = df[idx,]
  test = df[-idx,]
  ctrl <- trainControl(method = "repeatedcv",repeats =5,summaryFunction = twoClassSummary, classProbs = TRUE)
  pls.grid = expand.grid(ncomp = c(1,3,5,10,20))
  plsFit <- train( labs ~ ., data = df, method = "pls",trControl = ctrl, metric= "ROC",tuneGrid = pls.grid)
  return(plsFit)
}

plsPrediction <-function(plsmod, ndata){
  plsPred <- predict(plsFit, newdata = test[,-ncol(test)], type = "prob")
  plsPred[,3] = ifelse(plsPred[,1] > 0.5, "yes", "no")
  return(plsPred)
}
