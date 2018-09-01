load("~/Documents/codel/tcga_aneuploidy_glmnetmodels.RData")
hg19cyto = read.delim("~/Documents/codel/hg19.cytoBand.txt", stringsAsFactors = F)

library(glmnet)
st = NULL
for (m in 2:length(glmnetmodels)){
  md = glmnetmodels[[m]]
  xx <- coef(md, s = "lambda.min") 
  
  cytoInfo <-function(df0){
    df0$genename = unlist(lapply(df0$name,FUN = function(x) return(strsplit(x,"\\|")[[1]][1])))
    df0$entrez = unlist(lapply(df0$name,FUN = function(x) return(strsplit(x,"\\|")[[1]][2])))
    df0 = merge(df0, hg19cyto,by = "entrez", all.x = T)
    df0$arm = substr(df0$cytoband,1,2)
    for (i in 1:nrow(df0)){
      if (!is.na(df0$arm[i])){
        if (!endsWith(df0$arm[i],"p") & !endsWith(df0$arm[i],"q")){
          df0$arm[i] = substr(df0$cytoband[i],1,3)
        }
      }
    }
    return(df0)
  }
  
  df0 <- data.frame(name = xx$`0`@Dimnames[[1]][xx$`0`@i + 1], coefficient = xx$`0`@x, stringsAsFactors = F)
  df1 <- data.frame(name = xx$`1`@Dimnames[[1]][xx$`1`@i + 1], coefficient = xx$`1`@x, stringsAsFactors = F)
  dfm1 <- data.frame(name = xx$`-1`@Dimnames[[1]][xx$`-1`@i + 1], coefficient = xx$`-1`@x, stringsAsFactors = F)
  
  df0 = cytoInfo(df0)
  df1 = cytoInfo(df1)
  dfm1 = cytoInfo(dfm1)
  df = rbind(df0,df1,dfm1)
  df = df[!duplicated(df$entrez),]
  armstat = as.data.frame(table(df$arm))
  names(armstat) = c("arm",gsub(pattern = "glmnet_","",names(glmnetmodels)[m]))
  st = merge(st,armstat,all.y  =T)
}

rownames(st )= st$arm
st = st[,-1]

st2plot = st[,1:39]
rownames(st2plot) = paste0("chr",rownames(st2plot))
st2plot = st2plot[names(st2plot),]
pheatmap::pheatmap(st2plot, cluster_rows = F, cluster_cols = F,border_color = NA,scale = "row",show_colnames = T,cex = 0.8,
                   color = colorRampPalette(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                                                  "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")))(50))

