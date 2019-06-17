load("~/Documents/codel/CNAPE/model/tcga_aneuploidy_glmnetmodels.RData")
#load("~/Documents/codel/tcga_aneuploidy_glmnetmodels.lasso.RData")
hg19cyto = read.delim("~/Documents/codel/CNAPE/scripts/hg19.cytoBand.txt", stringsAsFactors = F)

if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('pheatmap')) install.packages('pheatmap'); library('pheatmap')

par(las = 2)
pred.acc = as.data.frame(lapply(glmnetmodels,FUN = function(x) x$accuracy))
names(pred.acc) = gsub('glmnet_','',names(pred.acc))

df = pred.acc[,1:39]*100
x = melt(df)
arm.acc <- ggplot(x, aes(variable, value, color = variable)) + ggtitle("Arm-level Models")+
  geom_boxplot(show.legend = F) + ylim(50,100) +
  theme_classic() +theme(axis.text.x=element_text(angle = -90,hjust = 0,size = 7, colour = "black")) +
  labs(x = "", y = "Precition accuracy (%)")
#boxplot(pred.acc[,1:39]*100, ylim = c(50,100), ylab = "Prediction Accuracy (%)", cex.axis= 0.7)
df = pred.acc[,c(40:51,25:27,52:56,38:39)]*100
names(df) <- gsub('q','',names(df))
x = melt(df)
chrom.acc <- ggplot(x, aes(variable, value, color = variable)) +ggtitle("Chromosome-level Models")+
  geom_boxplot(show.legend = F) + ylim(50,100) +
  theme_classic() +theme(axis.text.x=element_text(angle = -90, hjust = 0,size = 7, colour = "black")) +
  labs(x = "", y = "Precition accuracy (%)")
#boxplot(df, ylim = c(50,100), ylab = "Prediction Accuracy (%)", cex.axis= 0.7)


st = NULL
for (m in 1:length(glmnetmodels)){
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
  df = df[!is.na(df$arm),]
  armstat = as.data.frame(table(df$arm))
  names(armstat) = c("arm",gsub(pattern = "glmnet_","",names(glmnetmodels)[m]))
  if (m==1){
    st = armstat
  } else{
    st = merge(st,armstat,all.y  =T)
  }

}

rownames(st )= st$arm
st = st[,-1]
st[is.na(st)]=0
rsums = apply(st,1,sum)
# for (i in 1:nrow(st)){
#   st[i,] = st[i,]/rsums[i]
# }
st2plot = st[,1:39] #arms
rownames(st2plot) = paste0("chr",rownames(st2plot))
st2plot = st2plot[names(st2plot),]
st2plot[is.na(st2plot)] =0 
for (i in 1:ncol(st2plot)){
  st2plot[,i] = st2plot[,i]/sum(st2plot[,i])
}
pheatmap::pheatmap(st2plot, cluster_rows = F, cluster_cols = F,border_color = NA,show_colnames = T,cex = 0.7,
                   #color = colorRampPalette(c("black","black","red","yellow","white"))(50))
                   legend_breaks = c(0,.05,0.1,0.15,0.2,0.25,0.3),
                   color = colorRampPalette(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                                                  "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")))(50))

st3plot = st[,c(40:51,25:27,52:56,38:39)] #chromosomes

st3plot2 = st3plot[1:ncol(st3plot),]
rownames(st3plot2) = names(st3plot2)
for(i in c(1:12,16:20)){
  a = paste0("chr",i)
  p = paste0(i,"p")
  q = paste0(i,"q")
  st3plot2[a,] = st3plot[p,] + st3plot[q,]
}

for(i in c(13:15,21:22)){
  a = paste0("chr",i,"q")
  #p = paste0(i,"p")
  q = paste0(i,"q")
  st3plot2[a,] = st3plot[q,]
}

names(st3plot2) <- gsub("q","",names(st3plot2))
rownames(st3plot2) <- gsub("q","",rownames(st3plot2))

for (i in 1:ncol(st3plot2)){
  st3plot2[,i] = st3plot2[,i]/sum(st3plot2[,i])
}
pheatmap::pheatmap(st3plot2, cluster_rows = F, cluster_cols = F,border_color = NA,show_colnames = T,cex = 0.8,
                   #color = colorRampPalette(c("black","black","red","yellow","white"))(50))
                   legend_breaks = c(0,.05,0.1,0.15,0.2,0.25,0.3),
  color = colorRampPalette(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                               "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")))(50))
