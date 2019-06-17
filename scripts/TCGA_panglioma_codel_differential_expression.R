setwd("~/Documents/glioma_genefusion/projects/data/CNAPE/")
tcgaexp = read.delim("GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
                     stringsAsFactors = F, row.names = 1, comment.char = "#")
names(tcgaexp) = substr(names(tcgaexp),1,15)

tcgaau =read.delim("LGGGBM.CNV.txt", stringsAsFactors = F, comment.char = "#", row.names = 1)[1:14]
tcgaau$Case = gsub(pattern = "-", replacement = ".", rownames(tcgaau))
rownames(tcgaau)=tcgaau$Case


codellab = read.delim("~/Documents/codel/tcga_panglioma_codel_withRNA_labels.txt", stringsAsFactors = F)

tcgaexp = tcgaexp[,names(tcgaexp) %in% codellab$Case]

tcgaau = tcgaau[tcgaau$Case %in% names(tcgaexp),]


##align the sample names
tcgaau = merge(tcgaau, codellab)
tcgaexp = tcgaexp[,tcgaau$Case]

## annotate the genes with their cytobands
library(dplyr)
tcgagn = data.frame(genename = strsplit(rownames(tcgaexp), "\\|" ) %>% lapply( "[[", 1 ) %>% unlist(), 
                    geneid = strsplit(rownames(tcgaexp), "\\|" ) %>% lapply( "[[", 2 ) %>% unlist(),
                    stringsAsFactors = F)


idx = duplicated(tcgagn$genename)
tcgagn = tcgagn[-idx,]
tcgaexp = tcgaexp[-idx,]
hg19 = read.delim("~/Documents/pancanfusions/hg19.v19.mymap.txt", stringsAsFactors = F)
hg19 = hg19[,c(3,5,6,7)]
hg19 = hg19[hg19$chr != "chrM",]
hg19 = hg19[!duplicated(hg19$genename),]
cm = read.delim("~/Documents/dockerhome/tools/EXCAVATOR2_Package_v1.1.2/data/centromere/CentromerePosition_hg19.txt", stringsAsFactors = F)
for (i in 1:nrow(hg19)){
  si = cm$chromStart[cm$chrom==hg19$chr[i]]
  ei = cm$chromEnd[cm$chrom==hg19$chr[i]]
  if (hg19$end[i] <= si){
    hg19$arm[i] = paste0(hg19$chr[i],"p")
  } else if (hg19$start[i] >= ei){
    hg19$arm[i] = paste0(hg19$chr[i],"q")
  }else{
    hg19$arm[i] = "weired"
  }
}


tcgagn = merge(hg19, tcgagn, by = "genename")
tcgagn = tcgagn[!duplicated(tcgagn$genename),]
tcgagn$chr = factor(tcgagn$chr, levels = c(paste0("chr",1:22),"chrX","chrY"))
tcgagn = tcgagn[order(tcgagn$chr, tcgagn$start),]
tcgaexp = tcgaexp[rownames(tcgaexp) %in% paste(tcgagn$genename,tcgagn$geneid, sep = "|"),]
rownames(tcgaexp) = strsplit(rownames(tcgaexp), "\\|" ) %>% lapply( "[[", 1 ) %>% unlist()

#now compare codel with non codel
stopifnot(identical(tcgaau$Case, names(tcgaexp)))
expdata = as.data.frame(t(tcgaexp))

tcgacodeldge = data.frame(symbol = rownames(tcgaexp),
                          codelexp = 0,
                          ncdlexp = 0,
                          pvalue = 1,
                          stringsAsFactors = F)
for (i in 1:nrow(tcgacodeldge)){
  tcgacodeldge$codelexp[i] = mean(log2(expdata[tcgaau$codel == "codel",i] + 1))
  tcgacodeldge$ncdlexp[i] = mean(log2(expdata[tcgaau$codel != "codel",i] + 1))
  tcgacodeldge$pvalue[i] = t.test(log2(expdata[,i]+1) ~ tcgaau$codel)$p.value
}
tcgacodeldge$arm = tcgagn$arm[match(tcgacodeldge$symbol,tcgagn$genename)]
tcgacodeldge$diff = tcgacodeldge$codelexp - tcgacodeldge$ncdlexp
tcgacodeldge$qvalue = p.adjust(tcgacodeldge$pvalue, method = "bonferroni")
tcgacodeldge$color = ifelse(tcgacodeldge$qvalue >= 0.05 | abs(tcgacodeldge$diff) < 1, rgb(0,0,0,0.5),
                            ifelse(tcgacodeldge$diff < -1 , rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
tcgacodeldge$pch = ifelse(tcgacodeldge$arm =="chr1p",3, ifelse(tcgacodeldge$arm =="chr19q",4,20))
tcgacodeldge$cex = ifelse(tcgacodeldge$arm %in% c("chr1p","chr19q"),0.5,0.15)
tcgacodeldge = tcgacodeldge[order(tcgacodeldge$cex),]
plot(tcgacodeldge$diff, -log10(tcgacodeldge$qvalue), pch = tcgacodeldge$pch, cex =tcgacodeldge$cex , col = tcgacodeldge$color,
     xlab = "Fold Change (log2)", ylab = "Q-value (-log10)")
abline(h = -log10(0.05), lty = 2, lwd = 2)
abline(v = 1, lty = 2, lwd = 2, col = "red")
abline(v = -1, lty = 2, lwd = 2, col = "blue")
idx = which(tcgacodeldge$cex==0.5 & tcgacodeldge$color=="#FF000080")
points(tcgacodeldge$diff[idx],-log10(tcgacodeldge$qvalue[idx]), pch = tcgacodeldge$pch[idx], cex = 0.8, col = 'red')
table(tcgacodeldge$color,tcgacodeldge$cex)
legend(2,170,legend = c("genes on 1p","genes on 19q","others"), pch = c(3,4,20), cex = 0.7, bty='n')

library(ggplot2)
gn = "PTGFR"#show one gene as example
tmpdf = data.frame(exp = expdata[,gn], codel = tcgaau$codel)
ggplot(aes(x = codel, y=log2(exp + 1), fill = codel), data = tmpdf) + geom_violin() + geom_boxplot(width = 0.05) + 
  geom_jitter(width = 0.03, cex = 0.1)+theme_classic() + labs(x = "1p/19q codel status", y = paste0(gn," expression (log2)"))

###
###compare those with and without CDKN2A deep deletion
###
tcgaau$CDKN2A.1029 = ifelse(tcgaau$CDKN2A.1029 == -2,"dd","ndd")
tcgaau$CDKN2A.1029[is.na(tcgaau$CDKN2A.1029)] = "ndd"
tcgacdkn2adge = data.frame(symbol = rownames(tcgaexp),
                          cddexp = 0,
                          ncddexp = 0,
                          pvalue = 1,
                          stringsAsFactors = F)
for (i in 1:nrow(tcgacdkn2adge)){
  tcgacdkn2adge$cddexp[i] = mean(log2(expdata[tcgaau$CDKN2A.1029 == "dd",i] + 1))
  tcgacdkn2adge$ncddexp[i] = mean(log2(expdata[tcgaau$CDKN2A.1029 != "dd",i] + 1))
  tcgacdkn2adge$pvalue[i] = t.test(log2(expdata[,i]+1) ~ tcgaau$CDKN2A.1029)$p.value
}
tcgacdkn2adge$arm = tcgagn$arm[match(tcgacdkn2adge$symbol,tcgagn$genename)]
tcgacdkn2adge$diff = tcgacdkn2adge$cddexp - tcgacdkn2adge$ncddexp
tcgacdkn2adge$qvalue = p.adjust(tcgacdkn2adge$pvalue, method = "bonferroni")
tcgacdkn2adge$color = ifelse(tcgacdkn2adge$qvalue >= 0.05 | abs(tcgacdkn2adge$diff) < 1, rgb(0,0,0,0.5),
                             ifelse(tcgacdkn2adge$diff < -1 , rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
tcgacdkn2adge$pch = ifelse(tcgacdkn2adge$arm %in% c("chr9p","chr9q"),8,20)
tcgacdkn2adge$cex = ifelse(tcgacdkn2adge$arm %in% c("chr9p","chr9q"),0.7,0.5)
tcgacdkn2adge = tcgacdkn2adge[order(tcgacdkn2adge$cex),]
plot(tcgacdkn2adge$diff, -log10(tcgacdkn2adge$qvalue), pch = tcgacdkn2adge$pch, cex =tcgacdkn2adge$cex , col = tcgacdkn2adge$color,
     xlab = "Fold Change (log2)", ylab = "Q-value (-log10)")
abline(h = -log10(0.05), lty = 2, lwd = 2)
abline(v = 1, lty = 2, lwd = 2, col = "red")
abline(v = -1, lty = 2, lwd = 2, col = "blue")
idx = which(tcgacdkn2adge$cex==0.7 & tcgacdkn2adge$color=="#FF000080")

mygenestring = read.delim(paste0("~/Documents/codel/CNAPE/data/CDKN2A_string_interactions.tsv"), stringsAsFactors = F)
mystringnb = unique(c(mygenestring$X.node1, mygenestring$node2))
tcgacdkn2adge$strnb = ifelse(tcgacdkn2adge$symbol %in% mystringnb,1,0)
gndist = function(a,b,c,d){ #a,b: start, end of gene A; c,d: start, end of B
  if (c>b & d>b){ #B is downstram of A
    return(c-b)
  } else if (d<b & c<b){
    return(b-c) #B is upstram of A
  }else{ #overlap
    return(0)
  }
}

tcgagn$disttocdkn2a = 500000000
cs = tcgagn$start[tcgagn$genename == "CDKN2A"]#start of CDKN2A
ce = tcgagn$end[tcgagn$genename == "CDKN2A"]#end of CDKN2A
for (i in 1:nrow(tcgagn)){
  if (tcgagn$chr[i]=="chr9"){
    ts = tcgagn$start[i]
    te = tcgagn$end[i]
    tcgagn$disttocdkn2a[i] = gndist(cs,ce,ts,te)
  }
}
tcgacdkn2adge$dist2cdkn2a = tcgagn$disttocdkn2a[match( tcgacdkn2adge$symbol,tcgagn$genename)]
tcgacdkn2adge$physnb = ifelse(tcgacdkn2adge$dist2cdkn2a <20000000,1,0)#20M 
plot(tcgacdkn2adge$diff, -log10(tcgacdkn2adge$qvalue), pch = tcgacdkn2adge$pch, cex = 0.7*tcgacdkn2adge$physnb+0.1, col = tcgacdkn2adge$color,
     xlab = "Log2(fold change)", ylab = "Q-value (-log10)")
abline(h = -log10(0.05), lty = 2, lwd = 2)
abline(v = 1, lty = 2, lwd = 2, col = "red")
abline(v = -1, lty = 2, lwd = 2, col = "blue")
legend(2.5,60,legend = c("neighbours","others"), pch = c(8,20), cex = 0.7, bty='n')
table(tcgacdkn2adge$color,tcgacdkn2adge$physnb)

library(ggplot2)
gn = "PTGFR"#show one gene as example
tmpdf = data.frame(exp = expdata[,gn], codel = tcgaau$codel)
ggplot(aes(x = codel, y=log2(exp + 1), fill = codel), data = tmpdf) + geom_violin() + geom_boxplot(width = 0.05) + 
  geom_jitter(width = 0.03, cex = 0.1)+theme_classic() + labs(x = "1p/19q codel status", y = paste0(gn," expression (log2)"))


