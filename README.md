### Set libraries and variables

```{r}
library(ggplot2)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
nm = c("CLL", "CRC", "GLM", "MSC", "NPC", "PTC", "TSC", "MES")
```

### Combine p-values from multiple tests

Use data in Pvalues:

```{r}
svglite("p.svg",3,3)
print(ggplot(pv,aes(value, colour=Mark)) +  stat_ecdf(pad = FALSE) + 
    labs(x="P-value", y="ECDF") + theme(legend.position=c(0.75,0.25)))
```

### Correlations between epigenetic modifications and gene expression

Correlations between histone marks and gene expression (Exp)

```{r}
```
Correlations between histone marks differences and gene expression differences of two groups (dExp)

```{r}
```
Correlations between histone marks and time-coursed gene expression

```{r}
svglite('NE-Exp.svg',9,9)
par(mfrow = c(3,3), oma = c(2, 2, 0, 0), mar = c(2, 2, 2, 2), mgp = c(3, 1, 0), xpd = NA)
lapply(c(2:5,8,9,1),function(d){
plot(res[,d],type='l',xaxt = 'n', xlab=colnames(res)[d], ylab='Spearman r',cex.lab=1.5,cex.axis=1.5)
abline(v=2,col='lightgrey',lty=2)
axis(1, at=1:5, labels=rownames(res), cex.axis=1.5)
})
```
Correlations between dPCs and gene expression

```{r}
svglite("PC-Exp.svg",3,3)
ggplot(aes(y = value, x = variable, fill = variable), data = df) + geom_boxplot() + labs(x="", y="Spearman r") + theme(legend.position="none")
```

### Correlations between epigenetic marks

Correlations between epigenetic marks in each group

```{r}
```
### Interpreting dPCs

Variances explained by dPCs

```{r}
pcvar=cbind(melt(pvar),paste0("PC",1:4))
colnames(pcvar)=c("Type","value","PC")
svglite(paste0('PCvar.svg'),3,3)
print(ggplot(pcvar, aes(x=PC, y=value, colour=Type)) + geom_point(size=2) + labs(x="", y="Variance") + scale_y_continuous(labels = scales::percent) + theme(legend.position=c(0.85,0.7)))
```

Loadings of matrix D (observed differences)

```{r}
plotD = function (Dobs, labels, f='d') {
pc=prcomp(Dobs)$rotation
df=do.call("rbind",lapply(1:ncol(Dobs),function(i) data.frame(Loadings=pc[,i], Mark = labels, PC=i)))
svglite(paste0(f,'.svg'),3,3)
print(ggplot(df, aes(x=PC, y=Loadings, colour=Mark)) + geom_point(size=2))
dev.off()
}
```

Relationships of histone mark signals and dPC1

```{r}
```

### Interpreting ranks


### Genomic view of interested loci

Differential epigenetic modified enhancers

```{r}
ggplot(gr, aes(Position, value, fill=Mark)) + geom_area(stat="identity", position="identity") + facet_grid(Track~.) + scale_fill_manual(values=c("#e41a1c","#377eb8","#4daf4a")) + scale_x_continuous(label=gpos) + scale_y_continuous(breaks=uq) + geom_vline(xintercept = c(172,372,1848,1850)) + xlab("chr8") + theme(panel.background=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
```

Heatmap view of oncogenes, tumor suppressor genes, and housekeeping genes

```{r}
nm=c("H3K27ac","H3K4me1","H3K4me3","H3K9me3","H3K27me3","H3K36me3")
lbl=c(rep("OG",82),rep("TSG",63),rep("HKG",11))
map=lapply(1:6,function(i)
EnrichedHeatmap(normalizeToMatrix(gr[[i]], tss, value_column = "value", 
     extend = 5000, mean_mode = "w0", w = 100), col = c("white", col[i]), name = nm[i], split = lbl, 
     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4)),height=unit(2, "cm")), row_title_rot = 0,cluster_rows=FALSE))
map[[1]]+map[[2]]+map[[3]]+map[[4]]+map[[5]]+map[[6]]
```

Empirical cumulative distribution by ranking dPC1

```{r}
lapply(nm,function(f){
svglite(paste0(f,'.svg'),3,3)
x=read.csv(paste0(f,'rank.csv'))
y=read.csv(paste0(f,'marker.csv'))
df=data.frame(value=sort(match(y[,1],x[,1])),type='PromEnh')
df=rbind(df,data.frame(value=sort(match(y[,1],x[,2])),type='PromOnly'))
print(ggplot(df,aes(value, colour=type))+stat_ecdf(pad = FALSE)+ 
	labs(x="Rank", y="ECDF") + theme(legend.position=c(0.8,0.2)))
dev.off()
})
```
