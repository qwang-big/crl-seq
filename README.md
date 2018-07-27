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
ggplot(df, aes(Var1, value)) + geom_boxplot(aes(fill=type),outlier.size = NA) + labs(x="", y="Spearman r") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
```
Correlations between histone marks differences and gene expression differences of two groups (dExp)

```{r}
ggplot(df, aes(Var1, value)) + geom_boxplot(aes(fill=Var1)) + labs(x="", y="Spearman r") + 
    theme(legend.position="none",axis.text.x = element_text(angle = 30, hjust = 1))
```
Correlations between histone marks and time-coursed gene expression

```{r}
x=x[,c(2:5,8,9,1)]
x=melt(as.matrix(x))
svglite("NE-Exp.svg",5,5)
ggplot(x, aes(Var1, value, group=1))+geom_line()+facet_wrap(~Var2, scales="free")+geom_vline(xintercept=2, linetype="dashed", col="red")+labs(x="",y="Spearman r")+theme(axis.text.x = element_text(angle = 30, hjust = 1))
```
Correlations between dPCs and gene expression

```{r}
svglite("PC-Exp.svg",3,3)
ggplot(aes(y = value, x = variable, fill = variable), data = df) + geom_boxplot() + labs(x="", y="Spearman r") + theme(legend.position="none")
```

### Correlations between epigenetic marks

Correlations between epigenetic marks in each group

```{r}
svglite('cormarks.svg',10,2)
ggplot(m, aes(cell, Var1, fill = value)) + geom_tile(colour = "white") + 
  facet_grid(.~Var2) + scale_fill_gradient2(low="blue", high="red") +
  theme(axis.text.x = element_text(angle = 90, size=7, vjust = 0.5))+labs(x="",y="")
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

Relationships of histone mark signals and PC1

```{r}
svglite('d.svg',5,4)
ggplot(rint(df2), aes(x,y))+geom_point(aes(col=PC1),size=.1) + scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", space ="Lab" )+ scale_x_continuous(breaks=seq(0,40,20))+ scale_y_continuous(breaks=seq(0,40,20)) +facet_grid(type~mark)+labs(x="",y="")
```

### Weighting promoter-enhancer interactions

Probability density functions of promoter-enhancer interactions
```{r}
m=read.table('CD34.pd')
m=m[m[,1]>=100000,]
m=m[m[,1]<=2000000,]
m=m[m[,2]>=10,]
x=m[,1]
y=m[,2]
model=lm(log(y) ~ x)
m[,2]=log(m[,2])
svglite('CD34.svg',3,3)
ggplot(m,aes(V1, V2))+geom_point(size=1,col='#3288bd', alpha=0.2)+geom_abline(intercept=model$coefficients[1], slope=model$coefficients[2], color="red")+labs(title="CD34",x="distance (bp)",y="count (log)")
```

Determining distance decay coefficients
```{r}
dpos=function(x) -round(exp(x))
ggplot(df, aes(LogCo, AUC, col=type)) + geom_smooth(method="loess", se=F) + theme(legend.position=c(0.15,0.25))  + scale_x_continuous(label=dpos) + xlab("Decay coefficients")
```

### Interpreting ranks

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

Empirical cumulative distribution by ranking all dPCs

```{r}
df=melt(data.frame(lapply(getPCPromId(cll),function(d) sort(match(x,d)))))
colnames(df)[1]="PC"
df[,1]=as.character(df[,1])
df1=df[grep("PC",df[,1]),]
df2=df[grep("Prom",df[,1]),]
df2[,1]=str_replace(df2[,1],"Prom","PC")
svglite('auc.svg',3,3)
ggplot()+stat_ecdf(data=df1,aes(value, colour=PC),pad = TRUE)+ 
stat_ecdf(data=df2,aes(value, colour=PC),pad = FALSE, linetype="dotted")+
	labs(x="Rank", y="ECDF") + theme(legend.position=c(0.8,0.2))
```

AUCs of randomized promoter-enhancer interaction network using PC1

```{r}
dname=function(i) nm[i]
ggplot(df, aes(type))+geom_ribbon(aes(ymin = q1, ymax = q3), fill = "grey70") + geom_line(aes(y = enh), color="red")+ geom_line(aes(y = prom), color="red", linetype="dashed")+ scale_x_continuous(breaks=1:7,label=dname) + labs(x="",y="AUC")
```

AUCs of all PCs from PromOnly and PromEnh

```{r}
ggplot(df, aes(type))+geom_ribbon(aes(ymin = minPC, ymax = maxPC), fill = "lightblue", alpha=0.6) + geom_line(aes(y = PC1), color="red")+
geom_ribbon(aes(ymin = minProm, ymax = maxProm), fill = "grey70", alpha=0.5) + geom_line(aes(y = Prom1), color="red", linetype="dashed")+ scale_x_continuous(breaks=1:7,label=dname) + labs(x="",y="AUC")
```

AUCs of different PC transformation functions

```{r}
m=data.frame(id=rep(1:1000,6),value=unlist(fun(1000)),fun=rep(c("ReLU","Sigmoid","Logit","Inverse\nexponential", "Exponential", "Identity"),each=1000))
ggplot(m, aes(x=id, y=value, colour=fun)) + geom_line()+theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position = "none")+labs(x="dPCs",y="Weights")
levels(df$Functions) = levels(m$fun)
svglite("funAUC.svg",4,3)
ggplot(df, aes(x=Var2, y=value, colour=Functions)) + geom_point(size=2, alpha=0.6) + labs(x="", y="AUC")
```

Comparing positions of well known oncogenes with respect to PromEnh and PromOnly ranks

```{r}
x=read.table('~/rnk.xls',header=F)
colnames(x)=c('gene','pe','po','cell')
x[,2]=17888-x[,2]
x[,3]=17888-x[,3]
svglite("rnk.svg",10,3)
ggplot(x, aes(x=pe, y=po)) +geom_point(aes(col=gene))+facet_grid(.~cell) +geom_abline(slope=1,linetype="dashed")+ xlim(c(0, 17888))+ylim(c(0, 17888))+labs(x="PromEnh", y="PromOnly")+geom_text(aes(label=gene), size=2.5,vjust=2)+theme(axis.ticks=element_blank(),legend.position="none")
```

Wilcoxon-Mann-Whitney test of OG/TSG/HKG ranks comparing to uniformly distributed ranks

```{r}
svglite("oncoRank.svg",3,3)
ggplot(df, aes(Rank, type, fill = Pvalue)) + geom_tile(colour = "white") + 
  facet_grid(Class~name,scales='free') + scale_fill_gradient2(low="white", high="red", midpoint = 0.05) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+labs(x="",y="")
```
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

