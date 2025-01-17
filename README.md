
# Plot figures for the publication
### Setup theme for publication

```{r}
tsize <- 9; tsize2 <- 8; pchsize <- 1; lkey <- 3.5 # lsize <- 0.3
theme_set(theme_bw())
theme_update(panel.grid.minor = element_line(colour = NA),
panel.grid.major = element_line(colour = NA))
```
### Set libraries and variables

```{r}
library(ggplot2)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
bocx=function(x,lam=2/11) (x^lam - 1)/lam
fisherp=function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
data=list(all, aml, cll, crc, glm, mcl, mcll, mes, mm, msc, npc, ptc, tsc)
nm = c("ALL", "AML", "CLL", "CRC", "LGG", "MCL", "mCLL", "MES", "MM", "MSC", "NPC", "PTC", "TSC")
cl = c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Cancer", "Development", "Cancer", "Development", "Development", "Cancer", "Development")
```
### PC composition
```{r}
lbl = list("NPC"=c('DNAMeth','H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9ac','H3K9me3'),
"MSC"=c('DNAMeth','H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9ac','H3K9me3'),
"MES"=c('DNAMeth','H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9ac','H3K9me3'),
"TSC"=c('DNAMeth','H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9ac','H3K9me3'),
"CLL"=c('H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input','DNAMeth'),
"PTC"=c('H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input','DNAMeth'),
"CRC"=c('H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input','DNAMeth'),
"AML"=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","DNAMeth"),
"ALL"=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","DNAMeth"),
"mCLL"=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","DNAMeth"),
"MM"=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","DNAMeth"),
"MCL"=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","DNAMeth"),
"LGG"=c('H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input','DNAMeth','H3K9ac'))
j = list("NPC"=c(1,2,3,5,6),
"MSC"=c(1,2,3,5,6),
"MES"=c(1,2,3,5,6),
"TSC"=c(1,2,3,5,6),
"CLL"=c(1,2,6,8),
"PTC"=c(1,2,4,6,8),
"CRC"=c(1,2,4,6,8),
"AML"=1:6,
"ALL"=1:6,
"mCLL"=1:6,
"MM"=1:6,
"MCL"=1:6,
"LGG"=c(1,2,4,6))
df = do.call("rbind",lapply(1:13, function(i) data.frame(melt(prcomp(data[[i]]$Dobs)$rotation),rep(lbl[[nm[i]]][j[[nm[i]]]],ncol(data[[i]]$Dobs)),nm[i])))
colnames(df)=c("id","PC","Loadings","Mark","type")
ggplot(df, aes(x=Mark, y=Loadings, fill=Mark)) + geom_bar(stat='identity')+facet_grid(PC~type) +coord_flip()+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5))
```
### PC variance
```{r}
df = do.call("rbind",lapply(1:13, function(i) cbind(melt(data[[i]]$proj),nm[i])))
df = df[!is.na(df[,3]),]
df = df[df[,1]=="percent_var",2:4]
colnames(df)=c("PC","var","type")
ggplot(df, aes(x=PC, y=var, colour=type)) + geom_point(size=2) + labs(x="PC", y="Variance") + scale_y_continuous(labels = scales::percent) + theme(legend.position=c(0.85,0.56))
```
### ECDF
```{r}
markers$mCLL = markers$CLL
df=do.call("rbind",lapply(c(3,4,5,7,10,11,12,13), function(i){
df=data.frame(value=sort(match(markers[[nm[i]]],data[[i]]$pg$PC1)),type='Irene')
df=rbind(df,data.frame(value=sort(match(markers[[nm[i]]],data[[i]]$pg$prom)),type='Promoter'))
cbind(df,case=nm[i])
}))
ggplot(df,aes(value, colour=type))+stat_ecdf(pad = FALSE)+labs(x="Rank", y="ECDF") +facet_wrap(.~case)+ theme(legend.position=c(0.8,0.2))
```
### Oncogenes
```{r}
wtest = function(m,n,l=17888) -log(wilcox.test(sort(match(n,m)), seq(1,l,ceiling(l/length(n))))$p.value)
og =markers$OG
tsg=markers$TSG
hkg=markers$HKG
tn=c("OG","TSG","HKG")
df=do.call("rbind",lapply(1:13, function(i) {
rbind(data.frame(Pvalue=c(wtest(data[[i]]$pg$PC1,og),wtest(data[[i]]$pg$PC1,tsg),wtest(data[[i]]$pg$PC1,hkg)), name=tn, type=nm[i], Class=cl[i], Rank="Irene"),
data.frame(Pvalue=c(wtest(data[[i]]$pg$prom,og),wtest(data[[i]]$pg$prom,tsg),wtest(data[[i]]$pg$prom,hkg)), name=tn, type=nm[i], Class=cl[i], Rank="Promoter"))
}))
ggplot(df, aes(Rank, type, fill = Pvalue)) + geom_tile(colour = "white") + 
  facet_grid(Class~name,scales='free') + scale_fill_gradient2(low="white", high="red",name="-logP") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+labs(x="",y="")
```
### Frequency of oncogenes
```{r}
vg=function(gs) unlist(lapply(gs,function(g) V(g)$name))
lf=function(x) length(which(x))/length(x)
sp=function(x) {x=melt(x);colnames(x)=c("type","perc","value");x$type=nm[cl=="Cancer"][x$type];x}
dname=function(i) seq(10,50,5)[i]
genes=lapply(which(cl=="Cancer"), function(i){
lapply(seq(0.9,0.5,-0.05),function(w) vg(exportMultinets(edgeRank(data[[i]]$pg$PC1,hprd, cutoff=w), 15)))
})
data=data[which(cl=="Cancer")]
df=do.call("rbind",lapply(1:length(data), function(i){
c(unlist(lapply(1:9,function(j) lf(genes[[i]][[j]] %in% og))),
unlist(lapply(1:9,function(j) lf(data[[i]]$pg$PC1[1:length(genes[[i]][[j]])] %in% og))))
}))
data=list(cll, crc, glm, msc, npc, ptc, tsc, mcll)
nm = c("CLL", "CRC", "LGG", "MSC", "NPC", "PTC", "TSC", "mCLL")
genes=lapply(data, function(d){
lapply(seq(0.9,0.5,-0.05),function(w) vg(exportMultinets(edgeRank(d$pg$PC1,hprd, cutoff=w), 15)))
})
df=do.call("rbind",lapply(1:length(data), function(i){
c(unlist(lapply(1:9,function(j) lf(genes[[i]][[j]] %in% markers[[nm[i]]]))),
unlist(lapply(1:9,function(j) lf(data[[i]]$pg$PC1[1:length(genes[[i]][[j]])] %in% markers[[nm[i]]]))))
}))
ggplot() + geom_line(data=sp(df[,1:9]), aes(perc, value, col=type), size=1) +
  geom_line(data=sp(df[,10:18]), aes(perc, value, col=type), linetype="dotted", size=1) +
  labs(x="Top rank %",y="Frequency") +scale_x_continuous(breaks=seq(1,9,2),label=dname)
```
### Tissue specificity
```{r}
v=lapply(strsplit(readLines("ARCHS4_Tissues.txt"),"\t"), function(d) d[3:length(d)])
names(v)=lapply(strsplit(readLines("ARCHS4_Tissues.txt"),"\t"), function(d) d[1])
auc=lapply(data, function(d){
unlist(lapply(v, function(a) getAUC(d$pg$PC1,a)-getAUC(d$pg$prom,a)))
})
names(auc)=nm
auc=data.frame(auc)
auc=melt(as.matrix(auc))
names(auc)=c("Tissues","Cases","AUC")
auc$Cases = factor(as.character(auc$Cases), levels = c("ALL", "AML", "CLL", "mCLL", "MCL", "MM", "CRC", "LGG", "PTC", "MES", "MSC", "NPC", "TSC"))
ggplot(data = auc, aes(x=Cases, y=Tissues, fill=AUC)) +theme(axis.text.x = element_text(angle = 30, hjust = 1))+ geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
```
### Rewiring
```{r}
k=c(3,4,5,7,10,11,12,13)
df=do.call("rbind",lapply(k, function(i){
data.frame(auc=unlist(lapply(1:100, function(d) getAUC(pageRank(data[[i]]$gr[,c("id","PC1")], H1, rewire=TRUE)$PC1, markers[[nm[i]]]))), type=nm[i])
}))
df1=data.frame(auc=unlist(lapply(k, function(i){getAUC(data[[i]]$pg$PC1, markers[[nm[i]]])})), type=nm[k])
df2=data.frame(auc=unlist(lapply(k, function(i){getAUC(data[[i]]$pg$prom,markers[[nm[i]]])})), type=nm[k])
ggplot() + geom_boxplot(data=df, aes(type, auc),outlier.size = NA) + geom_point(data=df1, aes(type, auc),size=4,shape=19,col="red") + geom_point(data=df2, aes(type, auc),size=4,shape=18,col="red") + labs(x="Test cases",y="AUC") + theme(axis.text.x = element_text(angle = 30, hjust = 1))
write.table(df,file="irene-supp/AUC/rnd.txt",sep="\t",quote=F)
```
### Well-known oncogenes

```{r}
genes=c("BRAF","EGFR","KRAS","MYC","PAX5","PTEN","TP53")
df=do.call("rbind",lapply(c(1:7,9,12), function(i){
df=data.frame(Irene=match(genes,rev(data[[i]]$pg$PC1)),Promoter=match(genes,rev(data[[i]]$pg$prom)),gene=genes,type=nm[i])
}))
pos1=function(d) as.integer(d/17888)
ggplot(df, aes(x=Irene, y=Promoter)) +geom_point(aes(col=gene))+facet_wrap(~type) +geom_abline(slope=1,linetype="dashed")+ xlim(c(0, 17888))+ylim(c(0, 17888))+geom_text(aes(label=gene), size=2.5,vjust=1)+theme(axis.ticks=element_blank(),legend.position="none")+scale_x_continuous(label=pos1,breaks=c(0, 17888))+scale_y_continuous(label=pos1,breaks=c(0, 17888))
```

# Plot figures for the dissertation
### Alternative theme

```{r}
tsize <- 9; tsize2 <- 8; pchsize <- 1; lkey <- 3.5 # lsize <- 0.3
theme_set(theme_bw())
theme_update(panel.grid.minor = element_line(colour = NA),
panel.grid.major = element_line(colour = NA))

#theme_set(theme_bw(base_size=22))
theme_set(theme_bw())
theme_update(
plot.margin = unit(c(0.4,0.5,0.1,0), "lines"),
panel.spacing = unit(0.25, "lines"),
panel.grid.minor = element_line(colour = NA),
panel.grid.major = element_line(colour = NA),
panel.background=element_rect(fill = NA, colour = "black"),
panel.border = element_rect(colour="black", size=0),
plot.title = element_text(size = tsize, vjust = 0.5, hjust=0.5),
axis.title.x = element_text(size = tsize, vjust = 0.35),
axis.title.y = element_text(size = tsize, hjust = 0.5, vjust = 0.4, angle = 90),
axis.text.x = element_text(size = tsize2, margin=margin(0,0,0.1,0,"mm")),
axis.text.y = element_text(size = tsize2, margin=margin(0,0,0,0.1,"mm")),
axis.ticks.length=unit(1, units="mm"),
#axis.line = element_segment(colour = ‘black’, size = 1),
legend.key=element_rect(colour = NA),
legend.title=element_text(size = tsize2-1, hjust = 0),
legend.text=element_text(size = tsize2-1, hjust = 0),
legend.background=element_rect(colour=NA, size=0),
legend.spacing = unit(0, "mm"),
legend.key.size=unit(lkey,"mm"),
legend.key.width = unit(lkey*1.5, "mm"),
strip.background = element_rect(fill = NA, linetype=NULL, size=0, colour="white"),
strip.text.x = element_text(size=tsize, vjust=0.7, hjust= 0.5)
)
```

### Data validity

Distribution densities

```{r}
library(quantro)
load(paste0(s,'.hg19.rda'))
data=data*1000/abs(bed[,3]-bed[,2])
data=log(data+1)
mark=str_extract(meta$file,"H3K\\w+|WGB|frac|Input")
mark[mark=='WGB' | mark=='frac']='WGBS'
mark=as.factor(mark)
svglite(paste0(s,'.svg'),5,5)
matdensity(data,mark,main=s,xlab="",ylab="")
legend('topright', levels(mark), col = RColorBrewer::brewer.pal(8,"Dark2"), lty = 1, lwd = 3)
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
svglite("dExp.svg",3,3)
ggplot(df, aes(Mark, value)) + geom_hline(yintercept=0, col='red',linetype='dashed') + geom_point(aes(col=Differences), size=2, alpha=0.6) + labs(x="", y="Spearman r") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
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
svglite('d2.svg',4,3)
ggplot(df, aes(x=Mark, y=Loadings, fill=Mark)) + geom_bar(stat='identity')+facet_grid(PC~type) +coord_flip()+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5))
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
ggplot(x, aes(x=pe, y=po)) +geom_point(aes(col=gene))+facet_wrap(~cell) +geom_abline(slope=1,linetype="dashed")+ xlim(c(0, 17888))+ylim(c(0, 17888))+labs(x="PromEnh", y="PromOnly")+geom_text(aes(label=gene), size=2.5,vjust=2)+theme(axis.ticks=element_blank(),legend.position="none")
```

Wilcoxon-Mann-Whitney test of OG/TSG/HKG ranks comparing to uniformly distributed ranks

```{r}
svglite("oncoRank.svg",3,3)
ggplot(df, aes(Rank, type, fill = Pvalue)) + geom_tile(colour = "white") + 
  facet_grid(Class~name,scales='free') + scale_fill_gradient2(low="white", high="red", midpoint = 0.05) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+labs(x="",y="")
```

Rank and p-values of tissue specificity of PromEnh and PromOnly ranks

```{r}
svglite("tissueRank.svg",6,6)
ggplot(m, aes(variable, ts, fill = score)) + geom_tile(colour = "white") + geom_text(aes(label=rank)) + facet_grid(fs~.,scales='free') + scale_fill_gradient2(low="white", high="red") +
  labs(x="",y="")
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

