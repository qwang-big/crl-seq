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


