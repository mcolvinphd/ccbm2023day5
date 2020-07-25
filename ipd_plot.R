pops <- read.table("ipd_pop.txt",header=T)
name.list <- names(pops)
pops[1,]
ymax=max(pops)
ymin=min(pops)
plot(pops[,1],type="b",col=1,ylim=c(ymin,ymax),pch=1,lwd=2,ylab="Population",xlab="Generation")
for (i in seq(2,length(name.list))) {
    lines(pops[,i],col=i,lwd=2)
    points(pops[,i],col=i,pch=i)
}
legend("topright",name.list,pch=c(1:length(name.list)),col=c(1:length(name.list)))
