bedL<-list()
#bdn <- c('hESC','mESC','hIMR90','mCX')
bedL[['hESC']]<-read.csv('topological_domains/human.ESC.norm.combined','\t',header=F)
bedL[['mESC']]<-read.csv('topological_domains/mouse.ESC.norm.combined','\t',header=F)
bedL[['hIMR90']]<-read.csv('topological_domains/human.IMR90.norm.combined','\t',header=F)
bedL[['mCX']]<-read.csv('topological_domains/mouse.CORTEX.norm.combined','\t',header=F)

diff<-list()
for(bd in 1:length(bedL)){
	bed = bedL[[bd]]
	diff[[names(bedL)[bd]]] <- bed[,3] - bed[,2]
}


pdf('topological.domian.size.distribution.pdf')
plot(density(diff$hESC),main='Distribution of Topological Domain Size',
	xlab='Size of individual topological domains (nucleotides/domain)')
lines(density(diff$mESC),col='red')
lines(density(diff$hIMR90),col='green')
lines(density(diff$mCX),col='blue')
legend(x=3.5e6,y=8e-7,legend=names(bedL),fill=c('black','red','green','blue'))
dev.off()
