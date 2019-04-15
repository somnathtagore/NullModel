#############################
#Script = nullmodel.R
#Author = Somnath Tagore
#Last Update = 03.20.2019
##############################

#Packages
library(viper)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(matrixStats)

#load regulon and prune
load('coad_regulon.rda')
pregul <- pruneRegulon(regul, cutoff = 50)

# tcga data
tcga.tumor<-read.table(file="../batch_corrected_data/data/batch-correct/tcga-tumor/coad-tcga-tumor.txt", header=TRUE)
tcga.tumor[1:3,1:3]

# gtex data
gtex.tissue <- read.table("../batch_corrected_data/data/batch-correct/gtex/colon-gtex.txt", header=TRUE)
gtex.tissue[1:3,1:3]

# tpm normalization
for(i in 1:ncol(tcga.tumor)){
  tcga.tumor[,i] <- 1E6*tcga.tumor[,i]/sum(tcga.tumor[,i])
}

for(i in 1:ncol(gtex.tissue)){
  gtex.tissue[,i] <- 1E6*gtex.tissue[,i]/sum(gtex.tissue[,i])
}

# identical genes

identical(rownames(tcga.tumor),rownames(gtex.tissue))
table(rownames(tcga.tumor)%in%rownames(gtex.tissue))
common_genes <- intersect(rownames(tcga.tumor),rownames(gtex.tissue))
length(common_genes)

tcga.tumor <- tcga.tumor[match(common_genes,rownames(tcga.tumor)),]
gtex.tissue <- gtex.tissue[match(common_genes,rownames(gtex.tissue)),]

# mix tissues
mixed.tcga.tumor.gtex.tissue <- cbind(tcga.tumor, gtex.tissue)
dim(mixed.tcga.tumor.gtex.tissue)

# viper signature
vipsig <- viperSignature(eset = as.matrix(mixed.tcga.tumor.gtex.tissue), 
                         ref = as.matrix(gtex.tissue), per = 1000, method = "ttest", verbose = TRUE)
# run viper
vipermat <- viper(vipsig, pregul,minsize = 5)

dim(vipermat)


gene.df <- tibble(#sample = colnames(vipermat.mix),
  sample = colnames(vipermat),
  #activity = vipermat.mix[match("1956", rownames(vipermat.mix)),],
  activity = vipermat[match("1956", rownames(vipermat)),],
  type = rep(c('tumor', 'normal'), c(285,339)))
ggplot(data = gene.df, mapping = aes(x = activity, fill = type)) +
  geom_density(alpha=0.4, position="identity")

#highval.regs <-c('1956','667','1756','6794','7157','7273')
highval.regs <-c('10690',
                 '6664',
                 '7157',
                 '285175'
                 )
highval.regs.names <-c(
  'FUT9',
  'SOX11',
  'TP53',
  'UNC80'
  )


install.packages("gridExtra")
library(gridExtra)

#par(mfrow=c(2,2))
mygrobs<-list()

saveindex<-list()
savesaveindex<-list()


for (i in 1:length(highval.regs)){
  if(!identical(vipermat[row.names(vipermat)==highval.regs[i]],numeric(0)))
  {
    saveindex[[i]]<-highval.regs[i] 
    #savesaveindex[[i]]<-highval.regs.names[i]
  }
  
}
for (i in 1:length(highval.regs)){
  if(!identical(vipermat[row.names(vipermat)==highval.regs[i]],numeric(0)))
  {
    #saveindex[[i]]<-highval.regs[i] 
    savesaveindex[[i]]<-highval.regs.names[i]
  }
  
}
saveindex<-saveindex[!sapply(saveindex,is.null)]
savesaveindex<-savesaveindex[!sapply(savesaveindex,is.null)]

#Density plot tumor+gtex
for (i in 1:length(saveindex)) {
  
  df <- tibble(#sample = colnames(vipermat.mix),
    sample = colnames(vipermat),
    #activity = vipermat.mix[match(highval.regs[i], rownames(vipermat.mix)),],
    #activity = vipermat[rownames(vipermat)==highval.regs[i]],
    activity = vipermat[match(saveindex[i], rownames(vipermat)),],
    type = rep(c('tumor', 'normal'), c(285,339)))
  
  #  mygrobs[[i]]<-ggplot(data = df, mapping = aes(x = activity, fill = type)) +
  
  mygrobs[[i]]<-ggplot(data = df, mapping = aes(x = activity, fill = type)) +  
    geom_density(alpha=0.4, position="identity") +
    ggtitle(savesaveindex[i])
  #return(identical(vipermat[row.names(vipermat)==highval.regs[i]],numeric(0)))
  
  
}
pdf(file = "density.pdf", width = 6.25, height = 4, 
    family = "Times", pointsize = 12)
#do.call(grid.arrange,mygrobs)
marrangeGrob(mygrobs,nrow=2,ncol=2)
#grid.arrange(mygrobs,ncol=2)
dev.off()

#Mixture model implementation
library(mixtools)
mygrobsem<-list()
gg.mixEM <- function(EM) {
  require(ggplot2)
  x       <- with(EM,seq(min(x),max(x),len=1000))
  pars    <- with(EM,data.frame(comp=colnames(posterior), mu, sigma,lambda))
  em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
  em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))
  ggplot(data.frame(x=EM$x),aes(x,y=..density..)) + 
    #type = rep(c('tumor', 'normal'), c(503,313)) +
    #geom_histogram(fill=NA,color="black")+
    geom_polygon(data=em.df,aes(x,y,fill=comp),color="grey50", alpha=0.5)+
    #geom_polygon(data=em.df,aes(x,y,fill=type),color="grey50", alpha=0.5)+
    scale_fill_discrete("Component\nMeans",labels=format(em.df$mu,digits=3))+
    theme_bw()
}

set.seed(1)    # for reproducible example
# b <- rnorm(2000000, mean=c(8,17), sd=2)
# c <- b[vipermat[match("1956", rownames(vipermat)),](length(b), 1000000) ]
#c2 <- normalmixEM(vipermat[match("1956", rownames(vipermat)),], k=3, lambda=NULL, mu=NULL, sigma=NULL) 
for (i in 1:length(saveindex)) {
  
  activity <- normalmixEM(vipermat[match(saveindex[i], rownames(vipermat)),], 
                          k=2, lambda=NULL, mu=NULL, sigma=NULL)
  mygrobsem[[i]]<-gg.mixEM(activity)     
}
pdf(file = "gmm.pdf", width = 6.25, height = 4, 
    family = "Times", pointsize = 12)
#do.call(grid.arrange,mygrobs)
marrangeGrob(mygrobsem,nrow=2,ncol=2)
#grid.arrange(mygrobs,ncol=2)
dev.off()

#Other features

