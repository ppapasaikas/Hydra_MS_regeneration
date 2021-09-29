library(gsubfn)
library(rsvd)
library(RColorBrewer)
library(mgcv)
library(parallel)
library(dimRed)
library(Rtsne)
library(gplots)
library(plyr)
library(matrixStats)

# @param m  a (potentially sparse) gene x cells count matrix
# @param f  a number between 0 and 1 the fraction of overdispersed genes to keep
# @return a vector of rownames of genes to keep.
select_variable_genes<-function(m,f) {
  zeroes=which(rowSums(m) <= max(1,min(rowSums(m))) )
  all.nz.genes   <- rownames(m[-zeroes,])
  m <- m[all.nz.genes,]
  df <- data.frame(mean = rowMeans(m + 1/ncol(m)), cv = apply(m,1,sd) / rowMeans(m + 1/ncol(m)), var = apply(m,1,var))
  df$dispersion <- with(df, var/mean)
  df$mean_bin <- with(df, cut(mean, breaks = c(-Inf, unique(quantile(mean, seq(0.1,1,0.02), na.rm = TRUE)), Inf)))
  var_by_bin <- data.frame(mean_bin = factor(levels(df$mean_bin), levels = levels(df$mean_bin)),
                           bin_median = as.numeric(tapply(df$dispersion, df$mean_bin, stats::median)),
                           bin_mad = as.numeric(tapply(df$dispersion, df$mean_bin, stats::mad)))[table(df$mean_bin) > 0,]
  df$bin_disp_median <- var_by_bin$bin_median[match(df$mean_bin, var_by_bin$mean_bin)]
  df$bin_disp_mad <- var_by_bin$bin_mad[match(df$mean_bin, var_by_bin$mean_bin)]
  df$dispersion_norm <- with(df, (dispersion - bin_disp_median)/(bin_disp_mad + 0.01) )
  
  n_genes_keep=ceiling(f*nrow(m) ) #In the end retain only the top 100*f% overdispersed genes
  disp_cut_off <- sort(df$dispersion_norm,decreasing=TRUE)[n_genes_keep]
  genes_keep <- which(df$dispersion_norm >= disp_cut_off)
  return(rownames(m)[genes_keep])
}



### Load the log cpm matrix:
# Log_SampleCounts <- read.delim(Data/SSpheres_Log_SampleCounts_CNTR.txt)
# NSampleCounts <- 2^(Log_SampleCounts+1)


#### Gene selection sets:
nzeroes <- which(rowSums(NSampleCounts > 0) > 4 ) # Genes that are non-zero in more than 3 samples.
expr <- which(rowSums(NSampleCounts) > 8 ) # Genes above arbitrary cpm threshold
nz.genes   <- rownames(NSampleCounts[nzeroes,])
exp.genes <- intersect( rownames(NSampleCounts[expr,]), nz.genes) 
var.genes   <- select_variable_genes(NSampleCounts[exp.genes,S1],0.5) # Select the top 50% variable genes according to the mean variance trend




#### Remove PC components that are not functionals of time:
times <- as.numeric(gsub("h_.*","",colnames(Log_SampleCounts) ))
use.genes <- var.genes

Data <- Log_SampleCounts[use.genes,] # Only for demonstration DO NOT USE

centering.vector <- apply(Data,1,mean)
fit <- prcomp(t(Data),center = TRUE, scale = FALSE)

test <- apply(t(fit$x),1, function(x) { gam( x ~ s(times), method = "GCV.Cp", gamma=1.0)})
gpvals <-  unlist(lapply(test,function(x)  summary(x)$s.pv  )  )
adj.gpvals<- p.adjust(gpvals,method="fdr")

pvals <- adj.gpvals

#pc.use <- which( pvals < 0.05 )
pc.var.prop <- (fit$sdev)**2/sum((fit$sdev)**2)
pc.use <- which( pvals < 0.05  )
pc.var.prop.use <-  sum(pc.var.prop[pc.use])
DataRe <- t(fit$x[,pc.use] %*% t(fit$rotation[,pc.use]))
DataRe <- DataRe + centering.vector # Reconstructed dataset after PC removal






################ PSEUDO TIME ORDERING:

####### 1. Identify genes with good fits  against time (time functionals) using Generalized Additive Models
use.genes=intersect(var.genes, rownames(DataRe) )

GSmooth=list()
GAM=list()
GAM.shuffl=list()
ADJ.PVALS=list()
interval.weights=c()

intervals=list( c(1:7), c(2:8), c(3:9), c(4:10), c(5:11), c(6:12), c(7:13), c(8:14), c(9:15),
                c(10:16), c(11:17), c(12:18), c(13:19), c(14:20), c(15:21), c(16:22)
)
interval.samples <- list()

cl <- makeCluster(getOption("cl.cores", 32)  )
clusterExport(cl, c("times","gam"),envir=environment())

for(i in 1:length(intervals)){
  pv.thr <- 0.001
  sg <- use.genes
  time <- intervals[[i]]
  
  use.samples <- grep  ( paste0( "^",time,"h",collapse="|" ), colnames(NSampleCounts) )
  use.samplesn <- colnames(NSampleCounts)[use.samples]
  interval.samples[[i]] <-  length(use.samples)
  
  Y <- t(apply(Data[,use.samplesn],1, function(x)  {-1+2*(x-min(x,na.rm=TRUE))/ diff(range(x,na.rm=TRUE))} ))
  Y <- Y[which( rowSums(is.na(Y)) < 1  ),]
  sg <- intersect(sg,rownames(Y))
  
  clusterExport(cl, c("use.samples"),envir=environment())
  GAM[[i]] <- parApply(cl,Y[sg,],1, function(x) { gam( x ~ s(times[use.samples], k=4), method = "GCV.Cp", drop.intercept=TRUE, gamma=1.25)})
  
  FITV  <-  unlist(lapply(GAM[[i]],function(x) diff(range(x$fitted.values ))  ))
  PVALS  <-  unlist(parLapply(cl,GAM[[i]],function(x)  summary(x)$s.pv  )  )
  ADJ.PVALS[[i]] <- p.adjust(PVALS,method="fdr")
  pv.thr <- max(pv.thr, quantile(ADJ.PVALS[[i]],0.01) )
  PVAL.sel <- names(which(ADJ.PVALS[[i]] < pv.thr ))
  FITV.sel <- names(which(FITV > 0.5))
  SEL <- intersect (PVAL.sel, FITV.sel)
  GSmooth[[i]] <- SEL 
  interval.weights[i] <- -mean( log10(ADJ.PVALS[[i]][ GSmooth[[i]]  ])) * sqrt(length(GSmooth[[i]])) /10  #evidence support for PTordering in the interval
  cat("Interval: ",i,"\n")
  
}
stopCluster(cl)








####### 2. PT ordering using isomap:

use.genes.list <- GSmooth

pseudo.time <- lapply(interval.samples, function(x){ rep(0, unlist(x) ) })

start.tps <- unlist(lapply(intervals,function(x) min(x) ))
end.tps <- unlist(lapply(intervals,function(x) max(x) ))

time.annot=gsub("h_.*","",colnames(NSampleCounts))

par(mfrow=c(4,1),xpd=TRUE,mar=c(4,4,2,8))
knns=sample( c(6:7),100,replace=TRUE ) # Perform multiple isomap iteration using a knn of 6 or 7
nexec=0
for(i in 1:(length(intervals))    ) {
  time <- intervals[[i]]
  start.tp=start.tps[i]
  end.tp=end.tps[i]
  use.samples <- grep  ( paste0( "^",time,"h",collapse="|" ), colnames(NSampleCounts) )
  use.samplesn <- colnames(NSampleCounts)[use.samples]
  
  use.genes <- unlist(use.genes.list[[i]])
  
  for (knn in knns) {
    is.error <- TRUE
    while(is.error){
      nexec<-nexec+1
      tryCatch({
        im <- dimRed::embed( cor(DataRe[sample(use.genes, ceiling(length(use.genes)/2)  ),use.samplesn] )  , "Isomap", knn=knn) 
        is.error <- FALSE
      },error=function(e){cat("Error\n")},
      finally={} 
      )
    }
    
    imX=im@data@data
    if (imX[1,1] > imX[nrow(imX),1] ) {imX[,1]=-imX[,1]}
    
    pt <- (start.tp-1.0)  + ((end.tp-start.tp+1) * (imX[,1]-min(imX[,1]))/diff(range(imX[,1])))
    pseudo.time[[i]] <-pseudo.time[[i]] + (1/length(knns)) * pt
  }
  names(pseudo.time[[i]] ) <- colnames(NSampleCounts)[use.samples]
}    

pseudo.time.matrix <- rbind.fill( lapply( pseudo.time,function(x) data.frame(t(x))  )   )

# Weighted consensus pseudo.time assignment:
cons.pseudo.time <- rep(0,ncol(DataRe)  )
names(cons.pseudo.time) <- colnames(DataRe)
for (smpl in colnames(DataRe) ){
  w_tot <- 0
  for (i in 1:length(pseudo.time)) {
    if ( smpl %in% names(pseudo.time[[i]])      ){
      w <- interval.weights[i]
      cons.pseudo.time[smpl] <- cons.pseudo.time[smpl] + w * pseudo.time[[i]][smpl]  
      w_tot <- w_tot+w
    }
  }
  cons.pseudo.time[smpl] <- cons.pseudo.time[smpl]  / w_tot
}



