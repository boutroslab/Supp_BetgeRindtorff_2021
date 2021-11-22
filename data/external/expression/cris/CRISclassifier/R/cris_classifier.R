cris_classifier<-function(
# file I/O
input.exp.filename,
output.name="CRIS",
# resampling to generate null dist
nresmpl=1000,

# Seed
rnd.seed=7392854
)
{
#  suppressWarnings()
  temp.nn.wt = "TRUE"
  dist.selection="cosine"
  GenePattern.output="TRUE"
  # Advanced setting

  norm.method="row.std" # "row.std.ref","ratio.ref"
  within.sig="FALSE"
  dchip.output="FALSE"
  signature.heatmap="TRUE"
  FDR.sample.bar=0.2 # NA if not needed
  plot.FDR="TRUE"
  col.range=3         # SD in heatmap
  heatmap.legend=signature.heatmap
  histgram.null.dist="FALSE" # histgram of null dist for the distance
  hist.br=30

  # for GenePattern
  
  nresmpl<-as.numeric(nresmpl)
  col.range<-as.numeric(col.range)
  hist.br<-as.numeric(hist.br)  # bin number for resampled dist histgram
  rnd.seed <- as.numeric(rnd.seed)
  if (FDR.sample.bar!="NA"){
    FDR.sample.bar <- as.numeric(FDR.sample.bar)
    if (is.numeric(FDR.sample.bar)==FALSE){
      stop("### Provide numerical value (0~1) for FDR.sample.bar! ###")
    }
  }

  # set random seed
  set.seed(rnd.seed)

  ### input ###

  # selected features used for prediction
  data("features")

  ## file format check
  if (length(features[1,])!=3 & length(features[1,])!=4){
    stop("### Please use features file format! ###")
  }
  if (length(features[1,])<4 & temp.nn.wt=="TRUE"){
    temp.nn.wt <- "FALSE"
  }
  third.col<-rownames(table(features[,3]))
  if (is.na(as.numeric(third.col[1]))){
    stop("### TRUEhe 3rd column of feature file should be numerical! ###")
  }

  feat.col.names<-colnames(features)
  feat.col.names[1:2]<-c("ProbeID","GeneName")
  colnames(features)<-feat.col.names

  num.features<-length(features[,1])
  num.cls<-length(table(features[,3]))
  feature.col.num <- length(features[1,])

  ord<-seq(1:num.features)
  features<-cbind(ord,features)  # add order column to "features"

  # expression data
  ## file format check
  if (length(grep("gct$", input.exp.filename))>0){
      exp.dataset<-read.delim(input.exp.filename,header=TRUE,skip=2,check.names=FALSE)
      colnames(exp.dataset)[1:2] <- c("ProbeID","GeneName")
  }else{
     exp.dataset<-read.delim(input.exp.filename,header=TRUE,check.names=FALSE)
     colnames(exp.dataset)[1] <- c("GeneName")
  }

  ## Other dataset's mean & SD for row normalization (optional)

  if (length(grep("gct$", input.exp.filename))>0){
    ProbeID<-exp.dataset[,1]
    gene.names<-exp.dataset[,2]
    num.samples<-(length(exp.dataset[1,])-2)
    exp.dataset<-exp.dataset[-c(1:2)]
    exp.for.sample.names<-read.delim(input.exp.filename,header=FALSE,skip=2)  # read sample names
    sample.names<-as.vector(as.matrix(exp.for.sample.names[1,3:length(exp.for.sample.names[1,])]))
  }else{
    ProbeID<-exp.dataset[,1]
    gene.names<-exp.dataset[,1]
    num.samples<-(length(exp.dataset[1,])-1)
    exp.dataset<-exp.dataset[-c(1)]
     exp.for.sample.names<-read.delim(input.exp.filename,header=TRUE,check.names=FALSE)
     sample.names<-as.vector(colnames(exp.for.sample.names)[-1])
  }
  print(num.samples)
  # row normalize

  normed.exp.dataset<-exp.dataset

  if (norm.method=="row.std"){
    exp.mean <- apply(exp.dataset,1,mean,na.rm=TRUE)
    exp.sd <- apply(exp.dataset,1,sd,na.rm=TRUE)
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd   
  }

  normed.exp.dataset<-cbind(ProbeID,normed.exp.dataset)

  # extract features from normed.exp.dataset
  out <- list()
  out[[1]] <- features
  out[[2]] <- normed.exp.dataset
  #return(out)
  exp.dataset.extract<-merge(features,normed.exp.dataset,sort=FALSE)
  if (length(exp.dataset.extract[,1])<1){
    stop("### No matched probes! ###")
  }
  #return(exp.dataset.extract)
  order.extract<-order(exp.dataset.extract[,2])
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  order.extract.after<-exp.dataset.extract[,2]
  exp.dataset.extract<-exp.dataset.extract[-2]

  if (temp.nn.wt=="FALSE"){
    features.extract<-exp.dataset.extract[,1:3]
    if (feature.col.num==4){
      exp.dataset.extract <- exp.dataset.extract[-4]
    }
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    num.features.extract<-length(features.extract[,1])

    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:3)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }

#  temp.nn.wt.vector <- rep(1,num.features)

  if (temp.nn.wt=="TRUE" & num.cls==2){
    features.extract<-exp.dataset.extract[,1:4]
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)

#    if (is.numeric(features[,4])){
      temp.nn.wt.vector <- as.numeric(as.vector(features.extract[,5]))
#    }else{
    if (is.numeric(temp.nn.wt.vector)==FALSE){
      stop("# Please use numeric values in 4th column!#")
    }

    num.features.extract<-length(features.extract[,1])

    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }

  # make template

  for (i in 1:num.cls){
    temp.temp<-as.numeric(as.vector(features.extract[,4]))
    temp.temp[temp.temp!=i]<-0
    temp.temp[temp.temp==i]<-1
    eval(parse(text=paste("temp.",i,"<-temp.temp",sep="")))
#    eval(parse(text=paste("temp\.",i,"<-temp\.temp",sep="")))  ### for < R-2.4.0
  }

  # weighted template (only for 2cls)

  if (temp.nn.wt=="TRUE" & num.cls==2){
    temp.1 <- temp.nn.wt.vector
    temp.2 <- -temp.nn.wt.vector
  }

  ### compute distance and p-value ###

  predict.label<-vector(length=num.samples,mode="numeric")
  dist.to.template<-vector(length=num.samples,mode="numeric")
  dist.to.cls1<-vector(length=num.samples,mode="numeric")

  rnd.feature.matrix<-matrix(0,nrow=num.features.extract,ncol=nresmpl)

  perm.dist.vector<-vector(length=nresmpl*num.cls,mode="numeric")
  nominal.p<-vector(length=num.samples,mode="numeric")
  BH.FDR<-vector(length=num.samples,mode="numeric")
  Bonferroni.p<-vector(length=num.samples,mode="numeric")
  #return(exp.dataset.extract)
  for (i in 1:num.samples){

    print(paste("sample # ",i,sep=""))

    current.sample <- as.vector(exp.dataset.extract[,i])

    # compute original distance

    orig.dist.to.all.temp <- vector(length=num.cls,mode="numeric")

    if (temp.nn.wt=="TRUE"){   # weight sample data
      current.sample <- current.sample*abs(temp.nn.wt.vector)
    }

    if (dist.selection=="cosine"){
      for (o in 1:num.cls){      # compute distance to all templates
        current.temp <- vector()
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        orig.dist.to.all.temp[o]<-sum(current.temp*current.sample)/
                  (sqrt(sum(current.temp^2))*sqrt(sum(current.sample^2)))
      }
    }
    if (dist.selection=="correlation"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
#        eval(parse(text=paste("current\.temp <- temp\.",o,sep="")))  ### for < R-2.4.0
        orig.dist.to.all.temp[o] <- cor(current.temp,current.sample,method="pearson",use="complete.obs")
      }
    }

    if (num.cls==2){           # find nearest neighbor (2 classes)
      if (orig.dist.to.all.temp[1]>=orig.dist.to.all.temp[2]){
        predict.label[i]<-1
        dist.to.template[i]<-1-orig.dist.to.all.temp[1]
        dist.to.cls1[i]<--(orig.dist.to.all.temp[1]+1)
      }
      if (orig.dist.to.all.temp[1]<orig.dist.to.all.temp[2]){
        predict.label[i]<-2
        dist.to.template[i]<-1-orig.dist.to.all.temp[2]
        dist.to.cls1[i]<-orig.dist.to.all.temp[2]+1
      }
    }

    if (num.cls>2){
      for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
        if (is.na(orig.dist.to.all.temp[o])!=TRUE){
          if (orig.dist.to.all.temp[o]==max(orig.dist.to.all.temp,na.rm=TRUE)){
            predict.label[i]<-o
            dist.to.template[i]<-1-orig.dist.to.all.temp[o]
            dist.to.cls1[i]<-(1-orig.dist.to.all.temp[o])+o
          }
        }
      }
    }
    #return(orig.dist.to.all.temp)
    # permutation test

    if (within.sig=="FALSE"){     # generate resampled features from all probes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(normed.exp.dataset[,(i+1)],num.features.extract,replace=FALSE)
      }
    }
    if (within.sig=="TRUE"){     # generate resampled features from only signature genes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(exp.dataset.extract[,i],num.features.extract,replace=FALSE)
      }
    }

    if (temp.nn.wt=="TRUE" & num.cls==2){
      rnd.feature.matrix <- rnd.feature.matrix*abs(temp.nn.wt.vector)
    }

    # compute distance to all templates
    if (dist.selection=="cosine"){          # cosine
      for (res in 1:num.cls){
        temp.resmpl <- vector()
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))

        prod.sum<-apply(t(t(rnd.feature.matrix)*temp.resmpl),2,sum)

        data.sq.sum<-apply(rnd.feature.matrix^2,2,sum)
        temp.sq.sum<-sum(temp.resmpl^2)

        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
               (1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
      }
    }
    #return(perm.dist.vector)
    if (dist.selection=="correlation"){          # correlation
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
                 (1-as.vector(cor(rnd.feature.matrix,temp.resmpl,method="pearson",use="complete.obs")))
      }
    }

    # compute nominal p-value

    combined.stats.rank<-rank(c(dist.to.template[i],perm.dist.vector))
    nominal.p[i]<-combined.stats.rank[1]/length(combined.stats.rank)

    # histgram of combined null distributions

    if (histgram.null.dist=="TRUE" & capabilities("png")==TRUE){
      png(paste("resampled_",dist.selection,"_dist_histgram_",sample.names[i],".png",sep=""))
      hist(c(dist.to.template[i],perm.dist.vector),br=hist.br,main=paste(sample.names[i],", # resampling: ",nresmpl,sep=""))
      dev.off()
    }

  } # main sample loop END

  # MCTRUE correction

  BH.FDR<-nominal.p*num.samples/rank(nominal.p)
  Bonferroni.p<-nominal.p*num.samples

  BH.FDR[BH.FDR>1]<-1
  Bonferroni.p[Bonferroni.p>1]<-1

  ### output ###

  # prediction results
  #return(dist.to.template)
  dist.to.cls1.rank <- rank(dist.to.cls1)
    predict.label2 <- vector()
  predict.label2[predict.label == 1] <- "CRIS-A"
  predict.label2[predict.label == 2] <- "CRIS-B"
  predict.label2[predict.label == 3] <- "CRIS-C"
  predict.label2[predict.label == 4] <- "CRIS-D"
  predict.label2[predict.label == 5] <- "CRIS-E"
  pred.summary <- cbind(sample.names,predict.label2,dist.to.template,dist.to.cls1.rank,
              nominal.p,BH.FDR,Bonferroni.p)
  #return(pred.summary)
  write.table(pred.summary,paste(output.name,"_prediction_result.xls",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE)

  # extracted features

  if (temp.nn.wt=="TRUE" & num.cls==2){
    write.table(features.extract[,2:5],paste(output.name,"_features.xls",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE)
  }
  if (temp.nn.wt=="FALSE"){
    write.table(features.extract[,2:4],paste(output.name,"_features.xls",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE)
  }

  # sorted exp dataset for heatmap (row normalized)

  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]

  sorted.dataset<-t(to.dataset.out)
  heatmap.dataset<-as.matrix(sorted.dataset)
#  if (.Platform$OS.type == "windows") {
#    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
#  }

  # sorted exp dataset for spreadsheets (not normalized)

  exp.dataset.gannot<-cbind(ProbeID,exp.dataset)
  exp.dataset.extract<-merge(features.extract,exp.dataset.gannot,sort=FALSE) # redefine exp.dataset.extract

  order.extract<-order(exp.dataset.extract[,2])   # sort genes
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  if (temp.nn.wt=="TRUE" & num.cls==2){
    exp.dataset.extract<-exp.dataset.extract[-c(1:5)]
  }
  if (temp.nn.wt=="FALSE"){
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
  }

  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]

  sorted.dataset<-t(to.dataset.out)
  sorted.dataset<-cbind(features.extract[,2:3],sorted.dataset)

  sorted.dataset.header<-c("ProbeID","GeneName",sample.names[order(t.dataset[,1])])
  sorted.dataset<-t(cbind(sorted.dataset.header,t(sorted.dataset)))

  if (.Platform$OS.type == "windows") {
    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
  }

  # output for GenePattern

  if (GenePattern.output=="TRUE"){

    # exp data
    write.table("#1.2",paste(output.name,"_sorted.dataset.gct",sep="")
                ,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(paste(num.features.extract,num.samples,sep="\t"),paste(output.name,"_sorted.dataset.gct",sep="")
                ,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
    write.table(sorted.dataset,paste(output.name,"_sorted.dataset.gct",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)

    # cls ffiles (unsorted, sorted)
    cls.out<-matrix(0,nrow=3,ncol=1)

    ## unsorted cls
    cls.out[1,]<-paste(num.samples," ",num.cls," 1",sep="")  # line 1
    cls.out[2,]<-paste("# ",paste(unique(predict.label),collapse=" ")) # line 2
    predict.label.out<-as.numeric(predict.label)-1                      # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")

    write.table(cls.out,paste(output.name,"_predicted_unsorted.cls",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

    ## sorted cls
    sorted.predict.label<-sort(as.numeric(predict.label))
    cls.out[2,]<-paste("# ",paste(unique(sorted.predict.label),collapse=" ")) # line 2
    predict.label.out<-sorted.predict.label-1    # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")

    write.table(cls.out,paste(output.name,"_predicted_sorted.cls",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

    # sample info
    sample.info<-cbind(sample.names,predict.label)
    write.table(sample.info,paste(output.name,"_sample_info.txt",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE)
  }

  # output for dChip


  # heatmap

  if (signature.heatmap=="TRUE" & capabilities("png")==TRUE){

    subclass.col.source <- c(rgb(255/255,70/255,0/255),rgb(180/255,0/255, 30/255),rgb(0/255,32/255,96/255),rgb(0/255,130/255, 50/255),rgb(0/255,200/255,170/255),"orange","lightblue","darkgreen")
    predict.col.vector <- unique(sort(predict.label))
    subclass.col <- subclass.col.source[predict.col.vector]

    heatmap.col <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")

    heatmap.dataset[heatmap.dataset>col.range] <- col.range
    heatmap.dataset[heatmap.dataset<-col.range] <- -col.range

    ncol.heat <- length(heatmap.dataset[1,])
    nrow.heat <- length(heatmap.dataset[,1])

    heatmap.dataset <- apply(heatmap.dataset,2,rev)

    num.pred <- as.vector(table(predict.label))
    num.pred.gene <- as.vector(table(features.extract[,4]))

    increment.sample <- cumsum(num.pred)
    increment.sample <- c(0,increment.sample)

    increment.gene <- cumsum(num.pred.gene)
    increment.gene <- c(0,increment.gene)

    png(paste(output.name,"_heatmp.png",sep=""))
      image(1:ncol.heat,1:nrow.heat,t(heatmap.dataset),axes=FALSE,col=heatmap.col,zlim=c(-col.range,col.range),xlim=c(-0.5,(ncol.heat+0.5+round(ncol.heat*0.05))),ylim=c(-0.5,(nrow.heat+0.5+round(nrow.heat*0.08))),xlab=NA,ylab=NA)

      for (c in 1:num.cls){                          # gene bar
        rect((ncol.heat+1),0.5,(ncol.heat+0.5+round(ncol.heat*0.05)),(nrow.heat+0.5-increment.gene[c]),col=subclass.col[c],xpd=TRUE,border=FALSE)
      }
      for (c in 1:length(num.pred)){               # sample bar
        rect((0.5+increment.sample[c]),(nrow.heat+2),(ncol.heat+0.5),(nrow.heat+0.5+round(nrow.heat*0.08)),col=subclass.col[c],xpd=TRUE,border=FALSE)
      }
    dev.off()

    # heatmap legend

    if (heatmap.legend=="TRUE"){
      png(paste(output.name,"_heatmap_legend.png",sep=""))
        par(plt=c(.1,.9,.45,.5))
        a=matrix(seq(1:15),ncol=1)
        image(a,col=heatmap.col,xlim=c(0,1),axes=FALSE,yaxt="n")
        box()
      dev.off()
    }

    # FDR sample bar

    if (FDR.sample.bar!="NA" & num.cls==2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- 3
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""))
        par(plt=c(.1,.9,.45,.5))
        a=matrix(fdr.bar.vector,ncol=1)
        image(a,col=c("red","blue","gray"),axes=FALSE,yaxt="n")
      dev.off()
    }

    if (FDR.sample.bar!="NA" & num.cls>2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- (num.cls+1)
      uni.fdr.bar.vector <- sort(unique(fdr.bar.vector))
      if (length(uni.fdr.bar.vector)>1){
        n.sig.cls <- length(uni.fdr.bar.vector)-1
        uni.fdr.bar.vector <- uni.fdr.bar.vector[1:n.sig.cls]
      }else{
        uni.fdr.bar.vector <- NULL
      }

      sig.subclass.col <- c(subclass.col.source[1:num.cls],"gray")
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""))
        par(plt=c(.1,.9,.45,.5))
        a=matrix(fdr.bar.vector,ncol=1)
        image(a,col=sig.subclass.col,axes=FALSE,yaxt="n")
      dev.off()
    }

  }

  # plot FDR

  if (plot.FDR=="TRUE" & capabilities("png")==TRUE){
    png(paste(output.name,"_FDR.png",sep=""))
      par(plt=c(0.1,0.95,0.4,0.6),las=2)
      plot(BH.FDR[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main="BH-FDR")
      box(lwd=2)
    dev.off()
  }


  # plot distance to template


}  # END main

