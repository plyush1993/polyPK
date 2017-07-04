#' @title Plot the correlation diagram of two datasets.
#'
#' @description A function to calculate the correlation coefficients and plot the correlation diagram (8 types) of two input datasets.
#'
#' @param dataset1 The first dataset (data frame with required format). The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected. The second row of the data frame should be the group information. The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). Please see the demo data for detailed format.This variable maybe the results of GetAbso/GetEndo/GetSecdAbso.
#' @param dataset2 The second dataset (data frame with required format). The form of dataset2 is the same as the form of dataset1.This variable maybe the results of GetAbso/GetEndo/GetSecdAbso
#' @param cor.method A character string indicating which correlation analysis ("pearson", "kendall", or "spearman") is to be used. Default: "pearson".
#' @param filepath A character string indicating the path where the results may be saved in.
#' @param fig.form The form of the correlation diagram.  figure.fig.form=c("heatmap","bubble","ordered.bubble","chord","square","ord.square","pie","ord.pie"). Default: "heatmap".
#' @param design A study design dataset(data frame with required format).Use data(StudyDesign) to see the detailed form.
#'
#' @return NULL
#'
#' @examples CorrPlot(B,C,"pearson",filepath=getwd(),fig.form="heatmap",design)
#'
#' @export
CorrPlot<-function(dataset1,dataset2,cor.method="pearson",filepath,fig.form="heatmap",design=FALSE){
  #library(corrplot)
  #library(circlize)
  # library(lattice)
  # library(mixOmics)
  ####---------correlation of all data-------------------
  tp1<-as.matrix(dataset1[3,-c(1:4)])
  tp2<-as.matrix(dataset2[3,-c(1:4)])
  tpall<-cbind(tp1,tp2)
  tpall<-as.numeric(tpall)
  timepoin<-tpall[!duplicated(tpall)]
  dataset1<-dataset1[-3,]
  dataset2<-dataset2[-3,]
  cor.data1<-as.data.frame(dataset1[-1,])
  cor.data2<-as.data.frame(dataset2[-1,])# two input data
  id1<-cor.data1[-c(1:2),1]
  id1<-as.matrix(id1)#id of the first data
  id2<-cor.data2[-c(1:2),1]
  id2<-as.matrix(id2)#id of the second data
  m1<-dim(cor.data1)[1]-1
  m2<-dim(cor.data2)[1]-1#the number of the metabolites

  mm1<-m1+1
  mm2<-m2+1
  all.m<-m1+m2+1
  cor.m<-matrix(0,nrow = all.m,ncol = all.m)# prepare a zero matrix with id in two data(both row and column)
  for(i in 1:m1){
    cor.m[i+1,1]<-id1[i]
    cor.m[1,i+1]<-id1[i]
  }

  for(i in 1:m2){
    cor.m[1,mm1+i]<-id2[i]
    cor.m[mm1+i,1]<-id2[i]
  }
  # set the matrix with IDs
  ######calculate the p value and r value (data1 correlation with data2,data1 correlation with data1,data2 correlation with data2)
  cor.m.r<-cor.m
  for(i in 2:mm1){#m1+1
    for(j in 2:mm2){#m2+1
      c.isoque<-cor.data1[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data2[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j+m1,i]<-cor.m[i,j+m1]<-cor.b1.c1$p.value
      cor.m.r[j+m1,i]<-cor.m.r[i,j+m1]<-cor.b1.c1$estimate

    }
  }

  for(i in 2:mm1){
    for(j in 2:mm1){
      c.isoque<-cor.data1[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data1[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j,i]<-cor.b1.c1$p.value
      cor.m.r[j,i]<-cor.b1.c1$estimate
    }
  }
  for(i in 2:mm2){
    for(j in 2:mm2){
      c.isoque<-cor.data2[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data2[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j+m1,i+m1]<-cor.b1.c1$p.value
      cor.m.r[j+m1,i+m1]<-cor.b1.c1$estimate
    }
  }
  cor.m[1,1]<-"P value"
  cor.m.r[1,1]<-"r value"
  cor.r<-cor.m.r
  cor.p.value<-cor.m
  # r p finished
  diroutall = paste(filepath, "/CorrelationResults", "/", sep = "")
  dir.create(diroutall)
  dirout = paste(diroutall, "/CorrelationResults(all)", "/", sep = "")
  dir.create(dirout)


  ###
  #r value(data without name)
  cor.rr<-cor.m.r[-1,-1]
  cor.rr<-as.matrix(cor.rr)
  cor.rr<-as.numeric(cor.rr)
  dim(cor.rr)<-dim(cor.m.r)-1
  #p value(data without name)
  cor.p<-cor.p.value[-1,-1]
  cor.p<-as.matrix(cor.p)
  cor.p<-as.numeric(cor.p)
  dim(cor.p)<-dim(cor.p.value)-1
  # add names
  id.m<-matrix(0,nrow =1,ncol = dim(cor.rr)[1])
  for(i in 1: m1){
    w<-as.character(cor.data1[i+1,2])
    id.m[i]<-w
  }
  for(i in 1: m2){
    w<-as.character(cor.data2[i+1,2])
    id.m[i+m1]<-w
  }
  colnames(cor.rr)<-id.m
 # rownames(cor.rr)<-id.m
  rownames(cor.rr)<-c(1:dim(cor.rr)[1])
  colnames(cor.p)<-id.m
  #rownames(cor.p)<-id.m
  rownames(cor.p)<-c(1:dim(cor.p)[1])
  pwdxdef = paste(dirout, "r-value.xlsx", sep = "")
  xlsx::write.xlsx(cor.rr, pwdxdef)

  pwdx = paste(dirout, "p-value.xlsx", sep = "")
  xlsx::write.xlsx(cor.p, pwdx)
  #p<0.05 r<-0
  rownames(cor.rr)<-id.m
  if(fig.form=="heatmap"){
    heatfigpath = paste(dirout, "correlation-matrix-HeatMap.pdf", sep = "")
    pdf(heatfigpath)
    gplots::heatmap.2(cor.rr,scale = "none",dendrogram ="none",cexRow = 0.5 ,cexCol = 0.5 ,main="HeatMap of correlation matrix",trace="none" )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }

  if(fig.form=="bubble"){
    rownames(cor.rr)<-id.m
    bfigpath = paste(dirout, "correlation-matrix-BubbleDiagram.pdf", sep = "")
    pdf(bfigpath)

    corrplot::corrplot(cor.rr,method = "circle",tl.cex =  0.5)
    title(main = list("BubbleDiagram of correlation matrix",cex = 1.5),outer = F)
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

    # dev.off()
  }
  if(fig.form=="ordered.bubble"){
    rownames(cor.rr)<-id.m
    obfigpath = paste(dirout, "correlation-matrix-ordered-BubbleDiagram.pdf", sep = "")
    pdf(obfigpath)
    corrplot::corrplot(cor.rr,order="AOE",method = "circle",tl.cex =  0.5)
    title(main = list("ordered BubbleDiagram of correlation matrix",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="chord"){
    cor.r.cir<-cor.rr
    cor.r.cir[cor.p>=0.05]<-0
    cir.m<-corrplot::corrplot(cor.r.cir)
    cfigpath = paste(dirout, "correlation-matrix-ChordDiagram.pdf", sep = "")
    pdf(cfigpath)
    circlize::chordDiagram(cir.m,symmetric = TRUE)
    title(main = "Chord Diagram of correlation matrix",outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ord.square"){

    ofigpath = paste(dirout, "ordered-correlogram-square.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "square",order = "AOE",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations(ordered)",cex = 1.5),outer = F )
    #dev.off()
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ord.pie"){

    ofigpath = paste(dirout, "ordered-correlogram-pie.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "pie",order = "AOE",tl.cex =  0.5)
    title(main =list ("Correlogram of data intercorrelations(ordered)",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="square"){

    ofigpath = paste(dirout, "correlogram-square.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "square",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="pie"){

    ofigpath = paste(dirout, "correlogram-pie.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "pie",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  dev.off()


  ####---------divided by genders/dataset1-------------
  allnumdata<-dataset1[,-c(1:4)]
  headall<-dataset1[,c(1:4)]
  samplenum<-dim(allnumdata)[2]
  maledata<-data.frame(0)
  femaledata<-data.frame(0)
  for(i in 1:samplenum){
    if(allnumdata[1,i]==1){
      male<-allnumdata[,i]
      maledata<-cbind(maledata,male)
    }
    if(allnumdata[1,i]==2){
      female<-allnumdata[,i]
      femaledata<-cbind(femaledata,female)
    }
  }
  maledata<-cbind(headall,maledata)
  femaledata<-cbind(headall,femaledata)
  maledata1<-maledata[,-5]
  femaledata1<-femaledata[,-5]
  ####---------divided by genders/dataset2-------------
  allnumdata<-dataset2[,-c(1:4)]
  headall<-dataset2[,c(1:4)]
  samplenum<-dim(allnumdata)[2]
  maledata<-data.frame(0)
  femaledata<-data.frame(0)
  for(i in 1:samplenum){
    if(allnumdata[1,i]==1){
      male<-allnumdata[,i]
      maledata<-cbind(maledata,male)
    }
    if(allnumdata[1,i]==2){
      female<-allnumdata[,i]
      femaledata<-cbind(femaledata,female)
    }
  }
  maledata<-cbind(headall,maledata)
  femaledata<-cbind(headall,femaledata)
  maledata2<-maledata[,-5]
  femaledata2<-femaledata[,-5]

  ####---------correlation of maledata-------------------
  cor.data1<-as.data.frame(maledata1[-1,])
  cor.data2<-as.data.frame(maledata2[-1,])# two input data
  id1<-cor.data1[-c(1:2),1]
  id1<-as.matrix(id1)#id of the first data
  id2<-cor.data2[-c(1:2),1]
  id2<-as.matrix(id2)#id of the second data
  m1<-dim(cor.data1)[1]-1
  m2<-dim(cor.data2)[1]-1#the number of the metabolites

  mm1<-m1+1
  mm2<-m2+1
  all.m<-m1+m2+1
  cor.m<-matrix(0,nrow = all.m,ncol = all.m)# prepare a zero matrix with id in two data(both row and column)
  for(i in 1:m1){
    cor.m[i+1,1]<-id1[i]
    cor.m[1,i+1]<-id1[i]
  }

  for(i in 1:m2){
    cor.m[1,mm1+i]<-id2[i]
    cor.m[mm1+i,1]<-id2[i]
  }
  # set the matrix with IDs
  ######calculate the p value and r value (data1 correlation with data2,data1 correlation with data1,data2 correlation with data2)
  cor.m.r<-cor.m
  for(i in 2:mm1){#m1+1
    for(j in 2:mm2){#m2+1
      c.isoque<-cor.data1[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data2[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j+m1,i]<-cor.m[i,j+m1]<-cor.b1.c1$p.value
      cor.m.r[j+m1,i]<-cor.m.r[i,j+m1]<-cor.b1.c1$estimate

    }
  }

  for(i in 2:mm1){
    for(j in 2:mm1){
      c.isoque<-cor.data1[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data1[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j,i]<-cor.b1.c1$p.value
      cor.m.r[j,i]<-cor.b1.c1$estimate
    }
  }
  for(i in 2:mm2){
    for(j in 2:mm2){
      c.isoque<-cor.data2[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data2[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j+m1,i+m1]<-cor.b1.c1$p.value
      cor.m.r[j+m1,i+m1]<-cor.b1.c1$estimate
    }
  }
  cor.m[1,1]<-"P value"
  cor.m.r[1,1]<-"r value"
  cor.r<-cor.m.r
  cor.p.value<-cor.m
  # r p finished
  dirout = paste(diroutall, "/CorrelationResults(male)", "/", sep = "")
  dir.create(dirout)
  ###
  #r value(data without name)
  cor.rr<-cor.m.r[-1,-1]
  cor.rr<-as.matrix(cor.rr)
  cor.rr<-as.numeric(cor.rr)
  dim(cor.rr)<-dim(cor.m.r)-1
  #p value(data without name)
  cor.p<-cor.p.value[-1,-1]
  cor.p<-as.matrix(cor.p)
  cor.p<-as.numeric(cor.p)
  dim(cor.p)<-dim(cor.p.value)-1
  # add names
  id.m<-matrix(0,nrow =1,ncol = dim(cor.rr)[1])
  for(i in 1: m1){
    w<-as.character(cor.data1[i+1,2])
    id.m[i]<-w
  }
  for(i in 1: m2){
    w<-as.character(cor.data2[i+1,2])
    id.m[i+m1]<-w
  }
  colnames(cor.rr)<-id.m
  #rownames(cor.rr)<-id.m
  rownames(cor.rr)<-c(1:dim(cor.rr)[1])
  colnames(cor.p)<-id.m
  #rownames(cor.p)<-id.m
  rownames(cor.p)<-c(1:dim(cor.p)[1])
  pwdxdef = paste(dirout, "r-value.xlsx", sep = "")
  xlsx::write.xlsx(cor.rr, pwdxdef)

  pwdx = paste(dirout, "p-value.xlsx", sep = "")
  xlsx::write.xlsx(cor.p, pwdx)
  #p<0.05 r<-0
  rownames(cor.rr)<-id.m
  if(fig.form=="heatmap"){
    heatfigpath = paste(dirout, "correlation-matrix-HeatMap.pdf", sep = "")
    pdf(heatfigpath)
    gplots::heatmap.2(cor.rr,scale = "none",dendrogram ="none",cexRow = 0.5 ,cexCol = 0.5 ,main="HeatMap of correlation matrix",trace="none" )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }

  if(fig.form=="bubble"){
    rownames(cor.rr)<-id.m
    bfigpath = paste(dirout, "correlation-matrix-BubbleDiagram.pdf", sep = "")
    pdf(bfigpath)

    corrplot::corrplot(cor.rr,method = "circle",tl.cex =  0.5)
    title(main = list("BubbleDiagram of correlation matrix",cex = 1.5),outer = F)
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ordered.bubble"){
    rownames(cor.rr)<-id.m
    obfigpath = paste(dirout, "correlation-matrix-ordered-BubbleDiagram.pdf", sep = "")
    pdf(obfigpath)
    corrplot::corrplot(cor.rr,order="AOE",method = "circle",tl.cex =  0.5)
    title(main = list("ordered BubbleDiagram of correlation matrix",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="chord"){
    cor.r.cir<-cor.rr
    cor.r.cir[cor.p>=0.05]<-0
    cir.m<-corrplot::corrplot(cor.r.cir)
    cfigpath = paste(dirout, "correlation-matrix-ChordDiagram.pdf", sep = "")
    pdf(cfigpath)
    circlize::chordDiagram(cir.m,symmetric = TRUE)
    title(main = "Chord Diagram of correlation matrix",outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ord.square"){

    ofigpath = paste(dirout, "ordered-correlogram-square.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "square",order = "AOE",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations(ordered)",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ord.pie"){

    ofigpath = paste(dirout, "ordered-correlogram-pie.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "pie",order = "AOE",tl.cex =  0.5)
    title(main =list ("Correlogram of data intercorrelations(ordered)",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="square"){

    ofigpath = paste(dirout, "correlogram-square.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "square",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="pie"){

    ofigpath = paste(dirout, "correlogram-pie.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "pie",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  dev.off()





  ####---------correlation of femaledata-------------------
  cor.data1<-as.data.frame(femaledata1[-1,])
  cor.data2<-as.data.frame(femaledata2[-1,])# two input data
  id1<-cor.data1[-c(1:2),1]
  id1<-as.matrix(id1)#id of the first data
  id2<-cor.data2[-c(1:2),1]
  id2<-as.matrix(id2)#id of the second data
  m1<-dim(cor.data1)[1]-1
  m2<-dim(cor.data2)[1]-1#the number of the metabolites

  mm1<-m1+1
  mm2<-m2+1
  all.m<-m1+m2+1
  cor.m<-matrix(0,nrow = all.m,ncol = all.m)# prepare a zero matrix with id in two data(both row and column)
  for(i in 1:m1){
    cor.m[i+1,1]<-id1[i]
    cor.m[1,i+1]<-id1[i]
  }

  for(i in 1:m2){
    cor.m[1,mm1+i]<-id2[i]
    cor.m[mm1+i,1]<-id2[i]
  }
  # set the matrix with IDs
  ######calculate the p value and r value (data1 correlation with data2,data1 correlation with data1,data2 correlation with data2)
  cor.m.r<-cor.m
  for(i in 2:mm1){#m1+1
    for(j in 2:mm2){#m2+1
      c.isoque<-cor.data1[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data2[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j+m1,i]<-cor.m[i,j+m1]<-cor.b1.c1$p.value
      cor.m.r[j+m1,i]<-cor.m.r[i,j+m1]<-cor.b1.c1$estimate

    }
  }

  for(i in 2:mm1){
    for(j in 2:mm1){
      c.isoque<-cor.data1[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data1[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j,i]<-cor.b1.c1$p.value
      cor.m.r[j,i]<-cor.b1.c1$estimate
    }
  }
  for(i in 2:mm2){
    for(j in 2:mm2){
      c.isoque<-cor.data2[i,-c(1:4)]
      c.isoque<-as.matrix(c.isoque)
      c.isoque<-as.numeric(c.isoque)
      c.isoque[is.na(c.isoque)]<-0
      b.octadeca<-cor.data2[j,-c(1:4)]
      b.octadeca<-as.matrix(b.octadeca)
      b.octadeca<-as.numeric(b.octadeca)
      b.octadeca[is.na(b.octadeca)]<-0
      cor.b1.c1<-cor.test(c.isoque,b.octadeca,method=cor.method)
      cor.m[j+m1,i+m1]<-cor.b1.c1$p.value
      cor.m.r[j+m1,i+m1]<-cor.b1.c1$estimate
    }
  }
  cor.m[1,1]<-"P value"
  cor.m.r[1,1]<-"r value"
  cor.r<-cor.m.r
  cor.p.value<-cor.m
  # r p finished
  dirout = paste(diroutall, "/CorrelationResults(female)", "/", sep = "")
  dir.create(dirout)
  ###
  #r value(data without name)
  cor.rr<-cor.m.r[-1,-1]
  cor.rr<-as.matrix(cor.rr)
  cor.rr<-as.numeric(cor.rr)
  dim(cor.rr)<-dim(cor.m.r)-1
  #p value(data without name)
  cor.p<-cor.p.value[-1,-1]
  cor.p<-as.matrix(cor.p)
  cor.p<-as.numeric(cor.p)
  dim(cor.p)<-dim(cor.p.value)-1
  # add names
  id.m<-matrix(0,nrow =1,ncol = dim(cor.rr)[1])
  for(i in 1: m1){
    w<-as.character(cor.data1[i+1,2])
    id.m[i]<-w
  }
  for(i in 1: m2){
    w<-as.character(cor.data2[i+1,2])
    id.m[i+m1]<-w
  }
  colnames(cor.rr)<-id.m
  #rownames(cor.rr)<-id.m
  rownames(cor.rr)<-c(1:dim(cor.rr)[1])
  colnames(cor.p)<-id.m
  #rownames(cor.p)<-id.m
  rownames(cor.p)<-c(1:dim(cor.p)[1])
  pwdxdef = paste(dirout, "r-value.xlsx", sep = "")
  xlsx::write.xlsx(cor.rr, pwdxdef)

  pwdx = paste(dirout, "p-value.xlsx", sep = "")
  xlsx::write.xlsx(cor.p, pwdx)
  #p<0.05 r<-0
  rownames(cor.rr)<-id.m
  if(fig.form=="heatmap"){
    heatfigpath = paste(dirout, "correlation-matrix-HeatMap.pdf", sep = "")
    pdf(heatfigpath)
    gplots::heatmap.2(cor.rr,scale = "none",dendrogram ="none",cexRow = 0.5 ,cexCol = 0.5 ,main="HeatMap of correlation matrix",trace="none" )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }

  if(fig.form=="bubble"){
    rownames(cor.rr)<-id.m
    bfigpath = paste(dirout, "correlation-matrix-BubbleDiagram.pdf", sep = "")
    pdf(bfigpath)

    corrplot::corrplot(cor.rr,method = "circle",tl.cex =  0.5)
    title(main = list("BubbleDiagram of correlation matrix",cex = 1.5),outer = F)
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ordered.bubble"){
    rownames(cor.rr)<-id.m
    obfigpath = paste(dirout, "correlation-matrix-ordered-BubbleDiagram.pdf", sep = "")
    pdf(obfigpath)
    corrplot::corrplot(cor.rr,order="AOE",method = "circle",tl.cex =  0.5)
    title(main = list("ordered BubbleDiagram of correlation matrix",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="chord"){
    cor.r.cir<-cor.rr
    cor.r.cir[cor.p>=0.05]<-0
    cir.m<-corrplot::corrplot(cor.r.cir)
    cfigpath = paste(dirout, "correlation-matrix-ChordDiagram.pdf", sep = "")
    pdf(cfigpath)
    circlize::chordDiagram(cir.m,symmetric = TRUE)
    title(main = "Chord Diagram of correlation matrix",outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ord.square"){

    ofigpath = paste(dirout, "ordered-correlogram-square.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "square",order = "AOE",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations(ordered)",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="ord.pie"){

    ofigpath = paste(dirout, "ordered-correlogram-pie.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "pie",order = "AOE",tl.cex =  0.5)
    title(main =list ("Correlogram of data intercorrelations(ordered)",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="square"){

    ofigpath = paste(dirout, "correlogram-square.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "square",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  if(fig.form=="pie"){

    ofigpath = paste(dirout, "correlogram-pie.pdf", sep = "")
    pdf(ofigpath)
    corrplot::corrplot(cor.rr,method = "pie",tl.cex =  0.5)
    title(main = list("Correlogram of data intercorrelations",cex = 1.5),outer = F )
    if(is.data.frame(design)){

      stusd<-design
      ####meal
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of meals(h)"){
          mealtime<-stusd[i,]
        }
      }
      mealtime<-mealtime[,-c(1:2)]
      mealp<-as.matrix(0)
      for(i in 1:length(mealtime)){
        if(!is.na(mealtime[i])){
          mealp<-cbind(mealp,mealtime[i])
        }
      }
      mealp<-mealp[,-1]
      mealp<-as.numeric(as.matrix(mealp))
      mealtp<-data.frame()
      for(lm in 1:length(mealp)){
        low<-mealp[lm]
        high<-mealp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            mealtp<-rbind(mealtp,dataxx)
        }

      }
      timecors<-mealtp[1,1]
      for(mp in 2:dim(mealtp)[1]){
        timecors<-paste(timecors,mealtp[mp,1])
      }
      ####------sleep
      for(i in 1:dim(stusd)[1]){
        if(stusd[i,1]=="Number of sleep(h)"){
          sleeptime<-stusd[i,]
        }
      }
      sleeptime<-sleeptime[,-c(1:2)]
      sleepp<-as.matrix(0)
      for(i in 1:length(sleeptime)){
        if(!is.na(sleeptime[i])){
          sleepp<-cbind(sleepp,sleeptime[i])
        }
      }
      sleepp<-sleepp[,-1]### the time of sleeps
      sleepp<-as.numeric(as.matrix(sleepp))
      sleeptp<-data.frame()
      for(lm in 1:length(sleepp)){
        low<-sleepp[lm]
        high<-sleepp[lm]+2
        for(j in 1:length(timepoin)){
          dataxx<- as.matrix(timepoin[j])
          dataxx<-as.numeric(dataxx)
          if(low<=dataxx&dataxx<=high)
            sleeptp<-rbind(sleeptp,dataxx)
        }
      }

      for(mps in 1:dim(sleeptp)[1]){
        timecors<-paste(timecors,sleeptp[mps,1])
      }

      title(sub=paste("Results of time windows",timecors,"may be impacted by diet and circadian"))
    }
    if(!is.data.frame(design)){
      title(sub="Results of time windows may be impacted by diet and circadian")}

  }
  dev.off()





  }
