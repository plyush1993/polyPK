#' @title Plot the PCA or PLSDA score figures and trajectories on input data.
#'
#' @description A function to plot the PCA or PLSDA figures of input data.
#'
#' @param scat.data The data under analysis (data frame with required format). The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected. The second row of the data frame should be the group information. The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). Please see the demo data for detailed format.
#' @param scform  The form of scat plot. scform=c ("PCA","PLSDA"). Default:"PCA"
#' @param num.of.cp The number of components to decompose. Default:2.
#' @param fold Integer: number of random permutations of response labels to estimate R2Y and Q2Y significance by permutation testing [default is 100 for single response models (without train/test partition), and 0 otherwise]
#' @param filepath A character string indicating the path where the results may be saved in.
#' @param design A study design dataset(data frame with required format).Use data(StudyDesign) to see the detailed form.
#'
#' @return NULL
#'
#' @examples ScatPlot(scat.data=A,scform="PCA",num.of.cp=3,filepath,design=design)
#'
#' @export
ScatPlot<-function(scat.data,scform="PCA",num.of.cp=2,fold=100,filepath,design=FALSE){
#  library(mixOmics)
  ####---------------------------preprocess the input data----------------------
  timep<-scat.data[3,]
  tpall<-timep<-as.matrix(timep[-c(1:4)])
  tpall<-as.numeric(tpall)
  timepoin<-tpall[!duplicated(tpall)]
  scat.data<-as.data.frame(scat.data[-3,])
  group.num<-as.matrix(scat.data[2,dim(scat.data)[2]])
  group.num<-as.numeric(group.num)     ###num 9
  group.begin<-as.matrix(scat.data[2,5])
  group.begin<-as.numeric(group.begin)  ###num 0
  data<-as.matrix(scat.data)
  input.matrix<-scat.data
  lg<-group.begin:group.num
  time.points<-length(lg)            #########10
  list.sl<-list()
  all.time.num<-dim(input.matrix)[2]-4###100
  tg<-matrix(0,nrow=1,ncol=length(lg))
  for(g in 1:length(lg)){
    for(i in 1:all.time.num){
      if(input.matrix[2,i+4]==lg[g])
        tg[g]<-tg[g]+1
    }
  }
  data<-data[-1,]
  dim(timep)<-c(tg[1],length(group.begin:group.num))
  legc<-legct<-timep[1,]

  ####--------------------------------PCA-------------------------------
  if(scform=="PCA"){

    data.ma<-data[-1,-c(1:4)]
    data.ma<-t(data.ma)
    data.ma<-as.matrix(data.ma)
    data.ma1<-as.numeric(data.ma)
    dim(data.ma1)<-dim(data.ma)

    for(i in 1:dim(data.ma1)[1]){
       if(sum(is.na(data.ma1[i,]))==dim(data.ma1)[2]){
      data.ma1[i,]<-0
       }
    }
    data.ma1<-scale(data.ma1,center = F, scale = TRUE)
    #data.ma1[is.na(data.ma1)]<-0
   # ma.pca<-mixOmics::pca(data.ma1,center = TRUE)
    wellpca<-pcaMethods::pca(data.ma1,method = "nipals",nPcs = num.of.cp)
    #r2x_pca<-R2VX(wellpca, direction = c("variables"), data = completeObs(wellpca), pcs = nP(wellpca))
    #r2y_pca<-R2VX(wellpca, direction = c("complete"), data = completeObs(wellpca), pcs = nP(wellpca))
    r2_pca<-wellpca@R2
    r2_pca<-round(r2_pca,2)
    q2_pca<-pcaMethods::Q2(wellpca, originalData = completeObs(wellpca),
       variables = 1:pcaMethods::nVar(wellpca))
    q2_pca<-round(q2_pca,2)
    ma.score<-wellpca@scores
    ma.loading<-wellpca@loadings
    diroutall = paste(filepath, "/PCAresults", "/", sep = "")
    dir.create(diroutall)
    dirout = paste(diroutall, "/PCA(all)", "/", sep = "")
    dir.create(dirout)
   # pca_path = paste(dirout, "Variance explained for PCA.pdf", sep = "")
   # sacurine.pca <- ropls::opls(data.ma1,printL=FALSE,plotL=FALSE)
   # typeC_pca <-"overview"
   # ropls::plot(sacurine.pca, typeVc = typeC_pca,file.pdfC=pca_path)

    pwdxdef = paste(dirout, "PCA-score.xlsx", sep = "")
    xlsx::write.xlsx(ma.score, pwdxdef)
    pwdxd = paste(dirout, "PCA-loading.xlsx", sep = "")
    xlsx::write.xlsx(ma.loading, pwdxd)

  #  q2figpath = paste(dirout, "PCA-scorePlot.pdf", sep = "")
  #  barplot(r2cum, main="Krzanowski CV", xlab="Number of PCs", ylab=expression(Q^2))

    tutticolors = matrix(c(4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink",
                           4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink"),
                         ncol = 1)
    col = c()
    for (i in 1:time.points) {
      col = c(col, tutticolors[i, ])
    }
    allcol<-c()
    for(i in 1:length(tg)){
      allcol<-c(allcol,rep(col[i],tg[i]))
    }



    #legc<-paste("time",c(group.begin:group.num))
    pcafigpath = paste(dirout, "PCA-scorePlot.pdf", sep = "")
    par(mai=c(1,0.5,0.5,0.3),cex.lab=1,lab=c(3,3,5),mgp=c(1,1,0))
    pdf(pcafigpath)

    plot(x=ma.score[,1],y=ma.score[,2],pch=16 ,xaxt="n", yaxt="n",col=allcol,main = "PCA score figure",xlab = "PC1",ylab = "PC2")
    xy<-par("usr")
    legend("bottomright",c(paste("R2Y(PC1)=",r2_pca[1]),paste("R2Y(PC2)=",r2_pca[2]),paste("Q2Y(PC1)=",q2_pca[1]),paste("Q2Y(PC2)=",q2_pca[2])), cex=0.7);#new add
    legend(x=xy[1]+xinch(4.5),y=xy[3]-yinch(0.1),legend=legc,xpd=TRUE,fill = col, col=col,cex = 0.6,ncol = 3,x.intersp=0.8,xjust=0,yjust=1)
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
    dev.off()
    #pca2
    pca2figpath = paste(dirout, "PCA-scorePlot(track).pdf", sep = "")
    pdf(pca2figpath)

    groupsh<-input.matrix[2,-c(1:4)]
    groupsh<-as.matrix(t(groupsh))
    groupsh<-as.data.frame(groupsh)

    g.ma.score<-cbind(groupsh,ma.score)
    PC1.mean<-tapply(g.ma.score[,2],g.ma.score[,1],mean)
    PC2.mean<-tapply(g.ma.score[,3],g.ma.score[,1],mean)
    plot(x=PC1.mean,y=PC2.mean,type = "b",pch=21,bg = "red")

    title("PCA Score(track)")
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
    text(PC1.mean, PC2.mean,legct,pos = 1)
    dev.off()
    ###

  }
  ####--------------------------------PLSDA-----------------------------
  if(scform=="PLSDA"){

    data.ma<-data[-1,-c(1:4)]
    data.ma<-t(data.ma)
    data.ma<-as.matrix(data.ma)
    data.ma1<-as.numeric(data.ma)
    dim(data.ma1)<-dim(data.ma)
    data.ma2<-scale(data.ma1,center = F, scale = TRUE)
    #data.ma1<-data.ma1[,order(data.ma1[dim(data.ma1)[1],])]
    y.group<-data[1,-c(1:4)]

    ma.plsda<-mixOmics::plsda(data.ma2,y.group,mode = "classic",ncomp = num.of.cp)
######add
    perf_plsda<-mixOmics::perf(ma.plsda,validation = "Mfold", folds = 5)
    overall_plsda<-perf_plsda$error.rate$overall
    #the Prediction error rate for each block of object$X and each dist
    y_group<-as.numeric(y.group)
    data.ma3<- data.ma2
    data.ma3[is.na(data.ma2)] <- 0
    #data.ma3<-as.numeric(data.ma3)
    ma.pls<-mixOmics::pls(data.ma3,y_group,mode = "classic",ncomp =  dim(data.ma3)[2])
    #X<-t(t(y_group))
    #X<-as.character(X)
  #  Y<-data.ma3
   # Y[is.na(Y)]<-0
  #  PT<-permKS(Y,X,alternative="two.sided",method="pclt")

    perf_pls<-mixOmics::perf(ma.pls,validation = "Mfold")
    r2<-round(perf_pls$R2,2)[c(1:2)]
    q2<-round(perf_pls$Q2.total,2)[c(1,3)]
   # perfigpath = paste(dirout, "PermutationPlot.pdf", sep = "")
   # pdf(perfigpath)
    data.ma4 = data.ma3
#    randtest_pls<-mdatools::randtest(data.ma4, y_group, ncomp = 4, center = T, scale = F, nperm = 1000,
    #         sig.level = 0.05, silent = TRUE, exclcols = NULL, exclrows = NULL)

 #   plotCorr(randtest_pls, 1)

    diroutall = paste(filepath, "/PLSDAresults", "/", sep = "")
    dir.create(diroutall)
    dirout = paste(diroutall, "/PLSDA(all)", "/", sep = "")
    dir.create(dirout)
    perfigpath = paste(dirout, "PermutationPlot.pdf", sep = "")
    sacurine.plsda <- ropls::opls(data.ma4, y_group,printL=FALSE,plotL=FALSE,permI=fold)
    typeC <-"permutation"
    ropls::plot(sacurine.plsda, typeVc = typeC,file.pdfC=perfigpath)
    pw_allover = paste(dirout, "ErrorRate.xlsx", sep = "")
    xlsx::write.xlsx(overall_plsda, pw_allover)

   # sacurine.pca <- ropls::opls(data.ma4)

    ####old
   # Q2_sam<-matrix(nrow=2,ncol=fold)
  #  R2_sam<-matrix(nrow=2,ncol=fold)
  #  for(i in 1:fold){
  #    X_sam<-sample(X)
   #   ma.pls<-mixOmics::pls(data.ma3,X_sam,mode = "classic",ncomp = 2)
   #   perf_pls_sam<-perf(ma.pls,validation = "Mfold")
   #   r2_sam<-round(perf_pls_sam$R2,2)
    #  q2_sam<-round(perf_pls_sam$Q2,2)
    #  Q2_sam[]

   # }

    #####

    loading.plsda<-ma.plsda$loadings

    pwdxdef = paste(dirout, "PLSDA-loading.xlsx", sep = "")
    xlsx::write.xlsx(loading.plsda$X, pwdxdef)
    score.plsda<-ma.plsda$variates$X
    pwdxd = paste(dirout, "PLSDA-score.xlsx", sep = "")
    xlsx::write.xlsx(score.plsda, pwdxd)

    tutticolors = matrix(c(4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink",
                           4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink"),
                         ncol = 1)
    col = c()
    for (i in 1:time.points) {
      col = c(col, tutticolors[i, ])
    }
    allcol<-c()
    for(i in 1:length(tg)){
      allcol<-c(allcol,rep(col[i],tg[i]))
    }


    pcafigpath = paste(dirout, "PLSDA-scorePlot.pdf", sep = "")
    par(mai=c(1,0.5,0.5,0.3),cex.lab=1,lab=c(3,3,5),mgp=c(1,1,0))
    pdf(pcafigpath)

    plot(x=score.plsda[,1],y=score.plsda[,2],pch=16 ,xaxt="n", yaxt="n",col=allcol,main = "PLSDA score figure",xlab = "PC1",ylab = "PC2")
    xy<-par("usr")
###add
    legend("bottomright",c(paste("R2Y(PC1)=",r2[1]),paste("R2Y(PC2)=",r2[2]),paste("Q2Y(PC1)=",q2[1]),paste("Q2Y(PC2)=",q2[2])), cex=0.7);#new add
    legend(x=xy[1]+xinch(4.5),y=xy[3]-yinch(0.1),legend=legc,xpd=TRUE,fill = col, col=col,cex = 0.6,ncol = 3,x.intersp=0.8,xjust=0,yjust=1)
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
    dev.off()

    plsda2figpath = paste(dirout, "PLSDA-scorePlot(track).pdf", sep = "")
    pdf(plsda2figpath)

    groupsh<-input.matrix[2,-c(1:4)]
    groupsh<-as.matrix(t(groupsh))
    groupsh<-as.data.frame(groupsh)
    g.ma.score<-cbind(groupsh,score.plsda)
    PC1.mean<-tapply(g.ma.score[,2],g.ma.score[,1],mean)
    PC2.mean<-tapply(g.ma.score[,3],g.ma.score[,1],mean)
    plot(x=PC1.mean,y=PC2.mean,type = "b",pch=21,bg = "red")

    title("PLSDA Score(track)")
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
    text(PC1.mean, PC2.mean,legct,pos = 1)
    dev.off()
  }
  ####--------------------------------divide----------------------
  allnumdata<-scat.data[,-c(1:4)]
  headall<-scat.data[,c(1:4)]
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
  maledata<-maledata[,-5]
  femaledata<-femaledata[,-5]
  #####---------------------male----------------------
  scat.data<-maledata
  scat.data<-as.data.frame(scat.data)
  group.num<-as.matrix(scat.data[2,dim(scat.data)[2]])
  group.num<-as.numeric(group.num)     ###num 9
  group.begin<-as.matrix(scat.data[2,5])
  group.begin<-as.numeric(group.begin)  ###num 0
  data<-as.matrix(scat.data)
  input.matrix<-scat.data
  lg<-group.begin:group.num
  time.points<-length(lg)            #########10
  list.sl<-list()
  all.time.num<-dim(input.matrix)[2]-4###100
  tg<-matrix(0,nrow=1,ncol=length(lg))
  for(g in 1:length(lg)){
    for(i in 1:all.time.num){
      if(input.matrix[2,i+4]==lg[g])
        tg[g]<-tg[g]+1
    }
  }
  data<-data[-1,]
  ####--------------------------------PCA-------------------------------
  if(scform=="PCA"){

    data.ma<-data[-1,-c(1:4)]
    data.ma<-t(data.ma)
    data.ma<-as.matrix(data.ma)
    data.ma1<-as.numeric(data.ma)
    dim(data.ma1)<-dim(data.ma)

    for(i in 1:dim(data.ma1)[1]){
      if(sum(is.na(data.ma1[i,]))==dim(data.ma1)[2]){
        data.ma1[i,]<-0
      }
    }
    data.ma1<-scale(data.ma1,center = F, scale = TRUE)
    #data.ma1[is.na(data.ma1)]<-0
    # ma.pca<-mixOmics::pca(data.ma1,center = TRUE)
    wellpca<-pcaMethods::pca(data.ma1,method = "nipals",nPcs = num.of.cp)
    ma.score<-wellpca@scores
    ma.loading<-wellpca@loadings
    r2_pca<-wellpca@R2
    r2_pca<-round(r2_pca,2)
    q2_pca<-pcaMethods::Q2(wellpca, originalData = completeObs(wellpca),
               variables = 1:pcaMethods::nVar(wellpca))
    q2_pca<-round(q2_pca,2)
    dirout = paste(diroutall, "/PCA(male)", "/", sep = "")
    dir.create(dirout)

    pwdxdef = paste(dirout, "PCA-score(male).xlsx", sep = "")
    xlsx::write.xlsx(ma.score, pwdxdef)
    pwdxd = paste(dirout, "PCA-loading(male).xlsx", sep = "")
    xlsx::write.xlsx(ma.loading, pwdxd)
    tutticolors = matrix(c(4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink",
                           4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink"),
                         ncol = 1)
    col = c()
    for (i in 1:time.points) {
      col = c(col, tutticolors[i, ])
    }
    allcol<-c()
    for(i in 1:length(tg)){
      allcol<-c(allcol,rep(col[i],tg[i]))
    }


    pcafigpath = paste(dirout, "PCA-scorePlot(male).pdf", sep = "")
    par(mai=c(1,0.5,0.5,0.3),cex.lab=1,lab=c(3,3,5),mgp=c(1,1,0))
    pdf(pcafigpath)

    plot(x=ma.score[,1],y=ma.score[,2],pch=16 ,xaxt="n", yaxt="n",col=allcol,main = "PCA score figure(male)",xlab = "PC1",ylab = "PC2")
    xy<-par("usr")
    legend("bottomright",c(paste("R2Y(PC1)=",r2_pca[1]),paste("R2Y(PC2)=",r2_pca[2]),paste("Q2Y(PC1)=",q2_pca[1]),paste("Q2Y(PC2)=",q2_pca[2])), cex=0.7);#new add
    legend(x=xy[1]+xinch(4.5),y=xy[3]-yinch(0.1),legend=legc,xpd=TRUE,fill = col, col=col,cex = 0.6,ncol = 3,x.intersp=0.8,xjust=0,yjust=1)
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
    dev.off()
    #pca2
    pca2figpath = paste(dirout, "PCA-scorePlot-track(male).pdf", sep = "")
    pdf(pca2figpath)

    groupsh<-input.matrix[2,-c(1:4)]
    groupsh<-as.matrix(t(groupsh))
    groupsh<-as.data.frame(groupsh)

    g.ma.score<-cbind(groupsh,ma.score)
    PC1.mean<-tapply(g.ma.score[,2],g.ma.score[,1],mean)
    PC2.mean<-tapply(g.ma.score[,3],g.ma.score[,1],mean)
    plot(x=PC1.mean,y=PC2.mean,type = "b",pch=21,bg = "red")
    title("PCA Score(track)-male")
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
    text(PC1.mean, PC2.mean,legct,pos = 1)
    dev.off()
    ###

  }
  ####--------------------------------PLSDA-----------------------------
  if(scform=="PLSDA"){

    data.ma<-data[-1,-c(1:4)]
    data.ma<-t(data.ma)
    data.ma<-as.matrix(data.ma)
    data.ma1<-as.numeric(data.ma)
    dim(data.ma1)<-dim(data.ma)
    data.ma2<-scale(data.ma1,center = F, scale = TRUE)
    #data.ma1<-data.ma1[,order(data.ma1[dim(data.ma1)[1],])]
    y.group<-data[1,-c(1:4)]

    ma.plsda<-mixOmics::plsda(data.ma2,y.group,mode = "classic",ncomp = num.of.cp)
    loading.plsda<-ma.plsda$loadings

    ######add
    perf_plsda<-mixOmics::perf(ma.plsda,validation = "Mfold", folds = 5)
    overall_plsda<-perf_plsda$error.rate$overall
    #the Prediction error rate for each block of object$X and each dist
    y_group<-as.numeric(y.group)
    data.ma3<- data.ma2
    data.ma3[is.na(data.ma2)] <- 0
    #data.ma3<-as.numeric(data.ma3)
    ma.pls<-mixOmics::pls(data.ma3,y_group,mode = "classic",ncomp = num.of.cp)
    perf_pls<-mixOmics::perf(ma.pls,validation = "Mfold")
    r2<-round(perf_pls$R2,2)[c(1:2)]
    q2<-round(perf_pls$Q2.total,2)[c(1:2)]
    data.ma4 = data.ma3


    dirout = paste(diroutall, "/PLSDA(male)", "/", sep = "")
    dir.create(dirout)
    perfigpath = paste(dirout, "PermutationPlot.pdf", sep = "")
    sacurine.plsda <- ropls::opls(data.ma4, y_group,printL=FALSE,plotL=FALSE,permI=fold)
    typeC <-"permutation"
    ropls::plot(sacurine.plsda, typeVc = typeC,file.pdfC=perfigpath)
    pw_allover = paste(dirout, "ErrorRate.xlsx", sep = "")
    xlsx::write.xlsx(overall_plsda, pw_allover)

    pwdxdef = paste(dirout, "PLSDA-loading(male).xlsx", sep = "")
    xlsx::write.xlsx(loading.plsda$X, pwdxdef)
    score.plsda<-ma.plsda$variates$X
    pwdxd = paste(dirout, "PLSDA-score(male).xlsx", sep = "")
    xlsx::write.xlsx(score.plsda, pwdxd)

    tutticolors = matrix(c(4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink",
                           4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink"),
                         ncol = 1)
    col = c()
    for (i in 1:time.points) {
      col = c(col, tutticolors[i, ])
    }
    allcol<-c()
    for(i in 1:length(tg)){
      allcol<-c(allcol,rep(col[i],tg[i]))
    }

    pcafigpath = paste(dirout, "PLSDA-scorePlot(male).pdf", sep = "")
    par(mai=c(1,0.5,0.5,0.3),cex.lab=1,lab=c(3,3,5),mgp=c(1,1,0))
    pdf(pcafigpath)

    plot(x=score.plsda[,1],y=score.plsda[,2],pch=16 ,xaxt="n", yaxt="n",col=allcol,main = "PLSDA score figure(male)",xlab = "PC1",ylab = "PC2")
    xy<-par("usr")
    legend("bottomright",c(paste("R2Y(PC1)=",r2[1]),paste("R2Y(PC2)=",r2[2]),paste("Q2Y(PC1)=",q2[1]),paste("Q2Y(PC2)=",q2[2])), cex=0.7);
    legend(x=xy[1]+xinch(4.5),y=xy[3]-yinch(0.1),legend=legc,xpd=TRUE,fill = col, col=col,cex = 0.6,ncol = 3,x.intersp=0.8,xjust=0,yjust=1)
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
    dev.off()

    plsda2figpath = paste(dirout, "PLSDA-scorePlot-track(male).pdf", sep = "")
    pdf(plsda2figpath)

    groupsh<-input.matrix[2,-c(1:4)]
    groupsh<-as.matrix(t(groupsh))
    groupsh<-as.data.frame(groupsh)
    g.ma.score<-cbind(groupsh,score.plsda)
    PC1.mean<-tapply(g.ma.score[,2],g.ma.score[,1],mean)
    PC2.mean<-tapply(g.ma.score[,3],g.ma.score[,1],mean)
    plot(x=PC1.mean,y=PC2.mean,type = "b",pch=21,bg = "red")
    title("PLSDA Score(track)-male")
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
    text(PC1.mean, PC2.mean,legct,pos = 1)
    dev.off()
  }



  #####---------------------female----------------------
  scat.data<-femaledata
  scat.data<-as.data.frame(scat.data)
  group.num<-as.matrix(scat.data[2,dim(scat.data)[2]])
  group.num<-as.numeric(group.num)     ###num 9
  group.begin<-as.matrix(scat.data[2,5])
  group.begin<-as.numeric(group.begin)  ###num 0
  data<-as.matrix(scat.data)
  input.matrix<-scat.data
  lg<-group.begin:group.num
  time.points<-length(lg)            #########10
  list.sl<-list()
  all.time.num<-dim(input.matrix)[2]-4###100
  tg<-matrix(0,nrow=1,ncol=length(lg))
  for(g in 1:length(lg)){
    for(i in 1:all.time.num){
      if(input.matrix[2,i+4]==lg[g])
        tg[g]<-tg[g]+1
    }
  }
  data<-data[-1,]
  ####--------------------------------PCA-------------------------------
  if(scform=="PCA"){

    data.ma<-data[-1,-c(1:4)]
    data.ma<-t(data.ma)
    data.ma<-as.matrix(data.ma)
    data.ma1<-as.numeric(data.ma)
    dim(data.ma1)<-dim(data.ma)

    for(i in 1:dim(data.ma1)[1]){
      if(sum(is.na(data.ma1[i,]))==dim(data.ma1)[2]){
        data.ma1[i,]<-0
      }
    }
    data.ma1<-scale(data.ma1,center = F, scale = TRUE)
    #data.ma1[is.na(data.ma1)]<-0
    # ma.pca<-mixOmics::pca(data.ma1,center = TRUE)
    wellpca<-pcaMethods::pca(data.ma1,method = "nipals",nPcs = num.of.cp)
    ma.score<-wellpca@scores
    ma.loading<-wellpca@loadings

    r2_pca<-wellpca@R2
    r2_pca<-round(r2_pca,2)
    q2_pca<-pcaMethods::Q2(wellpca, originalData = completeObs(wellpca),
               variables = 1:pcaMethods::nVar(wellpca))
    q2_pca<-round(q2_pca,2)
    dirout = paste(diroutall, "/PCA(female)", "/", sep = "")
    dir.create(dirout)

    pwdxdef = paste(dirout, "PCA-score(female).xlsx", sep = "")
    xlsx::write.xlsx(ma.score, pwdxdef)
    pwdxd = paste(dirout, "PCA-loading(female).xlsx", sep = "")
    xlsx::write.xlsx(ma.loading, pwdxd)
    tutticolors = matrix(c(4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink",
                           4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink"),
                         ncol = 1)
    col = c()
    for (i in 1:time.points) {
      col = c(col, tutticolors[i, ])
    }
    allcol<-c()
    for(i in 1:length(tg)){
      allcol<-c(allcol,rep(col[i],tg[i]))
    }

    pcafigpath = paste(dirout, "PCA-scorePlot(female).pdf", sep = "")
    par(mai=c(1,0.5,0.5,0.3),cex.lab=1,lab=c(3,3,5),mgp=c(1,1,0))
    pdf(pcafigpath)

    plot(x=ma.score[,1],y=ma.score[,2],pch=16 ,xaxt="n", yaxt="n",col=allcol,main = "PCA score figure(female)",xlab = "PC1",ylab = "PC2")
    xy<-par("usr")
    legend("bottomright",c(paste("R2Y(PC1)=",r2_pca[1]),paste("R2Y(PC2)=",r2_pca[2]),paste("Q2Y(PC1)=",q2_pca[1]),paste("Q2Y(PC2)=",q2_pca[2])), cex=0.7);#new add
    legend(x=xy[1]+xinch(4.5),y=xy[3]-yinch(0.1),legend=legc,xpd=TRUE,fill = col, col=col,cex = 0.6,ncol = 3,x.intersp=0.8,xjust=0,yjust=1)
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
    dev.off()
    #pca2
    pca2figpath = paste(dirout, "PCA-scorePlot-track(female).pdf", sep = "")
    pdf(pca2figpath)

    groupsh<-input.matrix[2,-c(1:4)]
    groupsh<-as.matrix(t(groupsh))
    groupsh<-as.data.frame(groupsh)

    g.ma.score<-cbind(groupsh,ma.score)
    PC1.mean<-tapply(g.ma.score[,2],g.ma.score[,1],mean)
    PC2.mean<-tapply(g.ma.score[,3],g.ma.score[,1],mean)
    plot(x=PC1.mean,y=PC2.mean,type = "b",pch=21,bg = "red")
    title("PCA Score(track)-female")
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
    text(PC1.mean, PC2.mean,legct,pos = 1)
    dev.off()
    ###

  }
  ####--------------------------------PLSDA-----------------------------
  if(scform=="PLSDA"){

    data.ma<-data[-1,-c(1:4)]
    data.ma<-t(data.ma)
    data.ma<-as.matrix(data.ma)
    data.ma1<-as.numeric(data.ma)
    dim(data.ma1)<-dim(data.ma)
    data.ma2<-scale(data.ma1,center = F, scale = TRUE)
    #data.ma1<-data.ma1[,order(data.ma1[dim(data.ma1)[1],])]
    y.group<-data[1,-c(1:4)]

    ma.plsda<-mixOmics::plsda(data.ma2,y.group,mode = "classic",ncomp = num.of.cp)
    loading.plsda<-ma.plsda$loadings

    ######add
    perf_plsda<-mixOmics::perf(ma.plsda,validation = "Mfold" ,folds = 5)
    overall_plsda<-perf_plsda$error.rate$overall
    #the Prediction error rate for each block of object$X and each dist
    y_group<-as.numeric(y.group)
    data.ma3<- data.ma2
    data.ma3[is.na(data.ma2)] <- 0
    #data.ma3<-as.numeric(data.ma3)
    ma.pls<-mixOmics::pls(data.ma3,y_group,mode = "classic",ncomp = num.of.cp)
    perf_pls<-mixOmics::perf(ma.pls,validation = "Mfold")
    r2<-round(perf_pls$R2,2)[c(1:2)]
    q2<-round(perf_pls$Q2.total,2)[c(1:2)]
    data.ma4 = data.ma3


    dirout = paste(diroutall, "/PLSDA(female)", "/", sep = "")
    dir.create(dirout)
    perfigpath = paste(dirout, "PermutationPlot.pdf", sep = "")
    sacurine.plsda <- ropls::opls(data.ma4, y_group,printL=FALSE,plotL=FALSE,permI=fold)
    typeC <-"permutation"
    ropls::plot(sacurine.plsda, typeVc = typeC,file.pdfC=perfigpath)
    pw_allover = paste(dirout, "ErrorRate.xlsx", sep = "")
    xlsx::write.xlsx(overall_plsda, pw_allover)


    pwdxdef = paste(dirout, "PLSDA-loading(female).xlsx", sep = "")
    xlsx::write.xlsx(loading.plsda$X, pwdxdef)
    score.plsda<-ma.plsda$variates$X
    pwdxd = paste(dirout, "PLSDA-score(male).xlsx", sep = "")
    xlsx::write.xlsx(score.plsda, pwdxd)

    tutticolors = matrix(c(4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink",
                           4,1, 2,3, 5, 6, 7, 8, "rosybrown4","orange",
                           "green4", "navy", "purple2",  "pink", "chocolate2",
                           "coral3", "khaki3", "thistle", "turquoise3", "palegreen1",
                           "moccasin", "olivedrab3", "azure4", "gold3", "deeppink"),
                         ncol = 1)
    col = c()
    for (i in 1:time.points) {
      col = c(col, tutticolors[i, ])
    }
    allcol<-c()
    for(i in 1:length(tg)){
      allcol<-c(allcol,rep(col[i],tg[i]))
    }

    pcafigpath = paste(dirout, "PLSDA-scorePlot(female).pdf", sep = "")
    par(mai=c(1,0.5,0.5,0.3),cex.lab=1,lab=c(3,3,5),mgp=c(1,1,0))
    pdf(pcafigpath)

    plot(x=score.plsda[,1],y=score.plsda[,2],pch=16 ,xaxt="n", yaxt="n",col=allcol,main = "PLSDA score figure(female)",xlab = "PC1",ylab = "PC2")
    xy<-par("usr")
    legend("bottomright",c(paste("R2Y(PC1)=",r2[1]),paste("R2Y(PC2)=",r2[2]),paste("Q2Y(PC1)=",q2[1]),paste("Q2Y(PC2)=",q2[2])), cex=0.7);
    legend(x=xy[1]+xinch(4.5),y=xy[3]-yinch(0.1),legend=legc,xpd=TRUE,fill = col, col=col,cex = 0.6,ncol = 3,x.intersp=0.8,xjust=0,yjust=1)
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
    dev.off()

    plsda2figpath = paste(dirout, "PLSDA-scorePlot-track(female).pdf", sep = "")
    pdf(plsda2figpath)

    groupsh<-input.matrix[2,-c(1:4)]
    groupsh<-as.matrix(t(groupsh))
    groupsh<-as.data.frame(groupsh)
    g.ma.score<-cbind(groupsh,score.plsda)
    PC1.mean<-tapply(g.ma.score[,2],g.ma.score[,1],mean)
    PC2.mean<-tapply(g.ma.score[,3],g.ma.score[,1],mean)
    plot(x=PC1.mean,y=PC2.mean,type = "b",pch=21,bg = "red")
    title("PLSDA Score(track)-female")
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
    text(PC1.mean, PC2.mean,legct,pos = 1)
    dev.off()
  }











}
