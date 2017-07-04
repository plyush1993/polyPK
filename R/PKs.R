#' @title  Calculate the representative pharmacokinetics parameters and plot the time-intensity curves of specified compounds.
#'
#' @description A function to calculate the 7 pharmacokinetics parameters (Tmax, Cmax, AUC, CL, Tlast, Tfirst, Cmin) and plot the time-intensity curves for specified compounds.
#'
#' @param d.pk The data under analysis (data frame with required format).The first row is the variable names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected. The second row of the data frame should be the group information. The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose).
#' @param d.point The value of points in the time-intensity curve. d.point=c ("mean","median"). Default:"mean".
#' @param d.ebar The value of error bars. d.ebar=c ("SE","SD"). Defalut:"SE".
#' @param filepath A character string indicating the path where the results may be saved in.
#' @param design A study design dataset(data frame with required format).Use data(StudyDesign) to see the detailed form.
#'
#' @return list
#'
#' @examples PKs(d.pk=A,d.point="mean",d.ebar="SE",filepath=getwd(),design=design)
#'
#' @export
PKs<-function(d.pk,d.point="mean",d.ebar="SE",filepath=getwd(),design=FALSE){
 # library(PKNCA )
  #library(Hmisc)
  #--------------------prepare tha data--------------
  d.pk<-as.data.frame(d.pk)
  headall<-d.pk[,c(1:4)]
  meta.num<-dim(d.pk)[1]-3
  d.time<-d.pk[3,-c(1:4)]
  data_timepoints<-timec<-as.matrix(d.time)
  #timec<-round(timec,1)
  data_timepoints<-as.numeric(data_timepoints)
  timepoin<-data_timepoints[!duplicated(data_timepoints)]
  alltime<- length(timec)
  timecc<-matrix(0,nrow=dim(timec)[1],ncol=dim(timec)[2])

  timecc<-as.numeric(timec)
  d.group<-as.matrix(d.pk[2,-c(1:4)])
  d.group<-matrix(as.numeric(d.group), nrow = 1)
  d.gn<-as.matrix(d.group[length(d.group)])
  d.bn<-as.matrix(d.group[1])
  d.bn<-as.numeric(d.bn)
  d.gn<-as.numeric(d.gn)
  time<-c(d.bn:d.gn)
  lentime<-length(time)
  dim(timecc)<-c(alltime/lentime,lentime)
  timecourse<-timecc[1,]

  d.pk.data<-as.matrix(d.pk[,-c(1:4)])
  d.pk.data<-matrix(as.numeric(d.pk.data),nrow=nrow(d.pk.data))
  allnumdata<-d.pk.data
  d.pk.data<-d.pk.data[-c(1,3),]
  #------------choose the mean/median value to calculate the pharmacokinetics parameters
  if(d.point=="mean"){
    d.pk.m<-matrix(0,length(time),nrow=dim(d.pk.data)[1])
    for(i in 1: dim(d.pk.data)[1]){

      mean.data<-tapply(d.pk.data[i,],d.group,mean)
      for(j in 1: length(time)){
        d.pk.m[i,j]<-mean.data[j]
      }
    }
  }
  if(d.point=="median"){
    d.pk.m<-matrix(0,length(time),nrow=dim(d.pk.data)[1])
    for(i in 1: dim(d.pk.data)[1]){

      mean.data<-tapply(d.pk.data[i,],d.group,median)
      for(j in 1: length(time)){
        d.pk.m[i,j]<-mean.data[j]
      }
    }
  }
  d.pk.m[is.na(d.pk.m)]<-0
  aucdata<-matrix(0,nrow=dim(d.pk.data)[1]-1,ncol = 7)
  auc<-d.pk[-c(1:3),c(1:4)]# prepare a martrix for auc value
  auc<-cbind(auc,aucdata)
####----------------fill the 0 matrix(auc) with pharmacokinetics parameters----

  names(auc)<-c("Name","ID","m.z","R.T.min.","Tmax","Tlast","Tfirst","Cmax","Cmin","AUC","CL")
  for(i in 2:dim(d.pk.data)[1]){
    # id<-d.pk[1,2]
    myconc <- c(d.pk.m[i,])
    mytime<-c(time)
    auc[i-1,5]<-PKNCA::pk.calc.tmax(myconc,mytime)
    auc[i-1,6]<-PKNCA::pk.calc.tlast(myconc,mytime,check = TRUE)
    auc[i-1,7]<-PKNCA::pk.calc.tfirst(myconc,mytime,check = TRUE)
    auc[i-1,8]<-PKNCA::pk.calc.cmax(myconc, check = TRUE)
    auc[i-1,9]<-PKNCA::pk.calc.cmin(myconc,check = TRUE)
    auc[i-1,10]<-PKNCA::pk.calc.auc(myconc, mytime, interval=c(0, Inf))
    auc[i-1,11]<-PKNCA::pk.calc.cl(c(1:d.gn),auc[i-1,10])
  }
  pk.parameter<-auc #pk.parameter matrix of  id name and auc
  diroutall = paste(filepath, "/PKresults", "/", sep = "")
  dir.create(diroutall)
  dirout = paste(diroutall, "/PKresults(all)", "/", sep = "")
  dir.create(dirout)
  pwdxdef = paste(dirout, "PK-parameters.xlsx", sep = "")
  xlsx::write.xlsx(pk.parameter, pwdxdef,row.names = F)

  ##--------------fuction to calculate SD,SE value of data-----------
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
   # library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       median = median   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )

    # Rename the "mean" column
    datac <- plyr::rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
  }
  d.pk.data[is.na(d.pk.data)]<-0
  #pkfigure<-list()
  tgt<-t(d.pk.data)
  tgt<-as.data.frame(tgt)
  ##--------------plot the Time-Intensity-Curve figure with error bar--------------
 if(d.point=="mean"){
   if(d.ebar=="SE"){

    for(i in 1:meta.num){

      idnnm<-as.matrix(d.pk[i+2,1])
      obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
      pdf(obfigpath)
      tgt.da<-cbind(tgt[,1],tgt[,i+1])
      tgt.da<-as.data.frame(tgt.da)
      names(tgt.da)<-c("time.course","intensity")
      tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
      #"time point"<-tgc[,1]
      "time point"<-timecourse
      intensity<-tgc[,3]
      delta<-tgc[,6]
      Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
      lines(`time point`,intensity)
      #axis(1,tmcours,cex.axis=0.8)
      axis(1,round(timecourse,0),cex.axis=0.8)
      title(main =paste("Time-Intensity-Curve of",idnnm,"(mean-SE)"))
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
    }

  }

  if(d.ebar=="SD"){

    for(i in 1:meta.num){

      idnnm<-as.matrix(d.pk[i+2,1])
      obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
      pdf(obfigpath)
      tgt.da<-cbind(tgt[,1],tgt[,i+1])
      tgt.da<-as.data.frame(tgt.da)
      names(tgt.da)<-c("time.course","intensity")
      tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
      "time point"<-timecourse
      intensity<-tgc[,3]
      delta<-tgc[,5]
      #plot(tgc, aes(x=time.course, y=intensity))
      Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
      lines(`time point`,intensity)
      axis(1,round(timecourse,0),cex.axis=0.8)
      title(main =paste("Time-Intensity-Curve of",idnnm,"(mean-SD)"))
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
    }

  }
 }
  if(d.point=="median"){
    if(d.ebar=="SE"){

    for(i in 1:meta.num){

      idnnm<-as.matrix(d.pk[i+2,1])
      obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
      pdf(obfigpath)
      tgt.da<-cbind(tgt[,1],tgt[,i+1])
      tgt.da<-as.data.frame(tgt.da)
      names(tgt.da)<-c("time.course","intensity")
      tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
      "time point"<-timecourse
      intensity<-tgc[,4]
      delta<-tgc[,6]
      #plot(tgc, aes(x=time.course, y=intensity))
      Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
      lines(`time point`,intensity)

      axis(1,round(timecourse,0),cex.axis=0.8)
      title(main =paste("Time-Intensity-Curve of",idnnm,"(median-SE)"))
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
    }

  }

    if(d.ebar=="SD"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,4]
        delta<-tgc[,5]
        #plot(tgc, aes(x=time.course, y=intensity))
        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(median-SD)"))
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
      }

    }
  }
  ##-------------return the result------------------
 # pk<-list(PKparameter=pk.parameter,PKfigures=pkfigure)
  pk.parameterall<-pk.parameter


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
  ####---------------------------male-----------------------------####

  d.pk<-maledata
  meta.num<-dim(d.pk)[1]-3
  d.group<-as.matrix(d.pk[2,-c(1:4)])
  d.group<-matrix(as.numeric(d.group), nrow = 1)
  d.gn<-as.matrix(d.group[length(d.group)])
  d.bn<-as.matrix(d.group[1])
  d.bn<-as.numeric(d.bn)
  d.gn<-as.numeric(d.gn)
  time<-c(d.bn:d.gn)
  d.pk.data<-as.matrix(d.pk[,-c(1:4)])
  d.pk.data<-matrix(as.numeric(d.pk.data),nrow=nrow(d.pk.data))
  allnumdata<-d.pk.data
  d.pk.data<-d.pk.data[-c(1,3),]
  #####------------choose the mean/median value to calculate the pharmacokinetics parameters
  if(d.point=="mean"){
    d.pk.m<-matrix(0,length(time),nrow=dim(d.pk.data)[1])
    for(i in 1: dim(d.pk.data)[1]){

      mean.data<-tapply(d.pk.data[i,],d.group,mean)
      for(j in 1: length(time)){
        d.pk.m[i,j]<-mean.data[j]
      }
    }
  }
  if(d.point=="median"){
    d.pk.m<-matrix(0,length(time),nrow=dim(d.pk.data)[1])
    for(i in 1: dim(d.pk.data)[1]){

      mean.data<-tapply(d.pk.data[i,],d.group,median)
      for(j in 1: length(time)){
        d.pk.m[i,j]<-mean.data[j]
      }
    }
  }
  d.pk.m[is.na(d.pk.m)]<-0
  aucdata<-matrix(0,nrow=dim(d.pk.data)[1]-1,ncol = 7)
  auc<-d.pk[-c(1:3),c(1:4)]# prepare a martrix for auc value
  auc<-cbind(auc,aucdata)
  #----------------fill the 0 matrix(auc) with pharmacokinetics parameters----

  names(auc)<-c("Name","ID","m.z","R.T.min.","Tmax","Tlast","Tfirst","Cmax","Cmin","AUC","CL")
  for(i in 2:dim(d.pk.data)[1]){
    # id<-d.pk[1,2]
    myconc <- c(d.pk.m[i,])
    mytime<-c(time)
    auc[i-1,5]<-PKNCA::pk.calc.tmax(myconc,mytime)
    auc[i-1,6]<-PKNCA::pk.calc.tlast(myconc,mytime,check = TRUE)
    auc[i-1,7]<-PKNCA::pk.calc.tfirst(myconc,mytime,check = TRUE)
    auc[i-1,8]<-PKNCA::pk.calc.cmax(myconc, check = TRUE)
    auc[i-1,9]<-PKNCA::pk.calc.cmin(myconc,check = TRUE)
    auc[i-1,10]<-PKNCA::pk.calc.auc(myconc, mytime, interval=c(0, Inf))
    auc[i-1,11]<-PKNCA::pk.calc.cl(c(1:d.gn),auc[i-1,10])
  }
  pk.parameter<-auc #pk.parameter matrix of  id name and auc
  dirout = paste(diroutall, "/PKresults(male)", "/", sep = "")
  dir.create(dirout)
  pwdxdef = paste(dirout, "PK-parameters-male.xlsx", sep = "")
  xlsx::write.xlsx(pk.parameter, pwdxdef,row.names = F)

  ##--------------fuction to calculate SD,SE value of data-----------
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    # library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                         .fun = function(xx, col) {
                           c(N    = length2(xx[[col]], na.rm=na.rm),
                             mean = mean   (xx[[col]], na.rm=na.rm),
                             median = median   (xx[[col]], na.rm=na.rm),
                             sd   = sd     (xx[[col]], na.rm=na.rm)
                           )
                         },
                         measurevar
    )

    # Rename the "mean" column
    datac <- plyr::rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
  }
  d.pk.data[is.na(d.pk.data)]<-0
  #pkfigure<-list()
  tgt<-t(d.pk.data)
  tgt<-as.data.frame(tgt)
  ##--------------plot the Time-Intensity-Curve figure with error bar--------------
  if(d.point=="mean"){
    if(d.ebar=="SE"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,3]
        delta<-tgc[,6]

        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)

        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(mean-SE)"))
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
      }

    }

    if(d.ebar=="SD"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,3]
        delta<-tgc[,5]
        #plot(tgc, aes(x=time.course, y=intensity))
        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(mean-SD)"))
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
      }

    }
  }
  if(d.point=="median"){
    if(d.ebar=="SE"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,4]
        delta<-tgc[,6]
        #plot(tgc, aes(x=time.course, y=intensity))
        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(median-SE)"))
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
      }

    }

    if(d.ebar=="SD"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,4]
        delta<-tgc[,5]
        #plot(tgc, aes(x=time.course, y=intensity))
        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(median-SD)"))
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
      }

    }
  }

  ####---------------------------female-----------------------------####

  d.pk<-femaledata
  meta.num<-dim(d.pk)[1]-3
  d.group<-as.matrix(d.pk[2,-c(1:4)])
  d.group<-matrix(as.numeric(d.group), nrow = 1)
  d.gn<-as.matrix(d.group[length(d.group)])
  d.bn<-as.matrix(d.group[1])
  d.bn<-as.numeric(d.bn)
  d.gn<-as.numeric(d.gn)
  time<-c(d.bn:d.gn)
  d.pk.data<-as.matrix(d.pk[,-c(1:4)])
  d.pk.data<-matrix(as.numeric(d.pk.data),nrow=nrow(d.pk.data))
  allnumdata<-d.pk.data
  d.pk.data<-d.pk.data[-c(1,3),]
  #####------------choose the mean/median value to calculate the pharmacokinetics parameters
  if(d.point=="mean"){
    d.pk.m<-matrix(0,length(time),nrow=dim(d.pk.data)[1])
    for(i in 1: dim(d.pk.data)[1]){

      mean.data<-tapply(d.pk.data[i,],d.group,mean)
      for(j in 1: length(time)){
        d.pk.m[i,j]<-mean.data[j]
      }
    }
  }
  if(d.point=="median"){
    d.pk.m<-matrix(0,length(time),nrow=dim(d.pk.data)[1])
    for(i in 1: dim(d.pk.data)[1]){

      mean.data<-tapply(d.pk.data[i,],d.group,median)
      for(j in 1: length(time)){
        d.pk.m[i,j]<-mean.data[j]
      }
    }
  }
  d.pk.m[is.na(d.pk.m)]<-0
  aucdata<-matrix(0,nrow=dim(d.pk.data)[1]-1,ncol = 7)
  auc<-d.pk[-c(1:3),c(1:4)]# prepare a martrix for auc value
  auc<-cbind(auc,aucdata)
  #----------------fill the 0 matrix(auc) with pharmacokinetics parameters----

  names(auc)<-c("Name","ID","m.z","R.T.min.","Tmax","Tlast","Tfirst","Cmax","Cmin","AUC","CL")
  for(i in 2:dim(d.pk.data)[1]){
    # id<-d.pk[1,2]
    myconc <- c(d.pk.m[i,])
    mytime<-c(time)
    auc[i-1,5]<-PKNCA::pk.calc.tmax(myconc,mytime)
    auc[i-1,6]<-PKNCA::pk.calc.tlast(myconc,mytime,check = TRUE)
    auc[i-1,7]<-PKNCA::pk.calc.tfirst(myconc,mytime,check = TRUE)
    auc[i-1,8]<-PKNCA::pk.calc.cmax(myconc, check = TRUE)
    auc[i-1,9]<-PKNCA::pk.calc.cmin(myconc,check = TRUE)
    auc[i-1,10]<-PKNCA::pk.calc.auc(myconc, mytime, interval=c(0, Inf))
    auc[i-1,11]<-PKNCA::pk.calc.cl(c(1:d.gn),auc[i-1,10])
  }
  pk.parameter<-auc #pk.parameter matrix of  id name and auc
  dirout = paste(diroutall, "/PKresults(female)", "/", sep = "")
  dir.create(dirout)
  pwdxdef = paste(dirout, "PK-parameters-female.xlsx", sep = "")
  xlsx::write.xlsx(pk.parameter, pwdxdef,row.names = F)

  ##--------------fuction to calculate SD,SE value of data-----------
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    # library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                         .fun = function(xx, col) {
                           c(N    = length2(xx[[col]], na.rm=na.rm),
                             mean = mean   (xx[[col]], na.rm=na.rm),
                             median = median   (xx[[col]], na.rm=na.rm),
                             sd   = sd     (xx[[col]], na.rm=na.rm)
                           )
                         },
                         measurevar
    )

    # Rename the "mean" column
    datac <- plyr::rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
  }
  d.pk.data[is.na(d.pk.data)]<-0
  #pkfigure<-list()
  tgt<-t(d.pk.data)
  tgt<-as.data.frame(tgt)
  ##--------------plot the Time-Intensity-Curve figure with error bar--------------
  if(d.point=="mean"){
    if(d.ebar=="SE"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,3]
        delta<-tgc[,6]

        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(mean-SE)"))
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
      }

    }

    if(d.ebar=="SD"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,3]
        delta<-tgc[,5]
        #plot(tgc, aes(x=time.course, y=intensity))
        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(mean-SD)"))
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
      }

    }
  }
  if(d.point=="median"){
    if(d.ebar=="SE"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,4]
        delta<-tgc[,6]
        #plot(tgc, aes(x=time.course, y=intensity))
        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(median-SE)"))
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
      }

    }

    if(d.ebar=="SD"){

      for(i in 1:meta.num){

        idnnm<-as.matrix(d.pk[i+2,1])
        obfigpath = paste(dirout,"Time-Intensity-Curve of", idnnm,".pdf", sep = "")
        pdf(obfigpath)
        tgt.da<-cbind(tgt[,1],tgt[,i+1])
        tgt.da<-as.data.frame(tgt.da)
        names(tgt.da)<-c("time.course","intensity")
        tgc <- summarySE(tgt.da, measurevar="intensity", groupvars="time.course")
        "time point"<-timecourse
        intensity<-tgc[,4]
        delta<-tgc[,5]
        #plot(tgc, aes(x=time.course, y=intensity))
        Hmisc::errbar(`time point`,intensity,intensity+delta,intensity-delta,xaxt="n")
        lines(`time point`,intensity)
        axis(1,round(timecourse,0),cex.axis=0.8)
        title(main =paste("Time-Intensity-Curve of",idnnm,"(median-SD)"))
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
      }

    }
  }









  return(pk.parameterall)
  }

#pktest1<-PKs(d.pk,"median","SE",filepath=getwd())
