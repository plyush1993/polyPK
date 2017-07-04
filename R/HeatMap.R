#' @title Plot the heatmap of input data.
#'
#' @description A function to plot the heatmap and clusters of input data.
#'
#' @param data The data under analysis (data.frame with required format). The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected. The second row of the data frame should be the group information. The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). Please see the demo data for detailed format.
#' @param cluster A string indicating whether or in which direction the dendrograms should be drawn (“none”, “row”, “column” or “both”). Default:”both”.
#' @param scale  A character indicating whether the data should be centered and scaled before analysis and in which (“none”, “row” or “column”) direction. Default:"none".
#' @param filepath A character string indicating the path where the results may be saved in.
#' @param design A study design dataset(data frame with required format).Use data(StudyDesign) to see the detailed form.
#'
#' @return NULL
#'
#' @examples HeatMap(data=B,cluster="both",scale="row",filepath=getwd(),design=design)
#'
#' @export
HeatMap<-function(data,cluster="both",scale="row",filepath=getwd(),design=FALSE){
#library(gplots)
  ####----------prepare the data for plot---------------
  data<-as.data.frame(data)
  data_timepoints<-as.matrix(data[3,-c(1:4)])
  data_timepoints<-as.numeric(data_timepoints)
  timepoin<-data_timepoints[!duplicated(data_timepoints)]

  data.id<-data[-c(1:3),-c(1,3:4)]#choose the data without names
  id.name<-data.id[,1]#the IDs of the data
  rownames(data.id)<-id.name#rename the rownames with IDs
  data.d<-data.id[,-1]
  data.d[is.na(data.d)]<-0
  data.d<-data.or<-as.matrix(data.d)
 data.d<-as.numeric(data.d)
 dim(data.d)<-dim(data.or)
 rownames(data.d)<-id.name#rename the rownames with IDs
 ####------------plot the heatmap and save the reslut-------------------
 outpath=paste(filepath, "/HeatMap", "/", sep = "")
  dir.create(outpath)
   hfigpath = paste(outpath, "/HeatMap of data.pdf", sep = "")
   pdf(hfigpath)
   gplots::heatmap.2(data.d,scale = scale,dendrogram =cluster,cexRow = 0.5 ,cexCol = 0.3,trace="none" )
title("Heat map of data",outer = F )

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

