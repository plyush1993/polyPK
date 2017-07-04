#' @title preprocess the input data
#'
#' @description Preprocess the input data. Variables with a lot of zeros and outliers may be removed. Missing values may be imputed and filled by various methods. Data may be transformed by logarithm transformation.
#'
#' @param tes The data under pretreatment (data frame with required format). The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected. The second row of the data frame should be the group information. The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). Please see the demo data for detailed format.
#' @param rz  The percentage of zeros for variable elimination (Default:80). Variables with zero numbers higher than rz% will be removed.
#' @param mv The method of missing values imputation (Default:"mean"). mv=c ("mean", "groupmean", "median", "groupmedian", "groupmin", "min", "knn", "svd", "rf", "qrilc").
#' @param multiple The parameter for missing values imputation. Missing values will be replaced by multiple*mean/median/min (Default:0.1)
#' @param sv A logical value indicating whether to remove the outliers (Default:TRUE). The data which distance to the mean is bigger than 1.5 times of the difference value between lower quartile and upper quartile, should be identified as an outlier. And it will be replaced by the mean value of corresponding row.
#' @param log A logical value indicating whether to take the logarithm on the data (Default:FALSE)
#' @param filepath A character string indicating the path where the results may be saved in
#'
#' @return data.frame
#'
#' @examples DataPre(tes=postData,mv="mean",rz=80,multiple=0.1,sv=TRUE,log=FALSE,filepath=getwd())
#'
#' @export
DataPre <- function(tes,mv="mean",rz=80,multiple=0.1,sv=TRUE,log=FALSE,filepath=getwd()) {

  gend_g<-tes[1,]
  row_g<-tes[2,]#group information
  ti_g<-tes[3,]
  groupnumber<-tes[2,-c(1:4)]
  tes<-tes[-c(1:3),]
  headname<-tes[,c(1:4)]#headname information of the metabolites
  tntoz<-t3<-t1<-tes[,-c(1:4)]; #t1 all data

  #-------------------more than rz% zeros, delete the metabolites--------------------------------

  t2<-tesmv<-t1
  poz=rz*(length(tesmv)-4)/100
  f<-function(x){sum(x==0)}
  numof0<-apply(t2,1,f)
  numof0[is.na(numof0)] <- 0
  d2<-cbind(t2,numof0)
  d2<-cbind(headname,d2)
  tesrz<-subset(d2,subset=(numof0<poz))# subset num of zeros<80%(default)
  tes<-tesrz[,-dim(tesrz)[2]]
  #add row num of zero(95)
  tes[tes==0]<-NA;
  headname<-tes[,c(1:4)]
  tntoz<-t3<-t1<-tes[,-c(1:4)];

  #------------------function of replace the zeros or NA--------------------------------

  ###  method1:replace the zeros or NAs with the mean value of the group

  if(mv=="groupmean"){

    t_g<-rbind(groupnumber,t1)
    mean_gg<-data.frame()
    y<-as.numeric(groupnumber)
    for(i in 1:dim(t1)[1]){
      x<-as.numeric(t1[i,])
      mean_g<-t(tapply(x,y,mean,na.rm = TRUE))
      mean_g<-as.data.frame(mean_g)
      mean_gg<-rbind(mean_gg,mean_g)
    }
    mean_gg<-mean_gg*multiple



    if(y[1]!=0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,y[j]]
        }
      }
    }
    if(y[1]==0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,1]
        }
      }
    }

    t1[is.na(t1)]<-0
    tesmv<-t1
  }

  ###  method2:replace the zeros or NAs with the median value of the group

  if(mv=="groupmedian"){
    t_g<-rbind(groupnumber,t1)
    mean_gg<-data.frame()
    y<-as.numeric(groupnumber)
    for(i in 1:dim(t1)[1]){
      x<-as.numeric(t1[i,])
      mean_g<-t(tapply(x,y,median,na.rm = TRUE))
      mean_g<-as.data.frame(mean_g)
      mean_gg<-rbind(mean_gg,mean_g)
    }
    mean_gg<-mean_gg*multiple

    if(y[1]!=0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,y[j]]
        }
      }
    }
    if(y[1]==0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,1]
        }
      }
    }
    t1[is.na(t1)]<-0
    tesmv<-t1
  }
  ###  method3:replace the zeros or NAs with the minimum value of the group

  if(mv=="groupmin"){
    t_g<-rbind(groupnumber,t1)
    mean_gg<-data.frame()
    y<-as.numeric(groupnumber)
    for(i in 1:dim(t1)[1]){
      x<-as.numeric(t1[i,])
      mean_g<-t(tapply(x,y,min,na.rm = TRUE))
      mean_g<-as.data.frame(mean_g)
      mean_gg<-rbind(mean_gg,mean_g)
    }
    mean_gg<-mean_gg*multiple

    if(y[1]!=0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,y[j]]
        }
      }
    }
    if(y[1]==0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,1]
        }
      }
    }
    t1[is.na(t1)]<-0
    tesmv<-t1
  }
  ###  method4:replace the zeros or NAs with the mean value of the metabolite
  if(mv=="mean"){
    r<-apply(t1, 1, mean,na.rm = TRUE)
    r<-r*multiple
    for(i in 1:length(tes[,1])) tes[i,][is.na(tes[i,])]<-r[i]
    tesmv<-tes[,-c(1:4)]
  }
  ###  method5:replace the zeros or NAs with the median value of the metabolite
  if(mv=="median"){
    r<-apply(t1, 1, median,na.rm = TRUE)
    r<-r*multiple
    for(i in 1:length(tes[,1])) tes[i,][is.na(tes[i,])]<-r[i]
    tesmv<-tes[,-c(1:4)]
  }
  ###  method6:replace the zeros or NAs with the minimum value of the metabolite
  if(mv=="min"){
    r<-apply(t1, 1, min,na.rm = TRUE)
    r<-r*multiple
    for(i in 1:length(tes[,1])) tes[i,][is.na(tes[i,])]<-r[i]
    tesmv<-tes[,-c(1:4)]
  }
  #output tesmv(changed)
  ### method7: knn from tha package of "impute"
  if(mv=="knn"){
   # library(impute)
    t1knn <- impute::impute.knn(as.matrix(t1)) #t1knn data numeric
    tt<-as.data.frame(t1knn$data)#change the data form to dataframe-tt
    tes<-cbind(headname,tt)#combinate with headname-
    tesmv<-tes[,-c(1:4)]
  }
  ###method 8:rf-from the package of "missForest"
  if(mv=="rf"){
   # library(missForest)
    t1rf<-missForest::missForest(t(t1))#t-ed missforest
    ttrf<-t(t1rf$ximp)#data
    tr<-as.data.frame(ttrf)#dataframe
    tes<-cbind(headname,tr)#combinate
    tesmv<-tes[,-c(1:4)]
  }
  ####method 9:svd:from the package of "pcaMethods"
  if(mv=="svd"){
  #  library(pcaMethods)
    t1<-log(t1)
    ts <- pcaMethods::pca(t1, method="svdImpute", nPcs=5, center = TRUE)
    ttsvd <- pcaMethods::completeObs(ts) #data
    tsvd<-as.data.frame(ttsvd)
    tes<-cbind(headname,tsvd)
    tesmv<-tes[,-c(1:4)]
    tesmv<-exp(tesmv)
  }
  ####method 10 :qrilc, from teh package of "imputeLCMD"
  if(mv=="qrilc"){
   # library(imputeLCMD)
    t1<-log(t1)
    tq<-imputeLCMD::impute.QRILC(t1)
    ttqr<-as.data.frame(tq[1])
    tes<-cbind(headname,ttqr)
    tesmv<-tes[,-c(1:4)]
    tesmv<-exp(tesmv)
  }
  if(mv=="FASLE"){
    tesmv<-t1
  }
  #tesrz: name data nof0
  #output tesmv(changed)
  #--------------------------if the data take logarithms--------------------------------
  t3<-tesmv
  if(log=="TRUE"){
    tlog<-log(t3)

  }
  if(log=="e"){
    tlog<-log(t3)
  }
  if(log==2){
    tlog<-log2(t3)
  }
  if(log==10){
    tlog<-log10(t3)
  }

  if(log=="FALSE"){
    tlog<-t3

  }
  #output tlog-data
  ttlog<-cbind(headname,tlog) #name data
  #--------------------------------function to replace the outlier-----------------------
  #ttlog<-rbind(row_g,ttlog)
  if(sv=="TRUE"){
    #
    ndt01<-data.frame()
    for(i in 1: dim(tlog)[1]){

      value<-as.numeric(tlog[i,])
       QL <- quantile(value, probs = 0.25,na.rm=T)
      QU <- quantile(value, probs = 0.75,na.rm=T)
      QU_QL <- QU-QL
      # QL;QU;QU_QL
      which(value > QU + 1.5*QU_QL)
      value[which(value > QU + 1.5*QU_QL)]
      test01 <- value
      out_imp01 <- max(test01[which(test01 <= QU + 1.5*QU_QL)])
      test01[which(test01 > QU + 1.5*QU_QL)] <- out_imp01
      dt01<-as.data.frame(t(test01))
      ndt01<-rbind(ndt01,dt01)
    }
    ndt01<-t1<-as.matrix(ndt01)
    #  ndt01[is.na(ndt01)]<-
    groupnumber<-as.character(groupnumber)
    t_g<-rbind(groupnumber,t1)
    mean_gg<-data.frame()
    y<-as.numeric(groupnumber)

    for(i in 1:dim(t1)[1]){
      x<-as.numeric(t1[i,])
      mean_g<-t(tapply(x,y,mean,na.rm = TRUE))
      mean_g<-as.data.frame(mean_g)
      mean_gg<-rbind(mean_gg,mean_g)
    }
    mean_gg<-mean_gg*multiple

    if(y[1]!=0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,y[j]]
        }
      }
    }
    if(y[1]==0)
    {
      for(i in 1:dim(t1)[1]){
        for(j in 1:dim(t1)[2]){
          t1[i,j][is.na(t1[i,j])]<-mean_gg[i,1]
        }
      }
    }
    t1[is.na(t1)]<-0
    ndt01<-t1

    ###
    ttlog<-cbind(headname,ndt01)
    ttlog<-rbind(as.matrix(ti_g) ,as.matrix(ttlog) )
    ttlog<-as.data.frame(ttlog)
  }
  if(sv=="FALSE"){
    ttlog<-rbind(ti_g,ttlog)
  }
  ttlog<-rbind(row_g,ttlog)
  ttlog<-rbind(gend_g,ttlog)

#-------------------------------------saving the result-----------------------------------------------
  dirout = paste(filepath, "/preprocessed-data", "/", sep = "")
  dir.create(dirout)
  pwdxdef = paste(dirout, "preprocessed-data.xlsx", sep = "")
  xlsx::write.xlsx(ttlog,pwdxdef,row.names = F)


  return(ttlog)
}
