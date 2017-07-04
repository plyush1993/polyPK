#' @title Get the differential compounds across all the time points
#'
#' @description A function to get all the differential compounds between the pre-dose and every post-dose datasets
#'
#' @param pre_p  the datasets of pre-dose (data frame) with an indicator of time points (grouping variable) at the second row.
#' @param pos_p The post-dose dataset (data frame) with an indicator of time points (grouping variable) at the second row.
#' @param mv The method of missing values imputation (Default:"mean"). mv=c ("mean", "groupmean", "median", "groupmedian", "groupmin", "min", "knn", "svd", "rf", "qrilc")
#' @param rz  The percentage of zeros for variable elimination (Default: 80).
#' @param multiple The parameter for missing values imputation. Missing values will be replaced by multiple*mean/median/min (Default: 0.1).
#' @param sv A logical value indicating whether to remove the outliers (Default:TRUE). The data which distance to the mean is bigger than 1.5 times of the difference value between lower quartile and upper quartile, should be identified as an outlier. And it will be replaced by the mean value of corresponding row.
#' @param log A logical value indicating whether to take the logarithm on the datasets (Default:TRUE)
#' @param t The method for differential compounds identification. C ("Ttest", "MWtest"). Default: "Ttest". Compounds with p values less than 0.05 were taken as differential ones.
#' @param r.adj The methods for p values adjustment. r.adj=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Default: "fdr".
#' @param filepath A character string indicating the path where the results may be saved in.
#' @param simidata The similar compounds of drug and pre-dose metabolites,which is derived froem the Simi function.
#' @param design A study design dataset(data frame with required format).Use data(StudyDesign) to see the detailed form.
#'
#' @return NULL
#'
#' @examples GetDiffData(preData,postData,simidata)
#'
#' @export
GetDiffData<-function(preData,postData,simidata,mv="mean",rz=80,multiple=0.1,sv=TRUE,log=FALSE,t="Ttest",r.adj="fdr",filepath=getwd(),design=F){
  #library(plyr)
  #library(sqldf)
  pre_p<-as.data.frame(preData)
  pos_p<-as.data.frame(postData)
  ####--------------the function of prepocssing the input data-----------------####
  DataPre_gdd <- function(tes,mv="mean",rz=80,multiple=0.1,sv=TRUE,log=FALSE) {

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


    return(ttlog)
  }

  ####------------------get two prepocessed datasets---------------------------####
  pre<-DataPre_gdd(pre_p,mv,rz,multiple,sv,log)
  pos<-DataPre_gdd(pos_p,mv,rz,multiple,sv,log)
  pre.gender<-pre[1,]
  pos.gender<-pos[1,-c(1:4)]
  row.gender<-cbind(pre.gender,pos.gender)
  pre.group<-pre[2,]
  pos.group<-pos[2,-c(1:4)]
  row.group<-cbind(pre.group,pos.group)
  pre.tp<-pre[3,]
  pos.tp<-pos[3,-c(1:4)]
  row.tp<-cbind(pre.tp,pos.tp)
####--------------------find the three kinds of IDs-------------------------####
  a<-data.frame(pre[-c(1:3),2])
  b<-data.frame(pos[-c(1:3),2])
  pret<-sqldf::sqldf('SELECT * FROM [a] EXCEPT SELECT * FROM [b]')#pre have;pos not have
  post<-sqldf::sqldf('SELECT * FROM [b] EXCEPT SELECT * FROM [a]')# only post have
  pboth<-sqldf::sqldf('SELECT * FROM [a] intersect SELECT * FROM [b]')#both have
  post<-as.matrix(post)
  pret<-as.matrix(pret)
  pboth<-as.matrix(pboth)
  dr<-data.frame()
  for(j in 1:length(pret)){

    for(i in 1:dim(pre)[1]){
      if(pre[i,2]==pret[j]){
        dr<-rbind(pre[i,],dr)
      }
    }
  }
  drpre<-dr #pre have/ pos not
  dr_p<-data.frame()
  for(j in 1:length(pret)){

    for(i in 1:dim(pre_p)[1]){
      if(pre_p[i,2]==pret[j]){
        dr_p<-rbind(pre_p[i,],dr_p)
      }
    }
  }
  drpre_p<-dr_p #pre have/ pos not(original)
  dr<-data.frame()
  for(j in 1:length(post)){

    for(i in 1:dim(pos)[1]){
      if(pos[i,2]==post[j]){
        dr<-rbind(pos[i,],dr)
      }
    }
  }
  drpos<-dr   #pos have/ pre not
  dr_p<-data.frame()
  for(j in 1:length(post)){

    for(i in 1:dim(pos_p)[1]){
      if(pos_p[i,2]==post[j]){
        dr_p<-rbind(pos_p[i,],dr_p)
      }
    }
  }
  drpos_p<-dr_p#pos have/ pre not(orignal)
  dr<-data.frame()
  for(j in 1:length(pboth)){

    for(i in 1:dim(pre)[1]){
      if(pre[i,2]==pboth[j]){
        dr<-rbind(pre[i,],dr)
      }
    }
  }
  drbopr<-dr#both-pre
  dr<-data.frame()
  for(j in 1:length(pboth)){

    for(i in 1:dim(pos)[1]){
      if(pos[i,2]==pboth[j]){
        dr<-rbind(pos[i,],dr)
      }
    }
  }
  drbopo<-dr#both-pos 19(the last row is group)X1.1---
  drbopo1<-drbopo[,-c(1:4)]#x1.1--data
  drboth<-cbind(drbopr,drbopo1)# all time both data X0.1---

  group<-row.group
  testallgroup<-rbind(row.group,drboth)
  #testallgroup (both)
  groupnum<-as.numeric (apply(group[5:length(group)],1,max))# number of group(9)
  lg<-0:groupnum
  tg<-matrix(0,nrow=1,ncol=groupnum+1)
  for(g in 0:groupnum){
    for(i in 5:dim(testallgroup)[2]){
      if(testallgroup[1,i]==lg[g+1])
        tg[g+1]<-tg[g+1]+1
    }}
  li<-list()
  tb<-testallgroup
  for(g in 1:length(tg)){
    n<-g
    li[[g]]<-tb[-1,5:(4+tg[g])]
    tb<-tb[,-c(5:(4+tg[g]))]

  }
  lipre<-as.vector(li[[1]])

  p_group_2<-p_group<-matrix(0,nrow = dim(lipre)[1],ncol=length(tg)-1)
  ####---------------------------two method of test(t-test /Wilcox)---------------------------####
  if(t=="Ttest"){
    orglipre<-lipre
    lipre<-as.matrix(lipre)
    lipre<-as.numeric(lipre)
    dim(lipre)<-dim(orglipre)

    for(g in 2:length(tg)) {
      lipos<-as.vector(li[[g]])
      lipos<-as.matrix(lipos)
      lipos<-as.numeric(lipos)
      dim(lipos)<-dim(orglipre)
      for(i in 1:dim(lipre)[1]){
        # t_prepos<-wilcox.test(lipre[i,],lipos[i,])
        t_prepos<-t.test(lipre[i,],lipos[i,])
        p_group_2[i,g-1]<-t_prepos$p.value
      }
    }
  }
  if(t=="MWtest"){
    for(g in 2:length(tg)) {
      lipos<-as.vector(li[[g]])
      for(i in 1:dim(lipre)[1]){
        t_prepos<-wilcox.test(as.numeric(lipre[i,]),as.numeric(lipos[i,]),exact = FALSE)
        p_group_2[i,g-1]<-t_prepos$p.value
      }
    }
  }
  p_group<-p_group_2
  dim(p_group_2)<-dim(p_group)[1]*dim(p_group)[2]
  p_group_adj<-p.adjust(p_group_2,method =r.adj)
  dim(p_group_adj)<-c(dim(lipre)[1],length(tg)-1)
  #### -----------get the significant different metabolites from "both"----------------####
  # if all number of a row >0.05 delete
  # if all row<0 ,back the row nunmber of testallgroup
  ff<-function(x){sum(x<0.05) }
  p_g_a_s<-apply(p_group_adj,1,ff)# each row number of <0.05
  TorF<-matrix(0,nrow=length(p_g_a_s),ncol=1)
  for(i in 1:length(p_g_a_s)){

    if(p_g_a_s[i]!=0)TorF[i]<-i

  }
  drboth_data<-testallgroup[-1,]# all both have data (without group)
  T_both<-cbind(TorF,drboth_data)
  c_t<-c(0)
  for(i in 1:dim(T_both)[1]){
    if(TorF[i]==0) c_t<-rbind(c_t,c(i))
  }
  if(length(c_t)==1){test_both<-drboth_data}else{
    test_both<-drboth_data[-c_t,] # both-tested metabolites 3part with drpos drpre
  }
  A<-plyr::rbind.fill(test_both,drpos,drpre)
  idboth<-as.matrix(test_both[,2])
  bo_p<-data.frame()
  for(j in 1:length(idboth)){

    for(i in 2:dim(pos_p)[1]){
      if(pos_p[i,2]==idboth[j]){
        bo_p<-rbind(pos_p[i,],bo_p)
      }
    }
  }
  bopos_p<-bo_p[,-c(1:4)]
  bo_p<-data.frame()
  for(j in 1:length(idboth)){

    for(i in 2:dim(pre_p)[1]){
      if(pre_p[i,2]==idboth[j]){
        bo_p<-rbind(pre_p[i,],bo_p)
      }
    }
  }
  bopre_p<-bo_p
  bop_p<-cbind(bopre_p,bopos_p)
  A_pre<-plyr::rbind.fill(bop_p,drpos_p,drpre_p)

  ###add group row
  A<-rbind(row.tp,A)
  A<-rbind(row.group,A)
  A<-rbind(row.gender,A)
  A_pre<-rbind(row.tp,A_pre)
  A_pre<-rbind(row.group,A_pre)
  A_pre<-rbind(row.gender,A_pre)

 # A_pre.t<-as.data.frame(A_pre.t,row.names = NULL)
  ####--------------save the result in a file---------------------####
  Aresult<-list(A=A,A_pre=A_pre,p=p_group,p_adj=p_group_adj)

  ####--------A-----------------####
  datax<-A
  simiid<-simidata
  numsimi<-length(simiid) # the number of input similar metabolites
  rowdatax<-dim(datax)[1]
  similist<-data.frame()# the rows number of the similar metabolites
  for(j in 1:numsimi){
    for(i in 1:rowdatax){
      if(datax[i,2]==simiid[j])
        similist<-rbind(similist,i)
    }
  }
  numSi<-dim(similist)[1]# the number of similar metabolites in datax
  similist<-as.matrix(similist)

  dirout = paste(filepath, "/DifferentialMetabolites", "/", sep = "")
  dir.create(dirout)
  pwdxdef = paste(dirout, "DifferentialMetabolites(preprocessed).xlsx", sep = "")
  xlsx::write.xlsx(A, pwdxdef,row.names = F)
  wb.rr <-xlsx::loadWorkbook(pwdxdef)
  # getSheets get sheet
  sheets.in.rr <- xlsx::getSheets(wb.rr)
  # sheet.rr -first sheet
  sheet.rr <- sheets.in.rr[[1]]

  fill <- xlsx::Fill(foregroundColor="lightblue", backgroundColor="lightblue",
                     pattern="SOLID_FOREGROUND")
  coldatax<-dim(datax)[2]
  if(numSi!=0){
    for(rumk in 1:numSi){
      j<-similist[rumk]
      rows.rr <- xlsx::readRows(sheet.rr, startRow = j+1,
                                endRow = j+1, startColumn = 1,  endColumn = coldatax)
      # read store and remove the data
      block <- xlsx::CellBlock(sheet.rr,startRow = j+1,
                               startColumn=1, noRows=1,
                               noColumns=coldatax,create=FALSE)
      for(i in 1:coldatax){
        # fill tha background
        xlsx::CB.setFill(block, fill, colIndex = i, rowIndex=1)
        # add the data
        xlsx::CB.setColData(block, rows.rr[,i], i, rowOffset=0,
                            showNA=F, colStyle=NULL)
      }
    }
    xlsx::saveWorkbook(wb.rr,pwdxdef);
  }
  if(numSi==0){

    xlsx::saveWorkbook(wb.rr,pwdxdef);
  }

  ####study design
  if(is.data.frame(design)){

    wb.rr1 <-xlsx::loadWorkbook(pwdxdef)
    # getSheets get sheet
    sheets.in.rr1 <- xlsx::getSheets(wb.rr1)
    # sheet.rr -first sheet
    sheet.rr1 <- sheets.in.rr1[[1]]
    fill1 <- xlsx::Fill(foregroundColor="yellow2", backgroundColor="yellow2",
                        pattern="SOLID_FOREGROUND")
    fill2 <- xlsx::Fill(foregroundColor="gray64", backgroundColor="gray64",
                        pattern="SOLID_FOREGROUND")

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
    mealp<-mealp[,-1]### the time of meals
    mealtimepoint<-length(mealp)###the number of meal times


    mealcol<-data.frame()
    for(i in 1:mealtimepoint){
      mealti<-mealp[i]
      mealti<-as.matrix(mealti)
      mealti<-as.numeric(mealti)
      lowmealti<-mealti#the time of meal
      highmealti<-mealti+2# 4 hours after meal
      for(j in 5:dim(datax)[2]){
        dataxx<- as.matrix(datax[3,j])
        dataxx<-as.numeric(dataxx)
        if(lowmealti<=dataxx&dataxx<=highmealti)
          mealcol<-rbind(mealcol,j)
      }
    }
    nummealcol<-dim(mealcol)[1]# the number of similar metabolites in datax
    mealcol<-as.matrix(mealcol)
    coldatax<-dim(datax)[1]#row numbers 48
    if(nummealcol!=0){
      for(rumkmeal in 1:nummealcol){

        j<-mealcol[rumkmeal]
        rows.rr1 <- xlsx::readColumns(sheet.rr1, startColumn = j,
                                      endColumn = j, startRow  = 1,  endRow = coldatax+1,as.data.frame=TRUE,header=T)
        nrowsrr1<-dim(rows.rr1)[1]

        'for(i in 1:nrowsrr1){
        vrowsr<-rows.rr1[i,1]
        vrowsr<-as.matrix(vrowsr)
        if(vrowsr=="ERROR")rows.rr1[i,1]<-"#N/A"
      }'

        # read store and remove the data
        block1 <- xlsx::CellBlock(sheet.rr1,startColumn = j, startRow  = 2, noRows=coldatax,
                                  noColumns=1,create=FALSE)
        for(i in 1:coldatax){
          # fill tha background
          xlsx::CB.setFill(block1, fill1, colIndex = 1, rowIndex=i)
          # add the data
          xlsx::CB.setRowData(block1, rows.rr1[i,], i, colOffset=0
                              ,showNA=T, rowStyle=NULL)

        }
    }
      xlsx::saveWorkbook(wb.rr1,pwdxdef)
  }
    if(nummealcol==0){
      xlsx::saveWorkbook(wb.rr1,pwdxdef)
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
    sleeptimepoint<-length(sleepp)###the number of sleep times
    sleepcol<-data.frame()
    for(i in 1:sleeptimepoint){
      sleepti<-sleepp[i]
      sleepti<-as.matrix(sleepti)
      sleepti<-as.numeric(sleepti)
      lowsleepti<-sleepti#the time of sleep
      highsleepti<-sleepti+2# 4 hours after sleep
      for(j in 5:dim(datax)[2]){
        dataxx<- as.matrix(datax[3,j])
        dataxx<-as.numeric(dataxx)
        if(lowsleepti<=dataxx&dataxx<=highsleepti)
          sleepcol<-rbind(sleepcol,j)
      }
    }
    numsleepcol<-dim(sleepcol)[1]# the number of similar metabolites in datax
    sleepcol<-as.matrix(sleepcol)
    wb.rr2 <-xlsx::loadWorkbook(pwdxdef)
    # getSheets get sheet
    sheets.in.rr2 <- xlsx::getSheets(wb.rr2)
    # sheet.rr -first sheet
    sheet.rr2 <- sheets.in.rr2[[1]]
    if(numsleepcol!=0){
      for(rumksleep in 1:numsleepcol){

        j<-sleepcol[rumksleep]
        rows.rr2 <- xlsx::readColumns(sheet.rr2, startColumn = j,
                                      endColumn = j, startRow  = 1,  endRow = coldatax+1,as.data.frame=TRUE,header=T)
        nrowsrr2<-dim(rows.rr2)[1]

        'for(i in 1:nrowsrr1){
        vrowsr<-rows.rr1[i,1]
        vrowsr<-as.matrix(vrowsr)
        if(vrowsr=="ERROR")rows.rr1[i,1]<-"#N/A"
      }'

        # read store and remove the data
        block2 <- xlsx::CellBlock(sheet.rr2,startColumn = j, startRow  = 2, noRows=coldatax,
                                  noColumns=1,create=FALSE)
        for(i in 1:coldatax){
          # fill tha background
          xlsx::CB.setFill(block2, fill2, colIndex = 1, rowIndex=i)
          # add the data
          xlsx::CB.setRowData(block2, rows.rr2[i,], i, colOffset=0
                              ,showNA=T, rowStyle=NULL)

        }
    }
      xlsx::saveWorkbook(wb.rr2,pwdxdef)
}
    if(numsleepcol==0){
      xlsx::saveWorkbook(wb.rr2,pwdxdef)
    }



}
  if(!is.data.frame(design)){
    xlsx::saveWorkbook(wb.rr,pwdxdef)
  }




  ####-----------------A-pre-----------

  datax<-A_pre
  simiid<-simidata
  numsimi<-length(simiid) # the number of input similar metabolites
  rowdatax<-dim(datax)[1]
  similist<-data.frame()# the rows number of the similar metabolites
  for(j in 1:numsimi){
    for(i in 1:rowdatax){
      if(datax[i,2]==simiid[j])
        similist<-rbind(similist,i)
    }
  }
  numSi<-dim(similist)[1]# the number of similar metabolites in datax
  similist<-as.matrix(similist)

  dxdef = paste(dirout, "DifferentialMetabolites(raw).xlsx", sep = "")

  xlsx::write.xlsx(A_pre, dxdef,row.names = F)
  wb.rr <-xlsx::loadWorkbook(dxdef)
  # getSheets get sheet
  sheets.in.rr <- xlsx::getSheets(wb.rr)
  # sheet.rr -first sheet
  sheet.rr <- sheets.in.rr[[1]]

  fill <- xlsx::Fill(foregroundColor="lightblue", backgroundColor="lightblue",
                     pattern="SOLID_FOREGROUND")
  coldatax<-dim(datax)[2]
  if(numSi!=0){
    for(rumk in 1:numSi){
      j<-similist[rumk]
      rows.rr <- xlsx::readRows(sheet.rr, startRow = j+1,
                                endRow = j+1, startColumn = 1,  endColumn = coldatax)
      # read store and remove the data
      block <- xlsx::CellBlock(sheet.rr,startRow = j+1,
                               startColumn=1, noRows=1,
                               noColumns=coldatax,create=FALSE)
      for(i in 1:coldatax){
        # fill tha background
        xlsx::CB.setFill(block, fill, colIndex = i, rowIndex=1)
        # add the data
        xlsx::CB.setColData(block, rows.rr[,i], i, rowOffset=0,
                            showNA=F, colStyle=NULL)
      }
    }
    xlsx::saveWorkbook(wb.rr,dxdef);
  }
  if(numSi==0){

    xlsx::saveWorkbook(wb.rr,dxdef);
  }
  ####study design
  if(is.data.frame(design)){

    wb.rr1 <-xlsx::loadWorkbook(dxdef)
    # getSheets get sheet
    sheets.in.rr1 <- xlsx::getSheets(wb.rr1)
    # sheet.rr -first sheet
    sheet.rr1 <- sheets.in.rr1[[1]]
    fill1 <- xlsx::Fill(foregroundColor="yellow2", backgroundColor="yellow2",
                        pattern="SOLID_FOREGROUND")
    fill2 <- xlsx::Fill(foregroundColor="gray64", backgroundColor="gray64",
                        pattern="SOLID_FOREGROUND")

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
    mealp<-mealp[,-1]### the time of meals
    mealtimepoint<-length(mealp)###the number of meal times


    mealcol<-data.frame()
    for(i in 1:mealtimepoint){
      mealti<-mealp[i]
      mealti<-as.matrix(mealti)
      mealti<-as.numeric(mealti)
      lowmealti<-mealti#the time of meal
      highmealti<-mealti+2# 4 hours after meal
      for(j in 5:dim(datax)[2]){
        dataxx<- as.matrix(datax[3,j])
        dataxx<-as.numeric(dataxx)
        if(lowmealti<=dataxx&dataxx<=highmealti)
          mealcol<-rbind(mealcol,j)
      }
    }
    nummealcol<-dim(mealcol)[1]# the number of similar metabolites in datax
    mealcol<-as.matrix(mealcol)
    coldatax<-dim(datax)[1]#row numbers 48
    if(nummealcol!=0){
      for(rumkmeal in 1:nummealcol){

        j<-mealcol[rumkmeal]
        rows.rr1 <- xlsx::readColumns(sheet.rr1, startColumn = j,
                                      endColumn = j, startRow  = 1,  endRow = coldatax+1,as.data.frame=TRUE,header=T)
        nrowsrr1<-dim(rows.rr1)[1]

        'for(i in 1:nrowsrr1){
        vrowsr<-rows.rr1[i,1]
        vrowsr<-as.matrix(vrowsr)
        if(vrowsr=="ERROR")rows.rr1[i,1]<-"#N/A"
      }'

        # read store and remove the data
        block1 <- xlsx::CellBlock(sheet.rr1,startColumn = j, startRow  = 2, noRows=coldatax,
                                  noColumns=1,create=FALSE)
        for(i in 1:coldatax){
          # fill tha background
          xlsx::CB.setFill(block1, fill1, colIndex = 1, rowIndex=i)
          # add the data
          xlsx::CB.setRowData(block1, rows.rr1[i,], i, colOffset=0
                              ,showNA=T, rowStyle=NULL)

        }
    }
      xlsx::saveWorkbook(wb.rr1,dxdef)
  }
    if(nummealcol==0){
      xlsx::saveWorkbook(wb.rr1,dxdef)
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
    sleeptimepoint<-length(sleepp)###the number of sleep times
    sleepcol<-data.frame()
    for(i in 1:sleeptimepoint){
      sleepti<-sleepp[i]
      sleepti<-as.matrix(sleepti)
      sleepti<-as.numeric(sleepti)
      lowsleepti<-sleepti#the time of sleep
      highsleepti<-sleepti+2# 4 hours after sleep
      for(j in 5:dim(datax)[2]){
        dataxx<- as.matrix(datax[3,j])
        dataxx<-as.numeric(dataxx)
        if(lowsleepti<=dataxx&dataxx<=highsleepti)
          sleepcol<-rbind(sleepcol,j)
      }
    }
    numsleepcol<-dim(sleepcol)[1]# the number of similar metabolites in datax
    sleepcol<-as.matrix(sleepcol)
    wb.rr2 <-xlsx::loadWorkbook(dxdef)
    # getSheets get sheet
    sheets.in.rr2 <- xlsx::getSheets(wb.rr2)
    # sheet.rr -first sheet
    sheet.rr2 <- sheets.in.rr2[[1]]
    if(numsleepcol!=0){
      for(rumksleep in 1:numsleepcol){

        j<-sleepcol[rumksleep]
        rows.rr2 <- xlsx::readColumns(sheet.rr2, startColumn = j,
                                      endColumn = j, startRow  = 1,  endRow = coldatax+1,as.data.frame=TRUE,header=T)
        nrowsrr2<-dim(rows.rr2)[1]

        'for(i in 1:nrowsrr1){
        vrowsr<-rows.rr1[i,1]
        vrowsr<-as.matrix(vrowsr)
        if(vrowsr=="ERROR")rows.rr1[i,1]<-"#N/A"
      }'

        # read store and remove the data
        block2 <- xlsx::CellBlock(sheet.rr2,startColumn = j, startRow  = 2, noRows=coldatax,
                                  noColumns=1,create=FALSE)
        for(i in 1:coldatax){
          # fill tha background
          xlsx::CB.setFill(block2, fill2, colIndex = 1, rowIndex=i)
          # add the data
          xlsx::CB.setRowData(block2, rows.rr2[i,], i, colOffset=0
                              ,showNA=T, rowStyle=NULL)

        }
    }
      xlsx::saveWorkbook(wb.rr2,dxdef)
}
    if(numsleepcol==0){
      xlsx::saveWorkbook(wb.rr2,dxdef)
    }



}
  if(!is.data.frame(design)){
    xlsx::saveWorkbook(wb.rr,dxdef)
  }

  ####-----------output P value------------
  xdef = paste(dirout, "p-value.xlsx", sep = "")
  xlsx::write.xlsx(p_group, xdef)
  xdef = paste(dirout, "p-value(adjusted).xlsx", sep = "")
  xlsx::write.xlsx(p_group_adj, xdef)
  return(Aresult)
}

