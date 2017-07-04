#' @title Get the altered endogenous metabolites from the differential compounds
#'
#' @description A function to get the altered endogenous metabolites by similarity analysis on the list of differential compounds and the list of pre-dose compounds.
#'
#' @param pre  The pre-dose dataset (data frame).
#' @param A The differential compounds which is derived from the GetDiffData function.
#' @param sim The parameter (percentage) for similarity analysis. Default: 80.
#' @param filepath A character string indicating the path where the results may be saved in.
#' @param simidata The similar compounds of drug and pre-dose metabolites,which is derived froem the Simi function.
#' @param design A study design dataset(data frame with required format).Use data(StudyDesign) to see the detailed form.
#'
#' @return data.frame
#'
#' @examples GetEndo(pre<-preData,A,simidata,sim=80,filepath=getwd(),design=design)
#'
#' @export
GetEndo<-function(pre,A,simidata,sim=80,filepath=getwd(),design=FALSE){

 # library(plyr)
 # library(sqldf)
  ####----------------leave the metabolites with more than sim% of zeros(both pre and A)------------###

  A<-as.data.frame(A)
  pre<-as.data.frame(pre)
  gender.ge<-A[1,]
  group.ge<-A[2,]
  tp.ge<-A[3,]
  A<-A[-c(1:3),]
  headname_ge<-pre[-c(1,3),c(1:4)]
  pre_data<-pre[-c(1,3),-c(1:4)]
  poz_ge=(100-sim)*(length(pre)-4)/100
  f_ge<-function(x){sum(x==0)}
  numof0_ge<-apply(pre_data,1,f_ge)
  numof0_ge[is.na(numof0_ge)] <- 0
  d2<-cbind(pre_data,numof0_ge)
  d2<-cbind(headname_ge,d2)
  pre_rz<-subset(d2,subset=(numof0_ge<poz_ge))# subset num of zeros<80%
  pre_rz<-pre_rz[,-length(pre_rz)]
  #A reduce 0

  A_m<-A[,-c(1:4)]
  poz_A=(100-sim)*(length(A_m)-4)/100
  A_m<-as.matrix(A_m)
  A_m[is.na(A_m)]<-1
  f<-function(x){sum(x==0)}
  numof0<-apply(A_m,1,f)
  d2_A<-cbind(A,numof0)
  A_rz<-subset(d2_A,subset=(numof0<poz_A))
  A_rz<-A_rz[,-length(A_rz)]
  ####--------------select the metabolites which both exist in "pre" and "A" datasets------------####
  a_ge<-data.frame(pre_rz[,2])
  b_ge<-data.frame(A_rz[,2])
  pboth_ge<-sqldf::sqldf('SELECT * FROM [a_ge] intersect SELECT * FROM [b_ge]')
  pboth_ge<-as.matrix(pboth_ge)
  dr_ge<-data.frame()
  for(j in 1:length(pboth_ge)){

    for(i in 1:dim(A_rz)[1]){
      if(A_rz[i,2]==pboth_ge[j]){
        dr_ge<-rbind(A_rz[i,],dr_ge)
      }
    }
  }
  drbopr_ge<-dr_ge#all data
  drbopr_ge<-rbind(tp.ge,drbopr_ge)
  drbopr_ge<-rbind(group.ge,drbopr_ge)
  drbopr_ge<-rbind(gender.ge,drbopr_ge)
  #drbopr_ge_name<-drbopr_ge[,c(1:4)]# name and ID
  #return(drbopr_ge)

  ####-----------output xlsx type-------------

  datax<-drbopr_ge
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

  dirout = paste(filepath, "/EndogenousMetabolites", "/", sep = "")
  dir.create(dirout)
  wdxdef = paste(dirout, "EndogenousMetabolites.xlsx", sep = "")
  xlsx::write.xlsx(datax, wdxdef,row.names = F)
  wb.rr <-xlsx::loadWorkbook(wdxdef)
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
                            showNA=T, colStyle=NULL)
   }
  }
    xlsx::saveWorkbook(wb.rr,wdxdef)
  }
  if(numSi==0){

    xlsx::saveWorkbook(wb.rr,wdxdef)
  }


if(is.data.frame(design)){

  wb.rr1 <-xlsx::loadWorkbook(wdxdef)
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
  highmealti<-mealti+2# 2 hours after meal
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
  xlsx::saveWorkbook(wb.rr1,wdxdef)
}
  if(nummealcol==0){
    xlsx::saveWorkbook(wb.rr1,wdxdef)
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
wb.rr2 <-xlsx::loadWorkbook(wdxdef)
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
  xlsx::saveWorkbook(wb.rr2,wdxdef)
}
if(numsleepcol==0){
  xlsx::saveWorkbook(wb.rr2,wdxdef)
}



}
if(!is.data.frame(design)){
  xlsx::saveWorkbook(wb.rr,wdxdef)
}

  return(drbopr_ge)
}
