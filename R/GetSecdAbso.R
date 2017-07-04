#' @title Get the secondary metabolites of the absorbed drug constitutes
#'
#' @description A function to get secondary metabolites of the absorbed drug constitutes.
#'
#' @param A The differential compounds dataset which is derived from the GetDiffData function.
#' @param B The altered endogenous metabolites dataset which is derived from the GetEndo function.
#' @param C The absorbed drug constitutes dataset which is derived from the GetAbso function.
#' @param sim The parameter (percentage) for similarity analysis. Default: 80.
#' @param filepath A character string indicating the path where the results may be saved in.
#' @param simidata The similar compounds of drug and pre-dose metabolites,which is derived froem the Simi function.
#' @param design A study design dataset(data frame with required format).Use data(StudyDesign) to see the detailed form.
#'
#' @return data.frame
#'
#' @examples GetSecdAbso(A,B,C,simidata,sim=80,filepath=getwd(),design=design)
#' @export
GetSecdAbso<-function(A,B,C,simidata,sim=80,filepath=getwd(),design=FALSE){
  #library(plyr)
  #library(sqldf)
  A<-as.data.frame(A)
  B<-as.data.frame(B)
  C<-as.data.frame(C)
  A_m<-A
  tp.ge<-A[3,]
  gender.ge<-A[1,]
  group.ge<-A[2,]#get the group information of data frame A
  A<-A[-c(1:3),]
  B<-B[-c(1:3),]
  C<-C[-c(1:3),]
  a<-data.frame(A[,2])
  b<-data.frame(B[,2])
  c<-data.frame(C[,2])
  a_b<-sqldf::sqldf('SELECT * FROM [a] EXCEPT SELECT * FROM [b]')#pre have;pos not have
  both_id<-sqldf::sqldf('SELECT * FROM [a_b] EXCEPT SELECT * FROM [c]')#pos have
  both_id<-as.matrix(both_id)
  abso<-data.frame()
  for(j in 1:length(both_id)){

    for(i in 1:dim(A)[1]){
      if(A[i,2]==both_id[j]){
        abso<-rbind(A[i,],abso)
      }
    }
  }
  poz_A=(100-sim)*(length(A_m)-4)/100
  abso_data<-as.matrix(abso[,-c(1:4)])
  abso_data[is.na(abso_data)]<-1
  f<-function(x){sum(x==0)}
  numof0<-apply(abso_data,1,f)
  d2_A<-cbind(abso,numof0)
  D_rz<-subset(d2_A,subset=(numof0<poz_A))
  D_rz<-D_rz[,-length(D_rz)]
  getabso<-D_rz#data
  getabso<-rbind(tp.ge,getabso)
  getabso<-rbind(group.ge,getabso)
  getabso<-rbind(gender.ge,getabso)


  ####-----------output xlsx type-------------

  datax<-getabso
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

  dirout = paste(filepath, "/SecondAbsorbedMetabolites", "/", sep = "")
  dir.create(dirout)
  dxdef = paste(dirout, "SecondAbsorbedMetabolites.xlsx", sep = "")
  xlsx::write.xlsx(datax, dxdef,row.names = F)
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





  return(getabso)
}
