#' @title Get the similar metabolites of two input datas
#'
#' @description A function which can get the similar metabolites of two datas.Especially the similar metabolites between drug and pre-dose metabolites.
#'
#' @param data1 The pre-dose dataset (data frame with required format).The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected.  The first row of the data frame should be the gender information."1"means male,and "2" means female.The second row of the data frame should be the group information.The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). The third row of the data frame should be the information of timepoints.Please see the demo data for detailed format.
#'@param data2 The drug constitutes dataset (data frame)
#' @param filepath A character string indicating the path where the results may be saved in.
#'
#' @return a list of repetitive rates
#'
#' @examples Simi(data1<-preData,data2<-drugData,filepath=getwd())
#'
#' @export
Simi<-function(data1,data2,filepath=getwd()){
  data1<-as.data.frame(data1)
  data2<-as.data.frame(data2)
  data1_id<-as.data.frame(data1[-c(1:3),2])
  data2_id<-as.data.frame(data2[,2])
  data1_num<-dim(data1_id)[1]#the number of metabolites in dataset1
  data2_num<-dim(data2_id)[1]#the number of metabolites in dataset2
  simi_id<-sqldf::sqldf('SELECT * FROM [data1_id] intersect SELECT * FROM [data2_id]')
  simi_id<-as.matrix(simi_id)
  simi_num<-length(simi_id)#the number of intersected metabolites
  simi<-data.frame()
  for(j in 1:simi_num){
    for(i in 1:dim(data1)[1]){
      if(data1[i,2]==simi_id[j]){
        simi<-rbind(data1[i,],simi)
      }
    }
  }
  ####-----------output the data in pathway
  dirout = paste(filepath, "/SimilarData", "/", sep = "")
  dir.create(dirout)
  pw = paste(dirout, "Similar-data.xlsx", sep = "")
  xlsx::write.xlsx(simi, pw,row.names = F)


  ###
  rer<-list()
  rer1<-simi_num/data1_num#repetitive rates in dataset1
  rer2<-simi_num/data2_num#repetitive rates in dataset2
  rer[["repetitive rates in data1"]]<-rer1
  rer[["repetitive rates in data2"]]<-rer2
  rer[["similar metabolites"]]<-simi_id
  return(rer)

}
