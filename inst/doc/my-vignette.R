## ----message = FALSE-----------------------------------------------------
library(polyPK)

## ----results = "hide",warning = FALSE,message = FALSE--------------------
data("postData")
pred_post<-polyPK::DataPre(tes=postData,mv="min",rz=80,sv=TRUE,log=FALSE,filepath=getwd())

## ------------------------------------------------------------------------
pred_post[c(1:10),c(1:10)]

## ----warning = FALSE,message = FALSE-------------------------------------
data("preData")
data("drugData")
Simi(data1<-preData,data2<-drugData,filepath=getwd())


## ----warning = FALSE,message = FALSE,results = "hide"--------------------

data("preData")
data("postData")
data("design")
data("simidata")
dif<-GetDiffData(preData,postData,simidata,mv="min",rz=80,sv=TRUE,log=FALSE,t="Ttest",r.adj="fdr",filepath=getwd(),design=FALSE)
#'


## ------------------------------------------------------------------------
prepoA<-dif$A
as.data.frame(prepoA)

## ------------------------------------------------------------------------
orgA<-dif$A_pre
as.data.frame(orgA)

## ------------------------------------------------------------------------
p<-dif$p
p

## ------------------------------------------------------------------------
padj<-dif$p_adj
padj

## ----warning = FALSE,message = FALSE-------------------------------------
 data("preData")
 data("A")
 data("design")
 data("simidata")
GetEndo(preData,A,simidata,sim=80,filepath=getwd(),design=FALSE)


## ----warning = FALSE,message = FALSE-------------------------------------
GetAbso(drugData, A, simidata,sim = 80, filepath=getwd(),design = FALSE)


## ----warning = FALSE,message = FALSE-------------------------------------
GetSecdAbso(A,B,C,simidata,sim=80,filepath=getwd(),design=FALSE)


## ----warning = FALSE,message = FALSE-------------------------------------
data("C")
data("design")
pks<-PKs(C,d.point="mean",d.ebar="SE",filepath=getwd(),design=FALSE)
knitr::kable(pks[c(1:10),],align = 'c')

## ----warning = FALSE,message = FALSE,results = "hide"--------------------
data("B")
data("C")
CorrPlot(dataset1=B,dataset2=C,cor.method="pearson",filepath=getwd(),fig.form="heatmap",design = FALSE)

## ----warning = FALSE,message = FALSE,results = "hide"--------------------
data("A")
ScatPlot(scat.data=A,scform="PCA",num.of.cp=2,filepath=getwd(),design=FALSE)

## ----warning = FALSE,message = FALSE,results = "hide"--------------------
data("A")
HeatMap(data=A,cluster="both",scale="row",filepath=getwd(),design=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  help(package = 'polyPK', help_type = 'html')
#  # or see a standalone list of vignettes
#  browseVignettes('polyPK')

