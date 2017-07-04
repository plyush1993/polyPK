\name{DataPre}
\alias{DataPre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocess the input data
}
\description{
Preprocess the input data. Variables with a lot of zeros and outliers may be removed. Missing values may be imputed and filled by various methods. Data may be transformed by logarithm transformation.
}
\usage{
DataPre(tes,mv="mean",rz=80,multiple=0.1,sv=TRUE,log=FALSE,filepath=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tes}{
  The data under pretreatment (data frame with required format).
The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected.  The first row of the data frame should be the gender information."1"means male,and "2" means female.The second row of the data frame should be the group information.The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). The third row of the data frame should be the information of timepoints.Please see the demo data for detailed format.
}

  \item{rz}{
  The percentage of zeros for variable elimination (Default:80). Variables with zero numbers higher than rz% will be removed.
}
  \item{mv}{
   The method of missing values imputation (Default:"mean"). mv=c ("mean", "groupmean", "median", "groupmedian", "groupmin", "min", "knn", "svd", "rf", "qrilc").
}
  \item{multiple}{
   The parameter for missing values imputation. Missing values will be replaced by multiple*mean/median/min (Default:0.1).
}
  \item{sv}{
  A logical value indicating whether to remove the outliers (Default:TRUE).
The data which distance to the mean is bigger than 1.5 times of the difference value between lower quartile and upper quartile, should be identified as an outlier. And it will be replaced by the mean value of corresponding row.
}
  \item{log}{
   A logical value indicating whether to take the logarithm on the data (Default:FALSE).
}
  \item{filepath}{
  A character string indicating the path where the results may be saved in.
}


}
\details{
nothing
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame of the prepocessed data

A folder named "preprocessed-data" containing a file of the prepocessed datasets will be created automatically.The file's name is"preprocessed-data.xlsx".
}
\references{
Hastie,Botstein,et al.Imputing Missing Data for Gene Expression Arrays, Stanford University Statistics Department Technical report (1999)
}
\author{
Mengci Li,Shouli Wang, Guoxiang Xie, Tianlu Chen, Wei Jia
}
\note{
nothing
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
nothing
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("preData")
DataPre(preData)
##The result will be saved at your current working directory of the R process.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ preprocess }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
