\name{Simi}
\alias{Simi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Get the same compounds in two datasets
}
\description{
A function which can get the same compounds in two datasets.Especially the same compounds of drug and pre-dose metabolome data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
%%     ~~Describe \code{x} here~~
The pre-dose dataset (data frame with required format).The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected.  The first row of the data frame should be the gender information."1"means male,and "2" means female.The second row of the data frame should be the group information.The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). The third row of the data frame should be the information of timepoints.Please see the demo data for detailed format.

}
  \item{data2}{
The drug constitutes dataset (data frame)
  }
   \item{filepath}{
A character string indicating the path where the results may be saved in.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The results can be an input argument "simidata" of \link{GetDiffData} ,\link{GetEndo},\link{GetAbso},\link{GetSecdAbso}
}
\value{
A list :
 \item{repetitive rates in data1}{The repetitive rates of same compounds in first metabolites dataset}
 \item{repetitive rates in data2}{The repetitive rates of same compounds in second metabolites dataset}
  \item{similar metabolites}{IDs of same compounds(matrix),which can be an input argument of \link{GetDiffData} ,\link{GetEndo},\link{GetAbso},\link{GetSecdAbso}}
%% ...
And a folder named "SimilarData" containing a file named "Similar-data.xlsx" will be created automatically, which is the same compounds dataset.
}

\author{
Mengci Li, Shouli Wang, Guoxiang Xie, Tianlu Chen, Wei Jia
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{GetDiffData}
\link{GetEndo}
\link{GetAbso}
\link{GetSecdAbso}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
\dontrun{data(preData)
data(drugData)
Simi(data1<-preData,data2<-drugData,filepath=getwd())}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
