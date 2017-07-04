\name{postData}
\alias{postData}
\docType{data}
\title{
The post-dose metabolites dataset (data frame)
%%   ~~ data name/kind ... ~~
}
\description{
The example data of post-dose metabolites,which can be an input argument of \link{GetDiffData},\link{GetEndo} and \link{DataPre}.The first row should be column names. The first row of the data frame should be the gender information.The second row of the data frame should be the group information.The format of group number should be"1","2","3","4"...(post-dose).. The third row of the data frame should be the information of timepoints.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("postData")}

\details{
nothing
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
nothing
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
nothing
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(postData)
## maybe str(postData) ; plot(postData) ...
}
\keyword{datasets}
