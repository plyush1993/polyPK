\name{ScatPlot}
\alias{ScatPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the PCA or PLSDA score figures and trajectories on input data
}
\description{
A function to plot the PCA or PLSDA figures of input data.
}
\usage{
ScatPlot(scat.data,scform="PCA",num.of.cp,filepath,design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scat.data}{
The data under analysis (data frame with required format).
The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected.  The first row of the data frame should be the gender information."1"means male,and "2" means female.The second row of the data frame should be the group information.The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). The third row of the data frame should be the information of timepoints.Please see the demo data for detailed format.
}
  \item{scform}{
  The form of scat plot. scform=c ("PCA","PLSDA"). Default:"PCA".
}
  \item{num.of.cp}{
The number of components to decompose. Default:2.
  }
    \item{filepath}{
A character string indicating the path where the results may be saved in.
}
  \item{design}{
(optional) a study design dataset(data frame with required format).Use data(design) to see the detailed format.Default:"FALSE"
}



}
\details{
nothing
}
\value{
A folder named "PCAresults" or "PLSDAresluts" with three folders :"PCA/PLSDA(all)","PCA/PLSDA(male)",and "PCA/PLSDA(female)". At each folder 4 files will be created automatically.

PCA(PLSDA)-loading.xlsx: The loading values of PCA (PLSDA) analysis.

PCA(PLSDA)-score.xlsx: The score values of PCA (PLSDA) analysis.

PCA(PLSDA)-scorePlot.pdf: A 2 dimensional scores plot of PCA (PLSDA) analysis.If the study design is given by right from, the time points of meals and sleeps will be described at the bottom of the picture.

PCA(PLSDA)-scorePlot(track).pdf: A trajectory plot derived from the PCA (PLSDA) scores plot in which samples of a group will be represented by one point (the center of the group) and will be connected by lines in time ascending order.If the study design is given by right fromat, the time points of meals and sleeps will be described at the bottom of the picture.
}
\references{
nothing
}
\author{
Mengci Li, Shouli Wang, Guoxiang Xie, Tianlu Chen, Wei Jia
}
\note{
nothing
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{CorrPlot}: plot the correlation diagram of two datasets.

\link{HeatMap}: plot heatmap of the input data.
}
\examples{
\dontrun{data("A")
ScatPlot(scat.data=A,scform="PCA",num.of.cp=2,filepath=getwd(),design=FALSE)}
##----the result is saved in your current working directory of the R process
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PCA}% use one of  RShowDoc("KEYWORDS")
\keyword{PLSDA}% __ONLY ONE__ keyword per line
