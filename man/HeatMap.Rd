 \name{HeatMap}
\alias{HeatMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the heatmap of input data
}
\description{
A function to plot the heatmap and clusters of input data.
}
\usage{
HeatMap(data,cluster="both",scale="row",filepath,design)
}

\arguments{
  \item{data}{
The data under analysis (data.frame with required format) The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected.  The first row of the data frame should be the gender information. "1" means male, and "2" means female. The second row of the data frame should be the group information. The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). The third row of the data frame should be the information of timepoints. Please see the demo data for detailed format.
}
  \item{cluster}{
A string indicating whether or in which direction the dendrograms should be drawn (“none”, “row”, “column” or “both”). Default:”both”.
}
  \item{scale}{
A character indicating whether the data should be centered and scaled before analysis and in which (“none”, “row” or “column”) direction. Default:"row".
}
  \item{filepath}{
A character string indicating the path where the results may be saved in.
}
\item{design}{
(optional) a study design dataset(data frame with required format).Use data(design) to see the detailed format. Default:"FALSE".
}
}
\details{
nothing
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A folder named “HeatMap” containing a PDF file (the heatmap figure) will be created automatically. If the study design is given by right format, the time points of meals and sleeps will be described at the bottom of the picture.
}
\references{
nothing
%% ~put references to the literature/web site here ~
}
\author{
Mengci Li, Shouli Wang, Guoxiang Xie, Tianlu Chen and Wei Jia
}
\note{
nothing
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{CorrPlot}: plot the correlation diagram of two datasets.

\link{ScatPlot}: plot the PCA or PLS score figures and trajectories on input data.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data("B")
HeatMap(data=B,cluster="both",scale="row",filepath=getwd(),design=FALSE)
##the result is saved in your current working directory of the R process
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ heatmap}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
