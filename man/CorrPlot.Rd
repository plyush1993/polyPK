\name{CorrPlot}
\alias{CorrPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the correlation diagram of two datasets
}
\description{
A function to calculate the correlation coefficients and plot the correlation diagram (8 types) of two input datasets.
}
\usage{
CorrPlot(dataset1,dataset2,cor.method="pearson",filepath,fig.form="heatmap",design)
}

\arguments{
  \item{dataset1}{
The first dataset (data frame with required format).
The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected.  The first row of the data frame should be the gender information."1"means male,and "2" means female.The second row of the data frame should be the group information.The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). The third row of the data frame should be the information of timepoints.Please see the demo data for detailed format.

This variable maybe the results of \link{GetEndo},\link{GetAbso},\link{GetSecdAbso}.
}
  \item{dataset2}{
The second dataset (data frame with required format).
The form of dataset2 is the same as the form of dataset1.
This variable maybe the results of \link{GetEndo},\link{GetAbso},\link{GetSecdAbso}.
}
  \item{cor.method}{
A character string indicating which correlation analysis ("pearson", "kendall", or "spearman") is to be used. Default: "pearson".
}
  \item{filepath}{
A character string indicating the path where the results may be saved in.
}
  \item{fig.form}{
The form of the correlation diagram.  figure.fig.form=c("heatmap","bubble","ordered.bubble","chord","square","ord.square","pie","ord.pie"). Default: "heatmap".
}
  \item{design}{
(optional) a study design dataset(data frame with required format).Use data(design) to see the detailed format.Default:"FALSE"
}
}
\details{
nothing
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A folder named “CorrelationResults” containing three folders:"CorrelationResults(all)","CorrelationResults(female)",and "CorrelationResults(male)". Each folder has three files will be created automatically.

p-value.xlsx: The p values of the correlation analysis.

r-value.xlsx: The r values of the correlation analysis.

correlation-matrix-HeatMap.pdf/correlation-matrix-ChordDiagram.pdf/ (ordered-)correlogram-square.pdf /(ordered-)correlogram-pie.pdf /correlation-matrix(-ordered)-BubbleDiagram.pdf: A PDF file which contains the selected form of correlation diagram.If the study design is given by right fromat, the time points of meals and sleeps will be described at the bottom of the picture.
}
\references{
nothing
%% ~put references to the literature/web site here ~
}
\author{
Mengci Li, Shouli Wang, Guoxiang Xie, Tianlu Chen, Wei Jia
}
\note{
nothing
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{ScatPlot}:plot the PCA or PLS score figures and trajectories on input data.

\link{HeatMap}: plot heatmap of the input data.

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
data("B")
data("C")
CorrPlot(B,C,filepath=getwd(),fig.form="heatmap",design=FALSE)
##----the result is saved in your current working directory of the R process
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ heatmap}% __ONLY ONE__ keyword per line
