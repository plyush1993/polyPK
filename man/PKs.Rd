\name{PKs}
\alias{PKs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the representative pharmacokinetics parameters and plot the time-intensity curves of specified compounds.
}
\description{
A function to calculate the 7 pharmacokinetics parameters (Tmax, Cmax, AUC, CL, Tlast, Tfirst, Cmin) and plot the time-intensity curves for specified compounds.
}
\usage{
PKs(d.pk,d.point="mean",d.ebar="SE",filepath=getwd(),design=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d.pk}{
The data under analysis (data frame with required format)
The first row should be column names. The first and the second column of the first row should be "Name" and "ID", and you can set 2 more tags at the third and the fourth column of the first row, such as "m.z" and "RT.min." or anything you like. From the fifth column till the end, sample indexes or names are expected.  The first row of the data frame should be the gender information."1"means male,and "2" means female.The second row of the data frame should be the group information.The first column of the second row should be "group", and you can add group indexes of the data from the fifth column at the second row. The format of group number should be "0"(pre-dose). "1","2","3","4"...(post-dose). The third row of the data frame should be the information of timepoints.Please see the demo data for detailed format.
}
  \item{d.point}{
The value to calculate the pharmacokinetics parameters,and the value of points in the time-intensity curve. d.point=c ("mean","median"). Default:"mean".
}
  \item{d.ebar}{
The value of error bars. d.ebar=c ("SE","SD"). Defalut:"SE".
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list of metabolites and 7 pharmacokinetics parameters (Tmax,Cmax,AUC,CL,Tlast,Tfirst,Cmin) of specified compound.

Cmax: The peak plasma concentration of a drug after administration.

AUC: Area under the Drug Concentration Curve(0- infinite).

CL: The rate of clear.

Tlast: The last time.

Tfirst: The first time.

Cmin: The least plasma concentration of a drug after administration.
%%  \item{comp2 }{Description of 'comp2'}

A folder named "PKresluts"will also be created automatically which contains three folders:"PKresluts(all)","PKresluts(male)","PKresluts(female)".Each folder has two kinds of files (.xlsx and .PDF).

The file named “PK-parameters.xlsx” contains the pharmacokinetics parameters and the one or more *.PDF files show the time-intensity curves of specified metabolites. Each metabolite has one PDF file named "Time-Intensity-Curve of (its own name).pdf".If the study design is given by right fromat, the time points of meals and sleeps will be described at the bottom of the picture.
}
\references{
1.  Veng Pedersen P. Mean time parameters in pharmacokinetics. Definition, computation and clinical implications.(1898)Clin Pharmacokinet.

2.  Krishnaswami S, Wang T, et al. Single- and multiple-dose pharmacokinetics of tofacitinib in healthy Chinese volunteers.(2015)Clin Pharmacol Drug Dev.
}
\author{
Mengci Li, Shouli Wang, Guoxiang Xie, Tianlu Chen, Wei Jia
}
\note{
%%  ~~further notes~~
nothing
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
nothing
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data("B")
PKs(B,d.point="mean",d.ebar="SE",filepath=getwd(),design=FALSE)
####the result is saved in your current working directory of the R process

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PK parameters}% use one of  RShowDoc("KEYWORDS")
\keyword{ time-intensity curve }% __ONLY ONE__ keyword per line
