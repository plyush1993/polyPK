\name{GetSecdAbso}
\alias{GetSecdAbso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the secondary metabolites of the absorbed drug constitutes
}
\description{
A function to get secondary metabolites of the absorbed drug constitutes.
}
\usage{
GetSecdAbso(A,B,C,simidata,sim=80,filepath,design)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{A}{
The differential compounds dataset which is derived from the \link{GetDiffData} function.
}
  \item{B}{
The altered endogenous metabolites dataset which is derived from the \link{GetEndo} function.
}
  \item{C}{
The absorbed drug constitutes dataset which is derived from the \link{GetAbso} function.
}
  \item{simidata}{
The same compounds of drug and pre-dose metabolome data, which is derived from \link{Simi}.
}
  \item{sim}{
The parameter (percentage) for similarity analysis. Default: 80.
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
A folder named "SecondAbsorbedMetabolites" containing a file named "SecondAbsorbedMetabolites.xlsx" will be created automatically which is the list and data of secondary metabolites of the absorbed drug constitutes. And the foreground color of the same compounds produced by \link{Simi} will be marked with light blue. If the study design is given by right from, the meal times and sleep times will be marked as yellow and grey.
}
\references{
1.	Ke Lan, Wei Jia, et al. An Integrated Metabolomics and Pharmacokinetics Strategy for Multi-Component Drug Evaluation. (2010) Current Drug Metabolism.

2.  Guoxiang Xie, Wei Jia, et al. Metabolic Fate of Tea Polyphenols in Humans. (2012) Journal of Proteome Research.

3.	Ke Lan, Wei Jia, et al. Towards Polypharmacokinetics: Pharmacokinetics of Multicomponent Drug and Herbal Medicines Using a Metabolomics Approach. (2013) Evidence-Based Complementary and Alternative Medicine.

4.	Wei Jia, Tai-ping Fan, et al. The polypharmacokinetics of herbal medicines. (2015) Science.

5.	Guoxiang Xie, Wei Jia, et al. Poly-pharmacokinetic study of a multicomponent herbal medicine in healthy Chinese volunteers. (2017) Clinical Pharmacology&Therapeutics.
%% ~put references to the literature/web site here ~
}
\author{
Mengci Li, Shouli Wang, Guoxiang Xie, Tianlu Chen and Wei Jia
}
\note{
The list of absorbed drug compounds was obtained by excluding compounds in \link{B} and \link{C} from that of \link{A}.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{GetDiffData}
\link{GetAbso}
\link{GetEndo}
\link{Simi}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data("A")
data("B")
data("C")
data("simidata")
data("design")
GetSecdAbso(A,B,C,simidata,sim=80,filepath=getwd(),design=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ second absorbed metabolites }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
