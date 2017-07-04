\name{GetDiffData}
\alias{GetDiffData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the differential compounds across all the time points
}
\description{
A function to get all the differential compounds between the pre-dose and every post-dose datasets.
}
\usage{
 GetDiffData(preData, postData, simidata, mv = "mean", rz = 80, multiple
                 = 0.1, sv = TRUE, log = FALSE, t = "Ttest", r.adj =
                 "fdr", filepath = getwd(), design = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{preData}{
The original pre-dose dataset (data frame) with an indicator of gender variable at the first row,grouping variable at the second row,and time points at the third row.
}
  \item{postData}{
The original post-dose dataset (data frame) with an indicator of gender variable at the first row,grouping variable at the second row,and time points at the third row.
}
  \item{simidata}{
The same compounds of drug and pre-dose metabolome data,which is derived from \link{Simi}.
}

  \item{rz}{
    The percentage of zeros for variable elimination(Default:80)
}
  \item{mv}{
   The method of missing values imputation(Default:"mean").mv=c("mean","groupmean","median","groupmedian","groupmin","min","knn","svd","rf","qrilc")
}
  \item{multiple}{
 The parameter for missing values imputation. Missing values will be replaced by multiple*mean/median/min (Default: 0.1).
}
  \item{sv}{
   A logical value indicating whether to remove the outliers (Default:TRUE).
The data which distance to the mean is bigger than 1.5 times of the difference value between lower quartile and upper quartile, should be identified as an outlier. And it will be replaced by the mean value of corresponding row.
}
  \item{log}{
A logical value indicating whether to take the logarithm on the datasets (Default:FALSE)
}
  \item{t}{
  The method for differential compounds identification. C ("Ttest", "MWtest"). Default: "Ttest". Compounds with p values less than 0.05 were taken as differential ones.
}
  \item{r.adj}{
The methods for p values adjustment. r.adj=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Default: "fdr".
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
  \item{A}{A list of all the differential compounds, with preprocessed data.}
  \item{A_pre}{A list of all the differential compounds, with original data.}
  \item{p}{The p values. The dimension of p is i*j, if there are i compounds exist in both pre-dose and post-dose data sets and there are j post-dose time points.}
  \item{p_adj}{The adjusted p values. The dimension of p_adj is the same as p.}

A folder named "DifferentialMetabolites" containing four files of the above 4 datasets will be created automatically.

The name of the file is: "p-value.xlsx", "p-value(adjusted).xlsx", "DifferentialMetabolites(preprocessed).xlsx" and "DifferentialMetabolites(raw).xlsx" respectively.

And the foreground color of the same compounds produced by \link{Simi} will be marked with light blue.If the study design is given by right fromat, the meal times and sleep times will be marked as yellow and grey.
}
\references{
1.	Ke Lan,Wei Jia,et al.An Integrated Metabolomics and Pharmacokinetics Strategy for Multi-Component Drug Evaluation.(2010)Current Drug Metabolism.

2.  Guoxiang Xie,Wei Jia,et al.Metabolic Fate of Tea Polyphenols in Humans.(2012)Journal of Proteome Research.

3.	Ke Lan,Wei Jia,et al.Towards Polypharmacokinetics:Pharmacokinetics of Multicomponent Drug and Herbal Medicines Using a Metabolomics Approach.(2013)Evidence-Based Complementary and Alternative Medicine.

4.	Wei Jia,Tai-ping Fan,et al.The polypharmacokinetics of herbal medicines.(2015) Science.

5.	Guoxiang Xie,Wei Jia,et al.Poly-pharmacokinetic study of a multicomponent herbal medicine in healthy Chinese volunteers.(2017)Clinical Pharmacology&Therapeutics.
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
\link{GetEndo}
\link{GetAbso}
\link{GetSecdAbso}
\link{Simi}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
\dontrun{data("preData")
data("postData")
data("simidata")
data("design")

GetDiffData(preData,postData,simidata,filepath=getwd(),design=FALSE)
}
##the result will be saved in your current working directory of the R process.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ differentian metabolites }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
