\name{cdcanthro}

\alias{cdcanthro}
\alias{age_in_months}
\alias{weight_kg}
\alias{height_cm}
\alias{ref_data}
\alias{nhanes}
\alias{data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GENERATE SEX- AND AGE-STANDARDIZED WEIGHT, HEIGHT, AND BMI METRICS FROM THE CDC GROWTH CHARTS
}
\description{
Generate z-scores, percentiles, and other metrics for weight, height, and BMI based on the 2000 CDC growth charts.
Has a single function, 'cdcanthro'.  Requires the package data.table to be
installed; library(cdcanthro) will also attach data.table.

The BMI metrics included z-scores and percentiles base on the growth charts, along with various newer metrics such as extended BMIz,
percent of the 50th and 95th percentiles.
}
\usage{
cdcanthro(data,
   age = age_in_months,
   wt = weight_kg,
   ht = height_cm,
   bmi = bmi,
   all = FALSE)
# Default for 'all' is FALSE - see Detailts
# OR
# cdcanthro(data, age_in_months, weight_kg, height_cm, bmi)

# Do NOT put arguments in quotation marks


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data.frame or data.table}
  \item{age}{age in months specified as accuately as possible.}
  \item{wt}{weight (kg).}
  \item{ht}{height (cm).}
  \item{bmi}{BMI, kg/m^2.}
}
\details{
Expects 'sex' to be a variable in the dataset. Can be coded as either 'boys/girls' or 'male/female' or '1/2'.  Character values can be upper or lower case; only the first character is considered.

Weight is in kg, and ht is in cm.  BMI is kg/m^2.

Age in months should be given as accurately as possible because the function linearly interpolates between ages.  If completed number of months is known (e.g., NHANES), add 0.5.

If age is in days, divide by 30.4375 so that a child who is 3672 days old would have an age in months of 120.641.

For additional information on age, see information on agemos at
\url{https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm}

If all=TRUE, all variables in Wei et al. paper will be output. Default is FALSE

}
\value{
Returns a data.table containing the original data and various weight, height, and BMI metrics.  Can convert this to a dataframe with 'setDF(output_data)'.

Variables in output:

waz, haz, bmiz: CDC --for-age z-scores for Weight, Height, and BMI. BMIz is based on 2000 CDC growth charts (non-obese children) and extended BMIz (obese children)

mod_waz, mod_haz, mod_bmiz: modified z-scores

ext_bmip and ext_bmiz: extended BMI percentile and z-score.  See note to BMIz

pre_2022_bmiz and pre_2022_bmip:  orignal calculations of BMIz and BMI percentile

bmip95: BMI expressed as percentage of 95th percentile, 120 percent is lower threshold for severe obeseity

if 'all = TRUE', then output other BMI metrics describe in Wei et al. paper.  Default is FALSE.  These express BMI as distance or percent distance from the median.  If percent of the median is desired, 100 can be added to the values.

}
\references{
Kuczmarski RJ, Ogden CL, Guo SS, Grummer-Strawn LM, Flegal KM, Mei Z, et al. 2000 CDC Growth Charts for the United States: methods and development. Vital and Health Statistics Series 11, Data from the National Health Survey 2002;11:1–190.

Wei R, Ogden CL, Parsons VL, Freedman DS, Hales CM. A method for calculating BMI z-scores and percentiles above the 95th percentile of the CDC growth charts. Annals of Human Biology 2020;47:514–21. \url{https://doi.org/10.1080/03014460.2020.1808065}.

Freedman DS, Woo JG, Ogden CL, Xu JH, Cole TJ. Distance and Percent Distance from Median BMI as Alternatives to BMI z-score. Br J Nutr 2019;124:1–8.
\url{https://doi.org/10.1017/S0007114519002046}.


}
\author{
David Freedman
}
\note{
Do NOT put arguments in quotation marks, such as cdcanthro(data,'age','wt','ht','bmi').  Use: cdcanthro(data, age, wt, ht, bmi)

Reference data are the merged LMS data files at
\url{https://www.cdc.gov/growthcharts/percentile_data_files.htm}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm}
}
\examples{
data = expand.grid(sex=1:2, agem=120.5, wtk=c(30,60), htc=c(135,144));
data$bmi = data$wtk / (data$htc/100)^2;
data = cdcanthro(data, age=agem, wt=wtk, ht=htc, bmi);
# OR data = cdcanthro(data, agem, wtk, htc, bmi);
round(data,2)
# setDF(data) to convert to a dataframe

nhanes   # NHANES data (2015/16 and 2017/18)
nhanes  = nhanes[!is.na(bmi),] # exclude subjects with missing wt/ht
nhanes$agemos = nhanes$agemos + 0.5   # because agemos is completed number of months
data = cdcanthro(nhanes, age=agemos, wt, ht, bmi, all=TRUE)
# OR data = cdcanthro(nhanes, agemos, wt, ht, bmi, all=TRUE)
round(data, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
