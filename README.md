<h2 align=center> cdcanthro: CDC ANTHROpometry values </h2>
GENERATE SEX- AND AGE-STANDARDIZED WEIGHT, HEIGHT, AND BMI METRICS FROM THE CDC GROWTH CHARTS (non-obese children) AND FROM THE ‘EXTENDED METHOD’ (children with obesity)


### Description

Generate z-scores, percentiles, and other metrics for weight, height, and BMI based on the 2000 CDC growth charts (Kuczmarski et al., 2002), BMI metrics proposed at a 2018 meeting (Freedman et al., 2019), and extended z-scores and percentiles for children with obesity (Wei et al., 2020). It has a single function, 'cdcanthro'. Requires the package data.table (≥ 1.13) to be installed; library(cdcanthro) will attach data.table.

The BMI metrics included z-scores and percentiles based on the growth charts and newer metrics that more accurately characterize BMIs above the CDC 97th percentile. Note that the output variables - bmiz and bmip - are based on a combination of the LMS-based z-scores (Cole and Green, 1992; Centers for Disease Control and Prevention (CDC), 2022) for children without obesity and extended bmiz and extended bmip for children with obesity. The LMS-based z-scores/percentiles are named 'original_bmiz' and 'original_bmip'.

### Installation
Run the following commands -

install.packages('data.table') # if not already installed

install.packages(
   'https://raw.github.com/CDC-DNPAO/CDCAnthro/master/cdcanthro_0.1.1.tar.gz', type='source', repos=NULL
 )

### Usage

cdcanthro(data, age = age_in_months, wt = weight_kg, ht = height_cm, bmi = bmi, all = FALSE)

The default for 'all' is FALSE - See Details

cdcanthro(data, age_in_months, weight_kg, height_cm, bmi)

#### Do NOT put arguments in quotation marks, such as cdcanthro(data, 'age', 'wt', 'ht', 'bmi'). Instead, use cdcanthro(data, age, wt, ht, bmi)

Arguments:

data: data frame or data.table

age: age in months specified as accurately as possible.

wt: weight (kg).

ht: height (cm).

bmi: BMI, kg/m^2.

### Details
Expects 'sex' to be a variable in the dataset. Can be coded as either 'boys/girls' or 'male/female' or '1/2'.  Character values can be in upper or lower case; only the first character is considered.

Age in months should be given as accurately as possible because the function linearly interpolates between ages. If only the completed number of months is known (e.g., NHANES), add 0.5. If age is in days, divide by 30.4375 so that a child who is 3672 days old would have an age in months of 120.641.

Weight is in kg, and ht is in cm. BMI is kg/m^2.

For additional information on age, see information on agemos at https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm (A SAS Program for the 2000 CDC Growth Charts (ages 0 to <20 years), 2022)

If all=TRUE, all variables in Freedman et al. paper (Freedman et al., 2019) will be output. Will also output the L, M, and S values for each child and the value of sigma for the half-normal distribution. Default is FALSE

The calculation of BMI z-scores for children without obesity is Z = (((BMI / M) ^ L) -1) / (L*S) where BMI is the child’s BMI, L is Box-Cox transformation for normality for the child’s sex and age, M is median, and S is coefficient of variation.  Reference data are the merged LMS data files at https://www.cdc.gov/growthcharts/percentile_data_files.htm (Centers for Disease Control and Prevention (CDC), 2022)

For children with obesity, BMI percentiles are calculated as  90 + 10*pnorm((BMI - p95) / sigma) where p95 is the sex-and age-specific 95th percentile, and sigma is the scale distribution of the half-normal distribution.


### Return Value

Returns a data.table containing the original data and various weight, height, and BMI metrics. Can convert this to a data frame with 'setDF(output_data)'.

#### Variables in output:

waz, haz, bmiz: CDC –for-age z-scores for Weight, Height, and BMI

mod_waz, mod_haz, mod_bmiz: modified z-scores

bmip and bmiz: These are based on the LMS method for children without obesity and the 'extended' method for children with obesity. See the Wei et al. (Wei et al., 2020) for the 'extended' method, which is based on modeling high BMIs as a half-normal distribution.

bmip95: BMI expressed as a percentage of 95th percentile, 120 percent is the lower threshold for severe obesity

If  'all = TRUE', the output contains other BMI metrics described in Freedman et al. paper. The default is FALSE. These express BMI as distance or percent distance from the median. If the percent of the median is desired, 100 can be added to the values.

Author(s): David Freedman

### References

A SAS Program for the 2000 CDC Growth Charts (ages 0 to <20 years) (2022). Available at: https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm (Accessed: 3 October 2022).

Centers for Disease Control and Prevention (CDC) (2022) Percentile data files with LMS values. Available at: http://www.cdc.gov/growthcharts/percentile_data_files.htm (Accessed: 3 October 2022).

Cole, T.J. and Green, P.J. (1992) ‘Smoothing reference centile curves: the LMS method and penalized likelihood.’, Statistics in medicine, 11(10), pp. 1305–19. Available at: https://doi.org/10.1002/sim.4780111005.

Freedman, D.S. et al. (2019) ‘Distance and Percent Distance from Median BMI as Alternatives to BMI z-score’, The British Journal of Nutrition, 124(5 (Special issue on New Anthropometric Cut-offs in Nutrition Research)), pp. 1–8. Available at: https://doi.org/10.1017/S0007114519002046.

Kuczmarski, R.J. et al. (2002) ‘2000 CDC Growth Charts for the United States: methods and development.’, Vital and health statistics. Series 11, Data from the National Health Survey, 11(246), pp. 1–190. Available at: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Citation&list_uids=12043359.

Wei, R. et al. (2020) ‘A method for calculating BMI z-scores and percentiles above the 95th percentile of the CDC growth charts’, Annals of Human Biology, 47(6), pp. 514–521. Available at: https://doi.org/10.1080/03014460.2020.1808065.

### Examples

data = expand.grid(sex=1:2, agem=120.5, wtk=c(30,60), htc=c(135,144));

data$bmi = data$wtk / (data$htc/100)^2;

data = cdcanthro(data, age=agem, wt=wtk, ht=htc, bmi, all=FALSE); # if default=TRUE then output all variables in Wei et al. paper

round(data,2)

setDF(data) # to convert to a data frame

OR data = cdcanthro(data, agem, wtk, htc, bmi);

---------------------

nhanes   # NHANES data (2015/16 and 2017/18)

nhanes  = nhanes[!is.na(bmi)]  # exclude subjects with missing wt/ht

nhanes$agemos = nhanes$agemos + 0.5   # because agemos is completed number of months

data = cdcanthro(nhanes, agemos, wt, ht, bmi, all=TRUE)

