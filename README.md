## cdcanthro - CDC ANTHROpometry values
GENERATE SEX- AND AGE-STANDARDIZED WEIGHT, HEIGHT, AND BMI METRICS FROM THE CDC GROWTH CHARTS

### Description

Generate z-scores, percentiles, and other metrics for weight, height, and BMI based on the 2000 CDC growth charts, BMI metrics proposed at a 2018 meeting, and extended z-scores and percentiles for children with obesity. Has a single function, 'cdcanthro'. Requires the package data.table to be installed; library(cdcanthro) will also attach data.table.

The BMI metrics included z-scores and percentiles base on the growth charts, along with various newer metrics that more accuately characterize BMIs above the CDC 97th percentile. Note that the output variables, bmiz and bmip, are based on a combination of the LMS-based z-scores for children without obesity and extend bmiz and extended bmip for children with obesity.  The LMS-based z-scores/percentiles are named 'original_bmiz' and 'original_bmip'.

### Installation
Run the following command in R:

install.packages(
   'https://raw.github.com/CDC-DNPAO/CDCAnthro/master/cdcanthro_0.1.1.tar.gz', repos=NULL
 )

### Usage

cdcanthro(data, age = age_in_months, wt = weight_kg, ht = height_cm, bmi = bmi, all = FALSE)

Default for 'all' is FALSE - see Detailts

cdcanthro(data, age_in_months, weight_kg, height_cm, bmi)

#### Do NOT put arguments in quotation marks, such as cdcanthro(data,'age','wt','ht','bmi'). Instead, use cdcanthro(data, age, wt, ht, bmi)

Arguments:

data: data.frame or data.table

age: age in months specified as accuately as possible.

wt: weight (kg).

ht: height (cm).

bmi: BMI, kg/m^2.

### Details
Expects 'sex' to be a variable in the dataset. Can be coded as either 'boys/girls' or 'male/female' or '1/2'.  Character values can be upper or lower case; only the first character is considered.

Age in months should be given as accurately as possible because the function linearly interpolates between ages. If completed number of months is known (e.g., NHANES), add 0.5. If age is in days, divide by 30.4375 so that a child who is 3672 days old would have an age in months of 120.641.

Weight is in kg, and ht is in cm. BMI is kg/m^2.

For additional information on age, see information on agemos at https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm

If all=TRUE, all variables in Wei et al. paper will be output. Default is FALSE

### Return Value

Returns a data.table containing the original data and various weight, height, and BMI metrics. Can convert this to a dataframe with 'setDF(output_data)'.

#### Variables in output:

waz, haz, bmiz: CDC –for-age z-scores for Weight, Height, and BMI

mod_waz, mod_haz, mod_bmiz: modified z-scores

bmip and bmiz: These are based on the LMS-method for children without obesity and the 'extended' method for children with obesity.  See the Wei et al. reference for the 'extended' method which is based on modeling high BMIs as a half-normal distribution.

bmip95: BMI expressed as percentage of 95th percentile, 120 percent is lower threshold for severe obeseity

if 'all = TRUE', then output other BMI metrics describe in Wei et al. paper.  Default is FALSE.  These express BMI as distance or percent distance from the median.  If percent of the median is desired, 100 can be added to the values.

Reference data are the merged LMS data files at https://www.cdc.gov/growthcharts/percentile_data_files.htm

Author(s): David Freedman

### References

Kuczmarski RJ, Ogden CL, Guo SS, Grummer-Strawn LM, Flegal KM, Mei Z, et al. 2000 CDC Growth Charts for the United States: methods and development. Vital and Health Statistics Series 11, Data from the National Health Survey 2002;11:1–190.

Wei R, Ogden CL, Parsons VL, Freedman DS, Hales CM. A method for calculating BMI z-scores and percentiles above the 95th percentile of the CDC growth charts. Annals of Human Biology 2020;47:514–21. https://doi.org/10.1080/03014460.2020.1808065.

Freedman DS, Woo JG, Ogden CL, Xu JH, Cole TJ. Distance and Percent Distance from Median BMI as Alternatives to BMI z-score. Br J Nutr 2019;124:1–8. https://doi.org/10.1017/S0007114519002046.

See Also
https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm

### Examples

data = expand.grid(sex=1:2, agem=120.5, wtk=c(30,60), htc=c(135,144));

data$bmi = data$wtk / (data$htc/100)^2;

data = cdcanthro(data, age=agem, wt=wtk, ht=htc, bmi, default=FALSE); # if default=TRUE then output all variables in Wei et al. paper

round(data,2)

setDF(data) # to convert to a dataframe

OR data = cdcanthro(data, agem, wtk, htc, bmi);

---------------------

nhanes   # NHANES data (2015/16 and 2017/18)

nhanes  = nhanes[!is.na(bmi)]  # exclude subjects with missing wt/ht

nhanes$agemos = nhanes$agemos + 0.5   # because agemos is completed number of months

data = cdcanthro(nhanes, agemos, wt, ht, bmi)
