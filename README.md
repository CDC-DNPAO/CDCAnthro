<h2 align=center> cdcanthro: CDC ANTHROpometry values </h2>
GENERATE SEX- AND AGE-STANDARDIZED WEIGHT, HEIGHT, AND BMI METRICS FROM THE CDC GROWTH CHARTS (non-obese children) AND FROM THE ‘EXTENDED METHOD’ (children with obesity)

<h3 align=center> This function is part of the cdcanthro package on CRAN: https://cran.r-project.org/package=cdcanthro </h3>

### News
August 31, 2026. See the function cdc_z in the cdcanthro package on CRAN for more information

Sept 23, 2023. The function now recognizes 'sex', 'SEX' or 'Sex' as names for the child's sex.  It also provides informative errors if one of these variables is not in your data.  Also, the output is limited to children whose age in months is >= 24 and < 240; the CDC growth charts apply only to these ages. 

### Description

Generate z-scores, percentiles, and other metrics for weight, height, and BMI based on the 2000 CDC growth charts (Kuczmarski et al., 2002), BMI metrics proposed at a 2018 meeting (Freedman et al., 2019), and extended z-scores and percentiles for children with obesity (Wei et al., 2020). It has a single function, 'cdcanthro'. Requires the package data.table (≥ 1.13) to be installed; library(cdcanthro) will attach data.table.

The BMI metrics included z-scores and percentiles based on the growth charts and newer metrics that more accurately characterize BMIs above the CDC 97th percentile. Note that the output variables - bmiz and bmip - are based on a combination of the LMS-based z-scores (Cole and Green, 1992; Centers for Disease Control and Prevention (CDC), 2022) for children without obesity and extended bmiz and extended bmip for children with obesity. The LMS-based z-scores/percentiles are named 'original_bmiz' and 'original_bmip'.

The calculations and output variables are similar to those in the SAS program at https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm.  However, using all = 'TRUE' in the function will output additional BMI metrics.


### Installation
The package is now on CRAN, so install.packages('cdcanthro')

