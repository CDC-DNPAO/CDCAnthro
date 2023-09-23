
# to go from bz to bmi: M*(1 +L*S*z))^(1/L);
# to go from extended bz to bmi: p95 + qnorm(((100*pnorm(ext.z) - 90)/10)) * sigma

cc <- function (...) as.character(sys.call()[-1])

set_cols_first <- function (DT, cols, intersection = TRUE) # thanks to hutils
   {
       if (intersection) {
           return(setcolorder(DT, c(intersect(cols, names(DT)),
               setdiff(names(DT), cols))))
       }
       else {
           return(setcolorder(DT, c(cols, setdiff(names(DT), cols))))
       }
   }

cz_score=function(var, l, m, s){ # LMS formula with modified (m) z-scores
      ls=l*s; invl=1/l
      z = (((var/m) ^ l) -1) / (ls) # z-score formula
      sdp2 = (m * (1 + 2*ls) ^ (invl)) - m; # BMI at modified z-score of 2 SDs
      sdm2 = m - (m * (1 - 2*ls) ^ (invl));
      mz=fifelse(var < m, (var - m)/(0.5*sdm2),
                   (var - m)/(sdp2*0.5) )
      list(z, mz)
   }

cdcanthro <- function(data, age=age_in_months,
                      wt=weight_kg, ht=height_cm, bmi=bmi,
                      all=FALSE)
{
   age_in_months <- weight <- height <- seq_ <- sex <- agey <- bz <-
      lwt2 <- mwt2 <- swt2 <- lbmi2 <- mbmi2 <- sbmi2 <- lht2 <- mht2 <- sht2 <-
      lwt1 <- mwt1 <- swt1 <- lbmi1 <- mbmi1 <- sbmi1 <- lht1 <- mht1 <- sht1 <-
      mbmi <- lbmi <- sbmi <- mref <- sref <- denom <- weight_kg <- height_cm <-
      bmiz <- l <- m <- s <- waz <- haz <- z1 <- z0 <- p95 <- bmip <-
      '_AGEMOS1' <- ebp <- ebz <- agemos <- agemos1 <- agemos2 <-
      sexn <- bmi_l <- bmi_s <- bmi_m <- wl <- wm <- wm <- hl <- hm <- hs <-
      bl <- bm <- bs <- ws <- nms <- mod_bmiz <- mod_waz <- mod_haz <-
      hap <- wap <- NULL
   # due to NSE notes in R CMD check

   if (class(data)[1]=='data.frame') setDT(data)

   data[, seq_ := seq_len(.N)]
   dorig <- copy(data)

   nms <- grep('^sex$',names(data),ignore.case = TRUE, value = TRUE)
   if (length(nms) != 1) {
      stop ("A child's sex can named 'sex', 'SEX', or 'Sex'. There can be only
              1 of these 3 variables in your data")
   }
   names(data)[which(names(data)==nms)] <- 'sex'

   data$age <- data[[deparse(substitute(age))]]
   data$wt <- data[[deparse(substitute(wt))]]
   data$ht <- data[[deparse(substitute(ht))]]

   # changes on Mar 14 2023
   if ('bmi' %in% names(data)){
      data$bmi <- data[[deparse(substitute(bmi))]]
   } else {
   data[,bmi:=wt/(ht/100)^2]
   }

   # changed on May 12 2022
   data[,sexn:=toupper(substr(sex,1,1))]
   data[,sexn:=fcase(
      sexn %in% c(1,'B','M'), 1L,
      sexn %in% c(2,'G','F'), 2L
   )]

   data <- data[between(age,24,240) & !(is.na(wt) & is.na(ht)),
                    .(seq_, sexn,age,wt,ht,bmi)];

   # 'ref_data' is CDCref_d.csv,  https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm
   dref <- ref_data[`_AGEMOS1`>23 & denom=='age']
   setnames(dref, tolower(names(dref)))
   setnames(dref, gsub('^_', '', names(dref)))
   setnames(dref,'sex','sexn')

   # values at 240.0 months: https://www.cdc.gov/growthcharts/percentile_data_files.htm
   d20 <- dref[agemos2==240,
               .(sexn,agemos2,lwt2,mwt2,swt2,lbmi2,mbmi2,sbmi2,lht2,mht2,sht2)]
   names(d20) <- gsub('2','',names(d20));

   dref <- dref[,.(sexn,agemos1,lwt1,mwt1,swt1,lbmi1,mbmi1,sbmi1,lht1,mht1,sht1)]
   names(dref) <- gsub('1','',names(dref));

   dref=rbindlist(list(dref,d20))
   dref[sexn==1, ':=' (mref=23.02029, sref=0.13454)] # checked on 7/9/22
   dref[sexn==2, ':=' (mref=21.71700, sref=0.15297)]

   # v=c('sexn','age','wl','wm','ws','bl','bm','bs','hl','hm','hs','mref','sref');
   v=cc(sexn,age,wl,wm,ws,bl,bm,bs,hl,hm,hs,mref,sref);
   setnames(dref,v)

   # interpolate reference data to match each agemos in input data
   if (length(setdiff(data$age,dref$age))>0) {
      uages=unique(data$age); uages
      db <- dref[sexn==1]
           fapp <- function(v,...)approx(db$age,v,xout=uages)$y
           db <- sapply(db[,..v],fapp)
      dg <- dref[sexn==2]
           fapp <- function(v,...)approx(dg$age,v,xout=uages)$y
           dg <- sapply(dg[,..v],fapp)
      dref <- as.data.table(data.frame(rbind(db,dg)))
   }

   du <- unique(data[,.(sexn,age)],by=c('sexn','age'))
   dref <- dref[du, on=c('sexn','age')]

   setkey(data,sexn,age); setkey(dref,sexn,age)
   dt <- dref[data];

   dt[,c('waz', 'mod_waz'):= cz_score(wt, wl, wm, ws)]
   dt[,c('haz', 'mod_haz'):= cz_score(ht, hl, hm, hs)]
   dt[,c('bz', 'mod_bmiz'):= cz_score(bmi, bl, bm, bs)]

   as.data.table(dt);  setnames(dt,c('bl','bm','bs'),c('bmi_l','bmi_m','bmi_s'))
   dt[,c('wl','wm','ws','hl','hm','hs'):=NULL]

   dt[,':=' (
      bmip=100*pnorm(bz),
      p50= bmi_m * (1 + bmi_l*bmi_s*qnorm(0.50))^(1 / bmi_l),
      p85= bmi_m * (1 + bmi_l*bmi_s*qnorm(0.85))^(1 / bmi_l),
      p95= bmi_m * (1 + bmi_l*bmi_s*qnorm(0.95))^(1 / bmi_l),
      p97= bmi_m * (1 + bmi_l*bmi_s*qnorm(0.97))^(1 / bmi_l),
      wap=100*pnorm(waz),  hap=100*pnorm(haz),

     # other BMI metrics -- PMID 31439056
      z1=((bmi/bmi_m) - 1) / bmi_s,  # LMS formula when L=1: ((BMI/M)-1)/S
      z0 = log(bmi/bmi_m) / bmi_s # LMS transformation with L=0, note these end in '0'
     )
     ][,':=' (
      dist_median = z1 * bmi_m * bmi_s, # un-adjusted distance from median with L=1
      adj_dist_median = z1 * sref * mref, # adjusted (to age 20.0 y) dist from median
      perc_median = z1 * 100 * bmi_s, # un-adjusted % from median
      adj_perc_median = z1 * 100*sref, # adjusted % from median
      log_perc_median = z0 * 100 * bmi_s, # un-adjusted % from median with L=0 (log scale)
      adj_log_perc_median = z0 * 100* sref,  # adjusted % from median w L=0 (log scale)
      bmip95=100*(bmi/p95)
   )]

   ## now create Extended z-score for BMI >=95th P
    dt[,':=' (ebz=bz, ebp=bmip, agey=age/12)]
    dt[, sigma:=fifelse(sexn==1, 0.3728 + 0.5196*agey - 0.0091*agey^2,
                               0.8334 + 0.3712*agey - 0.0011*agey^2)]
    dt[bmip>=95, ebp:=90 + 10*pnorm((bmi - p95) / round(sigma,8))]
    # sigma rounded to 8 to agree with NCHS, Craig Hales
    dt[bmip>=95 & ebp/100 < 1, ebz:=qnorm(ebp/100)]
    dt[ebp/100==1, ebz:=8.21] # highest possible value is 8.20945

   x <- cc(agey,mref,sref,sexn,wt,ht,bmi); dt[,(x):=NULL]

   setnames(dt,
      cc(bz,            bmip,          ebp,  ebz),
      cc(original_bmiz, original_bmip, bmip, bmiz)
   )

   v=cc(seq_, bmiz, bmip, waz, wap, haz, hap, p50,
       p95, bmip95, original_bmip, original_bmiz, perc_median,
       mod_bmiz, mod_waz, mod_haz)

   if(all == TRUE){
      v=c(v, 'bmi_l', 'bmi_m', 'bmi_s',  'sigma', 'adj_dist_median', 'dist_median',
          'adj_perc_median', 'log_perc_median', 'adj_log_perc_median')
   }

   dt <- dt[,..v]
   setkey(dt,seq_); setkey(dorig,seq_)
   dtot <- dt[dorig]
   set_cols_first(dtot,names(dorig))
   dtot[,seq_:=NULL]
   return(dtot[])
}

