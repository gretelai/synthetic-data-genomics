# This data structure provides information about each QTL
# analysis. Each analysis is specified by the following 3 parameters:
#
#   pheno     phenotype or trait to map
#   cov       covariates included in regression model
#   outliers  data points to remove according to their residuals
#             after removing linear effects of covariates.
#
# See comments below for additional details on the QTL analyses.
analyses <- list(
    
# MUSCLE AND BONE TRAITS
# ----------------------
# For all five muscle weights (TA, EDL, soleus, plantaris and
# gastrocnemius), we map QTLs conditioning on tibia length
# ("tibia"). For tibia length, we map QTLs conditioned on body weight.
#
# Tibia length explains 12-18% of variance in the muscle weights. The
# rationale for including tibia length as a covariate is bone length
# may somehow regulate muscle weight as well, and we would like to
# isolate the genetic factors that directly regulate development of
# the muscle tissues.
#  
# For bone-mineral density (BMD), we created a binary trait that
# signals "abnormal" BMD. We do not include any covariates when
# mapping QTLs for these traits. Note that body weight is also
# uncorrelated with BMD.
# 
# For all muscle and bone traits, we include a binary indicator for
# round SW16 as a covariate because the mice from this round showed
# substantial deviation in these traits compared to the rest of the
# mice.
TA      = list(pheno="TA",cov=c("SW16","tibia"),
               outliers=function (x) x < (-18)),
EDL     = list(pheno="EDL",cov=c("SW16","tibia"),
               outliers=function (x) x < (-5) | x > 4),
soleus  = list(pheno="soleus",cov=c("SW16","tibia"),
               outliers=function (x) x < (-4) | x > 4),
plant   = list(pheno="plantaris",cov=c("SW16","tibia"),
               outliers=function (x) x < (-9) | x > 8),
gastroc = list(pheno="gastroc",cov=c("SW16","tibia"),
               outliers=function (x) x < (-40) | x > 50),
tibia   = list(pheno="tibia",cov=c("SW6","SW16","sacweight"),
               outliers=function (x) x < (-1.5)),
BMD     = list(pheno="BMD",cov="SW16",outliers=function (x) x > 0.14),
abBMD   = list(pheno="abBMD",cov="SW16",outliers=NULL),

# OTHER PHYSIOLOGICAL TRAITS
# --------------------------
# Body weights bw1, bw2 and bw3 were measured on subsequent days of
# the methamphetamine sensitivity tests, and are highly correlated
# with each other (r^2 = 98%), so it is only necessary to map QTLs for
# one of them. The body weight measurements after sacrifice
# ("sacweight") show a considerable departure in Round SW17, so we
# include a binary indicator for this round as a covariate for
# sacweight. We include age as a covariate for the "bw0" body weight
# because it was measured while the mouse was still growing.
#
# Fasting glucose levels are explained partially by body weight (PVE =
# 6%), so we include body weight as a covariate for fasting glucose
# levels. Rounds SW1 and SW11 showed a considerable departure in
# fasting glucose levels from the other rounds, so we included binary
# indicators for these two rounds as covariates for fasting glucose
# levels.
bw0    = list(pheno="bw0",cov="glucoseage",
              outliers=function (x) x < (-8.5) | x > 8.5),
bw1    = list(pheno="bw1",cov=c("methage","SW17"),
              outliers=function (x) x < (-9) | x > 10),
ppiwt  = list(pheno="PPIweight",cov="SW17",outliers=NULL),
sacwt  = list(pheno="sacweight",cov="SW17",outliers=NULL),
fastgl = list(pheno="fastglucose",cov=c("SW1","SW11","bw0"),
              outliers=function (x) x < (-60) | x > 60),
tail   = list(pheno="taillength",outliers = NULL,
              cov=c("bw0","glucoseage","SW3","SW4","SW19","SW20","SW22",
                    "SW24")),
testis = list(pheno="testisweight",cov="sacweight",
              outliers=function (x) x < (-0.075)),

# FEAR CONDITIONING TRAITS
# ------------------------
# For all fear conditioning traits, the cage used for testing appears
# to have an effect on the phenotype, so we include binary indicators
# for cage as covariates for all FC phenotypes. Further, the FC
# phenotype measurements in Round SW17 show a noticeably different
# distribution in the FC phenotypes from the other rounds, so we
# include a binary indicator for round SW17 as a covariate in all FC
# traits.
#
# These analyses control for proportion of freezing on day 1 during
# exposure to the tone ("AvToneD1"). AvToneD1 explains 11-25% of the
# variance in the Day 2 and Day 3 freezing measures. Note that here we
# can map QTLs for freezing to the altered context on Day 3
# ("AvAltContextD3") as a quantitative trait after conditioning on
# AvToneD1 because the distribution for this trait is no longer quite
# so bimodal, and looks fairly "normal". So there is no need to map
# QTLs for the binary version of this trait.
#
# PreTrainD1 is a very ugly trait with massive box effects and a lot
# of low values, which might have to be removed as outliers. It is
# quite likely that these outliers represent the "deaf" mice that
# might be skewing the whole results. These outliers are present in
# every box, so not a box-specific effect.
d2ctxt     = list(pheno="AvContextD2",outliers=NULL,
                  cov=c("AvToneD1","FCbox1","FCbox2","FCbox3","SW17")),
d3altc     = list(pheno="AvAltContextD3",outliers=function (x) x > 1,
                  cov=c("AvToneD1","FCbox1","FCbox2","FCbox3","SW17")),
d3tone     = list(pheno="AvToneD3",outliers=NULL,
                  cov=c("AvToneD1","FCbox1","FCbox2","FCbox3","SW17")),
pretraind1 = list(pheno="PreTrainD1",outliers=NULL,
               cov=c("FCbox1","FCbox2","FCbox3","SW10","SW16","SW17","SW20")),
d1avtone   = list(pheno="AvToneD1",outliers=NULL,
               cov=c("FCbox1","FCbox2","FCbox3","SW10","SW7","SW14","SW20")),

# METHAMPHETAMINE SENSITIVITY, LOCOMOTOR ACTIVITY AND ANXIETY-LIKE BEHAVIOR 
# -------------------------------------------------------------------------
# We checked all the cages used in these tests to see whether the
# phenotypes measured using any given cage departed noticeably from
# the other cages. Cage #7 consistently has a large effect.
d1dist0to15 = list(pheno="D1totaldist0to15",cov="methcage7",
                   outliers=function (x) x < (-2000) | x > 2200),
d2dist0to15 = list(pheno="D2totaldist0to15",cov="methcage7",
                   outliers=function (x) x < (-2000) | x > 2500),
d3dist0to15 = list(pheno="D3totaldist0to15",cov="methcage7",
                   outliers=function (x) x > 8500),
d1dist0to30 = list(pheno="D1totaldist0to30",cov="methcage7",
                   outliers=function (x) x < (-3500) | x > 4000),
d2dist0to30 = list(pheno="D2totaldist0to30",cov="methcage7",
                   outliers=function (x) x > 4500),
d3dist0to30 = list(pheno="D3totaldist0to30",cov="methcage7",
                   outliers=function (x) x > 20000),
    
d1totdist5  = list(pheno="D1TOTDIST5",cov="methcage7",
                   outliers=function (x) x < (-1000) | x > 1000),
d1totdist10 = list(pheno="D1TOTDIST10",cov="methcage7",
                   outliers=function (x) x < (-750) | x > 750),
d1totdist15 = list(pheno="D1TOTDIST15",cov="methcage7",
                   outliers=function (x) x < (-750) | x > 750),
d1totdist20 = list(pheno="D1TOTDIST20",cov="methcage7",
                   outliers=function (x) x < (-750) | x > 750),
d1totdist25 = list(pheno="D1TOTDIST25",cov="methcage7",
                   outliers=function (x) x < (-750) | x > 750),
d1totdist30 = list(pheno="D1TOTDIST30",cov="methcage7",
                   outliers=function (x) x > 700),
    
d2totdist5  = list(pheno="D2TOTDIST5",cov="methcage7",
                   outliers=function (x) x < (-1250) | x > 1250),
d2totdist10 = list(pheno="D2TOTDIST10",cov="methcage7",
                   outliers=function (x) x > 1000),
d2totdist15 = list(pheno="D2TOTDIST15",cov="methcage7",
                   outliers=function (x) x > 850),
d2totdist20 = list(pheno="D2TOTDIST20",cov="methcage7",
                   outliers=function (x) x > 1000),
d2totdist25 = list(pheno="D2TOTDIST25",cov="methcage7",
                   outliers=NULL),
d2totdist30 = list(pheno="D2TOTDIST30",cov="methcage7",
                   outliers=function (x) x > 900),
    
d3totdist5  = list(pheno="D3TOTDIST5",cov="methcage7",
                   outliers=function (x) x > 2000),
d3totdist10 = list(pheno="D3TOTDIST10",cov="methcage7",
                   outliers=function (x) x > 4000),
d3totdist15 = list(pheno="D3TOTDIST15",cov="methcage7",
                   outliers=function (x) x > 4000),
d3totdist20 = list(pheno="D3TOTDIST20",cov="methcage7",
                   outliers=function (x) x > 4500),
d3totdist25 = list(pheno="D3TOTDIST25",cov="methcage7",
                   outliers=function (x) x > 4500),
d3totdist30 = list(pheno="D3TOTDIST30",cov="methcage7",
                   outliers=function (x) x > 3750),
    
D1ctrtime0to15 = list(pheno="D1ctrtime0to15",cov="methcage7",
                      outliers=function (x) x < (-0.5)),
D2ctrtime0to15 = list(pheno="D2ctrtime0to15",cov="methcage7",
                      outliers=function (x) x < (-0.75)),
D3ctrtime0to15 = list(pheno="D3ctrtime0to15",cov="methcage7",
                      outliers=function (x) x < (-0.75)),
D1ctrtime0to30 = list(pheno="D1ctrtime0to30",cov="methcage7",
                      outliers=function (x) x < (-0.6)),
D2ctrtime0to30 = list(pheno="D2ctrtime0to30",cov="methcage7",
                      outliers=function (x) x < (-0.75)),
D3ctrtime0to30 = list(pheno="D3ctrtime0to30",cov="methcage7",
                      outliers=function (x) x < (-0.85)),

D1hact0to15 = list(pheno="D1hact0to15",cov="methcage7",outliers=NULL),
D2hact0to15 = list(pheno="D2hact0to15",cov="methcage7",outliers=NULL),
D3hact0to15 = list(pheno="D3hact0to15",cov="methcage7",outliers=NULL),
D1hact0to30 = list(pheno="D1hact0to30",cov="methcage7",outliers=NULL),
D2hact0to30 = list(pheno="D2hact0to30",cov="methcage7",outliers=NULL),
D3hact0to30 = list(pheno="D3hact0to30",cov="methcage7",outliers=NULL),

D1vact0to15 = list(pheno="D1vact0to15",cov=c("methcage7","methcage8",
                   "methcage9","methcage10","methcage11","methcage12"),
                   outliers=function (x) x < (-0.85) | x > 0.85),
D2vact0to15 = list(pheno="D2vact0to15",cov=c("methcage7","methcage8",
                   "methcage9","methcage10","methcage11","methcage12"),
                   outliers=function (x) x < (-1) | x > 1),
D3vact0to15 = list(pheno="D3vact0to15",cov=c("methcage7","methcage8",
                   "methcage9","methcage10","methcage11","methcage12"),
                   outliers=function (x) x < (-1.25) | x > 1.25),
D1vact0to30 = list(pheno="D1vact0to30",cov=c("methcage7","methcage8",
                   "methcage9","methcage10","methcage11","methcage12"),
                   outliers=function (x) x < (-1) | x > 1),
D2vact0to30 = list(pheno="D2vact0to30",cov=c("methcage7","methcage8",
                   "methcage9","methcage10","methcage11","methcage12"),
                   outliers=function (x) x < (-1) | x > 1),
D3vact0to30 = list(pheno="D3vact0to30",cov=c("methcage7","methcage8",
                   "methcage9","methcage10","methcage11","methcage12"),
                   outliers=function (x) x > 1.5),

# PREPULSE INHIBITION (PPI) PHENOTYPES
# ------------------------------------
# All boxes appear to have some effect on some of the PPI phenotypes,
# with Box #3 having a particularly large effect on some phenotypes,
# so we include all PPI box indicators as covariates in analysis of the
# PPI phenotypes.
#
# We also map QTLs for habituation to pulses by analyzing the startle
# response during the fourth block of pulse-alone trials against the
# startle response during the first block of pulse-alone trials.
pp3PPIavg  = list(pheno="pp3PPIavg",outliers=function (x) x < (-0.9),
                   cov=c("PPIbox1","PPIbox2","PPIbox3","PPIbox4")),
pp6PPIavg  = list(pheno="pp6PPIavg",outliers=function (x) x < (-1.1),
                  cov=c("PPIbox1","PPIbox2","PPIbox3","PPIbox4")),
pp12PPIavg = list(pheno="pp12PPIavg",outliers=function (x) x < (-1),
                  cov=c("PPIbox1","PPIbox2","PPIbox3","PPIbox4")),
PPIavg     = list(pheno="PPIavg",outliers=function (x) x < (-1),
                  cov=c("PPIbox1","PPIbox2","PPIbox3","PPIbox4")),
PPIstartle = list(pheno="startle",outliers=function (x) x > 250,
                  cov=c("PPIbox1","PPIbox2","PPIbox3","PPIbox4")),
PPIhabit   = list(pheno="p120b4",outliers=function (x) x < (-140) | x > (140),
                  cov=c("p120b1","PPIbox1","PPIbox2","PPIbox3","PPIbox4")))
