###########################################################
###  Master file for the UK CVD microsimulation model   ###
###     for parameters input and condition setting      ###
###########################################################
# preparation -------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(parallel)
library(doSNOW)
library(data.table)

# workding directory
wd <- "Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/Final_Model_Code_Archive_Folder"

# run function file
source(file.path(wd, "scripts/simulation/fun_sim_boots.R"))

# input data directory
.input_data_dir <- file.path(wd, "input")

# output directory
output_dir <- file.path(wd, "output")

# scenario control --------------------------------------------------------

# tx (treatment effect): FALSE for no tx; otherwise the name of the tx file
# regimen: "none" for no treatment
## simva, atorva, rosuva, prava, fluva
## 5mg, 10mg, 20mg, 40mg, 80mg
## ex. regimen could be "atorva 40mg", "none", or "atorva 40mg + ezetimibe"
## note: if no tx is specified, regimen setting is ignored, 
# side_effects_flag: the general switch of side effect. TRUE = on; FALSE = off
## if side_effects_flag=FALSE, the following two side effect switches are ignored
# side_dm: side effect of statin on diabetes on/off; on is recommended 
# side_mus: side effect of statin on muscle diseases on/off; on is recommended
# adj_cancer: adjustment of cancer incidence after 80 years of age on/off; on is recommended
# pop: simulate primary "prim" or secondary "sec" prevention population, using different coefficients 
# n_sim: number of 1st order simulation; n_sim>500 is recommended
# id_list: "all", or input a particular ID list file name
# delay_age: start statin since age; 0 means no delay
# stop_age: stop statin after age; 0 means no stop
# delay_45: under 45 at entry delay statin for 5 years
# adh_ctrl: whether apply real-world (UK) adherence rate
# incr_tx: whether apply increasing tx

# define scenario
snr <- c(tx="tx_hpc", regimen="atorva 80mg", 
         side_effects_flag=T, side_dm=T, side_mus=T, 
         adj_cancer = T,
         pop=c("sec"), n_sim=2, id_list="all", 
         delay_age = 0, stop_age = 0, delay_45 = F, adh_ctrl = F, 
         incr_tx = F
)

# if your input include multiple sets of coefficients for probabilistic sensitivity analysis
# input the sequence number here to indicate which set of coefficients you want to use
# here 1 represent the deterministic analysis
rnum <- 1

# Input file control ------------------------------------------------------

# output file name
output_filename_prefix <- "Result" 

# input coefficients file
cf_filename <- "cf_hpc_v4_VDppUpdate"

# input baseline and time-varying data name prefixes
mx_b_filename_prefix <- "sample_b"
mx_t_filename_prefix <- "sample_t"

# all events
events_list <- list(
  events_nf = c ("mi", "stroke", "crv", "cancer_icd", "dm"),
  events_f = c("vd", "nvd"))

# input of the model parameter file
pf_filename <- "pf"

# simulation parameter control --------------------------------------------

# distribution - the final choice
# there are three options: Exponential, Weibull and Gompertz
# We recommend the following distributions
dist_list <- list( 
  prim = list(
    mi = "wei", 
    stroke = "wei", 
    crv = "wei",
    cancer_icd = "gom", 
    dm = "wei",  
    vd = "gom",
    nvd = "gom"), 
  sec = list(
    mi = "exp",
    stroke = "exp",
    crv = "gom", 
    cancer_icd = "gom", 
    dm = "wei",
    vd = "gom", 
    nvd = "gom")
)

# nonlinear age (TRUE/FALSE)
nonlinage <- TRUE # for QoL only

# use re-calibrated equation (TRUE/FALSE)
# TRUE means using UK Biobank calibrated coefficients; otherwise using CTT coefficents
# coefficients for incident diabetes are derived from UK Biobank only
calibrated_eqns <- TRUE

# adjust CRV-MI order, as majority of CRV happened after MI
adjust_crv = TRUE
# the three options above are recommended as default 

# just check if adjust cancer in elderly
adj_cancer = as.logical(snr["adj_cancer"])

# # Length of simulation
# lifetime # maximum 110 years old
stop_expr <- expression(floor(111 - (v_j["CurrAge_cent"]*10 + 60)))
# Or you can specify a length using the two lines below
# how_long = 10 # e.g. 10 years
# stop_expr <- expr(how_long + 1)

# Should saving be done overall (FALSE) or on patient-level (TRUE)
save_by_pat <- F
# Should patients be sampled
# Could be FALSE (for complete sample) or an integer
sample_pat <- F

# use pre-selected groups of individuals?
id_list <- as.character(snr["id_list"])

# number of simulations
n_sim <- as.numeric(snr["n_sim"])

# number of cores
n_cores <- 8 # round(detectCores() * 4/5)
# HPC environment input
# n_cores <- as.numeric(Sys.getenv('NSLOTS'))

# populations to be simulated 
prim_flag <- snr["pop"]

# treatment effects control -----------------------------------------------

# filename or FALSE
# if tx != FALSE, baseline LDL will be replaced by pre-treated LDL in simulation
tx <- snr["tx"]

if (tx == "FALSE") tx <- as.logical(tx)

# use the simplified version
regimen <- snr["regimen"]

# side effect flag: TRUE or FLASE
# if regimen is none, side_effects_flag is forced to be FALSE in simulation

# whether add side effects on incident diabetes
side_effects_flag <- as.logical(snr["side_effects_flag"])

side_dm <- as.logical(snr["side_dm"])

side_mus <- as.logical(snr["side_mus"])

delay_age <- as.numeric(snr["delay_age"])

stop_age <- as.numeric(snr["stop_age"])

delay_45 <- as.logical(snr["delay_45"])

adh_ctrl <- as.logical(snr["adh_ctrl"])

incr_tx <- as.logical(snr["incr_tx"])

# simulation command ------------------------------------------------------

# distributions
dist <- dist_list[[prim_flag]]

# output filename parameters

ptm <- proc.time()

# run the master function

alpha <- master(
  rnum = rnum,
  .input_data_dir = .input_data_dir, 
  cf_filename = cf_filename, 
  prim_flag = prim_flag,
  mx_t_filename_prefix = mx_t_filename_prefix,
  mx_b_filename_prefix = mx_b_filename_prefix,
  events_list = events_list,
  pf_filename = pf_filename,
  adjust_crv = adjust_crv,
  calibrated_eqns = calibrated_eqns,
  nonlinage = nonlinage,
  tx = tx,
  dist = dist,
  stop_expr = stop_expr, 
  sample_pat = sample_pat,
  save_by_pat = save_by_pat,
  n_sim = n_sim, 
  n_cores = n_cores, 
  output_dir = output_dir, 
  output_filename_prefix = output_filename_prefix, 
  regimen = regimen,
  side_effects_flag = side_effects_flag,
  side_dm = side_dm,
  side_mus = side_mus,
  id_list = id_list,
  adj_cancer = adj_cancer,
  delay_age = delay_age,
  stop_age = stop_age,
  delay_45 = delay_45, 
  adh_ctrl = adh_ctrl,
  incr_tx = incr_tx
)

print(proc.time() - ptm)
print(Sys.time())





