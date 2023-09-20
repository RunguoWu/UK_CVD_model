# UK_CVD_model
The folder contains the R codes for the UK CVD micro-simulation model and input files
"master_hpc_boots.R" is the master file to control input/output of the model
To run this model, you need to specify path, input and other parameters in the master file and run it
"fun_sim_boots.R" are the main model functions that the master file calls

In the "input" folder, "cf_hpc_v4_VDppUpdate.rds", "tx_hpc.rds" and "pf.rds" are equation coefficients and model parameters that our research project generated. They were derived from individual patient data from trials of Cholesterol Treatment Trialists' Collaboration (CTT), and calibrated into a UK population of 500,000 individuals, i.e. UK Biobank.  

Files in the input folder with the prefix "sample_" are sample datasets of five pseudo individuals' characteristics, for primary prevention (without cardiovascular disease history) and secondary prevention (with cardiovascular disease history) separately. The model reads in baseline (b) and time-varying (t) characteristics in two separate files, in order to speed up simulation. Please specify your input characteristics in the same format as the samples

## Baseline fixed variables

id_new: patient ID

Intercept: all equal 1

### 0/1 binary variables
male: 1 = male; 0 = female

race_afro, race_sa and race_otherna: African, South Asain and other ethnicities, respectively. Otherwise White ethnicity

unhealthy_diet: 1 = unhealthy diet; 0 = healthy diet. Please refer to our published paper for more details

h_mi_raw_only, cvd_only, h_pad_raw_only, othchd_only and number_events2: History of myocardial infarction (MI), history of cerebral vascular disease (stroke), history of peripheral artery disease (PAD), history of other coronary heart diseases, history of more than one of these diseases, respectively; In primary prevention population, all equal 0; in secondary prevention population, one of them equals to 1  

dmT1: history of Type 1 diabetes

underweight, overweight, obese, obese2 and obese3: BMI(kg/m2) < 18.5, 25-30, 30-35, 35-40 and 40+, respectively. Otherwise BMI 18.5-25

Txhypen: history of treated hypertension

smk_ex; smk_cur: ex-smoker and current smoker, respectively. Otherwise non-smoker

town1, town2, town4 and town5: indicate the corresponding socioeconomic deprivation quintile. town1=1 means the least deprived. 

severe_mental_illness: the history of severe mental illness

pa_low, pa_high and pa_mis: indicates the corresponding level of physical activity; pa_mis = 1 means missing value for physical activity level. 

### Continuous variables
lnbcreann: logarithm of creatinine (umol/L)

sys_bpS: systolic blood pressure (mmHg) centered at 140 divided by 20

dias_bpS: diastolic blood pressure (mmHg) centered at 80 divided by 10

lhdl: logarithm of high-density lipoprotein (mmol/L)

NewB_LDL-CL_cent: baseline measure of LDL cholesterol, centered at 3.6 (mmol/L)

LDL_nostatin_cent: baseline measure of LDL cholesterol, adjusted to no statin treatment status (if the patient is on statin treatment), centered at 3.6 (mmol/L)

hba1c: HbA1c (mmol/mol)

## Time-varying variables at the baseline
cycle: all equal 0, i.e. baseline

CurrAge_cent: age at baseline, centered at 60 divided by 10

mi_0_1, mi_1_2, mi_2_3 and mi_3_inf: MI in the same year, MI in previous year, MI in 2 years ago, and MI in at least 3 years ago. They all equal 0 at baseline, as they are used to record incident MI, rather than history. The same principle applies to stroke and CRV. 

dmPre_1, dmPre_3, dmPre_4, dm_0_10 and dm_10_inf: no diabetes & HbA1c<32, no diabetes & HbA1c 37-42, no diabetes & HbA1c 42-48, diabetes duration 0-10 years, and diabetes duration 10+ years, respectively. Please prespecify these variables using the baseline diabetes diagnosis and HbA1c. All equal 0 means no diabetes & HbA1c 32-37. They will be automatically updated in the model simulation.

cancer_bsl_*: baseline cancer in corresponding years of history.

cancer_icd_*: incident cancer that happens in a certain year ago in simulation; all equaling 0 at baseline; They will be automatically updated in the model simulation.

CurrAge_cent_int_*: interaction terms with age; please prespecify depending on the baseline status. They will be automatically updated in the model simulation.

mi_1_inf_qol_pp, stroke_1_inf_qol_pp, two_plus_qol, wstatinayr1, wstatinayr2And, ind6_1 and ind7_1, d_without_*, nvd and vd: all equal 0










