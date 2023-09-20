# UK_CVD_model
The folder contains the R codes for the UK CVD micro-simulation model and input files
"master_hpc_boots.R" is the master file to control input/output of the model
To run this model, you need to specify path, input and other parameters in the master file and run it
"fun_sim_boots.R" are the main model functions that the master file calls

In the "input" folder, "cf_hpc_v4_VDppUpdate.rds", "tx_hpc.rds" and "pf.rds" are equation coefficients and model parameters that our research project generated. They were derived from individual patient data from trials of Cholesterol Treatment Trialists' Collaboration (CTT), and calibrated into a UK population of 500,000 individuals, i.e. UK Biobank.  

Files in the input folder with the prefix "sample_" are sample datasets of five pseudo individuals' characteristics, for primary prevention (without cardiovascular disease history) and secondary prevention (with cardiovascular disease history) separately. The model reads in baseline (b) and time-varying (t) characteristics in two separate files, in order to speed up simulation. Please specify your input characteristics in the same format as the samples

The following are columns' names and their corresponding meaning
id_new: patient ID
Intercept: all equal to 1
male: 1 = male; 0 = female
race_afro; race_sa; race_otherna: equal to 1 means African, South Asain and other ethnicities, respectively; otherwise White ethnicity
unhealthy_diet: 1 = unhealthy diet; 0 = healthy diet. Please refer to our published paper for more details
h_mi_raw_only; cvd_only; h_pad_raw_only, othchd_only; number_events2: equal 1 means history of myocardial infarction (MI), history of cerebral vascular disease (stroke), history of peripheral artery disease (PAD), history of other coronary heart diseases, history of more than one of these diseases; In primary prevention population, all equal 0; in secondary prevention population, one of them equals to 1  
