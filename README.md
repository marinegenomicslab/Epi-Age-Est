Weber et al. CJFAS (in prep)

Data and scripts/code associated with "Novel epigenetic age estimation in wild-caught northern Gulf of Mexico reef fishes" D. Nick Weber, Andrew T. Fields, William F. Patterson III, Beverly K. Barnett, Christopher M. Hollenbeck, and David S. Portnoy. In prep.

Associated data can be found in the "Data" folder, incuding:
- "rsn_fulldata.csv" contains all raw CpG site loci and the associated methylated/unmethylated read counts for each site in each red snapper individual.
- "rdg_fulldata.csv" contains all raw CpG site loci and the associated methylated/unmethylated read counts for each site in each red grouper individual.
- "rsn_percent_wide.csv" contains the matrix of percent methylation values (post all filtering steps) necessary to recreate the elastic net regression results for red snapper.
- "rdg_percent_wide.csv" contains the matrix of percent methylation values (post all filtering steps) necessary to recreate the elastic net regression results for red grouper.

Associated code can be found in the "Scripts" folder, including:

- "BAYES_GLM.r" contains code associated with the Bayesian GLM run with 4,000 warm-up/sampling iterations.
- "BAYES_GLM_50k.r" contains code associated with the Bayesian GLM run with 50,000 warm-up/sampling iterations.
- "95% Confidence Intervals.Rmd" contains the code necessary to calculate 95% confidence intervals.
- "Elastic Net Regression" contains the code necessary to run the elastic net regression.
- "assem_config.file" contains the code and parameters used to create a de novo reference genome for red grouper.
