#!/bin/R 

#Changing when scientific numbers are written
options("scipen"=100, "digits"=4)
options(warn = 1)

#Loading libraries
library(tidyverse)
library(emmeans)
library(rstanarm)
library(coda)

#Defining Function
variable_importance <- function(fit, interval, variable_names){
  #A function to look at the posterior intervals and summarize if there is any variation as a result of the given variables
  posterior_interval(fit, interval) %>%
    as_tibble(rownames = 'parameter') %>%
    rename(lwr = 2, upr = 3) %>%
    filter(str_detect(parameter, str_c(variable_names, collapse = '|'))) %>%
    mutate(variable_testing = str_extract_all(parameter, str_c(variable_names, collapse = '|')),
           variable_testing = map_chr(variable_testing, str_c, collapse = ':')) %>%
    group_by(variable_testing) %>%
    summarise(important_var = if_else(all(lwr < 0 & upr > 0), 'no', 'yes'))
}

#Importing Data
tmp_out <- read.table("BAYES_split_data.txt", head=F)
tmp_col <- read.table("BAYES_data_names.txt", head=F)
tmp_col[,1] <- as.character(tmp_col[,1])
tmp_col[1,] <- "Loc"
colnames(tmp_out) <- t(tmp_col)
tmp_names <- as.factor(as.matrix(read.table("Sample_names.txt", head=F)))
tmp_strata <- read.table("BAYES_strata.txt", head=T)

#Preparing objects for data collection
GLM <- list()
GLM$Stats <- data.frame(matrix(nrow=length(tmp_out$Loc), ncol=2), row.names=tmp_out$Loc)
colnames(GLM$Stats) <- c("Max_Rhat", "Min_ESS")

#Loop over all loci
for(i in 1:nrow(tmp_out)){
print(as.character(tmp_out$Loc[i]))
tmp <- tmp_out[tmp_out$Loc[i],2:ncol(tmp_out)]
tmp2 <- matrix(tmp, byrow=T, ncol=2)
colnames(tmp2) <- c("Meth", "Total")
tmp3 <- cbind(tmp_strata,tmp2[match(tmp_strata$Sample,tmp_names),])
for(j in c("Meth", "Total")){tmp3[,j] <- unlist(tmp3[,j])}

#Bayesian model
bayes_model <- stan_glmer(cbind(Meth, Total - Meth) ~ Age + (1 | Sample), data = tmp3, family = binomial(link = "logit"), cores = 4, iter=50000)

#Output warnings
#write(c(as.character(tmp_out$Loc[i]), names(warnings())[1]),"Warnings.txt",append=T)

#Purge warning messages
#assign("last.warning", NULL, envir = baseenv())

#Gathering convergance stats
GLM$Stats[i, ] <- c(max(summary(bayes_model)[, "Rhat"]), min(summary(bayes_model)[, "n_eff"]))

#Significant posteriors
tmp_post <- posterior_interval(bayes_model, prob = 0.95)
VARS <- rownames(tmp_post)[abs(tmp_post[,2])>abs(tmp_post[,2]-tmp_post[,1])]

#Find significant factors
sig <- joint_tests(bayes_model)
sig_fact <- cbind(rep(as.character(tmp_out$Loc[i]),nrow(sig)), sig)

#Getting posterior ranges
test1 <- data.frame(emmeans(bayes_model, ~ Age, type = "response"))

test1 <- test1[!(is.na(test1$prob)),]

test1 <- cbind(rep(as.character(tmp_out$Loc[i]),nrow(test1)), test1)

#Jason test
jtest1 <- variable_importance(bayes_model, 0.95, c('Age'))
jtest2 <- as.data.frame(jtest1[jtest1$important_var=="yes",1])
if(nrow(jtest2)>0){jtest2 <- cbind(rep(as.character(tmp_out$Loc[i]),nrow(jtest2)), jtest2)}

#Writing outputs
write.table(GLM$Stats[i, ], "GLM_Stats.txt", col.names=F, row.names=T, quote=F, append=T)
write.table(t(as.vector(c(as.character(tmp_out$Loc[i]), VARS))), "GLM_Sig_post.txt", col.names=F, row.names=F, quote=F, append=T)
write.table(sig_fact, "GLM_Sig_jt.txt", col.names=F, row.names=F, quote=F, append=T)
write.table(jtest2, "GLM_Sig_js.txt", col.names=F, row.names=F, quote=F, append=T)

#Writing factor intervals
write.table(test1, "GLM_emmeans_1factor.txt", col.names=F, row.names=F, quote=F, append=T)
}
