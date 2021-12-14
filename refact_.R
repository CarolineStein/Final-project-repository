#Refactoring a piece of existing code
#Waning VE 

outfile <- "/ihme/homes/cstein87/vaceff/tmv/"

library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(formattable)

tsvacc <- read.csv("/ihme/homes/cstein87/vaceff/tsvacc.csv")

#keep only "0=inlcude" of "exclude this source from analysis" variable
ve1 <- subset(tsvacc, exclude == 0)

ve2 <- subset(ve1, select = c(study_id, author, X1st.dose.only., location_id, location_name, location_name2, 
                              vaccine_developer, variant, symptom_severity, severity,
                              sample_size, efficacy_mean, efficacy_lower, efficacy_upper, start_interval, end_interval))

#rename "X1st.dose.only." variable
ve2 <- rename(ve2, firstdose = X1st.dose.only.)

#rename "location_name2" variable
ve2 <- rename(ve2, variable_name = location_name2)

#weeks after second dose
ve2 <- mutate(ve2, mid_point = ((end_interval - start_interval)/2) + start_interval)

#CI to standard error
ve2$efficacy_lower <- as.numeric(as.character(ve2$efficacy_lower))
ve2 <- mutate(ve2, se = (efficacy_upper - efficacy_lower)/3.92)

#transform to logit
library(crosswalk, lib.loc = "/ihme/code/mscm/R/packages/")
logit <- delta_transform(mean = ve2$efficacy_mean, sd = ve2$se, transformation = "linear_to_logit")

names(logit) <- c("mean_logit", "sd_logit")
vacc_logit <- cbind(ve2, logit)
ve2 <- vacc_logit

#Meta-regression spline models with separate models by vaccine and outcome
library(mrbrt001, lib.loc = "/ihme/code/mscm/R/packages/")

vaccines <- unique(ve2$vaccine_developer)

dt_combined <- data.table()

#v <- as.character(vaccines[1])
#s <- as.character(symp[1])

for(v in vaccines) {
  vef_a1 <- subset(ve2, vaccine_developer == v)
  symp <- unique(vef_a1$severity)
  
  for(s in symp){
    vef <- subset(ve2, vaccine_developer == v & severity == s)
    
    #3.5.1 - Setting priors and shape constraints on splines
    ve6 <- MRData()
    ve6$load_df(
      data = vef,  
      col_obs = "mean_logit", 
      col_obs_se = "sd_logit",
      col_covs = list("mid_point"),
      col_study_id = "variable_name" )
    
    mod1 <- MRBRT(
      data = ve6,
      cov_models = list(
        LinearCovModel("intercept", use_re = TRUE),
        LinearCovModel(
          alt_cov = "mid_point",
          use_spline = TRUE,
          spline_knots = array(seq(0, 1, length.out = 4)),       # 1.0  9.5 18.0 26.5    #(seq(0, 1, by = 0.2)
          spline_degree = 2L,
          spline_knots_type = 'domain',
          spline_r_linear = TRUE,
          spline_l_linear = FALSE,
          prior_spline_monotonicity = 'decreasing'
          # prior_spline_convexity = "convex"
          # prior_spline_maxder_gaussian = array(c(0, 0.01))
          # prior_spline_maxder_gaussian = rbind(c(0,0,0,0,-1), c(Inf,Inf,Inf,Inf,0.0001))
        )
      )
    )
    
    mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
    
    df_pred3 <- data.frame(mid_point = seq(0, 52, by = 0.1))
    dat_pred3 <- MRData()
    dat_pred3$load_df(
      data = df_pred3, 
      col_covs=list('mid_point')
    )
    
    df_pred3$pred5 <- mod1$predict(data = dat_pred3)
    
    #spreadsheet 
    dtfinal <- copy(df_pred3)
    dtfinal <- data.table(dtfinal)
    
    #calculate the VE from logit scale
    dtfinal[, ve := plogis(pred5)]
    dtfinal[, vaccine := v]
    dtfinal[, severity := s]
    
    dt_combined <- rbind(dt_combined, dtfinal, fill=TRUE)
    
    
    #Plots by outcomes and vaccines
    #solid line for the mid_point period where there are data and dashed line for projection line
    dt_combined <- mutate(dt_combined, solid = ifelse(mid_point < 32.5, "dashed", "solid"))
    
    #to know the weight for each data point: inverse of se
    ve2 <- mutate(ve2, insesqua = (1/sqrt(se))/10)
    
    f1 <- "Times"
    
    dt_combined1 <- subset(dt_combined, vaccine == v & severity == s)
    dt_combined2 <- subset(ve2, vaccine_developer == v & severity == s)
    
    p <- ggplot() + 
      geom_line(data = dt_combined1, mapping = aes(x = mid_point, y = ve, linetype = solid), size = 1.5) +
      geom_point(data = dt_combined2, mapping = aes(x = mid_point, y = efficacy_mean, color = variable_name, size = insesqua)) +
      ylim(c(0, 1))+
      xlim(c(0, 52)) +
      theme_classic() +
      theme(text=element_text(size=14, family = f1),
            legend.position = "bottom") +
      guides(linetype = FALSE) +
      labs(y="Vaccine effectiveness", x = "Week after second dose", color = "Study", title = paste0("", v, " ", "at preventing", " ", s)) +
      geom_hline(aes(yintercept = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), linetype = "yellow"))+
      guides(size = FALSE)
    
    for (k in mod1$cov_models[[2]]$spline_knots) p <- p + geom_vline(xintercept = k)
    
    png(filename = paste0(outfile, "w_i", v, "_", s, "_.png"),
        height = 7, width = 12, units = "in", res = 300)
    print(p)
    dev.off()
    
  }
}
