rm(list=ls())
library(fanetc)

library(fanetc)
library(rstan)
library(truncnorm)
library(bayesplot)
library(tidybayes)
library(distortr)
load('./data_setup2.RData')

reg_d <- lapply(reg_d, function(dd){
  dd %>% 
    filter(age<=75) %>% 
    mutate(age_gp = cut(age, breaks = seq(0,75,5))) %>% 
    group_by(age_gp) %>% 
    mutate(age = median(age)) %>% 
    group_by(newpuma, year, age_gp,age) %>% 
    summarise(migrant = sum(migrant), 
              totalpop = sum(totalpop)) %>% 
    mutate(rate = migrant/totalpop) %>% 
    ungroup()
})

reg_d$`6`$statefip = 6
reg_d$`36`$statefip = 36
reg_d$`48`$statefip = 48

reg_d$`6`$state = 'California'
reg_d$`36`$state = 'New York'
reg_d$`48`$state = 'Texas'

aa <- paste0(reg_d$`36` %>% group_by(age) %>% slice(1)%$% age_gp, collapse=", ")

# reg_d$`36`
# reg_d$`36` <- reg_d$`36` %>% 
#   left_join(obs_counties, by = c("newpuma","statefip"))

reg_d$`36` %>% mutate(dd = totalpop-migrant) %>% arrange(-rate) 

set.seed(36)
# sim_pop <- 500
sim_pop <- 1000
reg_d$`36` <- reg_d$`36` %>% 
  rowwise() %>% 
  mutate(sim_dat = rpois(1, rate*sim_pop)) %>% 
  ungroup() %>% 
  mutate(sim_pop = sim_pop,
         sim_rate = sim_dat/sim_pop)

### Model data 
tmp <- reg_d$`36` %>% 
  filter(year <=2019)

tmp <- tmp %>% 
  mutate(age_gp_num = as.numeric(as.factor(age_gp)), 
         region = as.numeric(as.factor(newpuma)))

# tmp <- tmp %>%
# filter(region <= 6)

A = length(unique(tmp$age_gp))
R = length(unique(tmp$newpuma))
N = length(tmp$migrant)
T = length(unique(tmp$year))
age_unique = tmp$age_gp %>% unique()
region_unique = tmp$newpuma %>% unique()
pred_df = expand.grid(age_gp = age_unique, 
                      newpuma = region_unique) %>% 
  mutate(year = 2020)

### Models

summary_mdl <- function(model_sum){
  d <- model_sum$summary[,c("mean", "n_eff", "Rhat","sd","se_mean","2.5%","97.5%")]
  #extract age-specific rates (tau)
  est_tau <- model_sum$summary[,c("mean", "sd", "n_eff", "Rhat")][grep("tau\\[[0-9]*\\]",rownames(d)),1]
  
  #extract region-specific rates (gamma)
  est_gamma <- model_sum$summary[,c("mean", "sd", "n_eff", "Rhat")][grep("gamma\\[[0-9]*\\]",rownames(d)),1]
  
  #extract intercept  (alpha)
  est_alpha <- model_sum$summary[,c("mean", "sd", "n_eff", "Rhat")][grep("alpha",rownames(d)),1]
  
  #extract log rates
  r1 <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("log_rate",rownames(d)),1]
  r1_upper <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("log_rate",rownames(d)),3]
  r1_lower <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("log_rate",rownames(d)),2]
  
  #extract predictions intervals 
  # p1 <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("rate_p",rownames(d)),1]
  # p1_upper <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("rate_p",rownames(d)),3]
  # p1_lower <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("rate_p",rownames(d)),2]
  # p1_se <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("rate_p",rownames(d)),2]
  
  # extract poisson_rng outcomes
  yp_upper <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("y_p",rownames(d)),3]
  yp_lower <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("y_p",rownames(d)),2]
  
  #extract fixed effects 
  f1 <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("fixed_eff",rownames(d)),1][1:nrow(tmp)]
  
  #extract RW on age + fixed effects 
  f2 <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("deviation_fixed_eff",rownames(d)),1][1:nrow(tmp)]
  
  #extract RW on age + RW on time + fixed effects 
  f3 <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("both_deviation_fixed_eff",rownames(d)),1][1:nrow(tmp)]
  
  
  #extract predictions intervals year 2020
  p2020 <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("pred_rate",rownames(d)),1]
  p2020_lower <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("pred_rate",rownames(d)),2]
  p2020_upper <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("pred_rate",rownames(d)),3]
  
  #generated
  p1.1 <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("generated_rate",rownames(d)),1]
  p1.1_upper <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("generated_rate",rownames(d)),3]
  p1.1_lower <- model_sum$summary[,c("mean", "2.5%","97.5%")][grep("generated_rate",rownames(d)),2]
  
  
  
  est_age_rates <- c(est_alpha + mean(est_gamma), est_tau + est_alpha + mean(est_gamma))
  est_age_rates <- exp(est_age_rates)
  
  mig_sched = pred_df %>% dplyr::select(age_gp, newpuma)
  mig_sched$fit = model.matrix(~age_gp + as.factor(newpuma), data = mig_sched) %*% (matrix(c(est_alpha, est_tau,est_gamma)))
  
  tmp_results <- tmp
  
  tmp_results$est_log_rates = r1
  tmp_results$est_log_rates_upper = r1_upper
  tmp_results$est_log_rates_lower = r1_lower
  tmp_results$fixed_eff = f1
  tmp_results$deviation_fixed_eff = f2
  tmp_results$both_deviation_fixed_eff = f3
  tmp_results$gen_log_rates = p1.1
  tmp_results$gen_log_rates_lower = p1.1_lower
  tmp_results$gen_log_rates_upper = p1.1_upper
  
  
  pred_results <- pred_df %>% 
    left_join(reg_d$`36` %>% dplyr::select(newpuma, age_gp, year, rate))
  
  pred_results$pred_log_rates = p2020
  pred_results$pred_log_rates_lower = p2020_lower
  pred_results$pred_log_rates_upper = p2020_upper
  
  
  tmp_results <- tmp_results %>% 
    mutate(est_rates = exp(est_log_rates),
           fixed_eff_rates = exp(fixed_eff),
           deviation_fixed_eff = exp(deviation_fixed_eff),
           both_deviation_fixed_eff = exp(both_deviation_fixed_eff),
           fitted_rate_lower = exp(est_log_rates_lower),
           fitted_rate_upper = exp(est_log_rates_upper),
           pred_rates = gen_log_rates,
           pred_rate_upper = gen_log_rates_upper,
           pred_rate_lower = gen_log_rates_lower) %>%
    mutate(log_rate_obs = log(rate))
  
  pred_results <- pred_results %>% 
    mutate(pred_rates = exp(pred_log_rates),
           pred_rate_upper = exp(pred_log_rates_upper),
           pred_rate_lower = exp(pred_log_rates_lower)
    ) %>%
    mutate(log_rate_obs = log(rate))
  
  tmp_results <- full_join(tmp_results, pred_results)
  return(list(fit = tmp_results, mig_sched = mig_sched, est_age_rates=est_age_rates))
}

load("../Model results/finalmodel_rw_RI_moreits.RData")
model_sumrw <- model_sum
summary_rw <- summary_mdl(model_sumrw)
load("../Model results/finalmodel_spline_moreits.RData")
model_sumcs <- model_sum
summary_cs <- summary_mdl(model_sumcs)

#Check 


### Standard age schedule

mig_sched <- full_join(summary_rw$mig_sched %>% mutate(method = "RW"),
                       summary_cs$mig_sched %>% mutate(method = "Spline")) %>% 
  mutate(age = factor(age_gp, labels = unique(tmp$age)),
         age = as.numeric(as.vector(age)))

mdl_fit <- full_join(summary_rw$fit %>% mutate(method = "RW"),
                     summary_cs$fit %>% mutate(method = "Spline")) 

ggplot(mig_sched, aes(age, exp(fit), group = method, col = method))+ 
  geom_line()+
  facet_wrap(~newpuma)+ 
  theme_bw()+ 
  scale_y_continuous(name ="rate",
                     breaks = scales::pretty_breaks(n=5))+
  scale_x_continuous(name ="age",
                     breaks = scales::pretty_breaks(n=5))
  

in_mig_plot <- summary_rw$mig_sched %>% 
  group_by(age_gp) %>% 
  ggplot(aes(age_gp, exp(fit), col = factor(newpuma), group = newpuma)) + 
  geom_line(alpha=0.2, show.legend = FALSE)+ 
  theme_bw() +
  ggtitle("In-Migration Standard Age Schedule, New York") +
  scale_y_continuous(name ="rate",
                     breaks = scales::pretty_breaks(n=10))+
  xlab("age")+ 
  geom_line(data=data.frame(age = tmp$age_gp %>% unique(),
                            fit = summary_rw$est_age_rates,
                            newpuma = NA), aes(age, fit) ,inherit.aes = TRUE,col = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


ggsave(filename = "../Paper/Plots/in_mig_plot_RW.pdf",
       plot=in_mig_plot, 
       device = "pdf",
       dpi = 300,
       height = 4, 
       width = 6)



load("boundary_info.RData")

plot_ts_pred <- function(df,puma){
  df %>% 
    filter(newpuma == puma) %>% 
    ggplot(aes(year, rate))+
    facet_wrap(~age_gp, ncol = 5)+
    geom_line(aes(year, pred_rates, col = "pred", group= age_gp)) +
    geom_ribbon(aes(ymin = pred_rate_lower, ymax = pred_rate_upper, col = "pred", fill ="pred", group = age_gp), alpha = 0.2, show.legend = FALSE)+ 
    geom_point(aes(col = "data", group = age_gp),size=1)+
    scale_color_manual(name = "", values = c("pred" = "red", "data" = "black")) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}

#PUMA 2100
obs_counties %>% filter(statefip==36 & newpuma == 2100) # Doesn't correspond to an exact county. 
plot_fun(df= to_keep, newpuma = 2100)
# Corresponds to Greene county and columbia county. 
gg = plot_ts_pred(df = summary_rw$fit, 2100)+ 
  ggtitle("In-Migration Rates, Greene and Columbia County, New York")


ggsave(filename = "../Paper/Plots/puma2100_tsfit.pdf",
       plot=gg, 
       device = "pdf",
       dpi = 300,
       height = 5, 
       width = 6.5)

#PUMA 1300
obs_counties %>% filter(statefip==36 & newpuma == 1300) # Doesn't correspond to an exact county. 
plot_fun(df= to_keep, newpuma = 1300)
# Corresponds to Wyoming and Livingston county
gg = plot_ts_pred(df = summary_rw$fit, 1300)+ 
  ggtitle("In-Migration Rates, Wyoming and Livingston County, New York")


ggsave(filename = "../Paper/Plots/puma1300_tsfit.pdf",
       plot=gg, 
       device = "pdf",
       dpi = 300,
       height = 5, 
       width = 6.5)


#PUMA 1700
obs_counties %>% filter(statefip==36 & newpuma == 1700) # Corresponds to exactly 1 county
plot_fun(df= to_keep, newpuma = 1700)
# Corresponds to Schenectady
# check:
county_dic %>% filter(countyfip==93 & statefip==36)
gg = plot_ts_pred(df = summary_rw$fit, 1700)+ 
  ggtitle("In-Migration Rates, Schenectady County, New York")


ggsave(filename = "../Paper/Plots/puma1700_tsfit.pdf",
       plot=gg, 
       device = "pdf",
       dpi = 300,
       height = 5, 
       width = 6.5)


# PUMA 2300
county_dic %>% filter(county == 'Tompkins') #Home of Cornell university
obs_counties %>% filter(statefip==36 & countyfip==109)

obs_counties %>% filter(statefip==36 & newpuma == 2300) # Corresponds to exactly 1 county
# plot_fun(df= to_keep, newpuma = 1700)
# Corresponds to Tomkins County 
# check:
# county_dic %>% filter(countyfip==93 & statefip==36)
gg = plot_ts_pred(df = summary_rw$fit, 2300)+ 
  ggtitle("In-Migration Rates, Tompkins County, New York")


ggsave(filename = "../Paper/Plots/puma2300_tsfit.pdf",
       plot=gg, 
       device = "pdf",
       dpi = 300,
       height = 5, 
       width = 6.5)


# PUMA 2300
county_dic %>% filter(county == 'New York') #Manhattan
obs_counties %>% filter(statefip==36 & countyfip==61)

obs_counties %>% filter(statefip==36 & newpuma == 3800) # Corresponds to exactly 1 county
# plot_fun(df= to_keep, newpuma = 1700)
# Corresponds to Tomkins County 
# check:
# county_dic %>% filter(countyfip==93 & statefip==36)
gg = plot_ts_pred(df = summary_rw$fit, 3800)+ 
  ggtitle("In-Migration Rates, New York County, New York")


ggsave(filename = "../Paper/Plots/puma3800_tsfit.pdf",
       plot=gg, 
       device = "pdf",
       dpi = 300,
       height = 5, 
       width = 6.5)

year2020_results %>% 
  filter(capt == FALSE) %>% 
  filter(rate < fitted_rate_lower) %>% 
  arrange(desc(rate))


# PUMA 300
obs_counties %>% filter(statefip==36 & newpuma == 300) # Corresponds to more than one county
plot_fun(df= to_keep, newpuma = 300)
# Corresponds to Warren and Washington County 
# check:
# county_dic %>% filter(countyfip==93 & statefip==36)
gg = plot_ts_pred(df = summary_rw$fit, 300)+ 
  ggtitle("In-Migration Rates, Warren and Washington County, New York")


ggsave(filename = "../Paper/Plots/puma300_tsfit.pdf",
       plot=gg, 
       device = "pdf",
       dpi = 300,
       height = 5, 
       width = 6.5)


# 2020 analysis summary stats 
year2020_results <- summary_rw$fit %>% 
  filter(year==2020) %>% 
  dplyr::select(newpuma, year, age_gp, rate,fitted_rate_lower,fitted_rate_upper) %>% 
  mutate(capt = rate >= fitted_rate_lower & rate <= fitted_rate_upper)

summary_rw$fit %>% 
  filter(year < 2020) %>% 
  mutate(capt = rate >= fitted_rate_lower & rate <= fitted_rate_upper) %>% 
  filter(capt == FALSE)%$% rate %>% summary()



### Migration curves
plot_mig_curve <- function(df, puma){
df %>% 
  filter(year<2020) %>% 
  filter(newpuma == puma) %>% 
  ggplot(aes(age, rate))+
  geom_point(aes(col = "data"),size=0.7)+ 
  facet_wrap(~year,nrow=1)+ 
  geom_line(aes(age, deviation_fixed_eff, col = "expected + delta(a,r)"))+
  geom_line(aes(age, both_deviation_fixed_eff, col = "expected + delta(a,r) + delta(t,r)"))+
  # geom_line(aes(age, est_rates, col = "expected + delta(a,r) + delta(t,r) + delta(a,r,t)"))+
  # geom_ribbon(aes(ymin = gen_log_rates_lower, ymax = gen_log_rates_upper),
  #             alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(age, fixed_eff_rates, col = "expected"))+
  scale_color_manual(name = "", values = c("data" = "black",
                                           "expected" = "darkorange",
                                           "expected + delta(a,r)" = "red",
                                           "expected + delta(a,r) + delta(t,r)" = "blue"),
                     label = c("data" = "data",
                               "expected" = "expected",
                               "expected + delta(a,r)" = expression(expected~+~phi(a,r)),
                               "expected + delta(a,r) + delta(t,r)" = expression(expected~+~phi(a,r)~+~xi(t,r)))) +
  scale_x_continuous(name = "age",
                     breaks = c(3, 18, 33, 48, 63),
                     labels = c('(0,5]', '(15,20]','(30,35]', '(45,50]', '(60,65]'))+
  theme_bw()+
    theme(legend.position="bottom")+
    theme(axis.text.x = element_text(angle = 90, hjust=1))}

gg = plot_mig_curve(df = summary_rw$fit, 300)+ 
  ggtitle("In-Migration Rates, Warren and Washington County, New York")

gg2 = plot_mig_curve(df = summary_rw$fit, 2100)+ 
  ggtitle("In-Migration Rates, Greene and Columbia County, New York")

gg3 = plot_mig_curve(df = summary_rw$fit, 1300)+ 
  ggtitle("In-Migration Rates, Wyoming and Livingston County, New York")

library(ggpubr)
gg_tog <- ggarrange(gg2,gg,nrow=2, common.legend = TRUE, legend = "bottom")

ggsave(filename = "../Paper/Plots/components.pdf",
       plot=gg_tog, 
       device = "pdf",
       dpi = 300,
       height = 5, 
       width = 6.5)

summary_rw$fit %>% 
  filter(year<2020) %>% 
  dplyr::select(newpuma, year, age_gp, rate,pred_rate_lower,pred_rate_upper) %>%
  mutate(capt = rate >= pred_rate_lower & rate <= pred_rate_upper)%$% capt %>% mean()


summary_rw$fit %>% 
  filter(year<2020) %>% 
  dplyr::select(newpuma, year, age_gp, rate,pred_rate_lower,pred_rate_upper) %>%
  mutate(capt = rate >= pred_rate_lower & rate <= pred_rate_upper) %>% 
  filter(capt)

for_pred_plot <- summary_rw$fit %>% 
  filter(year==2020) %>% 
  dplyr::select(newpuma, year, age_gp, rate,pred_rate_lower,pred_rate_upper, pred_rates) %>%
  mutate(capt = rate >= pred_rate_lower & rate <= pred_rate_upper) %>% 
  mutate(lower_than_expected = capt == FALSE & rate < pred_rate_lower,
         higher_than_expected = capt == FALSE & rate > pred_rate_upper) %>% 
  mutate(exp = ifelse(capt == TRUE, "As expected", NA), 
         exp = ifelse(capt == FALSE & lower_than_expected == TRUE, "Lower than expected",exp),
         exp = ifelse(capt == FALSE & lower_than_expected == FALSE, "Higher than expected",exp))

gg1 = to_keep %>% 
  mutate(migpuma = as.numeric(as.vector(migpuma))) %>% 
  left_join(for_pred_plot, by = c("migpuma"="newpuma")) %>% 
  filter(age_gp %in% age_unique[1:9]) %>% 
  ggplot()+
  geom_polygon(aes( x = long, y = lat, group = group,fill= exp), color="black", show.legend = TRUE,linewidth=0.2) +
  theme_bw()+ 
  facet_wrap(~age_gp)+
  theme(legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank())+ 
  ggtitle("Migration rates for year 2020")

ggsave(filename = "../Paper/Plots/expected_g1.pdf",
       plot=gg1, 
       device = "pdf",
       dpi = 300,
       height = 8, 
       width = 8)

gg1 = to_keep %>% 
  mutate(migpuma = as.numeric(as.vector(migpuma))) %>% 
  left_join(for_pred_plot, by = c("migpuma"="newpuma")) %>% 
  filter(age_gp %in% age_unique[9:15]) %>% 
  ggplot()+
  geom_polygon(aes( x = long, y = lat, group = group,fill= exp), color="black", show.legend = TRUE,linewidth=0.2) +
  theme_bw()+ 
  facet_wrap(~age_gp)+
  theme(legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank())+ 
  ggtitle("Migration rates for year 2020")

ggsave(filename = "../Paper/Plots/expected_g2.pdf",
       plot=gg1, 
       device = "pdf",
       dpi = 300,
       height = 8, 
       width = 8)
