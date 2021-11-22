library(tidyverse)
library(magrittr)
library(brms)
options(buildtools.check = function(action) TRUE )

cats<-read_csv("csv.csv")

# fit model Wr by log-transformed discharge
wr_gamma_brm <- brm(Wr ~ log10(Discharge), 
                    family = Gamma(link = "log"),
                    data = cats,
                    prior = c(prior(normal(4,.25), class = "Intercept"),
                              prior(normal(0, 1), class = "b"),
                              prior(exponential(.05), class = "shape")),
                    cores = 4, chains = 4, iter = 1000)

print(summary(wr_gamma_brm), digits=4)
plot(wr_gamma_brm)

# Plot
conditional_effects(wr_gamma_brm)
plot(conditional_effects(wr_gamma_brm), points=T)

# Prediction interval
vals<-unique(cats$Discharge);vals
vals<-vals[order(vals)]; vals

discharge_seq <- tibble(Discharge = seq(from = min(vals), to = max(vals), by = 1))

mu_summary <-fitted(wr_gamma_brm, newdata = discharge_seq) %>%
  as_tibble() %>%
  bind_cols(discharge_seq) %>% 
  as.data.frame()
head(mu_summary)

pred_wr <-predict(wr_gamma_brm, newdata = discharge_seq,probs = c(.05,0.95)) %>%
  as_tibble() %>%
  bind_cols(discharge_seq) %>% 
  as.data.frame()
head(pred_wr)

## Final GGPLOT
cats %>%
  ggplot(aes(x = Discharge)) +
  coord_cartesian(ylim = c(NA,140),xlim=c(0,7000))+
  geom_ribbon(data = pred_wr, 
              aes(ymin = Q5, ymax = Q95),
              fill = "grey83") +
  geom_smooth(data = mu_summary,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              fill = "grey70", color = "black", alpha = 1, size = 1/2) +
  geom_point(aes(y = Wr),
             color = "black", shape = 1, size = 1.5, alpha = 2/3) +
  scale_x_continuous(name="Discharge (cfs)",breaks=seq(0,7000,1000))+
  scale_y_continuous(name="Relative Weight",breaks=seq(60,140,10))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill=NA))

##########################################################################
############   Add length by year effect   ################################
##########################################################################
# fit model Wr by log-transformed discharge
get_prior(Wr ~ log10(Discharge)+ ??, data=cats, family=Gamma(link="log"))

wr_tl_yr <- brm(Wr ~ log10(Discharge) + ??, 
                    family = Gamma(link = "log"),
                    data = cats,
                    prior = c(prior(normal(4,.25), class = "Intercept"),
                              prior(normal(0, 1), class = "b"),
                              prior(exponential(.05), class = "shape"),
                              prior(exponential(2), class = "sd")),
                    #sample_prior = "only",
                    cores = 4, chains = 1, iter = 1000)
print(summary(wr_tl_yr),digits = 4)
