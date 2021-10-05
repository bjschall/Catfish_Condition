library(FSA)
library(tidyverse)
library(magrittr)
library(brms)
cat<-read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Bayes Course/Applied Bayes/Scotland All Years.csv")
water<-read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Bayes Course/Applied Bayes/WaterLevels.csv")


water<-select(water, -Date)
cat<-filterD(cat, Spp=="CCF")
cats<-merge(cat, water, by="Year")        
(ccf.cuts<-psdVal("Channel Catfish"))
(wsCCF<-wsVal("Channel Catfish"))

cats%<>% mutate(LenCat=lencat(Length,w=10), 
              logtl=log10(Length),
              logwt=log10(Weight),
              gcat=lencat(Length, breaks=ccf.cuts, use.names = T,drop.levels = T),
              Ws=10^(wsCCF[["int"]]+wsCCF[["slope"]]*logtl), 
              Wr=Weight/Ws*100)
  
#cats$Discharge<-factor(cats$Discharge)

cats<-filterD(cats, !is.na(Wr) & Wr>=50 & Wr<=150)
Summarize(Wr~Discharge, cats)
boxplot(Wr~Discharge, cats)

############################################################################
#######  Basic linear model evaluating Discharge  ##########################
############################################################################
# y ~ dgamma(scale, shape)
# scale = mu/shape
# log(mu)  = a + bx
# a ~ normal(4.5,1)
# b ~ normal(-6,1)
# shape ~ ?

plot(nums-mean(nums))

N = 100  # number of simulations

# simulate priors
priors <- tibble(a = rnorm(N, 2.5, 1),
                 b = rnorm(N,0, 1),
                 shape = rexp(N, 2), # for brms
                 scale = rexp(N, 2), # for rethinking
                 sim = 1:N)

# data (only the x values, since we're simulating y and mu and pretending we don't have them yet)
x <- rnorm(100, 2000,1000)

# combine and simulate
prior_and_x <- priors %>% expand_grid(x = x) %>%    # combine priors and x's
  mutate(mu = exp(a + b*(x-mean(x))/sd(x)),                         # simulate regressions
         y = rgamma(nrow(.), scale = mu/shape, shape),   # simulate data (e.g., y_rep) - brms parameterization
         y_scale = rgamma(nrow(.), shape = mu/scale, scale))  # simulate data (e.g., y_rep) - rethinking parameterization
# plot
prior_and_x %>% 
  ggplot(aes(x = (x-mean(x))/sd(x), y = mu, group = sim)) + 
  geom_line() +
  theme_classic()+
  geom_point(aes(y = y)) +
  labs(y = "sim") +
  #scale_y_log10() +
  NULL

# fit with brms (NOTE: rethinking and brms use different parameterizations for gamma)
gamma_brm <- brm(Wr ~ Discharge_cent, 
                 family = Gamma(link = "log"),
                 data = cats %>% mutate(Discharge_cent = (Discharge - mean(Discharge))/sd(Discharge)),
                 prior = c(prior(normal(4.5,1), class = "Intercept"),
                           prior(normal(1, 1), class = "b"),
                           prior(exponential(2), class = "shape")),
                 cores = 4, chains = 4, iter = 2000, warmup=1000)

summary(gamma_brm)
plot(gamma_brm)

#sample from the posterior
ano_gamma_post<-posterior_samples(gamma_brm)

#derived quantities
derived<-ano_gamma_post %>% mutate(Wr_mean = exp(b_Intercept), 
                          Slope = exp(b_Discharge_cent))

plot(density(derived$Wr_mean))


# Plot
conditional_effects(gamma_brm)
plot(conditional_effects(gamma_brm), points=T)


############################################################################
######### Length-Wr Relationship Development ################################################
############################################################################

N = 100  # number of simulations

y = 0.001x^2 + 75

# simulate priors
priors2 <- tibble(a = rnorm(N, 70, 15),
                 b = rnorm(N,0.001, 0.0005), ### playing with the prior. 
                 sigma = rexp(N, 0.1),
                 #shape = rexp(N, 2), # for brms
                 #scale = rexp(N, 2), # for rethinking
                 sim = 1:N)

# data (only the x values, since we're simulating y and mu and pretending we don't have them yet)
x2 <- rnorm(100, 0, 100)

# combine and simulate
prior_and_x2 <- priors2 %>% expand_grid(x = x2) %>%    # combine priors and x's
  mutate(mu = (a + b*(x-mean(x)) + b*(x-mean(x))^2),   # simulate regressions
         y = rnorm(nrow(.), mu, sigma))   # simulate data
         
# plot
prior_and_x2 %>% 
  ggplot(aes(x = (x-mean(x)), y = mu, group = sim)) + 
  geom_line() +
  theme_classic()+
  #geom_point(aes(y = y)) +
  labs(y = "sim") +
  #scale_y_log10() +
  NULL

# check raw data plots for curvature using splines
cats %>% mutate(L_cent = Length - mean(Length), Year_F=factor(Year)) %>% 
  ggplot(aes(x = L_cent, y = Wr, color = Year_F)) + 
  geom_point() +
  geom_smooth() # adds a quick spline fit. There is curvature, but not necessarialy parabolic curvature to any of these groups


# BJ's original model
L_Wr_fit <- brm(Wr ~ (L_cent^2)*Year_F, 
                data = cats %>% mutate(L_cent = Length - mean(Length), Year_F=factor(Year)),
                family = gaussian(),
                prior = c(prior(normal(70, 25), class = "Intercept"),
                          prior(normal(0.001, 0.001), class = "b"), # this prior is fairly weak with respect to curvature. 
                          prior(exponential(.1), class = "sigma")), 
                cores = 4, chains = 1, iter = 1000, 
                #file = "L_Wr_fit.rds"
                )

summary(L_Wr_fit)
plot(L_Wr_fit)
plot(conditional_effects(L_Wr_fit), points=T)


# model with varying intercepts
L_Wr_fit_var <- brm(Wr ~ (L_cent^2)*Year_F + (1|Year_F), # (1|Year_F) says to fit the model with varying intercepts for year 
                data = cats %>% mutate(L_cent = Length - mean(Length), Year_F=factor(Year)),
                family = gaussian(),
                prior = c(prior(normal(70, 25), class = "Intercept"),
                          prior(normal(0.001, 0.001), class = "b"),
                          prior(exponential(.1), class = "sigma"),
                          prior(exponential(0.1), class = "sd")),    # new prior for the sd of varying intercepts (check McElreath book). I'm not sure if this is a good prior or not 
                cores = 4, chains = 1, iter = 1000, 
                #file = "L_Wr_fit.rds"
)


summary(L_Wr_fit_var)
plot(L_Wr_fit_var)
plot(conditional_effects(L_Wr_fit_var, re_formula = NULL), points=T)


# model with varying intercepts and much stronger b prior
L_Wr_fit_var2 <- brm(Wr ~ L_cent*Year_F + I(L_cent^2)*Year_F + (1|Year_F), # (1|Year_F) says to fit the model with varying intercepts for year 
                    data = cats %>% mutate(L_cent = Length - mean(Length), Year_F=factor(Year)),
                    family = gaussian(),
                    prior = c(prior(normal(70, 25), class = "Intercept"),
                              prior(normal(0.001, 0.0005), class = "b"), ### play with this prior...It's currently very strong
                              prior(exponential(.1), class = "sigma"),
                              prior(exponential(0.1), class = "sd")),    # new prior for the sd of varying intercepts (check McElreath book). I'm not sure if this is a good prior or not 
                    cores = 4, chains = 1, iter = 1000, 
                    #file = "L_Wr_fit.rds"
)


summary(L_Wr_fit_var2)
plot(L_Wr_fit_var2)
plot(conditional_effects(L_Wr_fit_var2, re_formula = NULL), points=T)


