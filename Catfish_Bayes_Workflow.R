library(FSA)
library(tidyverse)
library(magrittr)
library(brms)
options(buildtools.check = function(action) TRUE )

cat<-read_csv("Scotland All Years.csv")
water<-read_csv("WaterLevels.csv")


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
cats<-filterD(cats, !is.na(Wr) & Wr>=50 & Wr<=150)
 
cats$Discharge_cent<-(log10(cats$Discharge)-(mean(log10(cats$Discharge))))/sd(log10(cats$Discharge))

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

plot(Wr~(Discharge_cent), data=cats)

cats$dis<-log10(cats$Discharge)-mean(log10(cats$Discharge))
plot(Wr~dis, data=cats)


N = 100  # number of simulations

# simulate priors
priors <- tibble(a = rnorm(N, 4, .25),
                 b = rnorm(N,0, 1),
                 shape = rexp(N, .05),# for brms
                 scale = rexp(N, 0.5), # for Rethinking
                 sim = 1:N)

# data (only the x values, since we're simulating y and mu and pretending we don't have them yet)
x <- rnorm(N, 1500,500)

# combine and simulate
prior_and_x <- priors %>% expand_grid(x = x) %>%    # combine priors and x's
  mutate(mu = exp(a + b*(log10(x)-log10(mean(x)))/log10(sd(x))),                         # simulate regressions
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
gamma_brm <- brm(Wr ~ log10(Discharge), 
                 family = Gamma(link = "log"),
                 data = cats,
                 prior = c(prior(normal(4,.25), class = "Intercept"),
                           prior(normal(0, 1), class = "b"),
                           prior(exponential(.05), class = "shape")),
                 cores = 4, chains = 1, iter = 1000,
                 file="gamma_brm.rds")

summary(gamma_brm)
plot(gamma_brm)

#sample from the posterior
gamma_post<-posterior_samples(gamma_brm)

#derived quantities
derived<-gamma_post %>% mutate(Wr_mean = exp(b_Intercept), 
                          Slope = exp(b_Discharge_cent))
plot(density(derived$Wr_mean))

# Plot
conditional_effects(gamma_brm)
plot(conditional_effects(gamma_brm), points=T)

#############################################################################
####### Fit a gompertz model
#############################################################################
N = 100  # number of simulations

# simulate priors
priors <- tibble(A = rnorm(N, 120, 10),
                 k = rnorm(N,1, .5),
                 B = rnorm(N, 3,1.5),
                 sigma = exp(0.1),
                 sim = 1:N)

# data (only the x values, since we're simulating y and mu and pretending we don't have them yet)
xs <- seq(0, 5,.1)
# combine and simulate
prior_comb <- priors %>% expand_grid(x = xs) %>%    # combine priors and x's
  mutate(mu = A * exp(-exp(-(k*(x-B)))),                         # simulate regressions
         y = rnorm(nrow(.), mu, sigma))
         
# plot
prior_comb %>% 
  ggplot(aes(x = x, y = mu, group = sim)) + 
  geom_line() +
  theme_classic()+
  geom_point(aes(y = y)) +
  labs(y = "sim") +
  #scale_y_log10() +
  NULL

form_1 <- bf(Wr ~ A * exp( -exp( -(k * (Discharge - delay) ) ) ),
             A ~ 1 ,
             k ~ 1 ,
             delay ~ 1 ,
             nl = TRUE)


mod1 <- brm(form_1 , 
            data = cats,
            prior = c(prior(normal(120, 10), nlpar = "A", lb=0),
                      prior(normal(2, 1), nlpar = "k", lb=0),
                      prior(normal(4, 2), nlpar = "delay", lb=0),
                      prior(exponential(0.1), class="sigma")),
            family = gaussian(),
            chains = 1, cores=4)
summary(mod1)
plot(conditional_effects(mod1), points = TRUE)

##
Dc<-(cats$Discharge - mean(cats$Discharge))/sd(cats$Discharge)
Wrc<-(cats$Wr - mean(cats$Wr))/sd(cats$Wr)
plot(cats$Discharge, cats$Wr)
#y_r<-(min(Wrc)+(max(Wrc)-min(Wrc))*(1/(1+exp(-4*(Wrc-mean(Wrc))))))
gomp<-90*exp(-exp(-7*(Dc-3604)))
lines(cats$Discharge,gomp, col="blue")

library(FSA)
vbFuns()
gomp1 <- GompertzFuns(msg=TRUE)
(vbs <- GompertzStarts(Wr~Dc,type="Typical",plot=F))
vbf <- nls(mean~vb(AGE,Linf,K,t0),data=Mean)

age<-cats$Wr
plot(gomp1(Dc,Linf=140,gi=1,ti=-0.5)~Dc,type="b",pch=19)
fit1 <- nls(age~gomp1(Dc,Linf,gi,ti),start=list(Linf=120,gi=1,ti=-0.5))
summary(fit1,correlation=TRUE)
curve(gomp1(x,Linf=coef(fit1)),from=-3,to=3,col="red",lwd=10,add=TRUE)

############################################################################
######### Length-Wr Relationship Development ################################################
############################################################################

N = 100  # number of simulations

# simulate priors
priors2 <- tibble(a = rnorm(N, 70, 15),
                 b = rnorm(N,0.0005, 0.0002), ### playing with the prior. 
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
L_Wr_fit <- brm(Wr ~ Year_F+I(L_cent^2):Year_F, 
                data = cats %>% mutate(L_cent = Length - 500, Year_F=factor(Year)),
                family = gaussian(),
                prior = c(prior(normal(70, 20), class = "Intercept"),
                          prior(normal(0.0005, 0.0002), class = "b"), # this prior is fairly weak with respect to curvature. 
                          prior(exponential(.1), class = "sigma")), 
                cores = 4, chains = 1, iter = 1000, 
                #file = "L_Wr_fit.rds"
                )

summary(L_Wr_fit)
plot(L_Wr_fit)
plot(conditional_effects(L_Wr_fit), points=T)


# model with varying intercepts
L_Wr_fit_var <- brm(Wr ~ I(L_cent^2):Year_F + (1|Year_F), # (1|Year_F) says to fit the model with varying intercepts for year 
                data = cats %>% mutate(L_cent = Length - 500, Year_F=factor(Year)),
                family = gaussian(),
                prior = c(prior(normal(70, 20), class = "Intercept"),
                          prior(normal(0.0005, 0.0002), class = "b"),
                          prior(exponential(.1), class = "sigma"),
                          prior(exponential(0.1), class = "sd")),    # new prior for the sd of varying intercepts (check McElreath book). I'm not sure if this is a good prior or not 
                cores = 4, chains = 1, iter = 1000, 
                #file = "L_Wr_fit.rds"
)


summary(L_Wr_fit_var)
plot(L_Wr_fit_var)
plot(conditional_effects(L_Wr_fit_var), points=T)


# model with varying intercepts and much stronger b prior
L_Wr_fit_var2 <- brm(Wr ~ L_cent*Year_F + I(L_cent^2)*Year_F + Year_F, # (1|Year_F) says to fit the model with varying intercepts for year 
                    data = cats %>% mutate(L_cent = Length - mean(Length), Year_F=factor(Year)),
                    family = gaussian(),
                    prior = c(prior(normal(70, 25), class = "Intercept"),
                              prior(normal(0.001, 0.0005), class = "b"), ### play with this prior...It's currently very strong
                              prior(exponential(.1), class = "sigma")),    # new prior for the sd of varying intercepts (check McElreath book). I'm not sure if this is a good prior or not 
                    cores = 4, chains = 1, iter = 1000, 
                    #file = "L_Wr_fit.rds"
)


summary(L_Wr_fit_var2)
plot(L_Wr_fit_var2)
plot(conditional_effects(L_Wr_fit_var2, re_formula = NULL), points=T)


# simple linear model with varying intercepts 
L_Wr_fit2 <- brm(Wr ~ L_cent*Year_F+Year_F, 
                data = cats %>% mutate(L_cent = Length - 500, Year_F=factor(Year)),
                family = gaussian(),
                prior = c(prior(normal(90, 20), class = "Intercept"),
                          prior(normal(1, 0.5), class = "b"), # this prior is fairly weak with respect to curvature. 
                          prior(exponential(.1), class = "sigma")), 
                cores = 4, chains = 1, iter = 1000, 
                #file = "L_Wr_fit.rds"
)
plot(conditional_effects(L_Wr_fit2, re_formula = NULL), points=T)
