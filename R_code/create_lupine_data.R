
library(ipmr)
library(lme4)
library(dplyr)

# # Create a cleaner lupine data set for this tutorial

data <- read.csv("../learnr_ex/lupine_all.csv") %>%
  select(location, year, newid,
         log_area_t0,
         flow_t0,
         numrac_t0,
         surv_t1,
         log_area_t1,
         numrac_t0,
         numab_t0) %>%
  mutate(location = trimws(gsub("\\([0-9]\\)", "", location))) %>%
  filter(year > 2007)

seed <- read.csv("../learnr_ex/seedsperfruit.csv")
fruit <- read.csv("../learnr_ex/fruits_per_raceme.csv")
cons <- read.csv("../learnr_ex/consumption.csv") %>%
  group_by(Site) %>%
  summarise(mean_cons = mean(Mean_consumption, na.rm = TRUE)) %>%
  mutate(Site = toupper(Site))

abor <- data %>%
  group_by(location) %>%
  summarise(mean_abs = mean(numab_t0 / numrac_t0, na.rm = TRUE))

data <- left_join(data, cons, by = c("location" = "Site"))
data <- left_join(data, abor, by = "location")

mean_fruit <- mean(fruit$NumFruits)
mean_seed  <- mean(seed$SEEDSPERFRUIT)


data <- mutate(data,
               seeds = round((numrac_t0 * mean_cons * mean_abs) * mean_fruit * mean_seed)) %>%
  group_by(location) %>%
  slice_sample(prop = 0.6) %>%
  ungroup() %>%
  select(-c(numrac_t0, numab_t0:mean_abs))

write.csv(data,
          file = "ipmr_examples/data/lupine.csv",
          row.names = FALSE)

germ <- read.csv("../learnr_ex/seedbaskets.csv") %>%
  select(g0:g2) %>%
  summarise(g_0 = mean(g0),
            g_1 = mean(g1),
            g_2 = mean(g2)) %>%
  apply(2, FUN = function(x) x * 0.1)

data <- read.csv("ipmr_examples/data/lupine.csv")

## Models 

surv_mod <- glmer(surv_t1 ~ log_area_t0 + (1 | location) + (1 | year),
                  data = data,
                  family = binomial())

grow_mod <- lmer(log_area_t1 ~ log_area_t0 + (1 | location) + (1 | year),
                  data = data)

repr_mod <- glmer(flow_t0 ~ log_area_t0 + (1 | location) + (1 | year),
                  data = data,
                  family = binomial())

seed_mod <- glmer(seeds ~ log_area_t0 + (1 | location) + (1 | year),
                  data = data,
                  family = poisson())


# Coefficients into parameter list

s_pars <- as.list(fixef(surv_mod)) %>%
  setNames(paste("s_", c("int", "z"), sep = ""))
g_pars <- as.list(fixef(grow_mod)) %>%
  setNames(paste("g_", c("int", "z"), sep = "")) %>%
  c(., list(g_sd = sd(resid(grow_mod))))
r_p_pars <- as.list(fixef(repr_mod)) %>%
  setNames(paste("r_p_", c("int", "z"), sep = ""))
r_n_pars <- as.list(fixef(seed_mod)) %>%
  setNames(paste("r_n_", c("int", "z"), sep = ""))

s_yr_pars <- as.list(unlist(ranef(surv_mod)$year)) %>%
  setNames(paste("s_ran_", rownames(ranef(surv_mod)$year), sep = ""))
s_pl_pars <- as.list(unlist(ranef(surv_mod)$location)) %>%
  setNames(paste("s_ran_", rownames(ranef(surv_mod)$location), sep = ""))

g_yr_pars <- as.list(unlist(ranef(grow_mod)$year)) %>%
  setNames(paste("g_ran_", rownames(ranef(grow_mod)$year), sep = ""))
g_pl_pars <- as.list(unlist(ranef(grow_mod)$location)) %>%
  setNames(paste("g_ran_", rownames(ranef(grow_mod)$location), sep = ""))

r_p_yr_pars <- as.list(unlist(ranef(repr_mod)$year)) %>%
  setNames(paste("r_p_ran_", rownames(ranef(repr_mod)$year), sep = ""))
r_p_pl_pars <- as.list(unlist(ranef(repr_mod)$location)) %>%
  setNames(paste("r_p_ran_", rownames(ranef(repr_mod)$location), sep = ""))

r_n_yr_pars <- as.list(unlist(ranef(seed_mod)$year)) %>%
  setNames(paste("r_n_ran_", rownames(ranef(seed_mod)$year), sep = ""))
r_n_pl_pars <- as.list(unlist(ranef(seed_mod)$location)) %>%
  setNames(paste("r_n_ran_", rownames(ranef(seed_mod)$location), sep = ""))

## Seedling sizes 
sdl_size_pars <- list(sdl_mu = 2.7254,
                      sdl_sd = 0.9146)

germ          <- data.frame(
  g_0 = 0.02465088,
  g_1 = 0.01433191,
  g_2 = 0.0005732764
)

germ_pars     <- as.list(germ)

data_list <- c(
  s_pars, s_yr_pars, s_pl_pars,
  g_pars, g_yr_pars, g_pl_pars,
  r_p_pars, r_p_yr_pars, r_p_pl_pars,
  r_n_pars, r_n_yr_pars, r_n_pl_pars,
  germ_pars, sdl_size_pars
)

saveRDS(data_list, file = "ipmr_examples/data/lupine_pars.RDS")

par_set_inds <- list(yr = rownames(ranef(surv_mod)$year),
                     site = rownames(ranef(surv_mod)$location))

z_vec <- c(data$log_area_t0,
           data$log_area_t1)

L <- min(z_vec, na.rm = TRUE) * 1.2
U <- max(z_vec, na.rm = TRUE) * 1.2
m <- 200
## IPM 

lupine_ipm <- init_ipm("general", "di", "det") %>%
  define_kernel(
    name            = "P_yr_site",
    formula         = s_yr_site * G_yr_site * d_z,
    family          = "CC",
    s_yr_site       = plogis(s_int + s_ran_yr + s_ran_site + s_z * z_1),
    g_mu_yr_site    = g_int + g_ran_yr + g_ran_site + g_z * z_1,
    G_yr_site       = dnorm(z_2, g_mu_yr_site, g_sd),
    states          = list(c("z")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = TRUE,
    evict_fun       = truncated_distributions("norm", "G_yr_site")
  ) %>%
  define_kernel(
    name            = "F_yr_site",
    formula         = r_p_yr_site * r_n_yr_site * r_d * g_0,
    family          = "CC",
    r_p_yr_site     = plogis(r_p_int + r_p_ran_yr + r_p_ran_site + r_p_z * z_1),
    r_n_yr_site     = exp(r_n_int + r_n_ran_yr + r_n_ran_site + r_n_z * z_1),
    r_d             = dnorm(z_2, sdl_mu, sdl_sd),
    states          = list(c("z")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = TRUE,
    evict_fun       = truncated_distributions("norm", "r_d")
  ) %>%
  define_kernel(
    name            = "go_b1_yr_site",
    formula         = r_p_yr_site * r_n_yr_site * g_1,
    family          = "CD",
    r_p_yr_site     = plogis(r_p_int + r_p_ran_yr + r_p_ran_site + r_p_z * z_1),
    r_n_yr_site     = exp(r_n_int + r_n_ran_yr + r_n_ran_site + r_n_z * z_1),
    states          = list(c("z", "sb1")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = FALSE
  ) %>%
  define_kernel(
    name            = "go_b2_yr_site",
    formula         = r_p_yr_site * r_n_yr_site * g_2,
    family          = "CD",
    r_p_yr_site     = plogis(r_p_int + r_p_ran_yr + r_p_ran_site + r_p_z * z_1),
    r_n_yr_site     = exp(r_n_int + r_n_ran_yr + r_n_ran_site + r_n_z * z_1),
    states          = list(c("z", "sb2")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = FALSE
  ) %>%
  define_kernel(
    name      = "b2_to_b1",
    formula   = 1,
    family    = "DD",
    states    = list(c("sb1", "sb2")),
    data_list = list(),
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name      = "b1_to_plants",
    formula   = r_d,
    family    = "DC",
    r_d       = dnorm(z_2, sdl_mu, sdl_sd),
    states    = list(c("sb1", "z")),
    data_list = data_list,
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_yr_site","F_yr_site","go_b1_yr_site", "go_b2_yr_site",
                       "b2_to_b1", "b1_to_plants"),
      int_rule     = rep("midpoint", 6),
      state_start  = c("z", "z", "z", "z", "sb2", "sb1"),
      state_end    = c("z", "z", "b1", "sb2", "sb1", "z")
    )
  ) %>%
  define_domains(
    z = c(L, U, m)
  ) %>%
  define_pop_state(
    n_z_yr_site   = rep(1/202, 200),
    n_sb1_yr_site = 1/202,
    n_sb2_yr_site = 1/202
  ) %>%
  make_ipm(
    iterations      = 100,
    return_all_envs = TRUE
  )

range(lambda(lupine_ipm))
hist(lambda(lupine_ipm), breaks = 20)


lupine_stoch_ipm <- init_ipm("general", "di", "stoch", "kern") %>%
  define_kernel(
    name            = "P_yr_site",
    formula         = s_yr_site * G_yr_site * d_z,
    family          = "CC",
    s_yr_site       = plogis(s_int + s_ran_yr + s_ran_site + s_z * z_1),
    g_mu_yr_site    = g_int + g_ran_yr + g_ran_site + g_z * z_1,
    G_yr_site       = dnorm(z_2, g_mu_yr_site, g_sd),
    states          = list(c("z")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = TRUE,
    evict_fun       = truncated_distributions("norm", "G_yr_site")
  ) %>%
  define_kernel(
    name            = "F_yr_site",
    formula         = r_p_yr_site * r_n_yr_site * r_d * g_0,
    family          = "CC",
    r_p_yr_site     = plogis(r_p_int + r_p_ran_yr + r_p_ran_site + r_p_z * z_1),
    r_n_yr_site     = exp(r_n_int + r_n_ran_yr + r_n_ran_site + r_n_z * z_1),
    r_d             = dnorm(z_2, sdl_mu, sdl_sd),
    states          = list(c("z")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = TRUE,
    evict_fun       = truncated_distributions("norm", "r_d")
  ) %>%
  define_kernel(
    name            = "go_b1_yr_site",
    formula         = r_p_yr_site * r_n_yr_site * g_1,
    family          = "CD",
    r_p_yr_site     = plogis(r_p_int + r_p_ran_yr + r_p_ran_site + r_p_z * z_1),
    r_n_yr_site     = exp(r_n_int + r_n_ran_yr + r_n_ran_site + r_n_z * z_1),
    states          = list(c("z", "sb1")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = FALSE
  ) %>%
  define_kernel(
    name            = "go_b2_yr_site",
    formula         = r_p_yr_site * r_n_yr_site * g_2,
    family          = "CD",
    r_p_yr_site     = plogis(r_p_int + r_p_ran_yr + r_p_ran_site + r_p_z * z_1),
    r_n_yr_site     = exp(r_n_int + r_n_ran_yr + r_n_ran_site + r_n_z * z_1),
    states          = list(c("z", "sb2")),
    data_list       = data_list,
    uses_par_sets   = TRUE,
    par_set_indices = par_set_inds,
    evict_cor       = FALSE
  ) %>%
  define_kernel(
    name      = "b2_to_b1",
    formula   = 1,
    family    = "DD",
    states    = list(c("sb1", "sb2")),
    data_list = list(),
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name      = "b1_to_plants",
    formula   = r_d,
    family    = "DC",
    r_d       = dnorm(z_2, sdl_mu, sdl_sd),
    states    = list(c("sb1", "z")),
    data_list = data_list,
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_yr_site","F_yr_site","go_b1_yr_site", "go_b2_yr_site",
                       "b2_to_b1", "b1_to_plants"),
      int_rule     = rep("midpoint", 6),
      state_start  = c("z", "z", "z", "z", "sb2", "sb1"),
      state_end    = c("z", "z", "b1", "sb2", "sb1", "z")
    )
  ) %>%
  define_domains(
    z = c(L, U, m)
  ) %>%
  define_pop_state(
    n_z   = rep(1/202, 200),
    n_sb1 = 1/202,
    n_sb2 = 1/202
  ) %>%
  make_ipm(
    iterations      = 1000,
    return_all_envs = TRUE
  )

lambda(lupine_stoch_ipm)
