bf_census_f %>% pivot_longer(cols = -1) %>%  mutate(source = "UNPD", sex = "female", Age = age) %>%
  bind_rows(
    bf_census_m %>% pivot_longer(cols = -1) %>%  mutate(source = "UNPD", sex = "male", Age = age)
  ) %>%
bind_rows(
bf.pop.aggr.f %>% pivot_longer(cols = -1) %>%  mutate(source = "WPP", sex = "female", Age = as.numeric(gsub("-\\d+","",Age))) %>%
bind_rows(
  bf.pop.aggr.m %>% pivot_longer(cols = -1) %>%  mutate(source = "WPP", sex = "male", Age = as.numeric(gsub("-\\d+","",Age)))
  )
) %>%
bind_rows(
cbind(burkina.faso.females$baseline.pop.counts, burkina.faso.females$census.pop.counts) %>% as_tibble() %>%
  mutate(Age = seq(0, 80, by = 5), .before = 1) %>% pivot_longer(cols=-1) %>%
  mutate(source = "Mark", sex = "female")
) %>% select(-age) %>%
  group_by(source, name, sex) %>%
  summarise_at(vars(value),sum) %>% ungroup() %>%
  mutate(name = as.numeric(name)) %>%
  ggplot() + geom_line(aes(x = name, y = value, linetype = source, col = sex), lwd = 1.2)


haha <- fit.LQ.both.hmean.extrapolated
haha$par.full %>% split(names(.)) %>% .$h_params_m

get.h <- function(x) {
  x$par.full %>% split(names(.)) %>% .$h_params_f %>% as_tibble() %>% mutate(sex="female", year = bf.idx5$periods) %>%
    bind_rows(
      x$par.full %>% split(names(.)) %>% .$h_params_m %>% as_tibble() %>% mutate(sex="male", year = bf.idx5$periods)
    )
}

h.df <- get.h(fit.LQ.both.common.mean.extrapolated) %>% mutate(source = "AR  common mean") %>%
        bind_rows(get.h(fit.LQ.ARIMA) %>% mutate(source = "ARIMA"),
            get.h(fit.LQ.both.extrapolated) %>% mutate(source = "AR IGME mean"),
#            get.h(fit.LQ.both.hmean.extrapolated) %>% mutate(source = "hmean IGME"),
#            get.h(fit.LQ.both.vec.igme.MVN.prior) %>% mutate(source = "MVN IGME mean"),
#            get.h(fit.LQ.both.vec.igme.MVN.prior.noDHS) %>% mutate(source = "MVN IGME without DHS"),
            bind_rows(as_tibble(data.vec$h_mean_f) %>% mutate(sex = "female", year = bf.idx5$periods, source = "IGME Estimates"),
                      as_tibble(data.vec$h_mean_m) %>% mutate(sex = "male", year = bf.idx5$periods, source = "IGME Estimates")
            )
            )


ggplot(h.df) + geom_line(aes(x = year, y = value, col = sex, linetype = source), lwd = 1.2)



