
# Accounting for multiple individuals in a camera FOV by
# determining the availability of the closest individual
# using order statistics

library(tidyverse)
library(tidyr)
library(Rcpp)

sourceCpp("src/generate_detections.cpp")
source("R/group_size_funcs.r")


##----  Continuous data - plots of adjusted availability ----

grp_size<- c(1,3,6,9)
w<- 12
n_entries<- 1e6
grp_pars<- expand.grid(grp = grp_size, w=w)
n<- nrow(grp_pars)

gs_plots<- list()

for(i in 1:nrow(grp_pars)) {
  y_data<- generate_closest_dist(entries = n_entries,
                                 grp_size = grp_pars$grp[i],
                                 w = grp_pars$w[i])
  y_data<- data.frame(distance = y_data)
  x<- seq(0, grp_pars$w[i], length=100)
  Ax<- availability_cont(x, grp_pars$w[i], grp_pars$grp[i])
  Ax_data<- data.frame(x=x,Ax=Ax)

  gs_plots[[i]]<- ggplot(aes(distance), data=y_data) +
    geom_histogram(aes(y = after_stat(density)), binwidth=0.25, color="white",fill="grey50") +
    geom_line(aes(x=x, y=Ax), data=Ax_data, color="red", linewidth = 1) +
    scale_x_continuous(breaks = seq(0,grp_pars$w[i],1),limits = c(0,grp_pars$w[i])) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(x="Distance to closest individual (m)",y="PDF") +
    theme_bw() +
    theme(axis.title.x = element_text(face="bold", size=14),
          axis.title.y = element_text(face="bold", size=14),
          axis.text = element_text(size=12))
}

win.graph(10,10)
cowplot::plot_grid(plotlist = gs_plots, labels = paste("n = ",grp_pars$grp),
                   label_x=0.4, label_y=0.95, ncol = 2)

ggsave(filename = "outputs/group_size_PDF.png",
       width = 10, height = 10, units = "in")


##---- Conditional likelihood for continuous distance sampling data ----
##---- Bias when using all distances vs closest distance

sigma<- c(3,4,5)
w<- 12
grp_size<- 1:10
gr<- hn_func
entries<- 1e6
clustered = FALSE
any = FALSE

entries_per_grp<- rep(entries, length(grp_size))

res<- vector("list", length(sigma))

for(j in seq_along(sigma)) {

  cat("doing sigma, ",sigma[j],"\n")
  y<- generate_detections(entries_per_grp, gr=gr,
                          sigma=sigma[j],
                          w=w, closest=TRUE,
                          clustered = clustered,
                          any=any)

  ngrp<- length(y)

  ##  Closest distance and adjusted availability

  sigma_mle<- matrix(NA, ngrp, 2)

  for(i in 1:ngrp) {
    mle <- optim(
      par = 1,
      fn = nll.cond.point.hn,
      x = y[[i]]$distance,
      gs = grp_size[i],
      w=w,
      method="Brent",
      lower = -5,
      upper = 5,
      hessian = TRUE
    )
    H<- mle$hessian
    sigma_mle[i,1]<- exp(mle$par)
    sigma_mle[i,2]<- exp(mle$par) * sqrt(solve(H))
  }

  colnames(sigma_mle)<- c("Est","SE")
  sigma_adj<- as_tibble(sigma_mle) %>%
    mutate(Model = "Closest distance",
           group = grp_size)


  ## All distances and standard availability (i.e. gs=1)

  y<- generate_detections(entries_per_grp, gr=gr,
                          sigma=sigma[j],
                          w=w, closest=FALSE,
                          clustered = clustered,
                          any = any)

  sigma_mle<- matrix(NA, ngrp, 2)

  for(i in 1:ngrp) {
    mle <- optim(
      par = 1,
      fn = nll.cond.point.hn,
      x = y[[i]]$distance,
      gs = 1,
      w=w,
      method = "Brent",
      lower= -5,
      upper = 5,
      hessian = TRUE
    )

    H<- mle$hessian
    sigma_mle[i,1]<- exp(mle$par)
    sigma_mle[i,2]<- exp(mle$par) * sqrt(solve(H))
  }

  colnames(sigma_mle)<- c("Est","SE")
  sigma_unadj<- as_tibble(sigma_mle) %>%
    mutate(Model = "All distances",
           group = grp_size)

  sigma_est<- bind_rows(sigma_adj,sigma_unadj)
  sigma_est<- sigma_est %>% mutate(Est = (Est - sigma[j])/sigma[j],
                                   SE = SE/sigma[j],
                                   sigma=sigma[j])
  res[[j]]<- sigma_est

}

res<- list_rbind(res)
res_all<- res %>% mutate(sigma = paste("\u03C3 =",sigma))

win.graph(10,5)
res_all %>%
  ggplot(aes(group, Est, color=Model)) +
  geom_pointrange(aes(ymin = Est - 2*SE, ymax=Est + 2*SE),
                  position=position_dodge(width=0.2), size=0.7, linewidth = 1) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ sigma) +
  scale_x_continuous(breaks = grp_size,labels = grp_size) +
  scale_y_continuous(breaks = seq(-0.5, 5, 0.5),limits=c(-0.5,5)) +
  labs(x = "Group size", y = "Relative bias") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold"),
        legend.position = "bottom")

ggsave(filename = "outputs/LL_cont.png",
       width = 10, height = 5, units = "in")


##---- Binned distances and conditional likelihood ----

sigma<- c(3,4,5)
delta<- c(0.5,1,2,3)
dpars<- expand.grid(sigma=sigma,delta=delta)
w<- 12
grp_size<- 1:10
clustered<- FALSE
any<- FALSE

gr<- hn_func
entries<- 1e6
entries_per_grp<- rep(entries, length(grp_size))

res<- vector("list", nrow(dpars))

for(j in 1:nrow(dpars)) {
  cat(sprintf("\rParameter %d/%d   ", j, nrow(dpars)))
  flush.console()

bin_bp <- seq(0, w, by = dpars$delta[j])  # bin breakpoints
sigma_true<- dpars$sigma[j]

y<- generate_detections(entries_per_grp, gr=gr,
                        sigma=sigma_true,
                        w=w, closest=TRUE,
                        binned=TRUE, breaks=bin_bp,
                        clustered=clustered,
                        any=any)

ngrp<- length(y)

##  Do adjusted availability

sigma_mle<- matrix(NA, ngrp, 2)

for(i in 1:ngrp) {
  mle <- optim(
    par = 1,
    fn = nll.cond.binned.hn,
    counts = y[[i]],
    breaks=bin_bp,
    gs = grp_size[i],
    method="Brent",
    lower = -5,
    upper = 5,
    hessian = TRUE
  )

  H<- mle$hessian
  sigma_mle[i,1]<- exp(mle$par)
  sigma_mle[i,2]<- exp(mle$par) * sqrt(solve(H))
}

colnames(sigma_mle)<- c("Est","SE")
sigma_adj<- as_tibble(sigma_mle) %>%
  mutate(Model = "Closest distance",
         group = grp_size)

## Now unadjusted i.e. gs==1

y<- generate_detections(entries_per_grp, gr=gr,
                        sigma=sigma_true,
                        w=w, closest=FALSE,
                        binned=TRUE, breaks=bin_bp,
                        clustered=clustered,
                        any=any)

sigma_mle<- matrix(NA, ngrp, 2)

for(i in 1:ngrp) {
  mle <- optim(
    par = 1,
    fn = nll.cond.binned.hn,
    counts = y[[i]],
    gs=1,
    breaks = bin_bp,
    method="Brent",
    lower = -5,
    upper = 5,
    hessian = TRUE
  )
  H<- mle$hessian
  sigma_mle[i,1]<- exp(mle$par)
  sigma_mle[i,2]<- exp(mle$par) * sqrt(solve(H))
}

colnames(sigma_mle)<- c("Est","SE")
sigma_unadj<- as_tibble(sigma_mle) %>%
  mutate(Model = "All distances",
         group = grp_size)


sigma_est<- bind_rows(sigma_adj,sigma_unadj)
sigma_est<- sigma_est %>% mutate(Est = (Est - sigma_true)/sigma_true,
                                 SE = SE/sigma_true,
                                 sigma=sigma_true,
                                 delta=dpars$delta[j])
res[[j]]<- sigma_est
}

res<- list_rbind(res)
res<- res %>% mutate(sigma = paste("\u03C3 =",sigma))

win.graph(10,7)
res %>% mutate(delta = factor(delta)) %>%
  ggplot(aes(group, Est, color=delta)) +
  geom_pointrange(aes(ymin = Est - 2*SE, ymax=Est + 2*SE),
                  position=position_dodge(width=0.5), size=0.7, linewidth = 1) +
  geom_hline(yintercept = 0) +
  facet_grid(Model ~ sigma) +
  scale_x_continuous(breaks = grp_size,labels = grp_size) +
  scale_y_continuous(breaks = round(seq(-0.5, 5, 0.5),2),limits=c(-0.5, 5)) +
  labs(x = "Group size", y = "Relative bias", color="Bin width (m)") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold"),
        legend.position = "bottom")

ggsave(filename = "outputs/LL_bins.png",
       width = 10, height = 7, units = "in")


##---- Sims for variable group size frequencies ----

sigma<- c(3,4,5)
delta<- 2
w<- 12
entries<- 2e3
nsims<- 500
gr<- hn_func
clustered<- c(FALSE,TRUE)
any<- c(FALSE,TRUE)

dpars<- expand.grid(sigma=sigma, clustered = clustered, any=any)
npars<- nrow(dpars)

grp_freq<- list("Asocial" = c(0.97,0.02,0.01),
                "Fallow deer" = c(0.55,0.3,0.12,0.02,0.01),
                "Feral pig" = c(0.5,0.20,0.15,0.05,0.05,0.025,0.025))

bin_bp <- seq(0, w, by = delta)  # bin breakpoints

freq_list<- vector("list", length(grp_freq))

for(g in seq_along(grp_freq)) {

  entries_per_group<- entries * grp_freq[[g]]

  res<- vector("list", nsims)

  for(i in 1:nsims) {
    cat(sprintf("\rGroup %d/%d  |  Sim %d/%d   ", g, length(grp_freq), i, nsims))
    flush.console()

    ## Using closest distance and  adjusted availability

    sigma_mle<- matrix(NA, npars, 5)
    grp_size<- seq_along(entries_per_group)

    for(j in 1:npars) {

      y<- generate_detections(entries_per_group, gr=gr,
                              sigma=dpars$sigma[j],
                              w=w, closest=TRUE,
                              binned=TRUE, breaks=bin_bp,
                              any=dpars$any[j],
                              clustered=dpars$clustered[j])

      mle <- optim(
        par = 1,
        fn = nll.cond.binned.hn,
        counts = y,
        gs=grp_size,
        breaks = bin_bp,
        method = "Brent",
        lower= -5,
        upper= 5,
        hessian = TRUE
      )

      H<- mle$hessian
      sigma_mle[j,1]<- exp(mle$par)
      sigma_mle[j,2]<- exp(mle$par) * sqrt(solve(H))
      sigma_mle[j,3]<- dpars$sigma[j]
      sigma_mle[j,4]<- as.numeric(dpars$clustered[j])
      sigma_mle[j,5]<- as.numeric(dpars$any[j])

    }
    colnames(sigma_mle)<- c("Est","SE","sigma","clustered","any")
    sigma_adj<- as_tibble(sigma_mle) %>%
      mutate(Model = "Closest distance")


    ## using all distances and standard availability i.e., gs==1

    sigma_mle<- matrix(NA, npars, 5)
    grp_size<- rep(1, length(entries_per_group))

    for(j in 1:npars) {

      y<- generate_detections(entries_per_group, g=gr,
                              sigma=dpars$sigma[j],
                              w=w, closest=FALSE,
                              binned=TRUE, breaks=bin_bp,
                              any=dpars$any[j],
                              clustered=dpars$clustered[j])

      mle <- optim(
        par = 1,
        fn = nll.cond.binned.hn,
        counts = y,
        gs=grp_size,
        breaks = bin_bp,
        method = "Brent",
        lower= -5,
        upper= 5,
        hessian = TRUE
      )

      H<- mle$hessian
      sigma_mle[j,1]<- exp(mle$par)
      sigma_mle[j,2]<- exp(mle$par) * sqrt(solve(H))
      sigma_mle[j,3]<- dpars$sigma[j]
      sigma_mle[j,4]<- as.numeric(dpars$clustered[j])
      sigma_mle[j,5]<- as.numeric(dpars$any[j])

    }
    colnames(sigma_mle)<- c("Est","SE","sigma","clustered","any")
    sigma_unadj<- as_tibble(sigma_mle) %>%
      mutate(Model = "All distances")

    sigma_est<- bind_rows(sigma_adj,sigma_unadj)
    sigma_est<- sigma_est %>% mutate(Est = (Est - sigma)/sigma,
                                     SE = SE/sigma,
                                     sigma=sigma)
    res[[i]]<- sigma_est
  }

  res<- list_rbind(res) %>%
    mutate(Group = names(grp_freq)[g])

  freq_list[[g]]<- res

}

freq_res<- list_rbind(freq_list)

saveRDS(freq_res, "outputs/freq_distribtion.rds")

# Plot of distribution effects
#freq_res<- readRDS("outputs/freq_distribtion.rds")

win.graph(12,5)
freq_res %>% mutate(sigma = paste("\u03C3 =",sigma),
                    Model = factor(Model, levels = c("Closest distance","All distances")),
                    clustered = factor(clustered, labels=c("uniform","clustered")),
                    any = factor(any, labels = c("closest only","any encounter"))) %>%
  filter(clustered == "uniform" & any == "closest only") %>%
  ggplot(aes(Group, Est, fill = Model)) +
  geom_violin(trim=TRUE) +
  stat_summary(aes(fill=Model),position=position_dodge(width=0.9),fun = median,
               geom = "point", size=2) +
  geom_hline(yintercept = 0,linetype=2,color="red") +
  facet_wrap(~sigma) +
  scale_y_continuous(breaks = seq(-0.2, 1.8, 0.2), limits=c(-0.2, 1.8)) +
  labs(x = "Species", y = "Relative bias") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text = element_text(face="bold",size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold"),
        legend.position = "bottom")

ggsave(filename = "outputs/freq_closest_uniform.png",
       width = 12, height = 5, units = "in")
