# Plot of distribution effects - sensitivity to assumptions
# closest vs any individual triggering detection
# Uniform vs clustered distribution

freq_res<- readRDS("outputs/freq_distribtion.rds")

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

##-------------------------------------------
# Closest vs any individual triggering detection

win.graph(10,10)
freq_res %>% mutate(sigma = paste("\u03C3 =",sigma),
                    Model = factor(Model, levels = c("Closest distance","All distances")),
                    clustered = factor(clustered, labels=c("uniform","clustered")),
                    any = factor(any, labels = c("closest only","any encounter"))) %>%
  filter(clustered == "uniform") %>%
  ggplot(aes(Model, Est, fill = any)) +
  geom_violin(trim=TRUE) +
  stat_summary(aes(fill=any),position=position_dodge(width=0.9),fun = median,
               geom = "point", size=2) +
  geom_hline(yintercept = 0,linetype=2,color="red") +
  facet_grid(Group ~ sigma) +
  scale_y_continuous(breaks = seq(-0.2, 1.8, 0.2), limits=c(-0.2, 1.8)) +
  labs(x = "Model", y = "Relative bias",fill="Detection") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text = element_text(face="bold",size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold"),
        legend.position = "bottom")


# uniform vs clustered distribution
win.graph(10,10)
freq_res %>% mutate(sigma = paste("\u03C3 =",sigma),
                    Model = factor(Model, levels = c("Closest distance","All distances")),
                    clustered = factor(clustered, labels=c("uniform","clustered")),
                    any = factor(any, labels = c("closest only","any encounter"))) %>%
  filter(any == "closest only") %>%
  ggplot(aes(Model, Est, fill = clustered)) +
  geom_violin(trim=TRUE) +
  stat_summary(aes(fill=clustered),position=position_dodge(width=0.9),fun = median,
               geom = "point", size=2) +
  geom_hline(yintercept = 0,linetype=2,color="red") +
  facet_grid(Group ~ sigma) +
  scale_y_continuous(breaks = seq(-0.2, 1.8, 0.2), limits=c(-0.2, 1.8)) +
  labs(x = "Model", y = "Relative bias",fill="Distribution") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text = element_text(face="bold",size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold"),
        legend.position = "bottom")


