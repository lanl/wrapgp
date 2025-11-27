library(dplyr)
library(janitor)
library(ggplot2)
library(tidyr)

results <- read.csv("results.csv") %>% 
  mutate(method = factor(method, 
                         levels = c("wgp_ours", "wgp_jonas", "mallasto"),
                         labels = c("WGP (Ours)", "WGP (Jona-Lasinio)", "WGP (Mallasto-Feragen)"))) %>% 
  rename(Model = method) %>% 
  gather(metric, value, rmse, crps)

cols = c("WGP (Ours)" = "red", "WGP (Jona-Lasinio)" = "green", "WGP (Mallasto-Feragen)" = "blue")
pdf("log_results1.pdf", width = 5, height = 5)
results %>% 
  filter(metric == "rmse") %>% 
  ggplot(aes(x = factor(n), y = value, fill = Model)) +
  geom_boxplot() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 2) +
  xlab("Training Size") + ylab("RMSE") +
  #coord_cartesian(ylim = c(0, 2.2)) +
  theme_bw(base_size = 15) + 
  scale_fill_manual(values = cols) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.70, 0.6), text = element_text(size = 15),
        legend.background = element_rect(fill = "white", linetype = "solid", colour = "black"))
dev.off()

pdf("log_results2.pdf", width = 5, height = 5)
results %>% 
  filter(metric == "crps") %>% 
  ggplot(aes(x = factor(n), y = value, fill = Model)) +
  geom_boxplot() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 2) +
  xlab("Training Size") + ylab("CRPS") +
  #coord_cartesian(ylim = c(0.8, 1)) +
  theme_bw(base_size = 15) + 
  scale_fill_manual(values = cols, drop = F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none", text = element_text(size = 15))
dev.off()
