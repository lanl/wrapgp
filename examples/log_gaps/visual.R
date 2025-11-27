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
pdf("log_gap_results1.pdf", width = 5, height = 5)
results %>% 
  filter(metric == "rmse") %>% 
  ggplot(aes(x = factor(n), y = value, fill = Model)) +
  geom_boxplot() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 2) +
  xlab("Training Size") + ylab("RMSE") +
  #coord_cartesian(ylim = c(0, 1.3)) +
  theme_bw(base_size = 15) + 
  scale_fill_manual(values = cols) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(legend.position = "none", text = element_text(size = 15))
  theme(legend.position = c(0.70, 0.63), text = element_text(size = 15),
        legend.background = element_rect(fill = "white", linetype = "solid", colour = "black"))
dev.off()

pdf("log_gap_results2.pdf", width = 5, height = 5)
results %>% 
  filter(metric == "crps") %>% 
  ggplot(aes(x = factor(n), y = value, fill = Model)) +
  geom_boxplot() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 2) +
  xlab("Training Size") + ylab("CRPS") +
  theme_bw(base_size = 15) + 
  scale_fill_manual(values = cols) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")
dev.off()

dat <- read.csv("data/train_n250_seed4.csv")

pdf("log_gap_ex1.pdf", width = 5, height = 5)
plot(dat$X, dat$Y, pch = 20, xlab = "X", ylab = "Y")
polygon(c(0.15, 0.25, 0.25, 0.15), c(rep(10, 2), rep(-1, 2)),
        col = "grey", border = NA)
polygon(c(0.8, 0.9, 0.9, 0.8), c(rep(10, 2), rep(-1, 2)),
        col = "grey", border = NA)
legend("topright", legend = "No data", col = "grey", lwd = 10, bg = "white")
dev.off()

pred1 <- read.csv("pred/wgp_ours_n250_seed4.csv") %>% arrange(X)
pred2 <- read.csv("pred/wgp_jonas_n250_seed4.csv") %>% arrange(X)
pred3 <- read.csv("pred/mallasto_n250_seed4.csv") %>% arrange(X)

pdf("log_gap_ex2.pdf", width = 5, height = 5)
plot(dat$X, dat$Y, pch = 20, xlab = "X", ylab = "Y")
lines(pred1$X, pred1$m, col = "red", lwd = 2)
lines(pred2$X, pred2$m, col = "green", lwd = 2)
lines(pred3$X, pred3$m, col = "blue", lwd = 2)
legend("bottomright", legend = c("WGP (Ours)", "WGP (Jona-Lasinio)", "WGP (Mallasto-Feragen)"),
       col = c("red", "green", "blue"), lwd = 2, bg = "white", cex = 1)
dev.off()
