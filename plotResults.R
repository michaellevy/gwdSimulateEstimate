filter(tab, !is.na(std.error), estimate < 50) %>%
  ggplot(aes(x = theta_s, y = estimate, color = p.value)) + 
  geom_line() + geom_hline(aes(yintercept = 0), linetype = "dotted") +
  facet_grid(meanDegree ~ degPop, scales = "free_y")