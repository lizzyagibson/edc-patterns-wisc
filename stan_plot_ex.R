
  example("read_stan_csv")
  stan_plot(fit)
  stan_trace(fit)
  
  library(gridExtra)
  fit <- stan_demo("eight_schools")
  
  stan_plot(fit)
  stan_plot(fit, point_est = "mean", show_density = TRUE, fill_color = "maroon")
  
  
  # histograms
  stan_hist(fit)
  # suppress ggplot2 messages about default bindwidth
  quietgg(stan_hist(fit))
  quietgg(h <- stan_hist(fit, pars = "theta", binwidth = 5))
  
  # juxtapose histograms of tau and unconstrained tau 
  tau <- stan_hist(fit, pars = "tau")
  tau_unc <- stan_hist(fit, pars = "tau", unconstrain = TRUE) +
    xlab("tau unconstrained")
  grid.arrange(tau, tau_unc)
  
  # kernel density estimates
  stan_dens(fit)
  (dens <- stan_dens(fit, fill = "skyblue", ))
  dens <- dens + ggtitle("Kernel Density Estimates\n") + xlab("")
  dens
  
  (dens_sep <- stan_dens(fit, separate_chains = TRUE, alpha = 0.3))
  dens_sep + scale_fill_manual(values = c("red", "blue", "green", "black"))
  (dens_sep_stack <- stan_dens(fit, pars = "theta", alpha = 0.5,
                               separate_chains = TRUE, position = "stack"))
  
  # traceplot
  trace <- stan_trace(fit)
  trace +
    scale_color_manual(values = c("red", "blue", "green", "black"))
  trace +
    scale_color_brewer(type = "div") +
    theme(legend.position = "none")
  
  facet_style <- theme(strip.background = ggplot2::element_rect(fill = "white"),
                       strip.text = ggplot2::element_text(size = 13, color = "black"))
  (trace <- trace + facet_style)
  
  # scatterplot
  (mu_vs_tau <- stan_scat(fit, pars = c("mu", "tau"), color = "blue", size = 4))
  mu_vs_tau +
    ggplot2::coord_flip() +
    theme(panel.background = ggplot2::element_rect(fill = "black"))
