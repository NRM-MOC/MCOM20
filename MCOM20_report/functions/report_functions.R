mocis_heatmap <- function(data, trunc = 15, grp = "all"){
  if (grp != "all"){
    data <- filter(data, group == grp)
  }
  data %>% 
    ggplot(aes(x = var_name, y = station)) + 
    geom_tile(aes(fill = ifelse(abs(slope) > trunc, trunc*sign(slope), slope))) + 
    scale_fill_gradient2(limits = c(-trunc, trunc), low = muted("blue"), mid = "white", high = muted("red"), na.value = 'white') + 
    geom_text(aes(label = stars,  color = (abs(slope) > 10)), size = 3, show.legend = FALSE) + scale_color_manual(values = c("black", "white")) +  xlab("") + ylab("")+
    facet_col(gen_name~., scales = "free_y", space = "free") + 
    theme(legend.position = "top", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.title = element_blank()) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggtitle(ifelse(grp == "all", "", paste0(grp, ", yearly percent change (2010-2019)")))
}

common_intercept <- function(fit.lm){
  theta <- coef(fit.lm)
  d <- length(theta)
  D <- cbind(1, diag(d - 1))
  D[1, 2] <- 0
  intercepts <- t(D %*% theta)
  Sigma <- D %*% vcov(fit.lm) %*% t(D)
  w <- diag(Sigma)
  sum(intercepts/w)/sum(1/w)
}

print_ci <- function(estimate, lower, upper, dec = 2) {
  ifelse(is.numeric(estimate), paste0(round(estimate, dec), 
                                      " (", round(lower, dec), ", ", round(upper, 
                                                                           dec), ")"), "")
}
slope2rate <- function(x) {100*(exp(x) -1)}

basin_plot <- function(data,  grp = "all", cols = 2, unit = ""){
  if (grp != "all"){
    data <- filter(data, group == grp)
  }
  ggplot(data, aes(x = YEAR, y = CONC, color = basin, fill = basin)) +
    geom_ribbon(data = filter(data, p_vs_full > .05), aes(ymin = LOWER, ymax = UPPER), alpha = .15, colour = NA)+
    geom_line(data = filter(data, p_vs_full > .05)) +
    geom_line(data = filter(data, p_vs_full <= .05), linetype = "dashed") +
    facet_wrap(~var_name, scales = "free_y", ncol = cols) +
    scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1))) +
    scale_x_continuous(breaks=seq(2010, 2019, 2)) + #ALS Nov 2020
    # theme_bw() +
    theme(legend.position = "top", legend.title = element_blank(),           
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) + 
    ylab(unit) + xlab("") +
    ggtitle(ifelse(grp == "all", "", paste0(grp, ", fitted trends (2010-2019)")))
}

