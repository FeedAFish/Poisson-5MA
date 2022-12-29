library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(curry)

theme_set(theme_ipsum(base_family = "") + theme(
  axis.title.x = element_text(hjust = 0.5),
  axis.title.y = element_text(hjust = 0.5), plot.margin = margin(
    t = 0.5,
    r = 2, b = 0.5, l = 2, "cm"
  )
))

simulate_hpp_t <- function(lambda, t) {
  return(sort(runif(rpois(1, lambda * t), 0, t)))
}

simulate_claim_amount_process <- function(lambda, t, claim_amount_fn) {
  claim_t <- simulate_hpp_t(lambda, t)
  claim_amount_dist <- claim_amount_fn(length(claim_t))

  x <- claim_amount_dist$X
  e_x <- claim_amount_dist$E_X

  return(list(T = claim_t, S = cumsum(x), E_X = e_x))
}

simulate_risk_process <- function(u, c, t, lambda, claim_amount_fn,
                                  with_plot = TRUE, color = "black") {
  cap <- simulate_claim_amount_process(lambda, t, claim_amount_fn)
  claim_t <- c(0, cap$T)
  s <- c(0, cap$S)
  y_tilde <- s - claim_t * c
  risk <- u - y_tilde
  plot <- NULL

  if (with_plot) {
    claim_t_end <- c(cap$T, t)
    peaks <- u + c * claim_t_end - s
    plot <- ggplot() +
      geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        color = color,
        data = data.frame(
          x = claim_t, xend = claim_t_end, y = risk, yend = peaks
        )
      ) +
      geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        color = color,
        linetype = "dashed",
        data = data.frame(
          x = cap$T, xend = cap$T, y = head(peaks, -1), yend = tail(risk, -1)
        )
      )
  }

  return(list(
    u = u, c = c, t = t, S = s,
    lambda = lambda, p = lambda * cap$E_X / c,
    X = data.frame(T = claim_t, R = risk),
    Y_tilde = data.frame(T = claim_t, Y = y_tilde),
    ruin_time = suppressWarnings(min(which(risk < 0))),
    plot = plot
  ))
}

simulate_survival_prob <- function(u, c, t, lambda, claim_amount_fn, iter) {
  count <- 0
  for (i in 1:iter) {
    rp <- simulate_risk_process(u, c, t, lambda, claim_amount_fn, FALSE)
    if (rp$ruin_time == Inf) {
      count <- count + 1
    }
  }
  return(count / iter)
}

plot_surv_c <- function(u, c, t, lambda, fn, iter) {
  res <- NULL
  for (i in seq_along(c)) {
    res <- c(res, simulate_survival_prob(u, c[i], t, lambda, fn, iter))
  }
  return(ggplot() +
    labs(y = "Survival probability") +
    geom_point(aes(x = c, y = res)))
}

plot_surv_l <- function(u, c, t, lambda, fn, iter) {
  res <- NULL
  for (i in seq_along(lambda)) {
    res <- c(res, simulate_survival_prob(u, c, t, lambda[i], fn, iter))
  }
  return(ggplot() +
    labs(y = "Survival probability") +
    geom_point(aes(x = lambda, y = res)))
}

plot_surv_prob <- function(mode, u, c, t, lambda, fn, iter) {
  res <- NULL
  if (mode == 1) {
    res <- plot_surv_c(u, c, t, lambda, fn, iter)
  } else if (mode == 2) {
    res <- plot_surv_l(u, c, t, lambda, fn, iter)
  }
  return(res)
}

plot_c <- function(u, c, t, lambda, fn, iter = 1) {
  res <- NULL
  for (i in seq_along(c)) {
    test <- ggplot() +
      labs(
        title = paste(
          "Wallet process with c = ",
          toString(c[i])
        )
      )

    for (j in 1:iter) {
      test <- test + geom_line(
        data = simulate_risk_process(
          u, c[i], t, lambda, fn, FALSE
        )$X,
        aes(x = T, y = R),
        color = i
      )
    }

    res[[length(res) + 1]] <- test
  }
  return(res)
}

plot_l <- function(u, c, t, lambda, fn, iter = 1) {
  for (i in seq_along(lambda)) {
    test <- ggplot() +
      labs(
        title = paste(
          "Wallet process with lambda = ",
          toString(lambda[i])
        )
      )
    for (j in 1:iter) {
      test <- test + geom_line(
        data = simulate_risk_process(
          u, c, t, lambda[i], fn, FALSE
        )$X,
        aes(x = T, y = R),
        color = i
      )
    }

    res[[length(res) + 1]] <- test
  }
  return(res)
}

simu_risk_plot_mode <- function(mode, u, c, t, lambda, fn, iter = 1) {
  res <- list()
  if (mode == 1) {
    res <- plot_c(u, c, t, lambda, fn, iter)
  } else if (mode == 2) {
    res <- plot_l(u, c, t, lambda, fn, iter)
  }
  return(res)
}

plot_fn <- function(u, c, t, lambda, fn, iter = 1, color = 1, prefix = NULL) {
  test <- ggplot() +
    labs(
      title = paste(
        "Wallet process with mu = ",
        toString(prefix)
      )
    )
  for (j in 1:iter) {
    test <- test + geom_line(
      data = simulate_risk_process(
        u, c, t, lambda, fn, FALSE
      )$X,
      aes(x = T, y = R),
      color = color
    )
  }
  return(test)
}
