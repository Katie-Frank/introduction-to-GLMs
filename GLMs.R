library(ggplot2); theme_set(theme_bw())

# Exponential example
la1 <- expression(paste(lambda, " = 1"))
la2 <- expression(paste(lambda, " = 2"))

ggplot(data.frame(time = c(0, 5)), aes(time)) +
  stat_function(fun = dexp, args = list(rate = 1), show.legend = TRUE) +
  stat_function(fun = dexp, args = list(rate = 2), color = "red") +
  annotate("text", x = 2, y = .5, label = la1, size = 4) +
  annotate("text", x = .7, y = 1.25, label = la2, size = 4, color = "red") +
  xlab("Time to failure") +
  ylab("PDF")

# GLM example with survival data
dat <- read.table("leukemia.txt", header = TRUE, sep = "\t")

ggplot(dat, aes(log_wbc, surv_weeks)) +
  geom_point() +
  xlab(expression(log[10]("initial white blood cell count"))) +
  ylab("Survival time (weeks)")

# IRLS procedure
irls_exp <- function(x, y, be1, be2, max_iter = 100, tol = 1e-5){
  X <- matrix(c(rep(1, length(x)), x), ncol = 2)
  be <- matrix(c(be1, be2), ncol = 1)
  be_prev <- be
  
  for(iter in 1:max_iter){
    z <- X %*% be + y / (exp(X %*% be)) - 1
    cov_mat <- solve(t(X) %*% X) 
    be <- cov_mat %*% t(X) %*% z
    
    if (abs(be[1] - be_prev[1]) < tol && abs(be[2] - be_prev[2]) < tol) {
      break
    } else {
      be_prev <- be
    }
  }
  list(be = be, cov_mat = cov_mat, num_iter = iter)
}

irls_exp(dat$log_wbc, dat$surv_weeks, 10, -1)

# glm fit
mod <- glm(surv_weeks ~ log_wbc, family = Gamma(link = "log"), data = dat)
summary(mod) # standard error estimates not right
summary(mod, dispersion = 1) # correction for standard errors

# glm fitted curve 
ggplot(dat, aes(log_wbc, surv_weeks)) +
  geom_point() + 
  geom_smooth(method = "glm", formula = y ~ x, se = FALSE, 
    method.args = list(family = Gamma(link = "log"))) + 
  xlab(expression(log[10]("initial white blood cell count"))) +
  ylab("Survival time (weeks)")

# standardized residuals
r <- (dat$surv_weeks - mod$fitted.values) / mod$fitted.values

# plot of standardized residuals vs. fitted values
ggplot(data.frame(fitted_values = mod$fitted.values, residual = r), 
  aes(fitted_values, residual)) +
  geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-2.75, 2.75)) + 
  xlab("Fitted value") +
  ylab("Residual")
