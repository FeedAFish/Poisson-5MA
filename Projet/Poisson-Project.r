library(ggplot2)
library(RColorBrewer)
################################################################################
#' \code{HomogeneousPoissonProcess} returns the value of the homogeneous Poisson
#' Process N_t, which corresponds to number of claims which arrived by time t.

#' @param lambda_Nt : parameter of the Poisson Process N_t.
#' @param t : time at which we count the number of claims that already happened.

HomogeneousPoissonProcess <- function(lambda_Nt,t){
  n <- 0 # initialization value of the homogeneous Poisson Process
  Wn <- 0 # corresponds to the interarrival time
  while(Wn <= t){
    Wn <- Wn + rexp(1,lambda_Nt)
    n <- n + 1
  }
  n
}
################################################################################
#' \code{risk process} returns the value of the risk process at time t, called 
#' R_t and defined according to the Cramer-Lundberg model.

#' @param u : the initial capital of the insurance.
#' @param c : premium rate.
#' @param t : time at which we return the value of the risk process.
#' @param lambda_Nt : parameter of the Poisson Process N_t.
#' @param lambda_Xn : parameter of the X_n variable (i.e. the amount of the n-th
#' claim, which follows an exponential law of parameter lambda_Xn)

risk_process <- function(u,c,t,lambda_Nt,lambda_Xn){
  if (t == 0){
    Rt = u # initialization of the risk process
    Rt
  }else{
    Xt <- 0
    Nt <- HomogeneousPoissonProcess(lambda_Nt,t)
    for(j in 1:Nt){
      Xt <- Xt + rexp(1,lambda_Xn)
    }
    Rt <- u + c*t - Xt # Cramer-Lundberg model
    Rt
  }
}
################################################################################
#' \code{simulation} returns the Risk Process R_t evaluated on each 
#' integer in the interval [0,t].

#' @param u : the initial capital of the insurance.
#' @param c : premium rate.
#' @param t : time at which we return the value of the risk process.
#' @param lambda_Nt : parameter of the Poisson Process N_t.
#' @param lambda_Xn : parameter of the X_n variable (i.e. the amount of the n-th
#' claim, which follows an exponential law of parameter lambda_Xn).
#' @param nb_simu : number of simulations of the risk process R_t.


simulation <- function(u,c,t,lambda_Nt,lambda_Xn,nb_simu){
  val2 <- NULL
  for(i in 1:nb_simu){
    val1 <- NULL
    for(k in 0:t){
      val1 <- c(val1,risk_process(u,c,k,lambda_Nt,lambda_Xn))
    }
    val2 <- c(val2,val1)
  }
  simulations <- data.frame(s = rep(c(1:nb_simu), each = t+1),
                             t = rep(c(0:t), nb_simu),
                             val = val2)
  simulations
}


# First simulation #############################################################
graph <- simulation(0,1,50,0.9,1,3) # simulation of 5 Risk Processes

ggplot(graph, aes(x=t, y=val, group=simu)) +
  geom_line(aes(color=as.character(simu))) + theme(legend.position="bottom") + 
  scale_x_continuous(name="Time") +
  scale_y_continuous(name="Rt: wallet at time t")

ggplot(graph, aes(x=t, y=val, group=s)) +
  geom_line(aes(color=as.character(s))) + theme(legend.position="bottom") +
  geom_abline(intercept = 0, slope = (1 - (0.9*(1/1))), color="red",size=1) +
  scale_x_continuous(name="Time") +
  scale_y_continuous(name="Rt: wallet at time t") 
# We plot in red the expectation of Rt : E[Rt] = u-(c-lambda*mu)t,
# where mu = E[X_1] = ... = E[X_n] = 1/lambda_Xn

################################################################################
# The previous simulations are generations of the R_t process, where the wallet
# is obtained for each s \in [0,t], for s an integer. i.e., the process restarts
# at each instant s+1; now we are going to write a code which substract the 
# amount of claims and takes into account what we had before. 
# Therefore, the wallet is a cumulative function of previous income and expenses.

#' \code{claims} returns the amount of the whole claims wallet.

#' @param t : time at which we return the value of the claims wallet.
#' @param lambda_Nt : parameter of the Poisson Process N_t.
#' @param lambda_Xn : parameter of the X_n variable (i.e. the amount of the n-th
#' claim, which follows an exponential law of parameter lambda_Xn).

claims <- function(t,lambda_Nt,lambda_Xn){
  if (t == 0){
    Xt = 0 # initialization
    Xt
  }else{
    Xt <- 0
    Nt <- HomogeneousPoissonProcess(lambda_Nt,t)
    for(j in 1:Nt){
      Xt <- Xt + rexp(1,lambda_Xn)
    }
    Xt
  }
}

################################################################################
#' \code{cumulative_simulation} returns the value of the cumulative wallet.

#' @param u : the initial capital of the insurance.
#' @param c : premium rate.
#' @param t : time at which we return the value of the wallet.
#' @param lambda_Nt : parameter of the Poisson Process N_t.
#' @param lambda_Xn : parameter of the X_n variable (i.e. the amount of the n-th
#' claim, which follows an exponential law of parameter lambda_Xn).
#' @param nb_simu : number of simulations of the risk process R_t.

cumulative_simulation <- function(u,c,t,lambda_Nt,lambda_Xn,nb_simu){
  val2 <- NULL
  for(i in 1:nb_simu){
    val1 <- NULL
    val1 <- u
    for(k in 1:t){
      val1 <- c(val1,val1[k] + c - claims(1,lambda_Nt,lambda_Xn))
    }
    val2 <- c(val2,val1)
  }
  simulations <- data.frame(sim = rep(c(1:nb_simu), each = t+1),
                             t = rep(c(0:t), nb_simu),
                             val = val2)
  simulations
}

graph <- cumulative_simulation(0,1,50,1,2,2)
ggplot(graph, aes(x=t, y=val, group=sim)) +
  geom_line(aes(color=sim)) + theme(legend.position="none") +
  scale_x_continuous(name="Time") +
  scale_y_continuous(name="Rt: Cumulative wallet at time t")

################################################################################
#' \code{ruin_probability} returns the value of Psi(u), i.e. the ruin 
#' probability of the associated ruin process R_t.

#' @param u : the initial capital of the insurance (u = R_0 > 0)
#' @param lambda : Poisson intensity of R_t.
#' @param alpha : parameter of the exponential distribution of the claim amounts.
#' @param c : premium rate.

psi_ruin_probability <- function(u,lambda,alpha,c){
  psi_ruin_probability <- (lambda/(alpha*c))*exp(-u*(alpha-(lambda/c)))
  psi_ruin_probability
}

psi_ruin_probability(u = 0,lambda = 0.9,alpha = 1,c = 1)
#psi_ruin_probability(u = 8,lambda = 0.9,alpha = 1,c = 1)
#psi_ruin_probability(u = 0,lambda = 1,alpha = 2,c = 1)

################################################################################
#' \code{graphic_ruin} plots the evolution of the ruin probability over time from the
#' values computed by the \code{psi_ruin_probability} function.

#' @param lambda : Poisson intensity of R_t.
#' @param alpha : parameter of the exponential distribution of the claim amounts.
#' @param c : premium rate.
#' @param n_max : maximum amount.

graphic_ruin <- function(lambda,alpha,c,n_max){ 
  amount <- c(NA)
  amount[1] <- 0
  proba <- c(NA)
  proba <- psi_ruin_probability(0,lambda,alpha,c)
  for(i in seq(from=0.1, to=n_max, by=0.05)){
    amount <- c(amount,i)
    proba <- c(proba,psi_ruin_probability(i,lambda,alpha,c))
  }
  df <- data.frame(x=amount, y=proba)
  print(qplot(x, y, data=df, geom="line", 
              xlab = "Initial amount", ylab = "Ruin probability") + 
          xlim(0, n_max) + ylim(0,1))
}

################################################################################
graphic_ruin(lambda = 0.9,alpha = 1,c = 1,n = 35) + 
  geom_point(aes(x=0, y=0.9), colour="blue") +
  geom_point(aes(x=8, y=0.4043961), colour="red")
################################################################################
graphic_ruin(lambda = 1,alpha = 2,c = 1,n = 8) + 
  geom_point(aes(x=0, y=0.5), colour="magenta")
################################################################################
#' We can also define the survival probability (which is the complementary of 
#' the ruin probability, i.e. theta(u) = 1 - psi(u)).
#' 
#' \code{theta_survival_probability} returns the value of theta(u), i.e. the survival
#' probability of the associated ruin process R_t.

#' @param u : the initial capital of the insurance (u = R_0 > 0)
#' @param lambda : Poisson intensity of R_t.
#' @param alpha : parameter of the exponential distribution of the claim amounts.
#' @param c : premium rate.
#' 

theta_survival_probability <- function(u,lambda,alpha,c){
theta_survival_probability <- 1 - psi_ruin_probability(u,lambda,alpha,c)
theta_survival_probability
}

################################################################################
#' \code{graphic_survival} plots the evolution of the survival probability over 
#' time from the values computed by the \code{theta_ruin_probability} function.

#' @param lambda : Poisson intensity of R_t.
#' @param alpha : parameter of the exponential distribution of the claim amounts.
#' @param c : premium rate.
#' @param n_max : maximum amount.

graphic_survival <- function(lambda,alpha,c,n_max){ 
  amount <- c(NA)
  amount[1] <- 0
  proba <- c(NA)
  proba <- theta_survival_probability(0,lambda,alpha,c)
  for(i in seq(from=0.1, to=n_max, by=0.05)){
    amount <- c(amount,i)
    proba <- c(proba,theta_survival_probability(i,lambda,alpha,c))
  }
  df <- data.frame(x=amount, y=proba)
  print(qplot(x, y, data=df, geom="line", 
              xlab = "Initial amount", ylab = "Survival probability") + 
          xlim(0, n_max) + ylim(0,1))
}

################################################################################
graphic_survival(lambda = 0.9,alpha = 1,c = 1,n = 35) + 
  geom_point(aes(x=0, y=1-0.9), colour="blue") +
  geom_point(aes(x=8, y=1-0.4043961), colour="red")
################################################################################
graphic_survival(lambda = 1,alpha = 2,c = 1,n = 8) + 
  geom_point(aes(x=0, y=0.5), colour="magenta")
