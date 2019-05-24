#' Frequency graph
#'
#' Drawing frequency graph (density plot) and analyzing the out put of
#' \code{\link{lalgp}} or \code{\link{manual_lalgp}}.
#'
#' @export
#' @param data list data of the result of functions; \code{\link{lalgp}} or \code{\link{manual_lalgp}}.
#'@examples
#' d <- Osayama
#' e <- lalgp(d)
#' lalgp_graph(e)

lalgp <- function(data){
  rho_max <- 100
  sigma_max <- 100
  m_waic <- 50
  m_graph <- 200
  sigma0 <- 1
  rho0 <- 10
  Samples <- 10000
  delta <-  1/25
  #LA_LGP() can compute WAIC and p(x) with fixed rho and sigma
  #------------------------------------------------------------------------------
  LA_LGP <- function(
    d,
    m,
    rho = 1,
    sigma = 1,
    f_0 = numeric(m),
    Loop = 5){

    n = length(d)
    max_d <- max(d)
    min_d <- min(d)
    #the definiton of various

    #     |<---------------  1 --------------->|
    #     delta                            delta            normalized_age
    #     |   |<------ 1 - 2*delta ------->|   |
    #------------------------------------------------------> age
    #     |   | <= min_d          max_d => |   |
    #     Delta                            Delta            real_age
    #     | <= Min_d                  Max_d => |

    Delta <- (max_d - min_d) / (1 - 2 * delta) * delta
    Min_d <- min_d - Delta
    Max_d <- max_d + Delta

    #nd is a normalized value of d
    nd <-  (d - Min_d) / (Max_d - Min_d)
    #X
    X <- numeric(m)
    for(i in 1:m){X[i] <- 1/(2*m)+(i-1)/m}
    X <- X*100
    #y
    y <- numeric(m)
    for(j in 1:m){
      for(i in 1:n){
        if((j-1) * (1/m) <= nd[i] & nd[i] <= j * (1/m)){
          y[j] = y[j]+1
        }
      }
    }


    #K
    K <- matrix(nrow = m,ncol = m,0)
    for(i in 1:m){
      for(j in 1:m){
        K[i,j] <- sigma^2 * exp(-(X[i] - X[j])^2 / (rho)^2)
        if(i==j) K[i,j] = K[i,j] + 0.01
      }
    }

    f <- matrix(f_0)

    #repeat
    #See Riihimaki and Vehtari (2014)
    #newton method for find estimate of ^f
    for(i in 0:Loop){
      #u
      u <- exp(f)/sum(exp(f))
      #W
      W <-n*( diag(u[,1])-u[,1] %*% t(u[,1]))
      #R
      R <- sqrt(n)*(diag(sqrt(u[,1])) - (u[,1] %*% t(u[,1])) %*% diag(u[,1]^(-1/2)))
      #L = (K-1 + W)-1
      I <- diag(m)
      L <- K %*% (I - R %*% solve(I + t(R) %*% K %*% R) %*% t(R) %*% K)
      #Wf
      Wf <- W %*% f
      #S <- nabla log(y|p)
      S <- (y - n * exp(f)/ sum(exp(f)))
      #fnew
      fnew <- L %*% (Wf + S)
      f <- fnew
    }

    #Random sampling from N(^f,(K-1 + W)-1)
    ignore_error_rmvnor <- function(Samples,fnew,L){
      f_result <- NULL
      f_result <- try(mvtnorm::rmvnorm(Samples,fnew,L), silent = T)
      if (class(f_result) == "try-error") {
        f_result <- matrix(ncol=m,nrow=2,0)
      }
      f_result
    }
    f_result <- ignore_error_rmvnor(Samples,fnew,L)



    #computation of WAIC
    p_result <- exp(f_result) / rowSums(exp(f_result))

    q1 <- colMeans(p_result)
    q2 <- log(p_result)
    q3 <- apply(q2, 2, var)

    lppd <- sum(y * log(q1))
    pwaic <- sum(y * q3)
    WAIC <- -2 * (lppd - pwaic)


    #pick up percentile


    answer <- list(WAIC = WAIC, p = p_result, data = d, rho = rho, sigma = sigma, f_new = colMeans(f_result))
    return(answer)
  }
  #-------------------------------------------------------------------------------


  # Find an initial f_0 that Newton's method converges
  ## Using sigma = 1 and rho = 10, Newton's method is less likely to be suboptimal
  #-------------------------------------------------------------------------------
  f_new <- LA_LGP(d = data, sigma = sigma0, rho = rho0, m = m_waic)$f_new
  #-------------------------------------------------------------------------------


  # Select appropriate rho and sigma using WAIC
  #-------------------------------------------------------------------------------
  #Example of WAIC search
  # n.> and n.< are Comparison calculation
  #ex) 4.> compare WAIC[1,4] > WAIC[1,5]
  #1., 2., ... are order of calculation

  #WAIC[rho,sigma]
  #         | 1 | 2 | 3 | 4 | 5 | 6 |...| max_rho |
  #1        | 1.> 2.> 3.> 4.<
  #2        |             5.> 7.<
  #3        |                 8.<
  #4        |                 9.>10.<
  #5        |
  #...      |
  #max_sigma|

  #above case, if 10. is larger than 8., stop computation. and
  #rho = 5 and sigma = 10 are adopted.

  WAIC <- matrix(ncol = rho_max, nrow = sigma_max, 10^6)
  cat("computing waic")
  r <- 1
  for(i in 1:(sigma_max - 1)){
    WAIC[i, r] <- LA_LGP(d = data, sigma = i, rho = r, f_0 = f_new, m = m_waic)$WAIC

    for(j in r:(rho_max - 1)){
      WAIC[i, j+1] <- LA_LGP(d = data, sigma = i, rho = (j+1), f_0 = f_new, m = m_waic)$WAIC
      cat(".")
      if(WAIC[i,j] < WAIC[i,(j+1)]){
        r <- j
        break
      }
    }
    if(i > 1){
      if(min(WAIC[(i-1),]) < min(WAIC[i,]))
        break
    }
  }
  cat("\n")
  cat("computing p(x)")
  for(i in 1:sigma_max){
    for(j in 1:rho_max){
      if(min(WAIC) == WAIC[i,j]){
        rho = j
        sigma = i
        break
      }
    }
  }
  #-------------------------------------------------------------------------------


  #Recalculate p (x) using appropriate rho and sigma
  #-------------------------------------------------------------------------------
  f_new <- LA_LGP(d = data, sigma = sigma0, rho = rho0, m = m_graph)$f_new
  ans <-  LA_LGP(d = data, sigma = sigma, rho = rho, f_0 = f_new, m = m_graph)

  p <- ans$p
  h <- 0.0001

  ##----
  n = length(data)
  max_d <- max(data)
  min_d <- min(data)
  Delta <- (max_d - min_d) / (1 - 2 * delta) * delta
  Min_d <- min_d - Delta
  Max_d <- max_d + Delta

  nd <-  (data - Min_d) / (Max_d - Min_d)
  X <- numeric(m_graph)
  for(i in 1:m_graph){X[i] <- 1/(2*m_graph)+(i-1)/m_graph}
  age_X <-  Min_d + (Max_d-Min_d) * X
  delta_X <- (max(age_X) - min(age_X))/m_graph
  f_mean <- ans$f_new
  p_mean <- exp(f_mean) / sum(exp(f_mean))
  ##---

  for(j in 1:10){

    qua <- matrix(0,ncol=2,nrow=m_graph)
    for(i in 1:m_graph){
      qua[i,] <- as.numeric(quantile(p[,i],c(h*j,1-h*j)))
    }

    count <- 0
    for(i in 1:Samples){
      if(sum(p[i,] < qua[,1]) > 0){
        count <- count + 1
      }
    }
    cat(".")
    if(count > 500){break}

  }
  ####
  for(j in 1:10){

    qua <- matrix(0,ncol=2,nrow=m_graph)
    for(i in 1:m_graph){
      qua[i,] <- as.numeric(quantile(p[,i],c(h*j,1-h*j)))
    }

    count <- 0
    for(i in 1:Samples){
      if(sum(p[i,] < qua[,2]) > 0){
        count <- count + 1
      }
    }
    cat(".")
    if(count > 500){break}

  }

  p_0025 <- qua[,1]
  p_0975 <- qua[,2]

  qua2 <- matrix(0,ncol=2,nrow=m_graph)
  for(i in 1:m_graph){
    qua2[i,] <- as.numeric(quantile(p[,i],c(0.025,0.975)))
  }

  p2_0025 <- qua2[,1]
  p2_0975 <- qua2[,2]

  p <- data.frame(age_X = age_X,
                  p_mean = p_mean/delta_X,
                  p_0025 = p_0025/delta_X,
                  p_0975 = p_0975/delta_X,
                  p2_0025 = p2_0025/delta_X,
                  p2_0975 = p2_0975/delta_X)


  anslist <- list(WAIC = ans$WAIC, p = p, data = data, rho = rho, sigma = sigma)
  return(anslist)
}
