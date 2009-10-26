mtsknn.neq = function(x,y,k,clevel=0.05)
{
# x, y are matrices or data frame with each row containing the coordinates of data points
# and with the last colume being the flags.

# preparation:

  x <- as.matrix(x)
  y <- as.matrix(y)
  if(ncol(x)!=ncol(y))return("The dimensions of two samples must match!!!")

  d <- ncol(x)
  n1 <- nrow(x)
  n2 <- nrow(y)
  
  if(n1 > n2){
     temp <- x
     x <- y
     y <- temp
     n1 <- nrow(x)
     n2 <- nrow(y)
  }

  b <- log(1/(1-clevel))

  q <- as.integer(n2/n1)
  p <- as.integer(log((q+1),2))
  m <- as.integer(n2/q)
  r <- n2-m*q
  starts <- seq(1,(q*m+1),by=m)

  if(r>0)starts <- c(starts[1:(q-r+1)],(starts[(q-r+2):(q+1)]+seq(1,r,by=1)))


      adjust.cl <- b/p

      y.permuted <- as.matrix(y[sample(c(1:n2)),])

      K <- 0
      y.temp <- NULL
      for(i in 1:p){
         x.sub <- rbind(x,y.temp)
         n1.sub <- nrow(x.sub)
         x.sub <- cbind(x,rep(1,n1))
         
         y.sub <- as.matrix(y.permuted[starts[2^(i-1)]:(starts[2^(i)]-1),])
         y.temp <- rbind(y.temp, y.sub)
         n2.sub <- nrow(y.sub)
         y.sub <- cbind(y.sub,rep(2,n2.sub))
         
         Set <- rbind(x.sub,y.sub)
         tSet <- t(Set)
         n.sub <- n1.sub+n2.sub
         output <- rep(0,n.sub)
         C.out <- .C("knn", as.double(tSet), as.integer(n.sub), as.integer(d), as.integer(k), as.integer(output))
         counts <- C.out[[5]]
         Tk <- sum(counts)/(n.sub*k)
         # Z scores and P values
         #V <- lam1*lam2+4*lam1^2*lam2^2
         V <- (n1.sub-1)*(n2.sub-1)/(n.sub-1)^{2}+4*((n1.sub-1)*(n1.sub-2)/((n.sub-1)*(n.sub-2)))*((n2.sub-1)*(n2.sub-2)/((n.sub-1)*(n.sub-2)))
         #Z <- (n*k)^(1/2)*(Tk-lam1^2-lam2^2)/sqrt(V)
         Z <- (n.sub*k)^(1/2)*(Tk-(n1.sub-1)*(n1.sub-2)/((n.sub-1)*(n.sub-2))-(n2.sub-1)*(n2.sub-2)/((n.sub-1)*(n.sub-2)))/sqrt(V)
         P <- pnorm(Z, lower.tail=F)
         if(P > adjust.cl) break
         K <- K+1
      }
      if(K==p) return("Null hypothesis is accepted")
      return("Null hypothesis is rejected")
  }








