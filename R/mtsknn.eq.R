mtsknn.eq = function(x,y,k,clevel=0.05)
{
# x, y are matrices or data frame with each row containing the coordinates of data points
# and with the last colume being the flags.

# preparation:

  x <- as.matrix(x)
  y <- as.matrix(y)
  if(ncol(x)!=ncol(y))return("The dimensions of two samples must match!!!")

  d <- ncol(x)-1
  n1 <- nrow(x)
  n2 <- nrow(y)
  
  b <- log(1/(1-clevel))
  
  if(n1 > n2){
     temp <- x
     x <- y
     y <- temp
     n1 <- nrow(x)
     n2 <- nrow(y)
  }
  

      q <- as.integer(n2/n1)
      m <- as.integer(n2/q)
      r <- n2-m*q
      starts <- seq(1,(q*m+1),by=m)
      if(r>0)starts <- c(starts[1:(q-r+1)],(starts[(q-r+2):(q+1)]+seq(1,r,by=1)))
      
      adjust.cl <- b/q
      
      y.permuted <- y[sample(c(1:n2)),]
      
      x <- cbind(x,rep(1,n1))
      
      K <- 0
      for(i in 1:q){
         y.sub <- y.permuted[starts[i]:(starts[i+1]-1),]
         n2.sub <- nrow(y.sub)
         y.sub <- cbind(y.sub,rep(2,n2.sub))
         Set <- rbind(x,y.sub)
         tSet <- t(Set)
         n <- n1+n2.sub
         output <- rep(0,n)
         C.out <- .C("knn", as.double(tSet), as.integer(n), as.integer(d), as.integer(k), as.integer(output))
         counts <- C.out[[5]]
         Tk <- sum(counts[(n1+1):n])/(n2.sub*k)
         # Z scores and P values
         #V <- 3/8  ?
         V <- (n1*(n2.sub-1))/((n-1)*(n-2))+((n2.sub-1)*n1*(n1-1))/((n)*(n-2)*(n-3))
         Z <- (n2.sub*k)^(1/2)*(Tk-(n2.sub-1)/(n-1))/sqrt(V)
         P <- pnorm(Z, lower.tail=F)
         if(P < adjust.cl) break
         K <- K+1
      }
      if(K==q) return("Null hypothesis is accepted")
      return("Null hypothesis is rejected")
}

