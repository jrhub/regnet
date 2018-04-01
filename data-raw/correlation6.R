
SimCorStru.AR <- function(p, r, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    a = (i-1)%/%n*n +1
    for(j in a:min(a+n-1,p)){
      sig[i,j] = r^abs(i-j)
    }
  }
  sig
}

SimCorStru.AR.C <- function(p, r){
  sig = matrix(0,p,p)
  for (i in 1: p){
    for(j in 1: p){
      sig[i,j] = r^abs(i-j)
    }
  }
  sig
}

SimCorStru.Clinic <- function(p, r){
  sig = matrix(r,p,p)
  # for (i in 1: p){
    # for(j in 1: p){
      # sig[i,j] = r
    # }
  # }
  diag(sig) = 1
  sig
}

SimCorStru.BLK <- function(p, r, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    a = (i-1)%/%n*n +1
    for(j in a:min(a+n-1,p)){
      sig[i,j] = r
    }
  }
  diag(sig) = 1
  sig
}

x.SimCorStru.BA.1 <- function(p, r, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    for(j in max(i-1,1):min(i+1,p)){
      sig[i,j] = r
    }
  }
  diag(sig) = 1
  sig
}

SimCorStru.BA.1 <- function(p, r, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    a = (i-1)%/%n*n +1
    for(j in a:min(a+n-1,p)){
      if(abs(j-i) <= 1) sig[i,j] = r^abs(i-j)
    }
  }
  sig
}

SimCorStru.BA.2 <- function(p, r, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    a = (i-1)%/%n*n +1
    for(j in a:min(a+n-1,p)){
      if(abs(j-i) <= 2) sig[i,j] = r^abs(i-j)
    }
  }
  sig
}

x.SimCorStru.BA.2 <- function(p, r, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    for(j in max(i-2,1):min(i+2,p)){
      sig[i,j] = r^abs(i-j)
    }
  }
  #diag(sig) = 1
  sig
}

x.SimCorStru.BA.2.2 <- function(p, r1, r2, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    for(j in max(i-2,1):min(i+2,p)){
      if(abs(i-j) == 1){
        sig[i,j] = r1
      }
      else if(abs(i-j) == 2){
        sig[i,j] = r2
      }
    }
  }
  diag(sig) = 1
  sig
}

x.SimCorStru.BA.3 <- function(p, r, n){
  sig = matrix(0,p,p)
  for (i in 1: p){
    for(j in max(i-3,1):min(i+3,p)){
      sig[i,j] = r^abs(i-j)
    }
  }
  #diag(sig) = 1
  sig
}

#round(sig, digits = 5)
#diag(sig)
#sig[1:5, 1:5]
#sig[1:5, 6:10]

# SimCorStru.BA.1(11, 0.9, 5)
# SimCorStru.BA.3(11, 0.6, 5)
# SimCorStru.BA.2.2(11, 0.4, 0.3, 5)
# SimCorStru.BLK(11, 0.9, 5)
# SimCorStru.AR(11, 0.9, 5)
