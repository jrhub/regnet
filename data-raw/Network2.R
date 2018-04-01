



# adjacency matrix
N.1 = function(X)
{
   n = nrow(X)
   p = ncol(X)
   r = cor(X)
   z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
   c0 = mean(sqrt(n-3)*z) + 2*sd(sqrt(n-3)*z)
   cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
   A = sign(r)*(abs(r) > cutoff)
   diag(A) = 0
   A
}

N.2 = function(X)
{
   n = nrow(X)
   p = ncol(X)
   r = cor(X)
   z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
   c0 = mean(sqrt(n-3)*z) + 2*sd(sqrt(n-3)*z)
   cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
   A = r*(abs(r) > cutoff)
   diag(A) = 0
   A
}

N.3 = function(X)
{
   n = nrow(X)
   p = ncol(X)
   r0 = cor(X)
   r = r0; r[which(r==1)] = 1 - 0.01
   z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
   c0 = mean(sqrt(n-3)*z) + 2*sd(sqrt(n-3)*z)
   cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
   r = r0
   A = (r)^5*(abs(r)>cutoff)
   diag(A) = 0
   A
}

N.4 = function(X)
{
   n = nrow(X)
   p = ncol(X)
   r0 = cor(X)
   r = r0; r[which(r==1)] = 1 - 0.01
   z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
   c0 = mean(sqrt(n-3)*z) + 2*sd(sqrt(n-3)*z)
   cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
   r = r0
   A = (abs(r))^5
   diag(A) = 0
   A
}

N.jr = function(X)
{
  n = nrow(X)
  p = ncol(X)
  r0 = cor(X)
  r = r0; r[which(r==1)] = 1 - 0.01
  z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
  c0 = mean(sqrt(n-3)*z) + 2*sd(sqrt(n-3)*z)
  cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
  r = r0
  A = (r)^5*(r>cutoff)
  diag(A) = 0
  A
}

N.jr2 = function(X)
{
  n = nrow(X)
  p = ncol(X)
  r = cor(X)
  #r = r0; r[which(r==1)] = 1 - 0.01
  #z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
  #c0 = mean(sqrt(n-3)*z) + 2*sd(sqrt(n-3)*z)
  #cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
  #r = r0
  A = (r)^5*(r>0)
  diag(A) = 0
  A
}
