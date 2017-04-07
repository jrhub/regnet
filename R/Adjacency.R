
Adjacency = function(x, alpha=5)
{
  n = nrow(x)
  p = ncol(x)
  r0 = stats::cor(x)
  r = r0; r[which(r==1)] = 1 - 0.01
  z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
  c0 = mean(sqrt(n-3)*z) + 2*stats::sd(sqrt(n-3)*z)
  cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
  r = r0
  A = (r)^alpha*(abs(r)>cutoff)
  diag(A) = 0
  A
}
