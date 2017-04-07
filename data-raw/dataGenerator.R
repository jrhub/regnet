rm(list=ls(all=TRUE))
library("Rlab")
library("MASS")

source("../LASSO/correlation6.R")

n = 200; p = 300;
n.sub = 5
corr = 0.9

sig = SimCorStru.AR(p, corr, n.sub)
#sig = SimCorStru.BA.2.2(p, 0.6, 0.33, n.sub)
#sig = SimCorStru.BA.1(p, corr, n.sub)
#sig = SimCorStru.BLK(p, corr, n.sub)

r = 4.5      # MCP
lam2 = 1
Pr = 1     # prob for bernulli

values.2 = c(0.1, 1, 10)
#-----------------------------------------------------------------------------
simulation <- function(){
  X0 = mvrnorm(n,rep(0,p),sig)
  #X0= GenSNP(X0)
  x = scale(X0, scale = apply(X0, 2, function(t) sd(t)*sqrt((n-1)/n)))        #normalize the data
  x = cbind(rep(1,n), x)
  y1 = x %*% b.true
  y2=exp(y1)/(1+exp(y1))

  y=matrix(0,n,1);
  for(i in 1:n) {
    y[i,1]=rbinom(1,1,y2[i]) ;
  }

  list(x,y)
}

GenSNP <- function(x){
  for(i in 1:p){
    v = x[,i]
    Q = quantile(v)
    x[,i] = 1
    x[v<Q[2],i] = 0
    x[v>Q[4],i] = 2
  }
  x
}

balance <- function(){
  e = n/50
  ID = 0
  while((ID <(n/2-e))||(ID >(n/2+e))){ # control the ratio of case/control
    data = simulation();
    ID = sum(data[[2]]);                # check the ratio of case/control;
  }
  data
}

sparsity <- function(n.sub, Pr){
  k = rep(0,p);  g = p/n.sub;
  index = sample(seq(1, g, 1), round(g*0.1/Pr), replace = FALSE)
  begining = (index-1)*n.sub + 1
  for(i in begining){
    k[i:(i+n.sub-1)] = rbern(n.sub, Pr)
  }
  k
}
#sum(sparsity(5))

TruePos <- function(b){
  tp = 0; fp = 0;
  for(i in which(b != 0)){
    if(i %in% index) tp = tp + 1
    else fp = fp +1
  }
  list(tp, fp)
}

##############################################################################

b.true = c(0.3, runif(p, 0.25,0.75) * sparsity(n.sub, Pr))   # 5 per cluster
index = which(b.true != 0)
T = length(index); F = p + 1 - T

data = balance()
fae= data[[1]]
tab= data[[2]]

test = list(X=fae[,-1], Y=drop(tab), beta=b.true)

ptm <- proc.time()
#result = CV.NetLogistic(fae[,-1], tab, r=4.5, folds=5)
#result = CV.McpLogistic(fae[,-1], tab)
result = CV.ElasLogistic(fae[,-1], tab, alpha=1)
print(proc.time() - ptm)

result$lambda
result$mcr
result$MCR


#b = NetLogistic(fae[,-1], tab, result$lambda[1,2], result$lambda[1,1])
#b = McpLogistic(fae[,-1], tab, result$lambda[1])
b = ElasLogistic(fae[,-1], tab, result$lambda[1], alpha=1)
#TruePos(b)
TruePositive(b,b.true)
