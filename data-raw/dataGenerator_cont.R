library("glmnet")
library("MASS")
library("Rlab")

source("./data-raw/Network2.R")
source("./data-raw/correlation6.R")

n = 200; p = 300;
n.sub = 5
corr = 0.8

sig = SimCorStru.AR(p, corr, n.sub)
#sig = SimCorStru.BA.1(p, corr, n.sub)
#sig = SimCorStru.BLK(p, corr, n.sub)

#eigen(sig)$values

sim = 1
Pr = 1     # prob for bernulli
#------------------------------------------- function -------------------------------------------------------
initiation <- function(x, y){
  lasso.cv = cv.glmnet(x,y, alpha=1, nfolds=5)       #compute lambda by crossvalidation
  alpha <<- lasso.cv$lambda.min
  lasso.fit = glmnet(x,y,family="gaussian", alpha=1, nlambda=100)
  coef0 <<- as.vector(predict(lasso.fit, s=alpha, type="coefficients"))[-2]    #initialize the coefficients vector
}

TruePos <- function(b){
  index = which(b.true != 0)
  pos = which(b != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
}

#--------------------------------------- simulate data set ------------------------
simulation <- function(){
  X0 = mvrnorm(n,rep(0,p),sig)
  x = scale(X0, scale = apply(X0, 2, function(t) sd(t)*sqrt((n-1)/n)))        #normalize the data
  y = x %*% b.true + rnorm(n,0,1)
  list(x,y)
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
MSE<- function(b){
  sum((b-b.true)^2)
}
#---------------------------------------------------------------------------------------------
TP = rep(0,sim); FP = rep(0,sim);
TP.o = rep(0,sim); FP.o = rep(0,sim);
TP.g = rep(0,sim); FP.g = rep(0,sim);
#TP.e = rep(0,sim); FP.e = rep(0,sim);
alpha = 0; coef0 = rep(0, p); #coef.e = 0
e = n/50##**************$^%$#^#
b0 = rep(0, p)

b.true = runif(p, 0.25,0.75) * sparsity(n.sub, Pr)   # 5 per cluster
#sum(b.true^2)/4    # mean(b.true) / (sum(b.true^2)/4)^0.5
index = which(b.true != 0)
T = length(index); F = p - T

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   Main Loop ++++++++++++++++++++++++++++++++++++++++++++++++

data = simulation();
x = data[[1]]; y = data[[2]]

initiation(x, y)
# b = coef0
b = rep(0, ncol(x))

lambda.n = rev(exp(seq(1,45,1)/5 -7))
# lambda.n = NULL
penalty="lasso"
r=1.5
init = NULL; #"zero"
lamb.2=1

ptm <- proc.time()
res = cv.regnet(x, y, response=c("conti"), penalty=penalty, lamb.1=lambda.n, lamb.2=lamb.2,
                folds=5, r = r, alpha.i = 0.5, initiation = init)
print(proc.time() - ptm)
res$lambda

b.rgn = regnet(x, y, response=c("continuous"), penalty=penalty, res$lambda[1], res$lambda[2], r = r, initiation = init)
# b.rgn = fit$b
# plot(y, fit$residual); abline(a=0, b=0, lty=2)

# b.net = regnet(x, y, response=c("continuous"), penalty="network", 0.3290141 , 1, r = 1.5, initiation = "zero")
counts = TruePos(b.rgn[-1]);
TP.r = counts[[1]];  FP.r = counts[[2]]

counts.g = TruePos(coef0);
TP.g = counts.g[[1]];  FP.g = counts.g[[2]];

#------------------------------------------------------------------------
res$lambda
cat("regnet  TP:   ", TP.r, "   FP: ", FP.r, " MSE: ", MSE(b.rgn[-1]), " intercept", "\n")
# cat("MPC      TP.o: ", TP.o, " FP.o: ", FP.o,"\n")
cat("Lasso    TP.g: ", TP.g, " FP.o: ", FP.g, " MSE: ", MSE(coef0), "\n","\n")

######################################################################################
Net.MCP <- function(lam1, lam2, b, r){
  for( k in 1: p){
    b.old = b[k]
    m = min((k+1),p)
    t = y - x[,-k] %*% b[-k]
    z = t(x[,k]) %*% t /n + lam2 * (t(a[k,m:p]) %*% b[m:p])    #compute Zk
    u = 1 + lam2 * sum(abs(a[k,m:p]))
    if(abs(z) > (r*lam1*u)) b[k] = z / u
    else b[k] = Soft(z, lam1) / (u- 1/r)
    # t = t - x[,k] * (b[k] - b.old)
  }
  b
}


Soft <- function(z, lambda){
  if (z > lambda) z - lambda                 #soft thresholding operate
  else if (z < -lambda) z + lambda
  else 0
}
run.net <- function(lam1,lam2, b, r){
  count = 0
  while(count < 30){
    b.new = Net.MCP(lam1,lam2, b, r)
    dif = sum(abs(b - b.new))/n
    #cat("L1 Diff: ", dif, ",  lamb1: ", lam1, ",  lamb2: ", lam2, ",  sim: ", s,"\n")      # print the test MSE for iterations
    if( dif < 0.001) break                              # stop the iteration if diverge
    else{
      b = b.new; count = count+1
    }
  }
  b.new

}
validation <- function(b){
  mse = sum((y2 - x2 %*% b)^2)/length(y2)
}
r = 1.5      # MCP
lam2 = 1
for(s in 1:sim)
{
  print(s)
  # train = simulation();
  test = simulation();
  # x = train[[1]]; y = train[[2]]
  x2 = test[[1]]; y2 = test[[2]]
  a = N.3(x)

  # initiation()
  # b = coef0
  b = rep(0, ncol(x))

  values = rev(c(seq(alpha/5, alpha, by= alpha/5), seq(alpha*1.5, alpha*3, by= alpha/2)))
  mse.seq = rep(999, length(values))
  #------------- the loop iterates over different lambda values ---------------

  for(i in 1:length(values)){  # Network
    b = run.net(values[i], lam2, b, r)
    mse.seq[i] = validation(b)                                                          # run lasso
    if(mse.seq[i] <= min(mse.seq[-i])) { b.net.R = b; lambda.net = values[i] }      # update lambda.opt & b.opt
  }
  counts = TruePos(b.net.R);
  TP[s] = counts[[1]];  FP[s] = counts[[2]]
  # graph(values, mse.seq, paste("Network ", s, "/",sim))

  counts.g = TruePos(coef0);
  TP.g = counts.g[[1]];  FP.g = counts.g[[2]];
  #------------------------------------------------------------------------
  cat("Network  TP:   ", TP[s], "   FP: ", FP[s], " lamb1: ", lambda.net, " MSE: ", MSE(b.net.R), "\n")
  cat("Lasso    TP.g: ", TP.g, " FP.o: ", FP.g,"\n","\n")

} ##### end of simulation replicates

cbind(b.true, b.net[-1], b.net.R)[b.true!=0, ]
cbind(b.true, b.net[-1], b.net.R)[b.true==0, ]


test = simulation();
x2 = test[[1]]; y2 = test[[2]]
validation <- function(b){
  sum((y2 - x2 %*% b)^2)/length(y2)
}
validation(b.net[-1])
validation(b.net.R)

###################################### SKCM #################################################
load(file ="../TCGA/SKCM_G567_s200r3.RData")
load(file ="../TCGA/SKCM_clinic.RData")
head(clc)
X=as.matrix(cbind(clc$age, as.numeric(as.character(clc$gender)),G)); X[1:5, 1:9]
# X=as.matrix(cbind(clc$age, G)); X[1:5, 1:9]
inds = which(apply(X, 2, function(t) sum(t==0)) >300); inds
# X = X[, -inds]
Y = (log(clc$breslow)); hist(Y)

lambda.n = rev(exp(seq(1,45,1)/5 -7))
response = "continuous"; standardize=TRUE; verbo = FALSE
penalty = "network"
folds = 5
r=5
clv = (1:2)
init =  "elnet"
lamb.1= lambda.n;
lamb.2=10
res=cv.regnet(X, Y, response = response, penalty = penalty, folds = folds, lamb.1=lamb.1, lamb.2=lamb.2, r=r,
              clv =clv, initiation =init)
res$lambda

fit = regnet(X, Y, response=c("continuous"), penalty=penalty, res$lambda[1], res$lambda[2], r = r,
             clv =clv, initiation = init)
# fit = regnet(X, Y, response=c("continuous"), penalty=penalty, 0.05, 1, r = r, clv =clv, initiation = init)
b.rgn = fit$b
plot(Y, fit$residual); abline(a=0, b=0, lty=2)
sum(b.rgn!=0)
as.numeric(b.rgn[b.rgn!=0])

data = as.data.frame(cbind(Y, X[, 1:150]))
fit = glm(Y~., data=data)
plot(Y, fit$residuals); abline(a=0, b=0, lty=2)

plot(X[,5], Y)

lasso.cv = cv.glmnet(X,Y, alpha=1, nfolds=5)       #compute lambda by crossvalidation
alpha = lasso.cv$lambda.min
lasso.fit = glmnet(X,Y,family="gaussian", alpha=1, nlambda=100, intercept = TRUE)
coef0 = as.vector(predict(lasso.fit, s=alpha, type="coefficients"))[-2]    #init
as.numeric(coef0[coef0!=0])
residual = Y - cbind(X) %*% coef0
plot(Y, residual); abline(a=0, b=0, lty=2)
