lambda.n = rev(exp(seq(1,45,1)/5 -7))
lambda.m = rev(exp(seq(1,45,1)/5 -7))
#lambda.n = rev(exp(seq(1,45,1)/4 -9))
#lambda.m = rev(exp(seq(1,45,1)/4 -9))
lambda.e = rev(exp(seq(1,45,1)/4 -9))
lambda.l = rev(exp(seq(1,45,1)/4 -9))

initiation <- function(x, y, alpha){
  lasso.cv <- glmnet::cv.glmnet(x,y, family="binomial", alpha=alpha, nfolds=5)       #compute lambda by crossvalidation
  lambda <- lasso.cv$lambda.min
  #cat("initiation lambda: ", lambda, "\n")
  lasso.fit <- glmnet::glmnet(x,y,"binomial", alpha=alpha, nlambda=50)
  coef0 <- as.vector(stats::predict(lasso.fit, s=lambda, type="coefficients"))[-2]    #initialize the coefficients vector
}

validation <- function(b, x2, y2, n){
  yi = x2 %*% b
  yi=1/(1+exp(-yi))
  y = ifelse(yi>0.5, 1, 0)
  sum(abs(y2 - y))/n        # misclassification rate
}

Soft <- function(z, lambda){
  if (z > lambda) z - lambda                 #soft thresholding operate
  else if (z < -lambda) z + lambda
  else 0
}

TruePositive <- function(b, b.true){
  index = which(b.true != 0)
  pos = which(b != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
}
