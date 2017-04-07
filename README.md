
<!-- README.md is generated from README.Rmd. Please edit that file -->
regnet
======

This package provide procedures for fitting network-based regularization, minimax concave penalty (MCP) and lasso penalty for generalized linear model. This first version includes procedures for logistic regression only. We plan to add functions for survival response and gene expression as soon as we could.

Examples
--------

### Example.1 (CV.NetLogistic)

    result = CV.NetLogistic(regnet$X, regnet$Y, lamb.2 = 1, r = 4.5)  
    result$lambda  
    b = NetLogistic(regnet$X, regnet$Y, result$lambda[1,2], result$lambda[1,1])  
    index = which(regnet$beta != 0)  
    pos = which(b != 0)  
    tp = length(intersect(index, pos))  
    fp = length(pos) - tp  
    list(tp=tp, fp=fp)  

### Example.2 (CV.McpLogistic)

    result = CV.McpLogistic(regnet$X, regnet$Y, r = 4.5)  
    result$lambda  
    b = McpLogistic(regnet$X, regnet$Y, result$lambda[1])  
    index = which(regnet$beta != 0)  
    pos = which(b != 0)  
    tp = length(intersect(index, pos))  
    fp = length(pos) - tp  
    list(tp=tp, fp=fp)
