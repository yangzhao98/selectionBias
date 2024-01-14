
## --- PerformSelectionBias2SMR()
## - Ref. IJE (2018); 48(3): 691-701. doi: https://doi.org/10.1093/ije/dyy202
## - Perform simulation studies to investigate the impact of selection bias on
## -   the 2-sample MR analysis, in which
## --- Params:
## -  (1) nTotal = the total sample size
## -  (2) nSample = the sample size of the biased samplings
## -  (3) betaGX = SNP-Exposure association
## -  (4) betaUS = Confounder-Exposure association
## -  (5) betaXY = Exposure-Outcome causal effect
## -  (6) betaUY = Confounder-Outcome association
## -  (7) SelectionBiasbeta0 = Pr(Selection cases)
## -  (8) SelectionBiasbetaX = the effect of exposure on the magnitude of selection bias
## -  (9) SelectionBiasbetaU = the effect of confounder on the magnitude of selection bias
## -  (10)SelectionBiasbetaY = the effect of outcome on the magnitude of selection bias
## -  (11)TrimCutoff = the cutoff of the Inverse-Probability-Weighting, in which
## -         the top 1% IPW was replaced with the 99th percentile of IPW
## -  (12)Seed = the value of random seed
## -------------------------------------------------------------------------- ##
## -- Date: November 6, 2019
## -- 
## -- Version 2: Adding a new method using residual weighted (ResW) method
## --
## -------------------------------------------------------------------------- ##

PerformSelectionBias2SMR <- function(nTotal, nSample, 
                                     betaGX, betaUX, 
                                     betaXY, betaUY,
                                     SelectionBiasbeta0, 
                                     SelectionBiasbetaX,
                                     SelectionBiasbetaU,
                                     SelectionBiasbetaY,
                                     TrimCutoff,
                                     nSimulation) {
    suppressPackageStartupMessages(library(TwoSampleMR))
    suppressPackageStartupMessages(library(MRInstruments))
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(ggpubr))
    suppressPackageStartupMessages(library(ggsci))
    
    SelectionBias2SMR <- function(nTotal = nTotal, nSample = nSample, 
                                  betaGX = betaGX, betaUX = betaUX, 
                                  betaXY = betaXY, betaUY = betaUY,
                                  SelectionBiasbeta0 = SelectionBiasbeta0, 
                                  SelectionBiasbetaX = SelectionBiasbetaX,
                                  SelectionBiasbetaU = SelectionBiasbetaU,
                                  SelectionBiasbetaY = SelectionBiasbetaY,
                                  TrimCutoff = TrimCutoff,
                                  Seed = 6) {
        set.seed(Seed)
        ## Random error term
        epsilonX <- rnorm(nTotal)
        epsilonY <- rnorm(nTotal)
        ## Gene expression
        G <- rnorm(nTotal)
        ## Confounder info
        U <- rnorm(nTotal)
        X <- betaGX * G + betaUX * U + sqrt(1 - betaGX^2 - betaUX^2) * epsilonX
        Y <- betaXY * X + betaUY * U + sqrt(1 - betaXY^2 - betaUY^2) * epsilonY
        
        ## Selection info
        logitpi <- SelectionBiasbeta0 + SelectionBiasbetaX * X + 
            SelectionBiasbetaU * U + SelectionBiasbetaY * Y
        
        SBP <- function(x) { exp(x)/(exp(x) + 1) }
        SelectionProb <- SBP(logitpi)
        #SProb <- SelectionProb/sum(SelectionProb)
        idx <- runif(nTotal) < SelectionProb
        
        fit <- glm(idx ~ X + U + Y, family = binomial)
        SProb <- predict(fit, type = "response")
        #SProb <- predProb/sum(predProb)
        
        ## Inverse probability-of-selection weighting method
        IPSW <- rep(NA, rep = nTotal)
        IPSW[idx == 1] <- 1/SProb[idx == 1]
        IPSW[idx == 0] <- 1/(1 - SProb[idx == 0])
        ## Trimming of weights by setting the largest 1% of weights to be equal to the 99th perdentile
        cutProb <- quantile(IPSW, probs = TrimCutoff)
        IPSW[IPSW > cutProb] <- cutProb
        
        ## Residual weighting method
        ResW <- as.numeric(fit$residuals)
        ResW[ResW == 0] <- 1
        cutResWU <- quantile(ResW, probs = TrimCutoff)
        ResW[ResW > cutResWU] <- cutResWU
        cutResWL <- quantile(ResW, probs = 1 - TrimCutoff)
        ResW[ResW < cutResWL] <- cutResWL
        
        datO <- data.frame(Xs = X, Ys = Y,
                           Us = U, Gs = G,
                           IPSW = IPSW,
                           ResW = ResW,
                           SProb = SProb)
        dat <- datO[idx,] %>%
            arrange(desc(SProb)) %>%
            sample_n(nSample)  ## Random selection 
        #   slice(1:nSample)   ## Select with the highest probability
        #ndat <- nrow(datO)
        #idx2 <- sample(ndat, 
        #               nSample,
        #               replace = FALSE)
        
        ## Selected participants
        #dat <- datO[idx2,]
        medianSelectionProb <- median(dat$SProb)
        
        idxT <- sample(nTotal, nSample, replace = FALSE)
        datT <- datO[idxT,]
        ## Without IPSW method
        iv2slsT <- AER::ivreg(Ys ~ Us + Xs | Xs + Gs,
                             data = datT)
        betaXYTrue <- summary(iv2slsT)$coefficients[3,1:2]
        pValueTrue <- 2*pnorm(-abs(betaXYTrue[1]/betaXYTrue[2]))
        biasTrue <- betaXYTrue[1] - betaXY        
        
        ## Without IPSW method
        iv2sls <- AER::ivreg(Ys ~ Us + Xs | Xs + Gs,
                             data = dat)
        betaXYbar <- summary(iv2sls)$coefficients[3,1:2]
        pValue <- 2*pnorm(-abs(betaXYbar[1]/betaXYbar[2]))
        bias <- betaXYbar[1] - betaXY
        
        ## With IPSW method
        iv2slsIPSW <- AER::ivreg(Ys ~ Us + Xs | Xs + Gs,
                                 weights = IPSW,
                                 data = dat)
        betaXYbarIPSW <- summary(iv2slsIPSW)$coefficients[3,1:2]
        pValueIPSW <- 2*pnorm(-abs(betaXYbarIPSW[1]/betaXYbarIPSW[2]))
        biasIPSW <- betaXYbarIPSW[1] - betaXY  
        
        ## With residual weighted method
        iv2slsResW <- AER::ivreg(Ys ~ Us + Xs | Xs + Gs,
                                 weights = ResW,
                                 data = dat)
        betaXYbarResW <- summary(iv2slsResW)$coefficients[3,1:2]
        pValueResW <- 2*pnorm(-abs(betaXYbarResW[1]/betaXYbarResW[2]))
        biasResW <- betaXYbarResW[1] - betaXY  
        
        Res <- data.frame(
            Method = c("2SLS(NoSB)", "2SLS(SB)", "IPSW-2SLS(SB)", "ResW-2SLS(SB)"),
            betaXY = c(betaXYTrue[1], betaXYbar[1], betaXYbarIPSW[1], betaXYbarResW[1]),
            sebetaXY = c(betaXYTrue[2] ,betaXYbar[2], betaXYbarIPSW[2], betaXYbarResW[2]),
            pNormal = c(pValueTrue, pValue, pValueIPSW, pValueResW),
            medianSelectionProb = c(medianSelectionProb,
                                    medianSelectionProb,
                                    medianSelectionProb,
                                    medianSelectionProb),
            Bias = c(biasTrue, bias, biasIPSW, biasResW)
        )
        return(Res)
    }
    
    
    ## Perform the simulation study
    ResSum <- data.frame()
    for (nSim in 1:nSimulation) {
        if( (nSim %% 50) == 0 ) cat("Simulation = ", nSim, "\n")
        tmp <- SelectionBias2SMR(
            nTotal = nTotal, nSample = nSample, 
            betaGX = alphaG, betaUX = alphaU, 
            betaXY = betaX, betaUY = betaU, 
            SelectionBiasbeta0 = SelectionBiasbeta0, 
            SelectionBiasbetaX = SelectionBiasbetaX, 
            SelectionBiasbetaU = SelectionBiasbetaU, 
            SelectionBiasbetaY = SelectionBiasbetaY,
            TrimCutoff = TrimCutoff,
            Seed = nSim
        )
        ResSum <- rbind(ResSum, tmp)
    }
    
    ## Type 1 error rate
    Type1Error <- ResSum %>%
        group_by(Method) %>%
        mutate(idx = (pNormal < 0.05),
               nTest = abs(betaXY/sebetaXY),
               dm = nTest > 1.96) %>%
        summarise(MeanbetaXY = mean(betaXY),
                  MedianbetaXY = median(betaXY),
                  MeansebetaXY = mean(sebetaXY),
                  MediansebetaXY = median(sebetaXY),
                  MedianSBProb = median(medianSelectionProb),
                  TypeI = mean(idx),
                  TypeINew = mean(dm))
    
    BiasCutoffU <- quantile(ResSum$Bias, probs = 0.99)
    ResSum$Bias[ResSum$Bias > BiasCutoffU] <- BiasCutoffU
    BiasCutoffL <- quantile(ResSum$Bias, probs = 0.01)
    ResSum$Bias[ResSum$Bias < BiasCutoffL] <- BiasCutoffL
    
    ## Histogram of selection bias
    p1 <- ggplot(data = ResSum, aes(Bias, fill = Method)) +
        geom_histogram(aes(y = ..density.., fill = Method), 
                       alpha = 0.5, position = "identity") + 
        geom_vline(aes(xintercept = 0), col = "grey") + 
        theme_classic() +
        theme(legend.position = c(0.2,0.8)) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        scale_color_jama() + scale_fill_jama()
    
    
    ## Boxplot of selection bias
    p2 <- ggplot(data = ResSum, 
                 aes(x = Method, y = Bias, col = Method)) +
        geom_boxplot() +
        stat_summary(fun.y = mean, geom = "point", 
                     shape = 23, size = 3) +
        geom_hline(aes(yintercept = 0), col = "grey") +
        theme_classic() +
        theme(legend.position = "none") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        scale_color_jama() + scale_fill_jama()
    
    
    ## Results
    return(list(SimulationResult = ResSum,
                Type1ErrorRate = Type1Error,
                SelectionBiasHistogram = p1,
                SelectionBiasBoxplot = p2))
    
}

nSimulation <- 300
nTotal <- 10e4
nSample <- 10e2
alphaG <- sqrt(0.02)
alphaU <- sqrt(0.2)
betaU <- sqrt(0.2)
betaX <- 0
gamma0 <- 1.0
TrimCutoff <- 0.99

## Scenario 1. Selection bias due to exposure
gammaX <- 1.0
gammaU <- 0.0
gammaY <- 1.0
S1Res <- PerformSelectionBias2SMR(nTotal = nTotal, nSample = nSample,
                                  betaGX = alphaG, betaUX = alphaU,
                                  betaXY = betaX, betaUY = betaU,
                                  SelectionBiasbeta0 = gamma0,
                                  SelectionBiasbetaX = gammaX,
                                  SelectionBiasbetaU = gammaU,
                                  SelectionBiasbetaY = gammaY,
                                  TrimCutoff = TrimCutoff,
                                  nSimulation = nSimulation)
# 
# ## Scenario 2. Selection bias due to confounder
# gammaX <- 0.0
# gammaU <- 2.0
# gammaY <- 0.0
# S2Res <- PerformSelectionBias2SMR(nTotal = nTotal, nSample = nSample, 
#                                   betaGX = alphaG, betaUX = alphaU, 
#                                   betaXY = betaX, betaUY = betaU,
#                                   SelectionBiasbeta0 = gamma0, 
#                                   SelectionBiasbetaX = gammaX,
#                                   SelectionBiasbetaU = gammaU,
#                                   SelectionBiasbetaY = gammaY,
#                                   TrimCutoff = TrimCutoff,
#                                   nSimulation = nSimulation)
# 
# ## Scenario 3. Selection bias due to outcome
# gammaX <- 0.0
# gammaU <- 0.0
# gammaY <- 2.0
# S3Res <- PerformSelectionBias2SMR(nTotal = nTotal, nSample = nSample, 
#                                   betaGX = alphaG, betaUX = alphaU, 
#                                   betaXY = betaX, betaUY = betaU,
#                                   SelectionBiasbeta0 = gamma0, 
#                                   SelectionBiasbetaX = gammaX,
#                                   SelectionBiasbetaU = gammaU,
#                                   SelectionBiasbetaY = gammaY,
#                                   TrimCutoff = TrimCutoff,
#                                   nSimulation = nSimulation)
# 
# ## Scenario 4. Selection bias due to both exposure and confounder
# gammaX <- 2.0
# gammaU <- 2.0
# gammaY <- 0.0
# S4Res <- PerformSelectionBias2SMR(nTotal = nTotal, nSample = nSample, 
#                                   betaGX = alphaG, betaUX = alphaU, 
#                                   betaXY = betaX, betaUY = betaU,
#                                   SelectionBiasbeta0 = gamma0, 
#                                   SelectionBiasbetaX = gammaX,
#                                   SelectionBiasbetaU = gammaU,
#                                   SelectionBiasbetaY = gammaY,
#                                   TrimCutoff = TrimCutoff,
#                                   nSimulation = nSimulation)
# 
