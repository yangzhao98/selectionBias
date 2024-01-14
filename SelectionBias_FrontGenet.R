SelectionBiasFG <- function(or_study) {
    # Reference: Use of multivariable Mendelian randomization to address biases 
    # due to competing risk before recruitment. (2020). C Mary Schooling, 
    # Priscilla M Lopez, Zhao Yang, J V Zhao, SL Au Yeung, Jian V Huang.
    # doi: https://doi.org/10.1101/716621.
    ## -----------------------------------------------------------------------------
    or_true <- 1.0
    Px  <- seq(0,0.95,0.05)
    Pcr <- seq(0,0.95,0.01)
    
    or_obs <- function(Px, Pcr, or_true) {
        return(or_true * (1 - (Px + Pcr))/((1 - Px)*(1 - Pcr)))
    }
    
    PxDat <- PcrDat <- orDat <- NULL
    for (i in 1:length(Px)) {
        for (j in 1:length(Pcr)) {
            if ((Px[i] + Pcr[j]) < 1){
                PxDat <- c(PxDat, Px[i])
                PcrDat <- c(PcrDat, Pcr[j])
                orDat <- c(orDat, or_obs(Px[i], Pcr[j], or_true))
            }
        }
    }
    
    PxLabel <- paste(Px*100, "%", sep = "")
    dat <- data.frame(Px = PxDat,
                      Pcr = PcrDat,
                      or_obs = orDat,
                      or_obs2 = 1/orDat)
    dat$Px <- factor(dat$Px, levels = Px, labels = PxLabel)
    
    p1 <- ggplot(data = dat, aes(x = Pcr, y = or_obs, col = Px)) +
        geom_line(aes(col = Px)) +
        geom_point(aes(col = Px), size = 0.3) +
        geom_hline(aes(yintercept = ifelse(or_study < 1, or_study, round(1/or_study,2))),
                   lty = 1, color = "gray") +
        geom_hline(aes(yintercept = 1), lty = 2, color = "gray") +
        theme_classic() +
        theme(legend.position = "right") +
        scale_x_continuous(name = "Proportion already dead from competing events",
                           breaks = seq(0,1,0.1),
                           labels = paste(seq(0,1,0.1), sep = ""),
                           limits = c(0,1), expand = c(0,0)) +
        # scale_y_continuous(name = "Observed odds ratio (Protective effect)",
        #                    breaks = seq(0,or_true,0.1),
        #                    labels = paste(seq(0,or_true,0.1), sep = ""),
        #                    limits = c(0,or_true+0.05), expand = c(0,0)) +
        scale_y_log10(name = "Observed odds ratio (Protective effect)",
                      breaks = c(0.02,0.05, seq(0.1,1.0,0.1)),
                      limits = c(0.02,or_true+0.05),
                      expand = c(0,0)) +
        guides(color = guide_legend(title = "Proportion already \ndead from exposure (%)")) +
        ggtitle(paste("The true odds ratio is ", or_true, sep = ""))
    p1
    p2 <- ggplot(data = dat, aes(x = Pcr, y = or_obs2, col = Px)) +
        geom_line(aes(col = Px)) +
        geom_point(aes(col = Px), size = 0.3) +
        geom_hline(aes(yintercept = 1), lty = 2, color = "gray") +
        geom_hline(aes(yintercept = ifelse(or_study > 1, or_study, round(1/or_study,2))),
                   lty = 1, color = "gray") +
        theme_classic() +
        theme(legend.position = "right") +
        scale_x_continuous(name = "Proportion already dead from competing events",
                           breaks = seq(0,1,0.1),
                           labels = paste(seq(0,1,0.1), sep = ""),
                           limits = c(0,1), expand = c(0,0)) +
        # scale_y_continuous(name = "Observed odds ratio (Protective effect)",
        #                    breaks = seq(0,or_true,0.1),
        #                    labels = paste(seq(0,or_true,0.1), sep = ""),
        #                    limits = c(0,or_true+0.05), expand = c(0,0)) +
        scale_y_log10(name = "Observed odds ratio (Hazard effect)",
                      breaks = seq(1,10,0.5),
                      limits = c(0.95, 10),
                      expand = c(0.,0)) +
        guides(color = guide_legend(title = "Proportion already \ndead from exposure (%)")) +
        ggtitle(paste("The true odds ratio is ", or_true, sep = ""))
    p2
    
    # p <- ggpubr::ggarrange(p1, p2, nrow = 1, labels = c("A", "B"),
    #                        common.legend = TRUE, legend = "bottom")
    # p
    if (or_study < 1) {
        datTable <- dat[round(dat$or_obs,2) == or_study,c("Px", "Pcr", "or_obs")]
        datTable$Pcr <- paste(datTable$Pcr*100, "%", sep = "")
        names(datTable) <- c("Pr(Deaths|Exposure)",
                             "Pr(Deaths|Competing events)",
                             "Observed OR")
        return(list(SBTable = datTable, plot = p1))
    } else {
        datTable <- dat[round(dat$or_obs2,2) == or_study, c("Px", "Pcr", "or_obs2")]
        datTable$Pcr <- paste(datTable$Pcr*100, "%", sep = "")
        names(datTable) <- c("Pr(Deaths|Exposure)",
                             "Pr(Deaths|Competing events)",
                             "Observed OR")
        return(list(SBTable = datTable, plot = p2))
    }
}


## betaGY = 0.97 for IL1RA-COPD
sbAnalysis1 <- SelectionBiasFG(or_study = 0.97)
sbAnalysis1$SBTable
sbAnalysis1 <- SelectionBiasFG(or_study = 0.80)
sbAnalysis1$SBTable
tiff(file = "SelectionBias-IL-1Ra.tiff", width = 550, height = 505.3)
sbAnalysis1$plot
dev.off()
## betaGY = 0.98 for IL18BP-COPD
sbAnalysis2 <- SelectionBiasFG(or_study = 0.98)
sbAnalysis2$SBTable
sbAnalysis2 <- SelectionBiasFG(or_study = 0.73)
sbAnalysis2$SBTable
tiff(file = "SelectionBias-IL-18BP.tiff", width = 550, height = 505.3)
sbAnalysis2$plot
dev.off()
