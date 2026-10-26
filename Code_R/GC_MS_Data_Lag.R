#This script aims at visualization GC-MS data from Lago
#Load necessary packages
{
  #Data reading
  library(here)
  
  #Text font
  library(showtext)
  font_add(family = "arial", regular = here('Fonts', 'arial.ttf'))
  showtext_auto()
  
  #Plot
  library(ggplot2)
  library(gridExtra)
  library(patchwork)
  
  #Tukey HSD
  library(multcomp)
  
  #Data oration 
  library(reshape2)
}

#Declare necessary functions
{
  lm_eqn <- function(data, xind, yind){
    m <- lm(data[, yind] ~ data[, xind]);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
}

#Read data - 9 dates
{
  #Read raw data
  rawGC <- read.csv(file = file.path('Data_R', '26_11_2025_lag_data_2025.csv'),
                    row.names = 1)
  
  #Remove NA with 0
  rawGC <- rawGC[-which(is.na(rawGC[, 8])),]
  
  #Read chronology
  chronLag <- read.csv(file = file.path('Output_R', 'Lagaccione_9_Dates.csv'),
                       row.names = 1)
  
  #Align the depth with age information
  rawGC <- cbind(rawGC, Age = approx(x = chronLag[, 1], y = chronLag[, 2],
                                     xout = rawGC[, 13] * 100)[[2]]/1000)
  
  tempVec <- approx(x = chronLag[, 1], y = chronLag[, 2],
                           xout = rawGC[, 13] * 100)[[2]]/1000
  
  #Remove rows without age information
  rawGC <- rawGC[-which(is.na(rawGC[, 14])),]
  
  #Loop to get the mean and variance of each sample
  #Create a matrix
  gcMeanMat9 <- matrix(nrow = length(unique(rawGC[, 13])),
                      ncol = 4,
                      dimnames = list(c(),
                                      c('Mean', 'sd', 
                                        colnames(rawGC)[13:14]
                                      )
                      )
  )
  
  #Fill in the matrix with batch information
  gcMeanMat9[,3:4] <- as.matrix(unique(rawGC[, 13:14]))
  
  #Loop to fill in the mean and sd for all samples
  for (loopI in 1:nrow(gcMeanMat9)) {
    
    if (length(which((rawGC[, 13] == gcMeanMat9[loopI, 3])
    )
    ) > 1
    ){
      
      gcMeanMat9[loopI, 1] <- mean(as.numeric(rawGC[which((rawGC[, 13] == gcMeanMat9[loopI, 3])
      ), 
      8]
      )
      )
      
      gcMeanMat9[loopI, 2] <- sd(as.numeric(rawGC[which((rawGC[, 13] == gcMeanMat9[loopI, 3])
      ), 
      8]
      )
      )
      
    } else {
      
      gcMeanMat9[loopI, 1] <- as.numeric(rawGC[which((rawGC[, 13] == gcMeanMat9[loopI, 3]) 
      ), 
      8]
      )
      
      gcMeanMat9[loopI, 2] <- 0
      
    }
    
  }
  
  rm(tempVec)
}

#Plot
{
  gcLagPlot9 <- ggplot(data = gcMeanMat9) +
    geom_point(aes(x = Age, y = Mean), size = 1.5, color = 'black', alpha = 0.5) +
    geom_errorbar(aes(x = Age, 
                      ymin = Mean - sd,
                      ymax = Mean + sd),
                  width = 0.1) +
    geom_line(aes(x = Age, y = Mean), linewidth = 0.5, color = 'red') +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 40, xmax = 41,
    #          ymin = -Inf, ymax = Inf) +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 37, xmax = 38.5,
    #          ymin = -Inf, ymax = Inf) +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 35, xmax = 36,
    #          ymin = -Inf, ymax = Inf) +
    scale_x_reverse(limit = c(46, 35), breaks = seq(46, 35, by = -1)) +
    labs(x = 'Age (ka)', y = 'pCA (ng/150gr)') +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(size = 16, family = 'arial'),
      axis.text.x = element_text(size = 14, family = 'arial'),
      
      axis.text.y = element_text(size = 14, family = 'arial'),
      axis.title.y = element_text(size = 16, family = 'arial'),
    )
  
  gcLagPlot9Dep <- ggplot(data = gcMeanMat9) +
    geom_point(aes(x = Depth, y = Mean), size = 1.5, color = 'black', alpha = 0.5) +
    geom_errorbar(aes(x = Depth, 
                      ymin = Mean - sd,
                      ymax = Mean + sd),
                  width = 0.1) +
    geom_line(aes(x = Depth, y = Mean), linewidth = 0.5, color = 'red') +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 40, xmax = 41,
    #          ymin = -Inf, ymax = Inf) +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 37, xmax = 38.5,
    #          ymin = -Inf, ymax = Inf) +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 35, xmax = 36,
    #          ymin = -Inf, ymax = Inf) +
    scale_x_reverse(limit = c(35, 19), breaks = seq(35, 19, by = -1)) +
    labs(x = 'Depth (m)', y = 'pCA (ng/150gr)') +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(size = 16, family = 'arial'),
      axis.text.x = element_text(size = 14, family = 'arial'),
      
      axis.text.y = element_text(size = 14, family = 'arial'),
      axis.title.y = element_text(size = 16, family = 'arial'),
    )
  
  #Output the age-depth sequence
  pdf(file = file.path('Figures_R', 'GC-MS_pCA_9_Dates.pdf'),   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 12)
  
  gcLagPlot9 /
    gcLagPlot9Dep
  
  dev.off()
}

#Read data - 8 dates
{
  #Read raw data
  rawGC <- read.csv(file = file.path('Data_R', '26_11_2025_lag_data_2025.csv'),
                    row.names = 1)
  
  #Remove NA with 0
  rawGC <- rawGC[-which(is.na(rawGC[, 8])),]
  
  #Read chronology
  chronLag <- read.csv(file = file.path('Output_R', 'Lagaccione_8_Dates.csv'),
                       row.names = 1)
  
  #Align the depth with age information
  rawGC <- cbind(rawGC, Age = approx(x = chronLag[, 1], y = chronLag[, 2],
                                     xout = rawGC[, 13] * 100)[[2]]/1000)
  
  tempVec <- approx(x = chronLag[, 1], y = chronLag[, 2],
                    xout = rawGC[, 13] * 100)[[2]]/1000
  
  #Remove rows without age information
  rawGC <- rawGC[-which(is.na(rawGC[, 14])),]
  
  #Loop to get the mean and variance of each sample
  #Create a matrix
  gcMeanMat8 <- matrix(nrow = length(unique(rawGC[, 13])),
                      ncol = 4,
                      dimnames = list(c(),
                                      c('Mean', 'sd', 
                                        colnames(rawGC)[13:14]
                                      )
                      )
  )
  
  #Fill in the matrix with batch information
  gcMeanMat8[,3:4] <- as.matrix(unique(rawGC[, 13:14]))
  
  #Loop to fill in the mean and sd for all samples
  for (loopI in 1:nrow(gcMeanMat8)) {
    
    if (length(which((rawGC[, 13] == gcMeanMat8[loopI, 3])
    )
    ) > 1
    ){
      
      gcMeanMat8[loopI, 1] <- mean(as.numeric(rawGC[which((rawGC[, 13] == gcMeanMat8[loopI, 3])
      ), 
      8]
      )
      )
      
      gcMeanMat8[loopI, 2] <- sd(as.numeric(rawGC[which((rawGC[, 13] == gcMeanMat8[loopI, 3])
      ), 
      8]
      )
      )
      
    } else {
      
      gcMeanMat8[loopI, 1] <- as.numeric(rawGC[which((rawGC[, 13] == gcMeanMat8[loopI, 3]) 
      ), 
      8]
      )
      
      gcMeanMat8[loopI, 2] <- 0
      
    }
    
  }
  
  rm(tempVec)
}

#Plot
{
  gcLagPlot8 <- ggplot(data = gcMeanMat8) +
    geom_point(aes(x = Age, y = Mean), size = 1.5, color = 'black', alpha = 0.5) +
    geom_errorbar(aes(x = Age, 
                      ymin = Mean - sd,
                      ymax = Mean + sd),
                  width = 0.1) +
    geom_line(aes(x = Age, y = Mean), linewidth = 0.5, color = 'red') +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 40, xmax = 41,
    #          ymin = -Inf, ymax = Inf) +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 37, xmax = 38.5,
    #          ymin = -Inf, ymax = Inf) +
    # annotate("rect", fill = "blue", alpha = 0.5,
    #          xmin = 35, xmax = 36,
    #          ymin = -Inf, ymax = Inf) +
    scale_x_reverse(limit = c(43, 35), breaks = seq(44, 35, by = -1)) +
    labs(x = 'Age (ka)', y = 'pCA (ng/150gr)') +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(size = 16, family = 'arial'),
      axis.text.x = element_text(size = 14, family = 'arial'),
      
      axis.text.y = element_text(size = 14, family = 'arial'),
      axis.title.y = element_text(size = 16, family = 'arial'),
    )
  
  #Output the age-depth sequence
  pdf(file = file.path('Figures_R', 'GC-MS_pCA_8_Dates.pdf'),   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6)
  
  gcLagPlot8
  
  dev.off()
  
  #Output the age-depth sequence
  pdf(file = file.path('Figures_R', 'GC-MS_pCA_Comp.pdf'),   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6)
  
  gcLagPlot8 / 
    gcLagPlot9
  
  dev.off()
}

#Read pollen concentration data
{
  #Setup parent directory for reading data
  setwd(here())
  setwd('..')
  tempDir <- paste(getwd(), '/Data_R/lgc_alistair.csv', sep = '', collapse = NULL)
  
  rawPollen <- read.csv(tempDir)
  
  #Combine pollen concentration with GC-MS signals
  tempVec <- approx(rawPollen[, 1], rawPollen[, 9],
                    xout = gcMeanMat[, 3])[[2]]
  
  tempVec2 <- approx(rawPollen[, 1], rawPollen[, 102],
                    xout = gcMeanMat[, 3])[[2]]
  
  tempVec3 <- log10(tempVec2)
  
  comGCPollen <- cbind(gcMeanMat, pollen = tempVec,
                       pollen_C = tempVec2,
                       pollen_C_Log = tempVec3)
  
}

#Plot scatters between GC and pollen
{
  scaPlotGCPollen <- ggplot() +
    geom_point(data = as.data.frame(comGCPollen[-47, ]),
               aes(y = Mean, x = pollen),
               size = 2, color = 'blue') + 
    geom_point(data = as.data.frame(comGCPollen[47,]),
               aes(y = comGCPollen[47, 1], x = comGCPollen[47, 5]),
               size = 2, color = 'red') + 
    geom_smooth(data = as.data.frame(comGCPollen[-47, ]),
                aes(y = Mean, x = pollen),
                method = "lm", color = 'blue') +
    geom_smooth(data = as.data.frame(comGCPollen),
                aes(y = Mean, x = pollen),
                method = "lm", color = 'red') +
    geom_text(
      aes(x = 20, y = 0.45,
          label = lm_eqn(comGCPollen[-47, ], 5, 1)), 
      size = 4,
      parse = TRUE) +
    labs(y = 'GC-MS (ng/150gr)',
         x = 'Pinus percentage (%)') + 
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 20, family = 'arial'),
          axis.text.y = element_text(size = 16, family = 'arial'),
          axis.title.x = element_text(size = 20, family = 'arial'),
          axis.text.x = element_text(size = 16, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  scaPlotGCPollenC <- ggplot() +
    geom_point(data = as.data.frame(comGCPollen[-46, ]),
               aes(y = Mean, x = pollen_C),
               size = 2, color = 'blue') + 
    geom_point(data = as.data.frame(comGCPollen[46, ]),
               aes(y = comGCPollen[46, 1], x = comGCPollen[46, 6]),
               size = 2, color = 'red') + 
    geom_smooth(data = as.data.frame(comGCPollen[-46, ]),
                aes(y = Mean, x = pollen_C),
                method = "lm", color = 'blue') +
    geom_smooth(data = as.data.frame(comGCPollen),
                aes(y = Mean, x = pollen_C),
                method = "lm", color = 'red') +
    geom_text(
      aes(x = 20000, y = 0.45,
          label = lm_eqn(comGCPollen[-46, ], 5, 1)), 
      size = 4,
      parse = TRUE) +
    labs(y = 'GC-MS (ng/150gr)',
         x = 'Pinus concentration (grains)') + 
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 20, family = 'arial'),
          axis.text.y = element_text(size = 16, family = 'arial'),
          axis.title.x = element_text(size = 20, family = 'arial'),
          axis.text.x = element_text(size = 16, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  scaPlotGCPollenCLog <- ggplot() +
    geom_point(data = as.data.frame(comGCPollen),
               aes(y = Mean, x = pollen_C_Log),
               size = 2, color = 'blue') + 
    geom_smooth(data = as.data.frame(comGCPollen),
                aes(y = Mean, x = pollen_C_Log),
                method = "lm", color = 'blue') +
    geom_text(
      aes(x = 3.5, y = 0.45,
          label = lm_eqn(comGCPollen, 5, 1)), 
      size = 4,
      parse = TRUE) +
    labs(y = 'GC-MS (ng/150gr)',
         x = 'Log-transformed concentration') + 
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 20, family = 'arial'),
          axis.text.y = element_text(size = 16, family = 'arial'),
          axis.title.x = element_text(size = 20, family = 'arial'),
          axis.text.x = element_text(size = 16, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )
  
  plotComp <- ggplot() +
    geom_bar(aes(x = comGCPollen[,4],
                 y = comGCPollen[,5]),
             fill = 'black',
             alpha = 1,
             stat = "identity",
             width = 0.1) +
    scale_x_reverse(limit = c(46, 35), breaks = seq(46, 35, by = -1)) +
    labs(y = 'Pollen percentage (%)',
         x = 'Age (ka)') + 
    theme_classic() + 
    theme(plot.title = element_blank(), 
          axis.title.y = element_text(size = 20, family = 'arial'),
          axis.text.y = element_text(size = 16, family = 'arial'),
          axis.title.x = element_text(size = 20, family = 'arial'),
          axis.text.x = element_text(size = 16, family = 'arial'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
          # axis.line.y = element_blank(),
          # axis.text.x = element_blank(), 
          # axis.line.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.ticks.x = element_blank()
    )

}

#Further investigation basing on the residual of linear model between GC and Pollen
{
  #Percentage
  {
    lm_GCPollen <- lm(comGCPollen[, 1] ~ comGCPollen[, 5])
    
    comGCPollenResi <- cbind(comGCPollen, residuals = lm_GCPollen$residuals)
    
    #Plot against ages
    gcLagPlotResi <- ggplot(data = as.data.frame(comGCPollenResi)) +
      geom_point(aes(x = Age, y = residuals), size = 1.5, color = 'black', alpha = 0.5) +
      geom_errorbar(aes(x = Age, 
                        ymin = residuals - sd,
                        ymax = residuals + sd),
                    width = 0.5) +
      geom_line(aes(x = Age, y = residuals), linewidth = 0.5, color = 'red') +
      annotate("rect", fill = "grey", alpha = 0.5, 
               xmin = 38, xmax = 41,
               ymin = -Inf, ymax = Inf) +
      annotate("rect", fill = "grey", alpha = 0.5, 
               xmin = 41.6, xmax = 43.3,
               ymin = -Inf, ymax = Inf) +
      annotate("rect", fill = "blue", alpha = 0.5, 
               xmin = 39, xmax = 40.3,
               ymin = -Inf, ymax = Inf) +
      scale_x_reverse(limit = c(46, 35), breaks = seq(46, 35, by = -1)) +
      labs(x = 'Age (ka)', y = 'pCA residuals (percentage)') +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_text(size = 16, family = 'arial'),
        axis.text.x = element_text(size = 14, family = 'arial'),
        
        axis.text.y = element_text(size = 14, family = 'arial'),
        axis.title.y = element_text(size = 16, family = 'arial'),
      )
  }
  
  #Concentration
  {
    lm_GCPollenC <- lm(comGCPollen[, 1] ~ comGCPollen[, 6])
    
    comGCPollenCResi <- cbind(comGCPollen, residuals = lm_GCPollenC$residuals)
    
    #Plot against ages
    gcLagPlotResiC <- ggplot(data = as.data.frame(comGCPollenCResi)) +
      geom_point(aes(x = Age, y = residuals), size = 1.5, color = 'black', alpha = 0.5) +
      geom_errorbar(aes(x = Age, 
                        ymin = residuals - sd,
                        ymax = residuals + sd),
                    width = 0.5) +
      geom_line(aes(x = Age, y = residuals), linewidth = 0.5, color = 'red') +
      annotate("rect", fill = "grey", alpha = 0.5, 
               xmin = 38, xmax = 41,
               ymin = -Inf, ymax = Inf) +
      annotate("rect", fill = "grey", alpha = 0.5, 
               xmin = 41.6, xmax = 43.3,
               ymin = -Inf, ymax = Inf) +
      annotate("rect", fill = "blue", alpha = 0.5, 
               xmin = 39, xmax = 40.3,
               ymin = -Inf, ymax = Inf) +
      scale_x_reverse(limit = c(46, 35), breaks = seq(46, 35, by = -1)) +
      labs(x = 'Age (ka)', y = 'pCA residuals (Concentration)') +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_text(size = 16, family = 'arial'),
        axis.text.x = element_text(size = 14, family = 'arial'),
        
        axis.text.y = element_text(size = 14, family = 'arial'),
        axis.title.y = element_text(size = 16, family = 'arial'),
      )
  }
  
  #Log_concentration
  {
    lm_GCPollenCLog <- lm(comGCPollen[, 1] ~ comGCPollen[, 7])
    
    comGCPollenCLogResi <- cbind(comGCPollen, residuals = lm_GCPollenCLog$residuals)
    
    #Plot against ages
    gcLagPlotResiCLog <- ggplot(data = as.data.frame(comGCPollenCLogResi)) +
      geom_point(aes(x = Age, y = residuals), size = 1.5, color = 'black', alpha = 0.5) +
      geom_errorbar(aes(x = Age, 
                        ymin = residuals - sd,
                        ymax = residuals + sd),
                    width = 0.5) +
      geom_line(aes(x = Age, y = residuals), linewidth = 0.5, color = 'red') +
      annotate("rect", fill = "grey", alpha = 0.5, 
               xmin = 38, xmax = 41,
               ymin = -Inf, ymax = Inf) +
      annotate("rect", fill = "grey", alpha = 0.5, 
               xmin = 41.6, xmax = 43.3,
               ymin = -Inf, ymax = Inf) +
      annotate("rect", fill = "blue", alpha = 0.5, 
               xmin = 39, xmax = 40.3,
               ymin = -Inf, ymax = Inf) +
      scale_x_reverse(limit = c(46, 35), breaks = seq(46, 35, by = -1)) +
      labs(x = 'Age (ka)', y = 'pCA residuals(log-concentration)') +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_text(size = 16, family = 'arial'),
        axis.text.x = element_text(size = 14, family = 'arial'),
        
        axis.text.y = element_text(size = 14, family = 'arial'),
        axis.title.y = element_text(size = 16, family = 'arial'),
      )
  }
}

#Output plots
{
  #Setup parent directory for reading data
  setwd(here())
  setwd('..')
  tempDir <- paste(getwd(), '/Figures_R/Lag_GC_Pollen.pdf', sep = '', collapse = NULL)
  
  #Combine plot
  pdf(file = tempDir,   # The directory you want to save the file in
      width = 18, # The width of the plot in inches
      height = 20)
  
  gcLagPlotResi / 
    gcLagPlotResiC /
    gcLagPlotResiCLog /
    (scaPlotGCPollen |  scaPlotGCPollenC | scaPlotGCPollenCLog)
  
  dev.off()
  
  #Setup parent directory for reading data
  setwd(here())
  setwd('..')
  tempDir <- paste(getwd(), '/Figures_R/Lag_GC_Pollen_Percentage.pdf', sep = '', collapse = NULL)
  
  #Combine plot
  pdf(file = tempDir,   # The directory you want to save the file in
      width = 18, # The width of the plot in inches
      height = 16)
  
  gcLagPlotResi / 
    plotComp /
    (scaPlotGCPollen |  scaPlotGCPollenC | scaPlotGCPollenCLog)
  
  dev.off()
}

