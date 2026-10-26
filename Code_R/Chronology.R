#This script aims at generating the chronology sequence of Lagaccione

#Load necessary packages
{
  #Plot
  library(ggplot2)
  library(gridExtra)
  
  #Bchron
  library(Bchron)
  
  #Font
  library(showtext)
  font_add(family = "arial", regular = file.path("Packages_R", "Fonts", "arial.ttf"))
  showtext_auto()
}

#Chronology
{
  #Read raw radiocarbon date for Bayesian approach
  rawRadData <- read.csv(file.path('Data_R', 'Lag_Chron_Bchron.csv'))
  
  #Run calibrate functions to get the calibrated dates at each radiocarbon dates
  dateCali <- BchronCalibrate(ages = rawRadData$C14_age,
                              ageSds = rawRadData$error,
                              calCurves = rawRadData$model, 
                              positions = rawRadData$depth,
                              ids = rawRadData$lab_ID)
  
  testa <- sampleAges(dateCali)
  
  testb <- apply(testa, 2, quantile, probs = c(0.05, 0.5, 0.95))
  
  testb <- t(testb)
  
  testc <- apply(testa, 2, sd)
  
  testb <- cbind(testb, testc)
  
  #Calibrated ages
  write.csv(testb, file = file.path('Output_R', 'Calibrated_Ages_24_Dates.csv'))
  
  pdf(file = file.path('Figures_R', 'Radiocarbon_Dates_24_Dates.pdf'),   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 12)
  
  #Plot the calibrated results
  plot(dateCali) +
    labs(x = 'Age (years BP)',
         y = 'Depth (cm)',
         title = 'Calibrated dates')
  
  dev.off()
  
  #Get the credible intervals for all radiocarbon dates
  age_Bey_20 <- sampleAges(dateCali)
  CI_Bey_20 <- apply(age_Bey_20, 2, quantile, prob = c(0.025, 0.975))
  
  #Run the Bayesian approach to get the age-depth model
  modBchron <- with(rawRadData, 
                    Bchronology(ages=C14_age,
                                ageSds=error, 
                                calCurves=model,
                                positions=depth, 
                                positionThicknesses=thickness,
                                ids=lab_ID, 
                                predictPositions=seq(1, 4200, by = 1)))
  #Output the age-depth sequence
  pdf(file = file.path('Figures_R', 'Chronology_24_Dates.pdf'),   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 12)
  
  #Plot the results
  plot(modBchron) +
    labs(x = 'Age (Cal yr BP)',
         y = 'Depth (cm)') +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, family = 'arial'), 
      axis.title.x = element_text(size = 12, family = 'arial'),  
      axis.title.y = element_text(size = 12, family = 'arial'), 
      axis.text = element_text(size = 12, family = 'arial'), 
      axis.text.y = element_text(size = 12, family = 'arial'), 
    )
  
  dev.off()
  
  #Get the data of predicted ages at given depth sequence
  predict <- predict(modBchron, 
                         newPositions = seq(1, 4200, by = 1))
  
  #Get the mean and 95% percentile values of the predicted ages at each 0.5 cm depth
  predictPercen <- matrix(data = NA, nrow = ncol(predict), 
                              ncol = 4,
                              dimnames = list(seq(1, 4200, by = 1),
                                              c('Depth_(cm)', 'Median', 'L95', 'U95')))
  
  #Fill in the matrix
  predictPercen[, 1] <- seq(1, 4200, by = 1)
  predictPercen[, 2] <- apply(predict, 2, median)
  loopI <- 1
  while (loopI <= nrow(predictPercen)) {
    
    predictPercen[loopI, 3:4] <- quantile(predict[, loopI], probs = c(0.05, 0.95))
    
    loopI <- loopI + 1
  }
  
  predictPercen <- as.data.frame(predictPercen)
  
  # #Output the Bayesian results
  write.csv(predictPercen,file = file.path('Output_R', 'Lagaccione_24_Dates.csv'))

}