rm(list=ls())

library(magick)
library(ggplot2)
library(netmeta)
library(tidyverse)
library(extrafont)
loadfonts()

# any issue or question? Let me know on Twitter @EGOstinelli

# data prep
setwd('~/Desktop/V plot') # set up your wd

outcome.names <- c('ES', 'EL', 'AS', 'AL', 'TS', 'TL', 'Sa') # a character vector with your outcomes
noutcome <- length(outcome.names)
reference = 'placebo'
db <- readRDS('db.rds') # load your dataset as data.frame (optional, you can provide an external list of interventions)
treat.names <- unique(db$your_interventions)
active.treat <- treat.names[!treat.names %in% reference]
active.treat <- active.treat[!active.treat %in% NA] # check, remove NAs

# df
data.ES <- readRDS('ES.rds') # load outcome-specific data.frames
data.EL <- readRDS('EL.rds')
data.AS <- readRDS('AS.rds')
data.AL <- readRDS('AL.rds')
data.TS <- readRDS('TS.rds')
data.TL <- readRDS('TL.rds')
data.Sa <- readRDS('Sa.rds')

# add estimates of placebo events for continuous outcomes
## several options are available, with different requirements. We opted for the 50% improvement from baseline, considering the direction of the scale.
## for this, we need information of baseline score and the best possible value of the scale.
## if you plan to use both endpoint and change scores, you will need to indicate this in a specific column (e.g. 'type' in our example).
## alternatively, these values can be #1 externally calculated or #2 provided by another study, or #3 be fictional.

data.ES$threshold <- data.ES$baseline+abs((data.ES$best_value_scale-data.ES$baseline)/2)
data.ES$prop <- 1-pnorm((data.ES$threshold-ifelse(data.ES$type=='endpoint', data.ES$mean, data.ES$baseline+data.ES$mean))/data.ES$sd)
data.ES$r <- round(data.ES$prop*data.ES$n,0)
data.ES$r[data.ES$intervention!=reference] <- NA
sum(data.ES$r[data.ES$intervention==reference], na.rm = T)/sum(data.ES$n[data.ES$intervention==reference & !is.na(data.ES$r)], na.rm = T)

data.EL$threshold <- data.EL$baseline+abs((data.EL$best_value_scale-data.EL$baseline)/2)
data.EL$prop <- 1-pnorm((data.EL$threshold-ifelse(data.EL$type=='endpoint', data.EL$mean, data.EL$baseline+data.EL$mean))/data.EL$sd)
data.EL$r <- round(data.EL$prop*data.EL$n,0)
data.EL$r[data.EL$intervention!=reference] <- NA
sum(data.EL$r[data.EL$intervention==reference], na.rm = T)/sum(data.EL$n[data.EL$intervention==reference & !is.na(data.EL$r)], na.rm = T)

store0 <- list(data.ES, data.EL, data.AS, data.AL, data.TS, data.TL, data.Sa) # data for each outcome

# netmeta objects
result.ES <- readRDS('ES_NMA.rds') # load previously saved netmeta objects (if already in your environment, you don't need this step)
result.EL <- readRDS('EL_NMA.rds')
result.AS <- readRDS('AS_NMA.rds')
result.AL <- readRDS('AL_NMA.rds')
result.TS <- readRDS('TS_NMA.rds')
result.TL <- readRDS('TL_NMA.rds')
result.Sa <- readRDS('Sa_NMA.rds')
store <- list(result.ES, result.EL, result.AS, result.AL, result.TS, result.TL, result.Sa) #NMA results

# this assumes you have used SMDs for your continuous values, using the Hasselblad & Hedges conversion
result.ES$TE.random <- (pi/sqrt(3))*result.ES$TE.random
result.ES$seTE <- (pi/sqrt(3))*result.ES$seTE

result.EL$TE.random <- (pi/sqrt(3))*result.EL$TE.random
result.EL$seTE <- (pi/sqrt(3))*result.EL$seTE

#####

final <- data.frame() # data frame to store results

for (i in 1:noutcome) {
  
  data <- store0[[i]]
  net1 <- store[[i]]
  
  # treatment estimate (odds ratio) from netmetas (different paths if MH or IV)
  OR.pla <- data.frame('drug' = colnames(net1$TE.random))
  OR.pla$logOR <- if (net1$method == 'MH') {net1$TE.fixed[,reference]} else {net1$TE.random[,reference]}
  OR.pla$seTE <- if (net1$method == 'MH') {net1$seTE.fixed[,reference]} else {net1$seTE.random[,reference]}
  OR.pla$OR <- exp(OR.pla$logOR)
  OR.pla <- OR.pla[-which(OR.pla$drug==reference),] # exclude comparison with placebo itself which is 0
  
  # meta analysis of event rates in placebo
  meta.pla = metaprop(event = round(data$r[data$intervention==reference]), n = data$n[data$intervention==reference], method = 'GLMM')
  rate.pla = exp(meta.pla$TE.fixed)/(1+exp(meta.pla$TE.fixed))
  odds.pla=rate.pla/(1-rate.pla)
  
  # calculate event rate for treatment
  OR.pla$event.rate <- round(OR.pla$OR*odds.pla/(1+OR.pla$OR*odds.pla),digits=3)
  
  # calculate Z-scores accounting for clinically important risk difference
  clinically.important.RD.0 <- 0.0
  risk.drugs.0 <- clinically.important.RD.0+rate.pla
  OR.import.0 <- risk.drugs.0/(1-risk.drugs.0)/((rate.pla)/(1-rate.pla))
  OR.pla$Zscore.0 <- (OR.pla$logOR-log(OR.import.0))/OR.pla$seTE
  
  outcome.result <- data.frame(outcome = paste0(outcome.names[i]), drug = OR.pla$drug, Zscore = OR.pla$Zscore.0,
                               event.rate = round(OR.pla$event.rate*100,1), rate.pla = rate.pla, 
                               logOR = OR.pla$logOR, seTE = OR.pla$seTE)
  
  for (j in 1:length(active.treat)) {
    if (!(active.treat[j] %in% outcome.result$drug)) {
      outcome.result[nrow(outcome.result)+1,] <- NA
      outcome.result$drug[nrow(outcome.result)] <- paste0(active.treat[j])
      outcome.result$outcome[nrow(outcome.result)] <- paste0(outcome.names[i])
    }
  }
  
  final <- rbind(final, outcome.result)
  
}

final_data <- final
final_data <- final_data[,c('outcome', 'drug', 'Zscore', 'event.rate', 'rate.pla')]

# add in control treatment (i.e. treatment = 1)
for (k in 1:noutcome) {
  final_data[nrow(final_data)+1,] <- NA
  final_data$drug[nrow(final_data)] <- reference
  final_data$outcome[nrow(final_data)] <- paste0(outcome.names[k])
  final_data$event.rate[nrow(final_data)] <- round(unique(final$rate.pla[final$outcome==outcome.names[k]])*100,1)[1]
  final_data$rate.pla[nrow(final_data)] <- unique(final$rate.pla[final$outcome==outcome.names[k]])[1]
  }

# you need to specify the interpretation of your outcomes (i.e. the more the better versus the less the better)
final_data$interpretation <- 'tmtb'
final_data$interpretation[final_data$outcome == 'AS'] <- 'tltb'
final_data$interpretation[final_data$outcome == 'AL'] <- 'tltb'
final_data$interpretation[final_data$outcome == 'TS'] <- 'tltb'
final_data$interpretation[final_data$outcome == 'TL'] <- 'tltb'
final_data$interpretation[final_data$outcome == 'Sa'] <- 'tltb'

final_data$Zscore2 <- final_data$Zscore #truncated z scores
final_data$Zscore2[final_data$Zscore2 < -3] = -3
final_data$Zscore2[final_data$Zscore2 > 3] = 3
final_data$Zscore2[final_data$interpretation == 'tltb'] <- -final_data$Zscore2[final_data$interpretation == 'tltb']

############## vitruvian plot ----
nrp <- final_data
nrp <- rename(nrp, value = event.rate)

# landmarks
eracle <- ceiling(max(nrp$value, na.rm = T)) # [ACTION] if you want to limit the plot until a specific value
piecount <- as.numeric(length(unique(nrp$outcome)))
linespace <- 10 # [ACTION] distance between lines, change as needed

# curved text ext (section title)
ctcount <- 12 # [ACTION] the higher this number, the smaller the relative space between the characters of section title curved text characters
ctangle <- (360/piecount)/ctcount
ctspace <- 1/ctcount
size.text.ext <- 6 # [ACTION] size of section title text
ctheight <- eracle*1.18 # [ACTION] distance of section title text from plot limits
radial.ext.col <- '#002147' # [ACTION] colour of section title text
title1.coord <- 2 # [ACTION] position of section title1 text
title2.coord <- 6 # [ACTION] position of section title2 text

# curved text int (outcome)
ctcount_end <- 27 # [ACTION] the higher this number, the smaller the relative space between the characters of outcome curved text characters
ctangle_end <- (360/piecount)/ctcount_end
ctspace_end <- 1/ctcount_end
radial.int.col <- '#002147' # [ACTION] colour of outcome text on light background
radial.int.col.dark <- '#ffffff' # [ACTION] colour of outcome text on dark background
radial.int.bg <- '#DDF1FB' # [ACTION] colour of outcome labels (background) #DDF1FB blue #E7E6E6 grey #F0ECFB purple #F1E0DD Lancet
radial.int.bg.2 <- '#C8E1F3' # [ACTION] colour of outcome labels (background)
radial.int.bg.3 <- '#A6D4EF' # [ACTION] colour of outcome labels (background)
radial.int.bg.4 <- '#4EA6DD' # [ACTION] colour of outcome labels (background)
radial.int.alpha <- 0.2 # [ACTION] alpha level of outcome text
size.text.int <- 6 # [ACTION] size of outcome text
outcome1.coord <- 1 # [ACTION] position of outcome1 label
outcome2.coord <- 2 # [ACTION] position of outcome2 label
outcome3.coord <- 3 # [ACTION] position of outcome3 label
outcome4.coord <- 4 # [ACTION] position of outcome4 label
outcome5.coord <- 5 # [ACTION] position of outcome5 label
outcome6.coord <- 6 # [ACTION] position of outcome6 label
outcome7.coord <- 7 # [ACTION] position of outcome7 label
outcome.height <- eracle*1.053 # [ACTION] change as needed
outcome.reverse.height <- eracle*1.05 # [ACTION] change as needed

# df var prep
nrp$outcome <- factor(nrp$outcome, levels = c('EL', 'AL', 'TL', 'Sa', 'TS', 'AS', 'ES'), # [ACTION] position from 12:00 clockwise
                      c('EL', 'AL', 'TL', 'Sa', 'TS', 'AS', 'ES'))

# folder prep
ifelse(dir.exists('vitruvian plots'), print('folder already existing'), dir.create('vitruvian plots'))

# function
nightingale <- function (drugname) {
  ggplot(drugname, aes(x = outcome, y = value)) +
    
    # graph
    geom_point(x = (piecount/2)+0.5, y = eracle*0.57, shape = 0, size = 201, colour = '#002147', fill = '#ffffff', stroke = 1.5) + # [ACTION] change y as needed, especially with higher values of eracle
    geom_rect(xmin = 0, xmax = piecount, ymin = 0, ymax = eracle*1.2, fill = '#ffffff') + # background
    
    # title [ACTION] change the following lines if you are plotting head-to-head comparisons
    geom_label(x = (piecount/2)+0.5, y = eracle*1.26, label = strrep(' ', 69), size = 7, fill = '#002147', 
               color = '#002147', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines')) +
    geom_text(x = (piecount/2)+0.5, y = eracle*1.26, label = (toupper(drugname$drug[1])), size = 7, colour = '#ffffff', family = 'Trebuchet MS', fontface = 'bold') +
    
    geom_rect(xmin = 0, xmax = piecount, ymin = 0, ymax = eracle, fill = '#f9fcff') +
    
    geom_col(aes(fill = round(Zscore2,2)), width = 1) +
    scale_fill_gradient2(low = '#EE7E7A', mid = '#F7EED6', high = '#2aad66', na.value = '#DDF1FB',
                         breaks = c(-2.575829, -1.959964, -1.644854, 0, 1.644854, 1.959964, 2.575829), limits = c(-3, 3), 
                         labels = c('p < 0.01', 'p = 0.05', 'p = 0.1', 'p = 1.00', 'p = 0.1', 'p = 0.05','p < 0.01')) +
    
    geom_rect(xmin = outcome1.coord-0.5, xmax = outcome1.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'EL']), 1, 0)) +
    geom_rect(xmin = outcome2.coord-0.5, xmax = outcome2.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AL']), 1, 0)) +
    geom_rect(xmin = outcome3.coord-0.5, xmax = outcome3.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TL']), 1, 0)) +
    geom_rect(xmin = outcome4.coord-0.5, xmax = outcome4.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'Sa']), 1, 0)) +
    geom_rect(xmin = outcome5.coord-0.5, xmax = outcome5.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TS']), 1, 0)) +
    geom_rect(xmin = outcome6.coord-0.5, xmax = outcome6.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AS']), 1, 0)) +
    geom_rect(xmin = outcome7.coord-0.5, xmax = outcome7.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'ES']), 1, 0)) +
    
    geom_hline(yintercept = seq(linespace*2, eracle, by = linespace*2), color = '#002147', size = 0.3, alpha = 0.2) +
    geom_hline(yintercept = seq(linespace, eracle, by = linespace*2), color = '#002147', size = 0.55, alpha = 0.2, linetype = 'dotted') +
    geom_vline(xintercept = seq(0.5, piecount+0.5, by = 1), color = '#002147', size = 0.3, alpha = 0.5) +
    geom_vline(xintercept = c(0.5, 3.5, 4.5), color = '#002147', linetype = 1, size = 0.7) +
    
    # outcome labels
    geom_rect(xmin = outcome1.coord-0.5, xmax = outcome1.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'EL']), 0, 1)) +
    geom_rect(xmin = outcome2.coord-0.5, xmax = outcome2.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AL']), 0, 1)) +
    geom_rect(xmin = outcome3.coord-0.5, xmax = outcome3.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TL']), 0, 1)) +
    geom_rect(xmin = outcome4.coord-0.5, xmax = outcome4.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'Sa']), 0, 1)) +
    geom_rect(xmin = outcome5.coord-0.5, xmax = outcome5.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TS']), 0, 1)) +
    geom_rect(xmin = outcome6.coord-0.5, xmax = outcome6.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AS']), 0, 1)) +
    geom_rect(xmin = outcome7.coord-0.5, xmax = outcome7.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'ES']), 0, 1)) +
    
    # columns and title
    geom_hline(yintercept = eracle, color = '#002147', linetype = 1, size = 1.2) +
    xlab(toupper(drugname$drug[1])) +
    
    # add active intervention values to the active intervention plots
    geom_label(x = outcome1.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'EL'],0)<10,
                                                                  paste0(' ', round(drugname$value[drugname$outcome == 'EL'],0), '%'), 
                                                                  paste0(round(drugname$value[drugname$outcome == 'EL'],0), '%')), 
              size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'EL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
  
    geom_label(x = outcome2.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'AL'])<10,
                                                                  paste0(' ', round(drugname$value[drugname$outcome == 'AL'],0), '%'), 
                                                                  paste0(round(drugname$value[drugname$outcome == 'AL'],0), '%')), 
              size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
  
    geom_label(x = outcome3.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'TL'])<10,
                                                                  paste0(' ', round(drugname$value[drugname$outcome == 'TL'],0), '%'), 
                                                                  paste0(round(drugname$value[drugname$outcome == 'TL'],0), '%')), 
              size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
  
    geom_label(x = outcome4.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'Sa'])<10,
                                                                  paste0(' ', round(drugname$value[drugname$outcome == 'Sa'],0), '%'), 
                                                                  paste0(round(drugname$value[drugname$outcome == 'Sa'],0), '%')), 
              size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'Sa']), 0, ifelse(drugname$drug == reference, 0, 1))) +
  
    geom_label(x = outcome5.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'TS'])<10,
                                                                  paste0(' ', round(drugname$value[drugname$outcome == 'TS'],0), '%'), 
                                                                  paste0(round(drugname$value[drugname$outcome == 'TS'],0), '%')), 
              size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TS']), 0, ifelse(drugname$drug == reference, 0, 1))) +
  
    geom_label(x = outcome6.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'AS'])<10,
                                                                  paste0(' ', round(drugname$value[drugname$outcome == 'AS'],0), '%'), 
                                                                  paste0(round(drugname$value[drugname$outcome == 'AS'],0), '%')), 
              size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AS']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    geom_label(x = outcome7.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'ES'])<10,
                                                                  paste0(' ', round(drugname$value[drugname$outcome == 'ES'],0), '%'), 
                                                                  paste0(round(drugname$value[drugname$outcome == 'ES'],0), '%')), 
              size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'ES']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    # add reference values to the reference plot
    geom_point(x = outcome1.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'EL']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    geom_text(x = outcome1.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'EL']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'EL']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    
    geom_point(x = outcome2.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AL']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    geom_text(x = outcome2.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'AL']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AL']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    
    geom_point(x = outcome3.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TL']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    geom_text(x = outcome3.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'TL']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TL']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    
    geom_point(x = outcome4.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'Sa']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    geom_text(x = outcome4.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'Sa']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'Sa']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    
    geom_point(x = outcome5.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TS']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    geom_text(x = outcome5.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'TS']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TS']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    
    geom_point(x = outcome6.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AS']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    geom_text(x = outcome6.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'AS']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AS']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    
    geom_point(x = outcome7.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'ES']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    geom_text(x = outcome7.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'ES']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'ES']), 0, ifelse(drugname$drug == reference, 1, 0))) +
    
    # add reference values in the active intervention plots
    geom_point(x = outcome1.coord, y = eracle*0.7, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'EL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    geom_text(x = outcome1.coord, y = eracle*0.7, label = paste0(round(drugname$rate.pla[drugname$outcome == 'EL']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'EL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    geom_point(x = outcome2.coord, y = eracle*0.7, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    geom_text(x = outcome2.coord, y = eracle*0.7, label = paste0(round(drugname$rate.pla[drugname$outcome == 'AL']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    geom_point(x = outcome3.coord, y = eracle*0.7, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    geom_text(x = outcome3.coord, y = eracle*0.7, label = paste0(round(drugname$rate.pla[drugname$outcome == 'TL']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TL']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    geom_point(x = outcome4.coord, y = eracle*0.7, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'Sa']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    geom_text(x = outcome4.coord, y = eracle*0.7, label = paste0(round(drugname$rate.pla[drugname$outcome == 'Sa']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'Sa']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    geom_point(x = outcome5.coord, y = eracle*0.7, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TS']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    geom_text(x = outcome5.coord, y = eracle*0.7, label = paste0(round(drugname$rate.pla[drugname$outcome == 'TS']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'TS']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    geom_point(x = outcome6.coord, y = eracle*0.7, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AS']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    geom_text(x = outcome6.coord, y = eracle*0.7, label = paste0(round(drugname$rate.pla[drugname$outcome == 'AS']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'AS']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    geom_point(x = outcome7.coord, y = eracle*0.7, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'ES']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    geom_text(x = outcome7.coord, y = eracle*0.7, label = paste0(round(drugname$rate.pla[drugname$outcome == 'ES']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'ES']), 0, ifelse(drugname$drug == reference, 0, 1))) +
    
    # polar coordinate system
    coord_polar(clip = 'off', direction = 1) +
    
    # title1 (section title)
    geom_text(x=(title1.coord-(ctspace*9)), y=ctheight, label='L', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*9) +
    geom_text(x=(title1.coord-(ctspace*8)), y=ctheight, label='O', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*8) +
    geom_text(x=(title1.coord-(ctspace*7)), y=ctheight, label='N', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*7) +
    geom_text(x=(title1.coord-(ctspace*6)), y=ctheight, label='G', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*6) +
    geom_text(x=(title1.coord-(ctspace*5)), y=ctheight, label='-', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*5) +
    geom_text(x=(title1.coord-(ctspace*4)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*4) +
    geom_text(x=(title1.coord-(ctspace*3)), y=ctheight, label='E', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*3) +
    geom_text(x=(title1.coord-(ctspace*2)), y=ctheight, label='R', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*2) +
    geom_text(x=(title1.coord-(ctspace*1)), y=ctheight, label='M', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)+ctangle*1) +
    geom_text(x=title1.coord, y=ctheight, label=' ', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)) +
    geom_text(x=(title1.coord+(ctspace*1)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*1) +
    geom_text(x=(title1.coord+(ctspace*2)), y=ctheight, label='R', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*2) +
    geom_text(x=(title1.coord+(ctspace*3)), y=ctheight, label='E', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*3) +
    geom_text(x=(title1.coord+(ctspace*4)), y=ctheight, label='A', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*4) +
    geom_text(x=(title1.coord+(ctspace*5)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*5) +
    geom_text(x=(title1.coord+(ctspace*6)), y=ctheight, label='M', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*6) +
    geom_text(x=(title1.coord+(ctspace*7)), y=ctheight, label='E', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*7) +
    geom_text(x=(title1.coord+(ctspace*8)), y=ctheight, label='N', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*8) +
    geom_text(x=(title1.coord+(ctspace*9)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title1.coord+0.5)-ctangle*9) +
    
    # title2 (section title)
    geom_text(x=(title2.coord-(ctspace*7)), y=ctheight, label='A', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)+ctangle*7) +
    geom_text(x=(title2.coord-(ctspace*6)), y=ctheight, label='C', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)+ctangle*6) +
    geom_text(x=(title2.coord-(ctspace*5)), y=ctheight, label='U', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)+ctangle*5) +
    geom_text(x=(title2.coord-(ctspace*4)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)+ctangle*4) +
    geom_text(x=(title2.coord-(ctspace*3)), y=ctheight, label='E', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)+ctangle*3) +
    geom_text(x=(title2.coord-(ctspace*2)), y=ctheight, label=' ', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)+ctangle*2) +
    geom_text(x=(title2.coord-(ctspace*1)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)+ctangle*1) +
    geom_text(x=title2.coord, y=ctheight, label='R', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)) +
    geom_text(x=(title2.coord+(ctspace*1)), y=ctheight, label='E', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)-ctangle*1) +
    geom_text(x=(title2.coord+(ctspace*2)), y=ctheight, label='A', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)-ctangle*2) +
    geom_text(x=(title2.coord+(ctspace*3)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)-ctangle*3) +
    geom_text(x=(title2.coord+(ctspace*4)), y=ctheight, label='M', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)-ctangle*4) +
    geom_text(x=(title2.coord+(ctspace*5)), y=ctheight, label='E', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)-ctangle*5) +
    geom_text(x=(title2.coord+(ctspace*6)), y=ctheight, label='N', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)-ctangle*6) +
    geom_text(x=(title2.coord+(ctspace*7)), y=ctheight, label='T', size = size.text.ext, colour = '#002147', angle = (360/piecount)*(piecount-title2.coord+0.5)-ctangle*7) +
    
    # symptom improvement (efficacy long)
    geom_text(x=outcome1.coord-(ctspace_end*10), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*10), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*9), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*8), y=outcome.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*7), y=outcome.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*6), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*5), y=outcome.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*4), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*3), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*2), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*1), y=outcome.height, label='n', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord, y=outcome.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*1), y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*2), y=outcome.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*3), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*4), y=outcome.height, label='m', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*5), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*6), y=outcome.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*7), y=outcome.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*8), y=outcome.height, label='v', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*9), y=outcome.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*10), y=outcome.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*10), alpha = radial.int.alpha, family='Andale Mono') +    
    
    # all-cause drop-outs (acceptability long)
    geom_text(x=outcome2.coord-(ctspace_end*9), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*8), y=outcome.height, label='l', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*7), y=outcome.height, label='l', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*6), y=outcome.height, label='-', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*5), y=outcome.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*4), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*3), y=outcome.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*2), y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*1), y=outcome.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*0), y=outcome.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*0), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*1), y=outcome.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*2), y=outcome.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*3), y=outcome.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*4), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*5), y=outcome.height, label='-', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*6), y=outcome.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*7), y=outcome.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*8), y=outcome.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*9), y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    
    # drop-outs due to AE (tolerability long)
    geom_text(x=outcome3.coord-(ctspace_end*9.5), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*9)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*8.5), y=outcome.reverse.height, label='E', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*9)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*7.5), y=outcome.reverse.height, label='A', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*8)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*6.5), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*7)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*5.5), y=outcome.reverse.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*6)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*4.5), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*3.5), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*4)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*2.5), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*3)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*1.5), y=outcome.reverse.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*2)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*0.5), y=outcome.reverse.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*1)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*0.5), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*0)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*1.5), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*1)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*2.5), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*2)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*3.5), y=outcome.reverse.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*3)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*4.5), y=outcome.reverse.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*4)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*5.5), y=outcome.reverse.height, label='-', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*6.5), y=outcome.reverse.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*6)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*7.5), y=outcome.reverse.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*7)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*8.5), y=outcome.reverse.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*8)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*9.5), y=outcome.reverse.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*9)+180, alpha = radial.int.alpha, family='Andale Mono') +
    
    # participants with AE (safety)
    geom_text(x=outcome4.coord-(ctspace_end*10), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*9.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*9), y=outcome.reverse.height, label='E', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*9.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*8), y=outcome.reverse.height, label='A', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*8.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*7), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*7.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*6), y=outcome.reverse.height, label='h', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*6.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*5), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*5.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*4), y=outcome.reverse.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*4.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*3), y=outcome.reverse.height, label='w', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*3.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*2), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*2.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*1), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*1.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*0), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*0.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*1), y=outcome.reverse.height, label='n', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*0.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*2), y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*1.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*3), y=outcome.reverse.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*2.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*4), y=outcome.reverse.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*3.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*5), y=outcome.reverse.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*4.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*6), y=outcome.reverse.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*5.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*7), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*6.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*8), y=outcome.reverse.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*7.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*9), y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*8.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*10), y=outcome.reverse.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*9.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    
    # drop-outs due to AE (tolerability short)
    geom_text(x=outcome5.coord-(ctspace_end*9.5), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*9)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*8.5), y=outcome.reverse.height, label='E', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*9)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*7.5), y=outcome.reverse.height, label='A', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*8)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*6.5), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*7)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*5.5), y=outcome.reverse.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*6)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*4.5), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*3.5), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*4)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*2.5), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*3)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*1.5), y=outcome.reverse.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*2)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*0.5), y=outcome.reverse.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*1)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*0.5), y=outcome.reverse.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*0)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*1.5), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*1)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*2.5), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*2)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*3.5), y=outcome.reverse.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*3)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*4.5), y=outcome.reverse.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*4)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*5.5), y=outcome.reverse.height, label='-', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*6.5), y=outcome.reverse.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*6)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*7.5), y=outcome.reverse.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*7)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*8.5), y=outcome.reverse.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*8)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*9.5), y=outcome.reverse.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*9)+180, alpha = radial.int.alpha, family='Andale Mono') +
    
    # all-cause drop-outs (acceptability short)
    geom_text(x=outcome6.coord-(ctspace_end*9), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*8), y=outcome.height, label='l', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*7), y=outcome.height, label='l', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*6), y=outcome.height, label='-', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*5), y=outcome.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*4), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*3), y=outcome.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*2), y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord-(ctspace_end*1), y=outcome.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)+ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*0), y=outcome.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*0), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*1), y=outcome.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*2), y=outcome.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*3), y=outcome.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*4), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*5), y=outcome.height, label='-', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*6), y=outcome.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*7), y=outcome.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*8), y=outcome.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome6.coord+(ctspace_end*9), y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome6.coord+0.5)-ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    
    # symptom improvement (efficacy short)
    geom_text(x=outcome7.coord-(ctspace_end*10), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*10), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*9), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*8), y=outcome.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*7), y=outcome.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*6), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*5), y=outcome.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*4), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*3), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*2), y=outcome.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord-(ctspace_end*1), y=outcome.height, label='n', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)+ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord, y=outcome.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*1), y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*2), y=outcome.height, label=' ', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*3), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*4), y=outcome.height, label='m', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*5), y=outcome.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*6), y=outcome.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*6), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*7), y=outcome.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*7), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*8), y=outcome.height, label='v', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*8), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*9), y=outcome.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*9), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome7.coord+(ctspace_end*10), y=outcome.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome7.coord+0.5)-ctangle_end*10), alpha = radial.int.alpha, family='Andale Mono') +
    
    # labels
    geom_text(x = 0.5, y = 10, label = ' 10%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-10)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 20, label = ' 20%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-20)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 30, label = ' 30%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-30)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 40, label = ' 40%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-40)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 50, label = ' 50%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-50)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 60, label = ' 60%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-60)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 70, label = ' 70%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-70)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 80, label = ' 80%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-80)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 90, label = ' 90%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-90)>=(linespace/2), ifelse(drugname$drug == reference, 0.15, 0), 0)) +
    
    # extra-lines
    geom_segment(x = 0.5, y = 0, xend = 0.5, yend = eracle*1.1, colour = "#002147", size = 0.7) +
    geom_segment(x = 3.5, y = 0, xend = 3.5, yend = eracle*1.1, colour = "#002147", size = 0.7) +
    geom_segment(x = 4.5, y = 0, xend = 4.5, yend = eracle*1.1, colour = "#002147", size = 0.7) +
    
    # theme
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = '#FFFFFF', color='#FFFFFF'),
          plot.margin = unit(c(0.2,0.1,0.7,0.1), 'in'), # it follows the trbl order (top, right, bottom, left)
          legend.position = 'none'
    )
}

plots <- lapply(split(nrp,nrp$drug), nightingale)

# loop - save
lapply(names(plots), function(x) ggsave(filename = paste('vitruvian plots/', x, '.png', sep=''), plot = plots[[x]], dpi = 300, width = 11, height = 11))

############## legend
legendplot <- ggplot(final_data, aes(outcome, drug)) + 
  geom_tile(aes(fill = round(Zscore2,2)), colour = 'white') + 
  geom_text(aes(label= event.rate), size = 6) +
  scale_fill_gradient2(low = '#EE7E7A', mid = '#F7EED6', high = '#2aad66', na.value = '#DDF1FB', breaks = c(-2.575829, -1.959964, -1.644854, 0, 1.644854, 1.959964, 2.575829), limits = c(-3, 3), labels = c('p < 0.01', 'p = 0.05', 'p = 0.1', 'p = 1.00', 'p = 0.1', 'p = 0.05','p < 0.01'))+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 15)) +
  labs(x = '',y = '') +
  theme(legend.title = element_blank(),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),legend.position = 'left', legend.text = element_text(size = 12)) +
  scale_x_discrete(position = 'top')

legend <- cowplot::get_legend(legendplot)
grid::grid.newpage()
png('vitruvian plots/legend.png', res = 600, width = 5.5, height = 5.5, units = 'in') 
grid::grid.draw(legend)
dev.off()

# synoptic charts ----
ifelse(dir.exists('synoptic chart'), print('folder already existing'), dir.create('synoptic chart'))
pattern <- paste0('.png$')
Lfile <- list.files('vitruvian plots/', pattern = pattern)

legend0 <- image_border(image_read('vitruvian plots/legend.png'), 'white', '70x120')
blank <- image_blank(3300, 3300, color = 'white', pseudo_image = '', defines = NULL)
BDZ.s <- image_border(image_read(paste('vitruvian plots/', Lfile[3], sep = '')), 'white', '70x120')
BDZ.i <- image_border(image_read(paste('vitruvian plots/', Lfile[1], sep = '')), 'white', '70x120')
BDZ.l <- image_border(image_read(paste('vitruvian plots/', Lfile[2], sep = '')), 'white', '70x120')
daridorexant <- image_border(image_read(paste('vitruvian plots/', Lfile[4], sep = '')), 'white', '70x120')
diphenydramine <- image_border(image_read(paste('vitruvian plots/', Lfile[5], sep = '')), 'white', '70x120')
doxepin <- image_border(image_read(paste('vitruvian plots/', Lfile[6], sep = '')), 'white', '70x120')
doxylamine <- image_border(image_read(paste('vitruvian plots/', Lfile[7], sep = '')), 'white', '70x120')
eszopiclone <- image_border(image_read(paste('vitruvian plots/', Lfile[8], sep = '')), 'white', '70x120')
lemborexant <- image_border(image_read(paste('vitruvian plots/', Lfile[10], sep = '')), 'white', '70x120')
melatonin <- image_border(image_read(paste('vitruvian plots/', Lfile[11], sep = '')), 'white', '70x120')
mirtazapine <- image_border(image_read(paste('vitruvian plots/', Lfile[12], sep = '')), 'white', '70x120')
placebo <- image_border(image_read(paste('vitruvian plots/', Lfile[13], sep = '')), 'white', '70x120')
propiomazine <- image_border(image_read(paste('vitruvian plots/', Lfile[14], sep = '')), 'white', '70x120')
quetiapine <- image_border(image_read(paste('vitruvian plots/', Lfile[15], sep = '')), 'white', '70x120')
ramelteon <- image_border(image_read(paste('vitruvian plots/', Lfile[16], sep = '')), 'white', '70x120')
seltorexant <- image_border(image_read(paste('vitruvian plots/', Lfile[17], sep = '')), 'white', '70x120')
suvorexant <- image_border(image_read(paste('vitruvian plots/', Lfile[18], sep = '')), 'white', '70x120')
trazodone <- image_border(image_read(paste('vitruvian plots/', Lfile[19], sep = '')), 'white', '70x120')
trimipramine <- image_border(image_read(paste('vitruvian plots/', Lfile[20], sep = '')), 'white', '70x120')
zaleplon <- image_border(image_read(paste('vitruvian plots/', Lfile[21], sep = '')), 'white', '70x120')
zolpidem <- image_border(image_read(paste('vitruvian plots/', Lfile[22], sep = '')), 'white', '70x120')
zopiclone <- image_border(image_read(paste('vitruvian plots/', Lfile[23], sep = '')), 'white', '70x120')

## multiple plots ----
row1 <- c(legend0, BDZ.s, BDZ.i, BDZ.l)
row2 <- c(blank, doxepin, eszopiclone, lemborexant)
row3 <- c(blank, melatonin, ramelteon, suvorexant)
row4 <- c(blank, trimipramine, zaleplon, zolpidem)
row1.img <- image_append(row1)
row2.img <- image_append(row2)
row3.img <- image_append(row3)
row4.img <- image_append(row4)
all.plots <- c(row1.img, row2.img, row3.img, row4.img)
all.plots.img <- image_append(all.plots, stack = TRUE)
image_write(all.plots.img, path = paste0('synoptic chart/synoptic chart.png'), format = 'png')

## all plots, horizontal ----
row1 <- c(placebo, BDZ.s, BDZ.i, BDZ.l, daridorexant, diphenydramine)
row2 <- c(doxepin, doxylamine, eszopiclone, lemborexant, melatonin, mirtazapine)
row3 <- c(propiomazine, quetiapine, ramelteon, seltorexant, suvorexant, trazodone)
row4 <- c(trimipramine, zaleplon, zolpidem, zopiclone, blank, legend0)
row1.img <- image_append(row1)
row2.img <- image_append(row2)
row3.img <- image_append(row3)
row4.img <- image_append(row4)
all.plots <- c(row1.img, row2.img, row3.img, row4.img)
all.plots.img <- image_append(all.plots, stack = TRUE)
image_write(all.plots.img, path = paste0('synoptic chart/all plots h.png'), format = 'png')

## all plots, vertical ----
row1 <- c(placebo, BDZ.s, BDZ.i, BDZ.l)
row2 <- c(daridorexant, diphenydramine, doxepin, doxylamine)
row3 <- c(eszopiclone, lemborexant, melatonin, mirtazapine)
row4 <- c(propiomazine, quetiapine, ramelteon, seltorexant)
row5 <- c(suvorexant, trazodone, trimipramine, zaleplon)
row6 <- c(zolpidem, zopiclone, blank, legend0)
row1.img <- image_append(row1)
row2.img <- image_append(row2)
row3.img <- image_append(row3)
row4.img <- image_append(row4)
row5.img <- image_append(row5)
row6.img <- image_append(row6)
all.plots <- c(row1.img, row2.img, row3.img, row4.img, row5.img, row6.img)
all.plots.img <- image_append(all.plots, stack = TRUE)
image_write(all.plots.img, path = paste0('synoptic chart/all plots v.png'), format = 'png')

