rm(list=ls())

library(magick)
library(ggplot2)
library(netmeta)
library(extrafont)
loadfonts()

noutcome <- 5 # number of outcomes
ntreat <- 7 # number of treatments

# generate dataset
generateData <- function (logOR = NULL, p.ref.1, p.ref.2) {
  
  t1 <- rep(combn(ntreat,2)[1,], each = 2) # arm 1 treatment
  t2 <- rep(combn(ntreat,2)[2,], each = 2) # arm 2 treatment
  nstudy <- length(t1) # number of studies
  
  OR <- exp(logOR)
  
  studlab <- seq(nstudy)
  n1 <- n2 <- round(runif(nstudy,300,600)) #sample size for each study
  
  p.ref <- runif(nstudy,p.ref.1, p.ref.2) # study-specific probability of an event in treatment 1
  odds.ref <- p.ref / (1 - p.ref)
  
  # define probabilities per treatment, per study arm
  odds.t1 <- odds.t2 <- r1 <- r2 <- vector()
  for(j in 1:nstudy){
    odds.t1[j] <- odds.ref[j] * OR[t1[j]]
    odds.t2[j] <- odds.ref[j] * OR[t2[j]]
  }
  p.t1 <- odds.t1/(1+odds.t1)
  p.t2 <- odds.t2/(1+odds.t2)
  
  for(j in 1:nstudy){
    r1[j] <- rbinom(1, n1[j], p.t1[j])
    r2[j] <- rbinom(1, n2[j], p.t2[j])
  }
  
  data01 <- data.frame(studlab = studlab, drug = t1, outcome = r1, n = n1)
  data02 <- data.frame(studlab = studlab, drug = t2, outcome = r2, n = n2)
  
  data <- rbind(data01, data02)
  data <- data[order(data$studlab, data$drug),]
  rownames(data) <- 1:dim(data)[1]
  
  return(data)
}

# last two parameters are range for probability of an event in treatment 1
set.seed(5)
data1 <- generateData(logOR = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 0.1, 0.15)
set.seed(6)
data2 <- generateData(logOR = c(0, -0.1, -0.3, -0.2, -0.1, -0.4, -0.3), 0.15, 0.20)
set.seed(7)
data3 <- generateData(logOR = c(0, 0.1, 0.3, -0.1, 0.1, 0.1, 0), 0.10, 0.12)
set.seed(8)
data4 <- generateData(logOR = c(0, 0.2, 0.05, 0, 0, 0.02, 0.1), 0.10, 0.17)
set.seed(9)
data5 <- generateData(logOR = c(0, 1.3, 0.5, 0.6, 0.8, 0.7, 1.2), 0.1, 0.15)
store0 <- list(data1, data2, data3, data4, data5) #data for each outcome

# Fit network meta analysis (for 5 outcomes)
p1 = pairwise(treat = drug, event = outcome, n = n, studlab = studlab, data = data1, sm = 'OR', allstudies=T)
result1 <- netmeta(p1)

p2 = pairwise(treat = drug, event = outcome, n = n, studlab = studlab, data = data2, sm = 'OR', allstudies=T)
result2 <- netmeta(p2)

p3 = pairwise(treat = drug, event = outcome, n = n, studlab = studlab, data = data3, sm = 'OR', allstudies=T)
result3 <- netmeta(p3)

p4 = pairwise(treat = drug, event = outcome, n = n, studlab = studlab, data = data4, sm = 'OR', allstudies=T)
result4 <- netmeta(p4)

p5 = pairwise(treat = drug, event = outcome, n = n, studlab = studlab, data = data5, sm = 'OR', allstudies=T)
result5 <- netmeta(p5)

store <- list(result1, result2, result3, result4, result5) # NMA results

#####

final <- data.frame() # data frame to store results

for(i in 1:noutcome){
  
  data <- store0[[i]]
  net1 <- store[[i]]
  
  # treatment estimate (odds ratio) from netmeta
  OR.pla <- data.frame('drug' = colnames(net1$TE.random))
  OR.pla$logOR <- net1$TE.random[,1]
  OR.pla$seTE <- net1$seTE.random[,1]
  OR.pla$OR <- exp(OR.pla$logOR)
  OR.pla <- OR.pla[-which(OR.pla$drug == 1),] #exclude comparison with placebo itself which is 0
  
  # meta analysis of event rates in placebo
  meta.pla = metaprop(event = round(data$outcome[data$drug==1]), n = data$n[data$drug==1], method = 'GLMM')
  rate.pla = exp(meta.pla$TE.fixed)/(1+exp(meta.pla$TE.fixed))
  odds.pla=rate.pla/(1-rate.pla)
  
  # calculate event rate for treatment
  OR.pla$value <- round(OR.pla$OR*odds.pla/(1+OR.pla$OR*odds.pla),digits=3)
  
  # Calculate Zscore accounting for clinically important risk difference
  clinically.important.RD.0 <- 0.0  
  risk.drugs.0 <- clinically.important.RD.0+rate.pla
  OR.import.0 <- risk.drugs.0/(1-risk.drugs.0)/((rate.pla)/(1-rate.pla))
  OR.pla$Zscore.0 <- (OR.pla$logOR-log(OR.import.0))/OR.pla$seTE
  
  outcome.result <- data.frame(outcome = paste('outcome', i), drug = OR.pla$drug, Zscore = OR.pla$Zscore.0,
                               value = round(OR.pla$value*100,1), rate.pla = rate.pla, 
                               logOR = OR.pla$logOR, seTE = OR.pla$seTE)
  final <- rbind(final, outcome.result)
  
}

final_data <- final
final_data <- final_data[,c('outcome', 'drug', 'Zscore', 'value', 'rate.pla')]

#add in control treatment (i.e. treatment = 1)
outcome.names <- c('outcome 1', 'outcome 2', 'outcome 3', 'outcome 4', 'outcome 5')

final_data$drug <- paste('treatment', final_data$drug)

for (k in 1:noutcome) {
  final_data[nrow(final_data)+1,] <- NA
  final_data$drug[nrow(final_data)] <- 'placebo'
  final_data$outcome[nrow(final_data)] <- paste0(outcome.names[k])
  final_data$value[nrow(final_data)] <- round(unique(final$rate.pla[final$outcome==outcome.names[k]])*100,1)[1]
  final_data$rate.pla[nrow(final_data)] <- unique(final$rate.pla[final$outcome==outcome.names[k]])[1]
}

final_data$interpretation <- 'tmtb' # the more, the better
final_data$interpretation[final_data$outcome == 'outcome 2'] <- 'tltb'
final_data$interpretation[final_data$outcome == 'outcome 3'] <- 'tltb'
final_data$interpretation[final_data$outcome == 'outcome 4'] <- 'tltb'

final_data$Zscore2 <- final_data$Zscore # truncated Z scores
final_data$Zscore2[final_data$Zscore2 < -3] = -3
final_data$Zscore2[final_data$Zscore2 > 3] = 3
final_data$Zscore2[final_data$interpretation == 'tltb'] <- -final_data$Zscore2[final_data$interpretation == 'tltb']

add_percent <- function(x){if(!is.na(x)){paste0(x, '%')} else{x}}
final_data$value2 <- sapply(final_data$value, add_percent)

############## Vitruvian plot ----
nrp <- final_data
nrp$drug[nrp$drug == 'treatment 2'] <- 'treatment 1'
nrp$drug[nrp$drug == 'treatment 3'] <- 'treatment 2'
simplify <- c('treatment 4', 'treatment 5', 'treatment 6', 'treatment 7')
nrp <- nrp[!nrp$drug %in% simplify,]

# landmarks
eracle <- ceiling(max(nrp$value)) # [ACTION] if you want to limit the plot until a specific value
piecount <- as.numeric(length(unique(nrp$outcome)))
linespace <- 10 # [ACTION] distance between lines, change as needed

# curved text ext (section title)
ctcount <- 8.5 # [ACTION] the higher this number, the smaller the relative space between the characters of section title curved text characters
ctangle <- (360/piecount)/ctcount
ctspace <- 1/ctcount
size.text.ext <- 7.5 # [ACTION] size of section title text
ctheight <- eracle*1.18 # [ACTION] distance of section title text from plot limits
radial.ext.col <- '#002147' # [ACTION] colour of section title text
title1.coord <- 5 # [ACTION] position of section title1 text
title2.coord <- 2 # [ACTION] position of section title2 text

# curved text int (outcome)
ctcount_end <- 27 # [ACTION] the higher this number, the smaller the relative space between the characters of outcome curved text characters
ctangle_end <- (360/piecount)/ctcount_end
ctspace_end <- 1/ctcount_end
radial.int.col <- '#002147' # [ACTION] colour of outcome text on light background
radial.int.col.dark <- '#ffffff' # [ACTION] colour of outcome text on dark background
radial.int.bg <- '#DDF1FB' # [ACTION] colour of outcome labels (background)
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
outcome.height <- eracle*1.053 # [ACTION] change as needed
outcome.reverse.height <- eracle*1.05 # [ACTION] change as needed

# df var prep
nrp$outcome <- factor(nrp$outcome, levels = c('outcome 1', 'outcome 2', 'outcome 3', 'outcome 4', 'outcome 5'), # [ACTION] position from 12:00 clockwise
                      labels = c(rep(paste('outcome', rep(1:piecount)))))

# folder prep
setwd('~/Desktop/V plot/graph8.0 - for V plot/figure2')
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
    geom_text(x = (piecount/2)+0.5, y = eracle*1.26, label = toupper(drugname$drug[1]), size = 7, colour = '#ffffff', family = 'Trebuchet MS', fontface = 'bold') +
    
    geom_rect(xmin = 0, xmax = piecount, ymin = 0, ymax = eracle, fill = '#f9fcff') +
    
    geom_col(aes(fill = round(Zscore2,2)), width = 1) +
    scale_fill_gradient2(low = '#fa3545', mid = '#F7EED6', high = '#00b050', na.value = '#DDF1FB', breaks = c(-2.575829, -1.959964, -1.644854, 0, 1.644854, 1.959964, 2.575829), limits = c(-3, 3), labels = c('p < 0.01', 'p = 0.05', 'p = 0.1', 'p = 1.00', 'p = 0.1', 'p = 0.05','p < 0.01')) +
    
    geom_rect(xmin = outcome1.coord-0.5, xmax = outcome1.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 1']), 1, 0)) +
    geom_rect(xmin = outcome2.coord-0.5, xmax = outcome2.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 2']), 1, 0)) +
    geom_rect(xmin = outcome3.coord-0.5, xmax = outcome3.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 3']), 1, 0)) +
    geom_rect(xmin = outcome4.coord-0.5, xmax = outcome4.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 4']), 1, 0)) +
    geom_rect(xmin = outcome5.coord-0.5, xmax = outcome5.coord+0.5, ymin = 0, ymax = eracle, fill = '#f9fcff', alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 5']), 1, 0)) +
    
    geom_hline(yintercept = seq(linespace*2, eracle, by = linespace*2), color = '#002147', size = 0.3, alpha = 0.2) +
    geom_hline(yintercept = seq(linespace, eracle, by = linespace*2), color = '#002147', size = 0.55, alpha = 0.2, linetype = 'dotted') +
    geom_vline(xintercept = seq(0.5, piecount+0.5, by = 1), color = '#002147', size = 0.3, alpha = 0.5) +
    
    # outcome labels
    geom_rect(xmin = outcome1.coord-0.5, xmax = outcome1.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 1']), 0, 1)) +
    geom_rect(xmin = outcome2.coord-0.5, xmax = outcome2.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 2']), 0, 1)) +
    geom_rect(xmin = outcome3.coord-0.5, xmax = outcome3.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 3']), 0, 1)) +
    geom_rect(xmin = outcome4.coord-0.5, xmax = outcome4.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 4']), 0, 1)) +
    geom_rect(xmin = outcome5.coord-0.5, xmax = outcome5.coord+0.5, ymin = eracle, ymax = eracle*1.1, fill = radial.int.bg, alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 5']), 0, 1)) +
    
    # columns and title
    geom_hline(yintercept = eracle, color = '#002147', linetype = 1, size = 1.2) +
    xlab(toupper(drugname$drug[1])) +
    
    # add intervention values to the intervention plots
    geom_label(x = outcome1.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'outcome 1'],0)<10,
                                                                   paste0(' ', round(drugname$value[drugname$outcome == 'outcome 1'],0), '%'), 
                                                                   paste0(round(drugname$value[drugname$outcome == 'outcome 1'],0), '%')), 
               size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 1']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_label(x = outcome2.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'outcome 2'],0)<10,
                                                                   paste0(' ', round(drugname$value[drugname$outcome == 'outcome 2'],0), '%'), 
                                                                   paste0(round(drugname$value[drugname$outcome == 'outcome 2'],0), '%')), 
               size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 2']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_label(x = outcome3.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'outcome 3'],0)<10,
                                                                   paste0(' ', round(drugname$value[drugname$outcome == 'outcome 3'],0), '%'), 
                                                                   paste0(round(drugname$value[drugname$outcome == 'outcome 3'],0), '%')), 
               size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 3']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_label(x = outcome4.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'outcome 4'],0)<10,
                                                                   paste0(' ', round(drugname$value[drugname$outcome == 'outcome 4'],0), '%'), 
                                                                   paste0(round(drugname$value[drugname$outcome == 'outcome 4'],0), '%')), 
               size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 4']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_label(x = outcome5.coord, y = eracle*0.84, label = ifelse(round(drugname$value[drugname$outcome == 'outcome 5'],0)<10,
                                                                   paste0(' ', round(drugname$value[drugname$outcome == 'outcome 5'],0), '%'), 
                                                                   paste0(round(drugname$value[drugname$outcome == 'outcome 5'],0), '%')), 
               size = 7, fill = '#787276', hjust = 0.5, color = '#f9fcff', family = 'Trebuchet MS', fontface = 'bold', label.padding = unit(0.6, 'lines'), 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 5']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
  
    # add placebo values to the placebo plot
    geom_point(x = outcome1.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 1']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    geom_text(x = outcome1.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 1']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 1']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    
    geom_point(x = outcome2.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 2']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    geom_text(x = outcome2.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 2']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 2']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    
    geom_point(x = outcome3.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 3']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    geom_text(x = outcome3.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 3']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 3']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    
    geom_point(x = outcome4.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 4']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    geom_text(x = outcome4.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 4']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 4']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    
    geom_point(x = outcome5.coord, y = eracle*0.84, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 22, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 5']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    geom_text(x = outcome5.coord, y = eracle*0.84, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 5']*100,0), '%'), 
              size = 7, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 5']), 0, ifelse(drugname$drug == 'placebo', 1, 0))) +
    
    # add placebo values in the intervention plots
    geom_point(x = outcome1.coord+(ctspace_end*2.8), y = eracle*0.855, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 1']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    geom_text(x = outcome1.coord+(ctspace_end*2.8), y = eracle*0.855, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 1']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 1']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_point(x = outcome2.coord+(ctspace_end*0.5), y = eracle*0.945, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 2']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    geom_text(x = outcome2.coord+(ctspace_end*0.5), y = eracle*0.945, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 2']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 2']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_point(x = outcome3.coord-(ctspace_end*2.5), y = eracle*0.9, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 3']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    geom_text(x = outcome3.coord-(ctspace_end*2.5), y = eracle*0.9, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 3']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 3']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_point(x = outcome4.coord-(ctspace_end*2.2), y = eracle*0.76, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 4']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    geom_text(x = outcome4.coord-(ctspace_end*2.2), y = eracle*0.76, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 4']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 4']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    geom_point(x = outcome5.coord+(ctspace_end*1.5), y = eracle*0.74, shape = 21, colour = '#787276', fill = '#DDF1FB', size = 11.7, 
               alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 5']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    geom_text(x = outcome5.coord+(ctspace_end*1.5), y = eracle*0.74, label = paste0(round(drugname$rate.pla[drugname$outcome == 'outcome 5']*100,0), '%'), 
              size = 3.5, colour = '#787276', family = 'Trebuchet MS', fontface = 'bold', 
              alpha = ifelse(is.na(drugname$value[drugname$outcome == 'outcome 5']), 0, ifelse(drugname$drug == 'placebo', 0, 1))) +
    
    # polar coordinate system
    coord_polar(clip = 'off', direction = 1) +
    
    # outcome1
    geom_text(x=outcome1.coord-(ctspace_end*4), y=outcome.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*3), y=outcome.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*2), y=outcome.height, label='m', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord-(ctspace_end*1), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)+ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord, y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*1), y=outcome.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*1), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*2), y=outcome.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*2), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*3), y=outcome.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*3), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome1.coord+(ctspace_end*4), y=outcome.height, label='n', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome1.coord+0.5)-ctangle_end*4), alpha = radial.int.alpha, family='Andale Mono') +
    
    # outcome2
    geom_text(x=outcome2.coord-(ctspace_end*6), y=outcome.reverse.height, label='y', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*6)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*5), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*4), y=outcome.reverse.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*4)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*3), y=outcome.reverse.height, label='l', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*3)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*2), y=outcome.reverse.height, label='i', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*2)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord-(ctspace_end*1), y=outcome.reverse.height, label='b', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)+ctangle_end*1)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord, y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5))+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*1), y=outcome.reverse.height, label='t', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*1)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*2), y=outcome.reverse.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*2)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*3), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*3)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*4), y=outcome.reverse.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*4)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*5), y=outcome.reverse.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome2.coord+(ctspace_end*6), y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome2.coord+0.5)-ctangle_end*6)+180, alpha = radial.int.alpha, family='Andale Mono') +
    
    # outcome3
    geom_text(x=outcome3.coord-(ctspace_end*3.5), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*3.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*2.5), y=outcome.reverse.height, label='h', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*2.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*1.5), y=outcome.reverse.height, label='c', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*1.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord-(ctspace_end*0.5), y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)+ctangle_end*0.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*0.5), y=outcome.reverse.height, label='d', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*0.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*1.5), y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*1.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*2.5), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*2.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome3.coord+(ctspace_end*3.5), y=outcome.reverse.height, label='h', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome3.coord+0.5)-ctangle_end*3.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    
    # outcome4
    geom_text(x=outcome4.coord-(ctspace_end*2.5), y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*2.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*1.5), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*1.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord-(ctspace_end*0.5), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)+ctangle_end*0.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*0.5), y=outcome.reverse.height, label='u', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*0.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*1.5), y=outcome.reverse.height, label='a', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*1.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome4.coord+(ctspace_end*2.5), y=outcome.reverse.height, label='n', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome4.coord+0.5)-ctangle_end*2.5)+180, alpha = radial.int.alpha, family='Andale Mono') +
    
    # outcome5
    geom_text(x=outcome5.coord-(ctspace_end*3.5), y=outcome.reverse.height, label='r', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*3.5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*2.5), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*2.5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*1.5), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*1.5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord-(ctspace_end*0.5), y=outcome.reverse.height, label='p', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)+ctangle_end*0.5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*0.5), y=outcome.reverse.height, label='o', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*0.5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*1.5), y=outcome.reverse.height, label='n', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*1.5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*2.5), y=outcome.reverse.height, label='s', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*2.5), alpha = radial.int.alpha, family='Andale Mono') +
    geom_text(x=outcome5.coord+(ctspace_end*3.5), y=outcome.reverse.height, label='e', size = size.text.int, colour = radial.int.col, angle = ((360/piecount)*(piecount-outcome5.coord+0.5)-ctangle_end*3.5), alpha = radial.int.alpha, family='Andale Mono') +
    
    # labels
    geom_text(x = 0.5, y = 10, label = ' 10%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-10)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 20, label = ' 20%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-20)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 30, label = ' 30%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-30)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 40, label = ' 40%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-40)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 50, label = ' 50%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-50)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 60, label = ' 60%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-60)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 70, label = ' 70%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-70)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 80, label = ' 80%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-80)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    geom_text(x = 0.5, y = 90, label = ' 90%', size = 5.5, colour = '#002147', vjust = '1', hjust = 'left', alpha = ifelse((eracle-90)>=(linespace/2), ifelse(drugname$drug == 'placebo', 0.15, 0), 0)) +
    
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
  geom_text(aes(label = value), size = 6) +
  scale_fill_gradient2(low = '#00b050', mid = '#F5CF47', high = '#fa3545', na.value = '#DDF1FB', breaks = c(-2.575829, -1.959964, -1.644854, 0, 1.644854, 1.959964, 2.575829), limits = c(-3, 3), labels = c('p < 0.01', 'p = 0.05', 'p = 0.1', 'p = 1.00', 'p = 0.1', 'p = 0.05','p < 0.01'))+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 15)) +
  labs(x = '',y = '') +
  theme(legend.title = element_blank(),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),legend.position = 'left', legend.text = element_text(size = 12)) +
  scale_x_discrete(position = 'top')

legend <- cowplot::get_legend(legendplot)
grid::grid.newpage()
png('legend.png', res = 600, width = 5.5, height = 5.5, units = 'in') 
grid::grid.draw(legend)
dev.off()

############## wrap
setwd('~/Desktop/V plot/graph8.0 - for V plot/figure2/')
multiple.plots <- function (x) {
  pattern <- paste0('.png$')
  Lfile <- list.files('vitruvian plots/', pattern = pattern)
  legend0 <- image_border(image_read('legend.png'), 'white', '70x120')
  image1 <- image_border(image_read(paste('vitruvian plots/', Lfile[1], sep = '')), 'white', '70x120')
  image2 <- image_border(image_read(paste('vitruvian plots/', Lfile[2], sep = '')), 'white', '70x120')
  image3 <- image_border(image_read(paste('vitruvian plots/', Lfile[3], sep = '')), 'white', '70x120')
  image4 <- image_border(image_read(paste('vitruvian plots/', Lfile[4], sep = '')), 'white', '70x120')
  image5 <- image_border(image_read(paste('vitruvian plots/', Lfile[5], sep = '')), 'white', '70x120')
  image6 <- image_border(image_read(paste('vitruvian plots/', Lfile[6], sep = '')), 'white', '70x120')
  image7 <- image_border(image_read(paste('vitruvian plots/', Lfile[7], sep = '')), 'white', '70x120')
  row1 <- c(legend0, image1)
  row2 <- c(image2, image3)
  row3 <- c(image4, image5)
  row4 <- c(image6, image7)
  row1.img <- image_append(row1)
  row2.img <- image_append(row2)
  row3.img <- image_append(row3)
  row4.img <- image_append(row4)
  all.plots <- c(row1.img, row2.img, row3.img, row4.img)
  all.plots.img <- image_append(all.plots, stack = TRUE)
  image_write(all.plots.img, path = paste0('synoptic chart.png'), format = 'png')
}
pangea <- lapply(unique(nrp$drug), multiple.plots)

some.plots <- function (x) {
  pattern <- paste0('.png$')
  Lfile <- list.files('vitruvian plots/', pattern = pattern)
  legend0 <- image_border(image_read('legend.png'), 'white', '70x120')
  image1 <- image_border(image_read(paste('vitruvian plots/', Lfile[1], sep = '')), 'white', '70x120')
  image2 <- image_border(image_read(paste('vitruvian plots/', Lfile[2], sep = '')), 'white', '70x120')
  image3 <- image_border(image_read(paste('vitruvian plots/', Lfile[3], sep = '')), 'white', '70x120')
  row1 <- c(legend0, image1)
  row2 <- c(image2, image3)
  row1.img <- image_append(row1)
  row2.img <- image_append(row2)
  all.plots <- c(row1.img, row2.img)
  all.plots.img <- image_append(all.plots, stack = TRUE)
  image_write(all.plots.img, path = paste0('figure2.jpeg'), format = 'jpeg', density = 1200, quality = 100)
}
pangea <- lapply(unique(nrp$drug), some.plots)
