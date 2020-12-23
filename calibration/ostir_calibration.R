#!/usr/bin/env Rscript

## This script takes an input CSV file that has delta G values **without a spacing term included**.
##
## First, it uses all sets of measurements with the same RBS distances to fit Beta and LogK.
## Then, it fits the coefficients for calculating the delta G of spacing from the RBS distance
## These values are then put into the OSTIR code to enable it to estimate translaton initiation rates.
##
## Parameterization of the model uses experimental data from:
##
## Salis HM, Mirsky EA, Voigt CA. 2009. Automated design of synthetic ribosome binding sites to control
## protein expression. Nat. Biotechnol. 27:946–950.
##
## The fitting follows the procedure from the paper, except:
##
## (1) The x-offset in the sigmoidal spacing term for the "compressed" equation is fit, rather than
##     being fixed at a value of 2. We also allow for a deltaG offset term added to all entries.
##     This fit constant deltaG offset gets moved into LogK so the offset is zero in the end
##
## (2) We use all measurements from their paper simultaneously to fit all parameters, rather than
##    fitting Beta and LogK from some measurements and fitting the spacing model using only set 1.
##

library(ggplot2)
library(tidyverse)

input.csv.path = "input.csv"
output.prefix = "output"
theme_set(theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

## Load data and set up the data frame
X = read_csv(input.csv.path,comment = "#")

## tidy data by method used for computating deltaG
X = X %>% gather("method", "dG", 7:ncol(X))

## add/reformat columns
X$Set = as.factor(X$Set)
X$spacing.factor = as.factor(X$Spacing)
X$Fluorescence.Log.Average = log(X$Fluorescence.Average)

## STEP 1: Fit Beta and LogK from sets of sequences with the same spacing

method.of.interest = "OSTIR.ViennaRNA.Turner2004.No.Spacing.dG"

X.filtered = X %>% filter(method==method.of.interest)

## How many measurements do we have with different spacing values?
#ggplot(X.filtered, aes(x=Spacing)) + geom_histogram(binwidth=1)

# Fit the overall best Beta and logK values by only comparing those that have the same spacing
spacing.factor.fit = lm(Fluorescence.Log.Average~dG+spacing.factor+0, data=X.filtered)

#Now we know the slope and the intercept assuming the minimum is at 5 bases spacing 
LogK = as.numeric(coef(spacing.factor.fit)["spacing.factor5"])
Beta = -as.numeric(coef(spacing.factor.fit)["dG"])

spacing.factor.fit = lm(Fluorescence.Log.Average~dG+spacing.factor, data=X.filtered)


## STEP 1: Fit the spring model

X.filtered = X.filtered %>% mutate(spacing.deviation.dG = -(1/Beta) * (Fluorescence.Log.Average - LogK) - dG)

# s is spacing
s.opt = 5

spring_model <- function(s, stretched.c1, stretched.c2, compressed.c1, compressed.c2, compressed.c3, y.zero) {
  compressed.c4 = 3
  s.diff = s-s.opt
  stretched = stretched.c1 * (s.diff)^2 + stretched.c2 * s.diff
  compressed = compressed.c1 / (1 + exp(compressed.c2 * (s.diff+compressed.c3)) )^compressed.c4
  r = s
  r[s >= s.opt] = stretched[s >= s.opt]
  r[s < s.opt] = compressed[s < s.opt]
  r = r + y.zero
  return(r)
}



data.to.spring.fit = X.filtered

## Alternative is to only use the spacing set 1
#data.to.spring.fit = X.filtered %>% filter(Set==1)
#data.to.spring.fit$spacing.deviation.dG = data.to.spring.fit$spacing.deviation.dG - min(data.to.spring.fit$spacing.deviation.dG)

spring.fit <- nls(spacing.deviation.dG~spring_model(Spacing,stretched.c1, stretched.c2, compressed.c1, compressed.c2, compressed.c3, y.zero),
               data = data.to.spring.fit, 
               start = list(
                 stretched.c1 = 0.048,
                 stretched.c2 = 0.24,
                 compressed.c1 = 12.2,
                 compressed.c2 = 2.5,
                 compressed.c3 = 2,
                 y.zero = 0
               ), 
               control=nls.control(maxiter=1000, tol = 1e-6, warnOnly=TRUE),
               trace = FALSE)

stretched.c1 = as.numeric(coef(spring.fit)["stretched.c1"])
stretched.c2 = as.numeric(coef(spring.fit)["stretched.c2"])
compressed.c1 = as.numeric(coef(spring.fit)["compressed.c1"])
compressed.c2 = as.numeric(coef(spring.fit)["compressed.c2"])
compressed.c3 = as.numeric(coef(spring.fit)["compressed.c3"])
compressed.c4 = 3
y.zero = as.numeric(coef(spring.fit)["y.zero"])

## Now we include the additional offset that we found in LogK!
LogK = LogK + y.zero*-Beta

# Graph comparison of fit (Note: this removes the y-offset present during fitting)
data.to.spring.fit = data.to.spring.fit %>% mutate(spacing.deviation.dG = -(1/Beta) * (Fluorescence.Log.Average - LogK) - dG)

X.filtered$spacing.penalty.dG = spring_model(X.filtered$Spacing, stretched.c1, stretched.c2, compressed.c1, compressed.c2, compressed.c3, 0 )

line.for.spring.fit = data.frame(Spacing=0:20, Set="0")
line.for.spring.fit$spacing.deviation.dG = spring_model(line.for.spring.fit$Spacing, stretched.c1, stretched.c2, compressed.c1, compressed.c2, compressed.c3, 0 )

ggplot(data.to.spring.fit, aes(x=Spacing, y=spacing.deviation.dG)) +  
  geom_point(aes(color=Set)) +
  geom_line(data=line.for.spring.fit) + 
  xlab("RBS distance (nt)") +
  ylab("Spacing deltaG")

#ggsave(paste0(output.prefix, ".spacing_dG.pdf"))
ggsave(paste0(output.prefix, ".spacing_dG.png"), width=4, height=3.6)



# Calculate the total deltaG
X.filtered$total.dG = X.filtered$dG + spring_model(X.filtered$Spacing, stretched.c1, stretched.c2, compressed.c1, compressed.c2, compressed.c3, 0 )

# Now calculate predicted values according to the equation

## We could look at only a subset of the data here
X.stats = X.filtered

X.stats = X.stats %>% 
  mutate(
    Log.Predicted.Translation.Initiation.Rate = -Beta * total.dG + LogK,
    Log10.Predicted.Translation.Initiation.Rate = Log.Predicted.Translation.Initiation.Rate / log(10),
    Log.Measured.Translation.Initiation.Rate = Fluorescence.Log.Average,
    Log10.Measured.Translation.Initiation.Rate = Log.Measured.Translation.Initiation.Rate / log(10),
    Log2.Translation.Initiation.Rate.Ratio =  (Log.Measured.Translation.Initiation.Rate - Log.Predicted.Translation.Initiation.Rate) / log(2),
    Abs.Log2.Translation.Initiation.Rate.Ratio =  abs(Log2.Translation.Initiation.Rate.Ratio)
  )

#Graph total dG versus measured rate
ggplot(X.stats, aes(x=total.dG, y=Log.Measured.Translation.Initiation.Rate, color=Set)) + 
  geom_point() + 
  geom_abline(slope=-Beta, intercept=LogK, color="black", linetype="solid")
#ggsave(paste0(output.prefix, ".predicted_versus_measured.pdf"))
#ggsave(paste0(output.prefix, ".total_dG_versus_measured_rate.png"), width=4, height=3.6)

# Graph comparison of predicted and measured values, coloring by set
ggplot(X.stats, aes(x=Log10.Predicted.Translation.Initiation.Rate, y=Log10.Measured.Translation.Initiation.Rate, color=Set)) +  
  geom_abline(slope=1, intercept=0, color="black", linetype="solid") +
  geom_point() +
  scale_x_continuous(limits=c(-0.6, 5.6), breaks=0:5) + xlab("Log10 Predicted Rate") +
  scale_y_continuous(limits=c(-0.6, 5.6), breaks=0:5) + ylab("Log10 Measured Rate")
#ggsave(paste0(output.prefix, ".predicted_versus_measured.pdf"))
ggsave(paste0(output.prefix, ".predicted_rate_versus_measured_rate.png"), width=4, height=3.6)


## Graph histogram of the log2 fold error to be sure it is symmetrical
ggplot(X.stats, aes(x=Log2.Translation.Initiation.Rate.Ratio)) + 
  geom_histogram(binwidth=0.5) + xlab("Log2FoldError")
#ggsave(paste0(output.prefix, ".log2_fold_error.pdf"))
ggsave(paste0(output.prefix, ".log2_fold_error.png"), width=4, height=3.6)

## Graph the absolute log2 fold error CDF
ggplot(X.stats, aes(x=Abs.Log2.Translation.Initiation.Rate.Ratio)) + 
  stat_ecdf(geom = "step") + 
  scale_x_continuous(expand = c(0,0), breaks=0:10) + xlab("Log2 Fold Error") +
  scale_y_continuous(expand = c(0,0), breaks=0.1*(0:10)) + ylab("Fraction < Log2 Fold Error") #+
  #geom_vline(xintercept=1, color="red") +
  #geom_vline(xintercept=2, color="red") +
  #geom_vline(xintercept=log2(10), color="red")
#ggsave(paste0(output.prefix, ".log2_fold_error.pdf"))
ggsave(paste0(output.prefix, ".log2_fold_error.png"), width=4, height=3)


## Call out some values at certain set points

cat("Fraction of measurements with < 2-fold error:",
  nrow(X.stats %>% filter(Abs.Log2.Translation.Initiation.Rate.Ratio < 1))/nrow(X.stats),
  "\n"
)
cat("Fraction of measurements with < 4-fold error:",
    nrow(X.stats %>% filter(Abs.Log2.Translation.Initiation.Rate.Ratio < 2))/nrow(X.stats),
    "\n"
)
cat("Fraction of measurements with <10-fold error:",
    nrow(X.stats %>% filter(Abs.Log2.Translation.Initiation.Rate.Ratio < log2(10)))/nrow(X.stats),
    "\n"
)

## Output the fit parameters to a file

params = data.frame(
  names = c(
    "Beta",
    "RTeff",
    "LogK",
    "K",
    "spacing.stretched.c1",
    "spacing.stretched.c2",
    "spacing.compressed.c1",
    "spacing.compressed.c2",
    "spacing.compressed.c3",
    "spacing.compressed.c4"
  ),
  values = c(
    Beta,
    1/Beta,
    LogK, 
    exp(LogK), 
    stretched.c1, 
    stretched.c2, 
    compressed.c1, 
    compressed.c2, 
    compressed.c3, 
    compressed.c4
  )
)

write_csv(params, paste0(output.prefix, ".fit_parameters.csv"))

## Output the fit information
X.stats$Predicted.Translation.Initiation.Rate = exp(X.stats$Log.Predicted.Translation.Initiation.Rate)
write_csv(X.stats, paste0(output.prefix, ".fit_values.csv"))
