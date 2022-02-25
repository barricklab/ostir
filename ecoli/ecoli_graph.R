library(ggplot2)
library(tidyverse)

X = read.csv("compiled_data.csv")

# Use exact matches to gene starts
X = X %>% filter(X$expression > 0.001)

X$gene_start = (!is.na(X$offset_base) & (X$offset_base==0))
X$log10expression = log10(X$expression)

X$gene_associated = (!is.na(X$offset_base))

my.theme = theme_bw()

#Plot distributions
my.bin.width = 0.1

#Plot on same y-scale
ggplot(X, aes(x=log10expression, fill=gene_start)) +
  coord_cartesian(xlim=c(-1, 5), ylim=c(0, 16000), expand=F) +
  geom_histogram(binwidth=my.bin.width, alpha=.5, position="identity") +
  my.theme

#Plot showing same scale bu
ggplot(X, aes(x=log10expression, fill=gene_start)) +
  coord_cartesian(xlim=c(-1, 5), ylim=c(0, 500)) +
  geom_histogram(binwidth=my.bin.width, alpha=.5, position="identity") +
  my.theme

#
ggplot(X, aes(x=log10expression, fill=gene_start)) +
  geom_histogram(binwidth=my.bin.width, alpha=.5, position="identity") +
  facet_grid(vars(gene_start), scales="free_y") +
  scale_x_continuous(limits=c(-2,5), expand=expansion(mult=c(0,0))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  my.theme

ggplot(X, aes(x=log10expression, fill=gene_associated)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  facet_grid(vars(gene_associated), scales="free_y")

ggplot(X, aes(x=log10expression, fill=gene_start)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  scale_y_log10() +
  facet_grid(vars(gene_start), scales="free_y")
