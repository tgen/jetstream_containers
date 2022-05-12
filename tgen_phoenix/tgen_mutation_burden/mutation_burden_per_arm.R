
## Copyright (c) 2020 Translational Genomics Research Institute
##
## This software may be modified and distributed under the terms
## of the MIT license.  See the LICENSE file for details.
##
## Major Contributors: Christophe Legendre
require(optparse)
require(ggplot2)

klc11 <<- c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")
klc10 <<- c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")
klc9 <<- c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")
klc8 <<- c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e")
klc7 <<- c("#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e")
klc6 <<- c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e")
klc5 <<- c("#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e")
klc4 <<- c("#a6611a", "#dfc27d", "#80cdc1", "#018571")
klc3 <<- c("#dfc27d","#80cdc1","#018571")
klc2 <<- c("#dfc27d","#018571")
klc1 <<- c("#dfc27d")


#@@@@@@@@@@@
# GET INPUTS
#@@@@@@@@@@@
option_list = list(
  make_option(c("", "--helpme"), action="store_true", default=NA, type='character', help="HELP"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="add more verbosity to the log file"),
  make_option(c("-o", "--outfilename"), action="store", default=NA, type='character', help="give a specific name tout the output files (both files, stats and plot"),
  make_option(c("-i", "--infile"), action="store", default=NA, type='character', help="file with tgen_mutation_burden results from different coordinates accross the genome; see dir examples in install directory"),
  make_option(c("-c", "--noPrefixChr"), action="store_true", default=FALSE, help="add more verbosity to the log file"),
  make_option(c("-s", "--addPctDbsnp"), action="store_true", default=FALSE, help="will add the dbsnp percentages to the plot if present in column #3 of INFILE")
)

opt = parse_args(OptionParser(option_list=option_list))
print(paste0("length opt is : ",length(opt)))

infile=opt$infile
outfilename=opt$outfilename
addPctDbsnp=opt$addPctDbsnp
contig_without_chr_prefix=opt$noPrefixChr
verbose=opt$verbose
help=opt$help

if(help){stop()}
if(is.na(outfilename)){ outfilename=infile }

if(verbose){
  print(paste("optparse version = ", packageVersion("optparse")) )
  print(paste("ggplot version = ", packageVersion("ggplot2")) )
  print("recap inputs ...")
  for ( P in c(infile, outfilename, addPctDbsnp) ){ print(P)}
  }


if(contig_without_chr_prefix){
  list_contig_ordered=c(seq(1,22), "X", "Y")
} else {
  list_contig_ordered=paste0("chr", c(seq(1,22), "X", "Y"))
}
## reading input file
d <- read.table(infile, sep="\t")
## forcing the type for some specific columns
d$V1 <- factor(d$V1, levels=list_contig_ordered)
d$V3 <- as.integer(d$V3)
d$V4 <- as.integer(d$V4)
d$V5 <- as.numeric(d$V5)
d$V7 <- as.numeric(d$V7)


if(addPctDbsnp){
  ## preProcess data if dbSNP values available
  colnames(d) <- c("chrom","region","mutationsCount","coverage","mutBurdenPerMB","region_name","pct_dbsnp")
  test <- is.na(d$pct_dbsnp)
  if(any(test)) {
  d[is.na(d$pct_dbsnp), ]$pct_dbsnp <- 0
  }
} else {
  d <- d[,c(1:6)]
  colnames(d) <- c("chrom","region","mutationsCount","coverage","mutBurdenPerMB","region_name")
}
if(verbose) { head(d) ; str(d) ; dim(d) }

## check if more than 11 regions:
nregions <- length(unique(d$region_name))
if( nregions > 11 ){
  ## creating a new variable for colors based on the exact number of regions
  assign(paste0("klc",nregions), rep(klc11, ceiling(nregions/11))[1:nregions])
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## calculate basic Stats accross regions of interest given one region per line in the input file
sumMutationCount <- sum(na.omit(d$mutationsCount))
sumCoverage <- sum(na.omit(d$coverage))
mutBurden <- sumMutationCount/sumCoverage*1000000
median <- median(na.omit(d$mutBurdenPerMB))
mean <- mean(na.omit(d$mutBurdenPerMB))
min <- min(na.omit(d$mutBurdenPerMB))
max <- max(na.omit(d$mutBurdenPerMB))

## print some basic stats
if(verbose) {print(paste("overall mutBurden per MB = ", mutBurden))
  print(paste("mean = ", mean))
  print(paste("median = ", median))
  print(paste("min = ", min))
  print(paste("max = ", max))
  print(paste("mutationsCount = ",sumMutationCount))
  print(paste("coverage = ",sumCoverage))
}
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


## save the basic stats
names=c("mutBurdenPerMB","mean","median","min","max","totMutationCount","totCoverage")
values=c(mutBurden,mean,median,min,max,sumMutationCount, sumCoverage)
df <- data.frame(calcul=names, value=values)

options(scipen = 100000000)
write.table(x = df, file = paste0(outfilename,".basic_stats.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
write.table(x = t(df), file = paste0(outfilename,".basic_stats.transposed.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## MAKE PLOTS
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if(addPctDbsnp && sum(d$pct_dbsnp)!= 0) {
  ## column should be named as follow:
  g <- ggplot(d, aes(x=region_name, y=mutBurdenPerMB)) +
    geom_bar(data = d, stat = "identity", aes(fill=d$pct_dbsnp) ) +
    geom_text(aes(label = d$mutationsCount), position = position_dodge(width=0.9), size = 2.5, vjust=+1 )  +
    facet_wrap(~chrom, ncol = 6) +
    ylab(label = "Mutations / MB") +
    xlab(label =  "Chromosome Region Name") +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(fill='% dbSNP') +
    theme(axis.title = element_text(size=13, face="bold"),
          axis.text = element_text(size=11),
          strip.text = element_text(size=11, face="bold"),
          panel.background = element_rect(fill = "grey96"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0))
    )

  ggsave(filename = paste(outfilename,".ggplot.with_pct_dbsnp.png",sep=""), plot = g, width = 10, height = 8)

  ## same plot as above but with values converted to log10
  g <- ggplot(d, aes(x=region_name, y=log10(mutBurdenPerMB+1))) +
    geom_bar(data = d, stat = "identity", aes(fill=d$pct_dbsnp) ) +
    geom_text(aes(label = d$mutationsCount), position = position_dodge(width=0.9), size = 2.5, vjust=+1 )  +
    facet_wrap(~chrom, ncol = 6) +
    ylab(label = "log10(Mutations) / MB") +
    xlab(label =  "Chromosome Region Name") +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(fill='% dbSNP') +
    theme(axis.title = element_text(size=13, face="bold"),
          axis.text = element_text(size=11),
          strip.text = element_text(size=11, face="bold"),
          panel.background = element_rect(fill = "grey96"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0))
    )

  ggsave(filename = paste(outfilename,".log10.ggplot.with_pct_dbsnp.png",sep=""), plot = g, width = 10, height = 8)

} else {

  ##color = c(get(paste0("klc",length(unique(d$region_name)))) ), fill=c(get(paste0("klc",length(unique(d$region_name)))) )
  ## without dbsnp information
  g <- ggplot(d, aes(x=region_name, y=mutBurdenPerMB, fill=region_name)) +
    geom_bar(data = d, stat = "identity"  ) +
    geom_text(aes(label = d$mutationsCount), position = position_dodge(width=0.9), size = 2.5, vjust=+1 )  +
    facet_wrap(~chrom, ncol = 6) +
    scale_fill_manual(values = c(get(paste0("klc",length(unique(d$region_name)))) ) ) +
    ylab(label = "Mutations / MB") +
    xlab(label =  "Chromosome Region Name") +
    theme(axis.title = element_text(size=13, face="bold"),
          axis.text = element_text(size=11),
          strip.text = element_text(size=11, face="bold"),
          panel.background = element_rect(fill = "grey96"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0))
    )

  ggsave(filename = paste(outfilename,".ggplot.png",sep=""), plot = g, width = 10, height = 8)

  ## same plot but with values converted to log10
  g <- ggplot(d, aes(x=region_name, y=log10(mutBurdenPerMB+1), fill=region_name)) +
    geom_bar(data = d, stat = "identity"  ) +
    geom_text(aes(label = d$mutationsCount), position = position_dodge(width=0.9), size = 2.5, vjust=+1 )  +
    facet_wrap(~chrom, ncol = 6) +
    scale_fill_manual(values = c(get(paste0("klc",length(unique(d$region_name)))) ) ) +
    ylab(label = "log10(Mutations) / MB") +
    xlab(label =  "Chromosome Region Name") +
    theme(axis.title = element_text(size=13, face="bold"),
          axis.text = element_text(size=11),
          strip.text = element_text(size=11, face="bold"),
          panel.background = element_rect(fill = "grey96"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0))
    )

  ggsave(filename = paste(outfilename,".log10.ggplot.png",sep=""), plot = g, width = 10, height = 8)

}

