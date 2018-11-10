# computeDiversity.R
# Norris Khoo 11/3/2018
# This script calculates per window pi per site values given
# --genotype: VCF file (matrix of genotypes consisting of individuals vs variant sites)
# --target: BED file of target sites (region to consider when calculating per window pi per site values)
# --name: File containing names of sample individuals
# --chromosome: Size of entire chromosome in bp (integer)
# --window: Size of each window in bp (integer)
# --output: Output filename (output is in the following format: windowStart windowEnd piWithinWindow numberTargetSitesWithinWindow piPerSite)
# Per window pi per site values provide insight into the distribution of genetic diversity within a chromosome

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option(c("--genotype"), type = "character", default = '', help = "VCF file"),
  make_option(c("--target"), type = "character", default = '', help = "BED file of target sites"),
  make_option(c("--name"), type = "character", default = '', help = "File containing names of sample individuals"),
  make_option(c("--chromosome"), type = "integer", default = 0, help = "Length (bp) of chromosome"),
  make_option(c("--window"), type = "integer", default = 0, help = "Length (bp) of window"),
  make_option(c("--output"), type = "character", default = '', help = "Output file"))
argument <- parse_args(OptionParser(option_list = options))

suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(IRanges))

options(scipen = 999)
options(digits = 13)

# Filters out variant sites containing missing genotype data for >=1 individual
# @param genotype: matrix of genotypes (individuals vs variant sites)
# @return genotype post filtering
filterCallable <- function(genotype){
  notCallable = rep(TRUE, length(genotype[,1]))

  for(i in seq(1, length(notCallable))){
    notCallable[i] = !any(genotype[i,] == '.')
  }

  return(genotype[notCallable,])
}

# Computes allele frequency amongst sample individuals for each variant site
# @param genotype: matrix of genotypes (individuals vs variant sites)
# @return vector of allele frequencies for each variant site
computeAF <- function(genotype){
  numberVariant = length(genotype[,1])
  numberIndividual = length(genotype[1,])

  numberAllele = 2.0 * numberIndividual
  vectorAF = rep(0, numberVariant)

  for(i in seq(1, numberVariant)){
    alternateAlleleCount = 0.0
    for(j in seq(1, numberIndividual)){
      if((genotype[i, j] == '0/1') | (genotype[i, j] == '1/0')){
        alternateAlleleCount = alternateAlleleCount + 1.0
      }
      if(genotype[i, j] == '1/1'){
        alternateAlleleCount = alternateAlleleCount + 2.0
      }
    }
    vectorAF[i] = alternateAlleleCount / numberAllele
  }

  return(vectorAF)
}

# Calculates per window pi per site (measure of genetic diversity)
# @param vectorAF: vector of allele frequencies for each variant site
# @param genotype: matrix of genotypes (individuals vs variant sites)
# @param windowSize: number of bp per window
# @param chromosomeLength: number of bp in entire chromosome
# @param outsideBound: IRange object containing regions outside of target sites
# @param variantLocation: vector of variant sites
# @return matrix of per window pi per site values in the following format:
# windowStart windowEnd piWithinWindow numberTargetSitesWithinWindow piPerSite
computePiPerSite <- function(vectorAF, genotype, windowSize, chromosomeLength, outsideBound, variantLocation){
  numberAllele = 2 * length(genotype[1,])
  numberWindow = ceiling(chromosomeLength / windowSize)
  column1 = rep(0, numberWindow)
  column2 = rep(0, numberWindow)
  column3 = rep(0, numberWindow)
  column4 = rep(0, numberWindow)
  column5 = rep(0, numberWindow)

  for(i in seq(1, numberWindow)){
    lowerBound = windowSize * (i - 1)
    upperBound = windowSize * i
    insideBound = seq(lowerBound, upperBound, by = 1)

    withinIndex = variantLocation %in% insideBound
    withinAF = vectorAF[withinIndex]

    column1[i] = lowerBound
    if(i != numberWindow){
      column2[i] = upperBound
    } else{
      column2[i] = chromosomeLength
    }

    windowIRange = IRanges(start = lowerBound, end = upperBound - 1)
    overlapIRange = intersect(windowIRange, outsideBound)
    overlapSize = sum(overlapIRange@width)
    column4[i] = overlapSize

    if(length(withinAF) == 0){
      column3[i] = NA
      column5[i] = NA
    } else{
      heterozygousEach = 2.0 * withinAF * (1.0 - withinAF)
      heterozygousAll = sum(heterozygousEach)
      pi = (numberAllele / (numberAllele - 1.0)) * heterozygousAll
      column3[i] = pi
      column5[i] = pi / overlapSize
    }
  }

  output = cbind(column1, column2, column3, column4, column5)
  return(output)
}

genotypeFilename = argument$genotype
genotypeInput = readGT(genotypeFilename)

variantRegex = as.numeric(gsub(pattern = "chr\\d+:(\\d+)_.*", replacement = "\\1", rownames(genotypeInput))) - 1
variantIRange = IRanges(start = variantRegex, end = variantRegex)

targetFilename = argument$target
targetInput = read.table(targetFilename, sep = '\t')
targetIRange = IRanges(start = targetInput[,2], end = targetInput[,3] - 1)

individualFilename = argument$name
individualInput = as.vector(read.table(individualFilename, sep = '\t')$V1)

usableVariant = intersect(variantIRange, targetIRange)
disjointUsableVariant = unlist(tile(usableVariant, width = 1))
filterInIndex = variantRegex %in% disjointUsableVariant@start
genotypeFilter = genotypeInput[filterInIndex, individualInput]

genotypeCallable = filterCallable(genotypeFilter)
AF = computeAF(genotypeCallable)

variantRegex = as.numeric(gsub(pattern = "chr\\d+:(\\d+)_.*", replacement = "\\1", rownames(genotypeCallable))) - 1

chromosomeLength = argument$chromosome
windowSize = argument$window

piPerSite = computePiPerSite(AF, genotypeCallable, windowSize, chromosomeLength, targetIRange, variantRegex)

outputFilename = argument$output
write.table(piPerSite, file = outputFilename, quote = FALSE, row.names = FALSE, col.names=FALSE, sep = '\t')

