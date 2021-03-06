# compute1DMinorAlleleFoldedSFS.R
# Norris Khoo 11/9/2018
# This script calculates the minor allele folded SFS for 1 population given
# --genotype: VCF file (matrix of genotypes consisting of individuals vs variant sites)
# --target: BED file of target sites (region to consider when counting number of alternate alleles per SNP)
# --name: File containing names of sample individuals
# --output: Output filename 
# The frequency of minor alleles provides insight into the demographic history of a population.
# The minor allele folded SFS can be input into programs such as fastsimcoal and dadi to generate a
# demographic model explaining a population's history

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option(c("--genotype"), type = "character", default = '', help = "VCF file"),
  make_option(c("--target"), type = "character", default = '', help = "BED file of target sites"),
  make_option(c("--name"), type = "character", default = '', help = "File containing names of sample individuals"),
  make_option(c("--output"), type = "character", default = '', help = "Output file"))
argument <- parse_args(OptionParser(option_list = options))

suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(plyr))

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

# Iterates through variant sites and determines which allele, 0 or 1, is the minor allele
# @param genotype: matrix of genotypes (individuals vs variant sites)
# @return vector of 0s and 1s denoting the minor allele
determineMinorAllele <- function(genotype){
  numberVariant = length(genotype[,1])
  numberIndividual = length(genotype[1,])

  numberAllele = 2 * numberIndividual
  minorAllele = rep(0, numberVariant)
  
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
    AF = alternateAlleleCount / numberAllele
    # If the allele frequency of both alleles is 0.5, denote 1 as minor allele
    if (AF > 0.5){
      minorAllele[i] = 0
    } else{
      minorAllele[i] = 1
    }
  }
  
  return(minorAllele)
}

# Iterates through variant sites and counts the number of minor alleles for each
# @param genotype: matrix of genotypes (individuals vs variant sites)
# @param minorAllele: vector of 0s and 1s denoting the minor allele
# @return vector of minor allele counts for each variant site
computeSFS <- function(genotype, minorAllele){
  numberVariant = length(genotype[,1])
  numberIndividual = length(genotype[1,])

  countVariant = rep(0, numberVariant)
  
  for (i in seq(1, numberVariant)){
    alleleCount = 0
    if (minorAllele[i] == 0){
      for (j in seq(1, numberIndividual)){
        if ((genotype[i, j] == '0/1') | (genotype[i, j] == '1/0')){
          alleleCount = alleleCount + 1
        }
        if (genotype[i, j] == '0/0'){
          alleleCount = alleleCount + 2
        }
      }
    } else{
      for (j in seq(1, numberIndividual)){
        if ((genotype[i, j] == '0/1') | (genotype[i, j] == '1/0')){
          alleleCount = alleleCount + 1
        }
        if (genotype[i, j] == '1/1'){
          alleleCount = alleleCount + 2
        }
      }
    }
    countVariant[i] = alleleCount
  }
  
  return(countVariant)
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
minorAllele = determineMinorAllele(genotypeCallable)

countVector = computeSFS(genotypeCallable, minorAllele)
numberIndividual = length(genotypeCallable[1,])

countVector = countVector[which(countVector != 0)]
countVector = countVector[which(countVector != 2 * numberIndividual)]

countSFS = count(countVector)$freq

if(numberIndividual %% 2 == 1){
  foldedSFS = rep(0, length(countSFS) / 2)
  for(i in seq(1, length(foldedSFS))){
    foldedSFS[i] = countSFS[i] + countSFS[length(countSFS) + 1 - i]
  }
} else{
  foldedSFS = rep(0, (length(countSFS) - 1) / 2)
  for(i in seq(1, length(foldedSFS))){
    foldedSFS[i] = countSFS[i] + countSFS[length(countSFS) + 1 - i]
  }
  foldedSFS[length(foldedSFS) + 1] = countSFS[length(foldedSFS) + 1]
}

outputFilename = argument$output
write.table(foldedSFS, file = outputFilename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
