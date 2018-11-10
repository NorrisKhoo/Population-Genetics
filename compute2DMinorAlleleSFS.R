# compute2DMinorAlleleSFS.R
# Norris Khoo 11/9/2018

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option(c("--genotype"), type = "character", default = '', help = "VCF file"),
  make_option(c("--target"), type = "character", default = '', help = "BED file of target sites"),
  make_option(c("--name1"), type = "character", default = '', help = "File containing names of sample individuals for species 1"),
  make_option(c("--name2"), type = "character", default = '', help = "File containing names of sample individuals for species 2"),
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

genotypeFilename = 'AW_selectVariant_chr38.vcf' # argument$genotype
genotypeInput = readGT(genotypeFilename)

variantRegex = as.numeric(gsub(pattern = "chr\\d+:(\\d+)_.*", replacement = "\\1", rownames(genotypeInput))) - 1
variantIRange = IRanges(start = variantRegex, end = variantRegex)

targetFilename = 'AW_callableNeutral_chr38.bed' # argument$target
targetInput = read.table(targetFilename, sep = '\t')
targetIRange = IRanges(start = targetInput[,2], end = targetInput[,3] - 1)
targetTotal = sum(targetIRange@width)

species1Filename = 'AW_name1.txt' # argument$name1
species1Input = as.vector(read.table(species1Filename, sep = '\t')$V1)

species2Filename = 'AW_name2.txt' # argument$name2
species2Input = as.vector(read.table(species2Filename, sep = '\t')$V1)

usableVariant = intersect(variantIRange, targetIRange)
disjointUsableVariant = unlist(tile(usableVariant, width = 1))
filterInIndex = variantRegex %in% disjointUsableVariant@start
genotypeFilter = genotypeInput[filterInIndex,]

genotypeCallable = filterCallable(genotypeFilter)
minorAllele = determineMinorAllele(genotypeCallable)

genotypeSpecies1 = genotypeCallable[, species1Input]
genotypeSpecies2 = genotypeCallable[, species2Input]

countVectorSpecies1 = computeSFS(genotypeSpecies1, minorAllele)
numberIndividualSpecies1 = length(genotypeSpecies1[1,])

countVectorSpecies2 = computeSFS(genotypeSpecies2, minorAllele)
numberIndividualSpecies2 = length(genotypeSpecies2[1,])

jointCountSFS = matrix(0, ncol = 2 * numberIndividualSpecies1 + 1, nrow = 2 * numberIndividualSpecies2 + 1)
for(i in seq(1, length(countVectorSpecies1))){
  column = countVectorSpecies1[i] + 1
  row = countVectorSpecies2[i] + 1
  jointCountSFS[row, column] = jointCountSFS[row, column] + 1
}
variantTotal = sum(jointCountSFS) - jointCountSFS[1, 1]
jointCountSFS[1, 1] = targetTotal - variantTotal

columnSpecies1 = rep('', 2 * numberIndividualSpecies1 - 1)
for(i in seq(0, 2 * numberIndividualSpecies1)){
  columnSpecies1[i + 1] = paste('d1_', toString(i), sep = '')
}

rowSpecies2 = rep('', 2 * numberIndividualSpecies2 - 1)
for(i in seq(0, 2 * numberIndividualSpecies2)){
  rowSpecies2[i + 1] = paste('d0_', toString(i), sep = '')
}

colnames(jointCountSFS) = columnSpecies1
rownames(jointCountSFS) = rowSpecies2

outputFilename = 'jointCountSFS.txt' # argument$output
write.table(jointCountSFS, file = outputFilename, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
