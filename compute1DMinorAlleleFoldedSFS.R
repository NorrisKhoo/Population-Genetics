suppressPackageStartupMessages(library(optparse))
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
    if (AF > 0.5){
      minorAllele[i] = 0
    } else{
      minorAllele[i] = 1
    }
  }
  
  return(minorAllele)
}

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

individualFilename = 'AW_name.txt' # argument$name
individualInput = read.table(individualFilename, sep = '\t')$V1

usableVariant = intersect(variantIRange, targetIRange)
disjointUsableVariant = unlist(tile(usableVariant, width = 1))
filterInIndex = variantRegex %in% disjointUsableVariant@start
genotypeFilter = genotypeInput[filterInIndex, individualInput]

genotypeCallable = filterCallable(genotypeFilter)
minorAllele = determineMinorAllele(genotypeCallable)

countVector = computeSFS(genotypeCallable, minorAllele)
numberIndividual = length(genotypeCallable[1,])

countVector = countVector[which(countVector != 0)]
countVector = countVector[which(countVector != numberIndividual)]

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

outputFilename = 'AW_countSFS_chr38.txt' # argument$output
write.table(foldedSFS, file = outputFilename, quote = FALSE, row.names = FALSE, col.names=FALSE, sep = '\t')

barplot(foldedSFS, names.arg = seq(1, length(foldedSFS)), ylab = 'Count', xlab = 'Frequency', main = 'Minor Allele Folded SFS', col = 'red')

