library(optparse)
library(VariantAnnotation)
library(data.table)
library(IRanges)

options(scipen = 999)
options(digits = 13)

filterCallable <- function(genotype){
  notCallable = rep(TRUE, length(genotype[,1]))
  
  for(i in seq(1, length(notCallable))){
    notCallable[i] = !any(genotype[i,] == '.')
  }
  
  return(genotype[notCallable,])
}

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

genotypeFilename = 'AW_selectVariant_chr38.vcf'
genotypeInput = readGT(genotypeFilename)

variantRegex = as.numeric(gsub(pattern = "chr\\d+:(\\d+)_.*", replacement = "\\1", rownames(genotypeInput))) - 1
variantIRange = IRanges(start = variantInput, end = variantInput)

targetFilename = 'AW_callableNeutral_chr38.bed'
targetInput = read.table(targetFilename, sep = '\t')
targetIRange = IRanges(start = targetInput[,2], end = targetInput[,3] - 1)

individualFilename = 'AW_name.txt'
individualInput = read.table(individualFilename, sep = '\t')$V1

usableVariant = intersect(variantIRange, targetIRange)
disjointUsableVariant = unlist(tile(usableVariant, width = 1))
filterInIndex = variantRegex %in% disjointUsableVariant@start
genotypeFilter = genotypeInput[filterInIndex, individualInput]

genotypeCallable = filterCallable(genotypeFilter)
AF = computeAF(genotypeCallable)

variantRegex = as.numeric(gsub(pattern = "chr\\d+:(\\d+)_.*", replacement = "\\1", rownames(genotypeCallable))) - 1

chromosomeFilename = 'canFam3.1Length_chr38.txt'
chromosomeLength = read.table(chromosomeFilename, sep = '\t')[1,2]

windowSize = 100000

piPerSite = computePiPerSite(AF, genotypeCallable, windowSize, chromosomeLength, targetIRange, variantRegex)

outputFilename = 'AW_selectVariant_diversity_chr38.txt'

write.table(piPerSite, file = outputFilename, quote = FALSE, row.names = FALSE, col.names=FALSE, sep = '\t')

