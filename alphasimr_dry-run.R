## R code to simulate DNA sequences for a breeding program resembling that of 
## a natural terrapin population.

## Primary life history characteristics to simulate:
# Overlapping generations
# High adult survivorship, long lifespan
# Low egg/hatchling survivorship, low adult recruitment
# Historic population bottle-neck, stabilized population size but without recovery

## Clear the environment
rm(list = ls())

## Set working directory
setwd("~/Dropbox/terp_parentage/")

library(AlphaSimR)

## The following is slightly modified from Martin Johnsson's code found here: 
# (https://github.com/EdinbR/edinbr-talks/blob/master/2019-03-20/alphasimr_demo.R)

## Generate founder chromosomes
 
FOUNDERPOP <- runMacs(nInd = 20,
                      nChr = 1,
                      segSites = 50,
                      inbred = FALSE,
                      species = "GENERIC")

## Simulation parameters
 
SP <- SimParam$new(FOUNDERPOP)
#SP$addTraitA(nQtlPerChr = 100,
#                   mean = 100,
#                   var = 10)
SP$addSnpChip(nSnpPerChr = 50)
SP$setSexes("yes_sys")

## Return genetic map of all segregation sites
getGenMap(SP)

## Founding population

pop <- newPop(FOUNDERPOP,
              sP = SP)

## Breeding
# Next steps, add overlapping generations, add low/high attrition rates for adults/offspring
 
print("Breeding")
breeding <- vector(length = 11, mode = "list")
breeding[[1]] <- pop
 
for (i in 2:11) {
    print(i)
    sires <- selectInd(pop = breeding[[i - 1]],
                       nInd = 10,
                       sex = "M",
                       use = "rand",
                       sP = SP)
 
    dams <- selectInd(pop = breeding[[i - 1]],
                      nInd = 10,
                      sex = "F",
                      use = "rand",
                      sP = SP)
 
    breeding[[i]] <- randCross2(males = sires,
                                females = dams,
                                nCrosses = 20,
                                nProgeny = 3,
                                simParam = SP)
}

## Return pedigree for generation 1
getPed(breeding[[1]])

## Return SNP genotypes for generation 1
pullSnpGeno(breeding[[1]], snpChip = 1, chr = NULL, asRaw = FALSE, simParam = SP)

## Return SNP haplotypes for generation 1
pullSnpHaplo(breeding[[1]], simParam=SP)

## Next steps: Combine genetic map with variants, then convert to vcf to act as input for WGS short-read sequencer?
# Look into simuG-- can potentially simulate terrapin-like vcfs at the front-end, and add AlphaSimR-simulated variants on the back-end