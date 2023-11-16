#!/usr/bin/env Rscript

library(AlphaSimR)

## The following is slightly modified from Martin Johnsson's code found here: 
# (https://github.com/EdinbR/edinbr-talks/blob/master/2019-03-20/alphasimr_demo.R)

## Generate founder chromosomes
 
FOUNDERPOP <- runMacs(nInd = 20,
                      nChr = 1,
                      segSites = 250,
                      inbred = FALSE,
                      species = "GENERIC")

## Simulation parameters
 
SP <- SimParam$new(FOUNDERPOP)
#SP$addTraitA(nQtlPerChr = 100,
#                   mean = 100,
#                   var = 10)
SP$addSnpChip(nSnpPerChr = 250)
SP$setSexes("yes_sys")

## Return genetic map of all segregation sites
#getGenMap(SP)

## Founding population

pop <- newPop(FOUNDERPOP,
              sP = SP)

## Breeding
# Next steps, add overlapping generations, add low/high attrition rates for adults/offspring
 
print("Breeding")
breeding <- vector(length = 3, mode = "list")
breeding[[1]] <- pop
 
for (i in 2:3) {
    print(i)
    sires <- selectInd(pop = breeding[[i - 1]],
                       nInd = 3,
                       sex = "M",
                       use = "rand",
                       sP = SP)
 
    dams <- selectInd(pop = breeding[[i - 1]],
                      nInd = 3,
                      sex = "F",
                      use = "rand",
                      sP = SP)
 
    breeding[[i]] <- randCross2(males = sires,
                                females = dams,
                                nCrosses = 10,
                                nProgeny = 3,
                                simParam = SP)
}

# Combine parent and offspring generation genotypes into one matrix
Gen2 <- pullSnpGeno(breeding[[2]], snpChip = 1, chr = NULL, asRaw = FALSE, simParam = SP)
Gen3 <- pullSnpGeno(breeding[[3]], snpChip = 1, chr = NULL, asRaw = FALSE, simParam = SP)
Gen <- rbind(Gen2,Gen3)
Gen <- cbind(as.numeric(rownames(Gen)),as.matrix(Gen)) #Add individual id to each genotype
write.table(Gen, "genotype.txt", col.names=FALSE, row.names = FALSE)

# Create potential sires file
# First column is offspring ids, followed by a list of potential sires
sires <- data.frame(as.numeric(rownames(Gen3)), matrix(rep(as.numeric(breeding[[2]]@id[which(breeding[[2]]@sex == "M")]),each=1,times=30), 30, 15, byrow = TRUE))
write.table(sires, "sires.txt", col.names=FALSE, row.names = FALSE)

write.table(getPed(breeding[[3]]), "Ped.txt")
cat("Ped.txt","\n")