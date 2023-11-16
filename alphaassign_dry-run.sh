# Downloading and testing AlphaAssign

# Download the AlphaAssign zip folder from https://github.com/AlphaGenes/AlphaAssign

# Install AlphaAssign using the included wheel file
pip install AlphaAssign-1.1.0-py3-none-any.whl

## Updating source code to fix error
## Open the following script
/usr/local/lib/python3.11/site-packages/tinyassign/tinyhouse/Pedigree.py

# Delete line 3
# Replace with the following:

#from numba import jit, float32, int32, int64, optional
#try:
#    from numba.experimental import jitclass
#except ModuleNotFoundError:
#    from numba import jitclass

# Test AlphaAssign by requesting man page
AlphaAssign

# Use Rscript_alphasimr.sh to extract simulated pedigree genotype data (parent-offspring, 2 gens total)
./Rscript_alphasimr.sh

# Make an output directory
mkdir outputs

# Run parentage analysis
AlphaAssign -genotypes genotypes.txt -potentialsires sires.txt -out output/out_genotypes

