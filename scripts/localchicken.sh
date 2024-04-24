SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42" # Random seed
MU="2e-9" # Mutation rate, 2e-9 from previous publication in chicken genome analysis
URTR="40" # upperRhoThetaRatio, note higher value cause slow simulation
GenTime="1" # Assumed generation time 1 year from previous pub
DIR="./local_out/" # Output need to be in the same directory
VCF="/home6/ycheng29/BIO624/data/local_chicken_autosomes_biallelic.vcf" # Can specify input vcf absolute path (can be different), default unphased GT
GENOME="/home6/ycheng29/BIO624/data/galGal4_filtered.bed" # Can specify input bed absolute path
MaxSites="4000" # default 1750 generates small window size 49kb; we want larger window size approx 500kb here
# MASK="./accessibility_mask.bed" # Optional, A BED-formatted accessibility mask (with non-overlapping ascending windows)
# In simulation module, add assumedGenTime, set to 1
# nTrain, nVali, nTest examples are all recommended to keep as default
# If any infered demographic history(stairwayplot_v1, SMC++, and MSMC) available, add demographicHistory with path
# It is possible to treat hemizygous chromosomes as "diploids with missing data" using the --forceDiploid option, however this is not recommended.

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --upperRhoThetaRatio ${URTR} \
    --assumedGenTime ${GenTime} \
    --maxSites ${MaxSites} \
    --nTrain 13000 \
    --nVali 2000 \
    --nTest 100 \
    --seed ${SEED}

# Train network
${TRAIN} \
    --projectDir ${DIR} \
    --nEpochs 20 \
    --nValSteps 10 \
    --seed ${SEED}

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED}

# Try basic prediction module first: took 11 mins for URTR=1

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${DIR} \
    --nSlice 10 \
    --nReps 10 \
    --seed ${SEED}
