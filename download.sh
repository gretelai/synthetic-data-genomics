# Following links are from https://datadryad.org/stash/dataset/doi:10.5061/dryad.2rs41
WORKING_DIR=mice_data_set

# Download experiment datasets
mkdir -p $WORKING_DIR/data $WORKING_DIR/out $WORKING_DIR/gemma/output
wget -O $WORKING_DIR/data/geno.txt.gz https://datadryad.org/stash/downloads/file_stream/4344
gunzip -f $WORKING_DIR/data/geno.txt.gz
wget -O $WORKING_DIR/data/map.txt https://datadryad.org/stash/downloads/file_stream/4342
wget -O $WORKING_DIR/data/pheno.csv https://datadryad.org/stash/downloads/file_stream/4340

# Download GEMMA
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz -P $WORKING_DIR/gemma/
gunzip -f $WORKING_DIR/gemma/gemma-0.98.4-linux-static-AMD64.gz
chmod u+x $WORKING_DIR/gemma/gemma-0.98.4-linux-static-AMD64
