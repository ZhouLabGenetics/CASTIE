#!/bin/bash
#
# SAIGE-QTL Dynamic Installation Script for Conda Environment
# Updated to match the current Dockerfile approach
#

### RUN THIS SCRIPT IN A DIRECTORY IN WHICH YOU'RE HAPPY TO INSTALL SAIGE-QTL DYNAMIC

echo "Installing SAIGE-QTL Dynamic in conda environment..."

# activate SAIGE env (use "source activate" if "conda activate" fails)
conda activate RSAIGE 2>/dev/null
if [ $? -ne 0 ]
then
   ACTPATH=`which activate 2>/dev/null` || ACTPATH=`whereis activate 2>/dev/null | awk '{print $2}'`
   DEACTPATH=`which deactivate 2>/dev/null` || DEACTPATH=`whereis deactivate 2>/dev/null | awk '{print $2}'`
   source $ACTPATH RSAIGE || { echo -e "Commands \"conda activate RSAIGE\" and \"source activate RSAIGE\" failed."; echo -e "Please check whether a) conda is installed and b) you have prepared the RSAIGE conda environment."; exit 1; }
fi

# set path for lib and compiler flags
FLAGPATH=`which python | sed 's|/bin/python$||'`

# set lib and compiler flags using python path
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"

# Clone the SAIGE-QTL Dynamic repo (this repo, not the original SAIGE)
src_branch=main
repo_src_url=https://github.com/weizhou0/SAIGE-QTL-doc.git
echo "Cloning SAIGE-QTL Dynamic repository..."
git clone --depth 1 -b $src_branch $repo_src_url SAIGE_QTL_Dynamic || {
    echo "Failed to clone repository. Trying alternative branch..."
    git clone $repo_src_url SAIGE_QTL_Dynamic
}

cd SAIGE_QTL_Dynamic

# Install GitHub dependencies using remotes
echo "Installing GitHub packages..."
R -e "remotes::install_github('leeshawn/MetaSKAT', dependencies=TRUE)"
R -e "remotes::install_github('cysouw/qlcMatrix', dependencies=TRUE)"

# Run configure to install C++ dependencies
echo "Configuring C++ dependencies..."
./configure

# Install the R package
echo "Installing SAIGE-QTL Dynamic R package..."
R CMD INSTALL .

# Test installation
echo "Testing installation..."
R -e "library(SAIGEQTL); cat('SAIGE-QTL Dynamic installed successfully\n')" || {
    echo "Installation test failed, but package may still work"
}

# Set up environment variables for script access
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
[ ! -s $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh ] && echo '#!/bin/sh' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
[ ! -s $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh ] && echo '#!/bin/sh' > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

# Add scripts to PATH
echo -e "export PATH_OLD=\$PATH\nexport PATH=\$PATH:$PWD/extdata" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo -e "export PATH=\$PATH_OLD\nunset PATH_OLD" >> $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

echo "Installation completed! Scripts are available at:"
echo "  - step1_fitNULLGLMM_qtl.R"
echo "  - step2_tests_qtl.R"
echo "  - step3_gene_pvalue_qtl.R"
echo "  - makeGroupFile.R"

# deactivate the env
[ -z $DEACTPATH ] && conda deactivate || source $DEACTPATH
