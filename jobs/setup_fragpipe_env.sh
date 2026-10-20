#!/bin/bash
# =============================================================================
# setup_fragpipe_env.sh -- record of how the FragPipe environment was built
# on the UCT HPC login node (14-15 Oct 2026). Read it; run steps one by one.
# Thesis: Methods 3.6 and Appendix A.
# =============================================================================
set -e
# 1. FragPipe 24.0 (Linux zip, 677,845,803 bytes)
mkdir -p ~/software && cd ~/software
wget https://github.com/Nesvilab/FragPipe/releases/download/24.0/FragPipe-24.0-linux.zip
unzip -q FragPipe-24.0-linux.zip
# 2. MSFragger 4.4.1, IonQuant 1.11.20, diaTracer 2.2.1: academic licence
#    forms in a browser (msfragger-upgrader.nesvilab.org/{upgrader,ionquant,diatracer}),
#    copied to the cluster with scp, unzipped, then copied INTO FragPipe's tools/:
#    cp -r ~/fragpipe_tools/{MSFragger-4.4.1,IonQuant-1.11.20,diatracer-2.2.1} ~/software/fragpipe-24.0/tools/
# 3. Python 3.11 (FragPipe 24 accepts 3.9-3.11; the cluster default is 3.12)
/opt/exp_soft/miniconda3-py3.12/bin/conda create -p ~/envs/fragpipe python=3.11 -y
# 4. Spectral-library module and the modules FragPipe checks for
~/envs/fragpipe/bin/pip install FragPipe-SpecLib==0.1.58
~/envs/fragpipe/bin/pip install lxml               # required by FragPipe's check; not a SpecLib dependency
~/envs/fragpipe/bin/pip install "pandas>=2.2,<3"   # SpecLib 0.1.58 fails with pandas 3 (itertuples error)
~/envs/fragpipe/bin/pip install "scikit-learn>=1.3,<1.6"   # SpecLib 0.1.58 fails with sklearn >=1.6 (__sklearn_tags__ error)
# Result used: numpy 1.26.4, pandas 2.3.3, scikit-learn 1.5.2, lxml 6.1.3, pyopenms 3.5.0, FragPipe-SpecLib 0.1.58
