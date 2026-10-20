#!/bin/bash
# =============================================================================
# run_intein_analysis.sh -- intein analysis of the finished FragPipe search
# (Dataset 3). Run on the login node from the repository's scripts/ folder
# copy, or set S to where the scripts are. Light work: minutes.
# Thesis: Methods 3.6-3.7; Results 4.4.
# =============================================================================

# Set PROJ to the directory holding the FragPipe run, the search
# database and the results, for example:
#   export PROJ=/scratch/$USER/TB_rif_DIA
: "${PROJ:?set PROJ to your project directory}"
set -euo pipefail
S=${S:-$HOME/intein-splicing-mtb/scripts}
W=$PROJ/fragpipe
DB=$PROJ/db/h37rv_uniprot_plus_spliced.fasta
OUT=$W/intein_analysis
mkdir -p $OUT
PY=/opt/exp_soft/miniconda3-py3.12/bin/python3       # any Python 3.8+; standard library only
$PY $S/fragpipe_intein_peptides.py --fasta $DB --psm $W/results/psm.tsv \
    --out-prefix $OUT/rifDIA --exclude-run 10989 --pepxml "$W/results/interact-*.pep.xml"
$PY $S/depletion_test.py --fasta $DB --counts $OUT/rifDIA_depletion_counts.tsv \
    --psm $W/results/psm.tsv --exclude-run 10989 > $OUT/rifDIA_depletion.txt
cat $OUT/rifDIA_depletion.txt
# keep a backed-up copy (scratch is not backed up)
mkdir -p ~/intein_results/TB_rif_DIA && cp -r $OUT ~/intein_results/TB_rif_DIA/
echo "copied to ~/intein_results/TB_rif_DIA/intein_analysis"
