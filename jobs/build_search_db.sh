#!/bin/bash
# =============================================================================
# build_search_db.sh -- the Dataset 3 search database (Methods 3.6)
# UniProt H37Rv (UP000001584, 3,997) + 3 spliced forms -> 4,000 entries;
# Philosopher adds 118 contaminants and reversed decoys -> 8,236 entries.
# =============================================================================

# Set PROJ to the directory holding the FragPipe run, the search
# database and the results, for example:
#   export PROJ=/scratch/$USER/TB_rif_DIA
: "${PROJ:?set PROJ to your project directory}"
set -e
cd $PROJ/db
cat ~/UP000001584_83332.fasta ~/intein_results/test_PXD009315/spliced_search/spliced_three.fasta \
    > h37rv_uniprot_plus_spliced.fasta
grep -c '>' h37rv_uniprot_plus_spliced.fasta                      # expect 4000
PHILO=~/software/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9
$PHILO workspace --init
$PHILO database --custom h37rv_uniprot_plus_spliced.fasta --contam
$PHILO workspace --clean
F=$(ls -t *-decoys-contam-h37rv_uniprot_plus_spliced.fasta.fas | head -1)
echo "total $(grep -c '>' $F)  decoys $(grep -c '>rev_' $F)"     # expect 8236 / 4118
