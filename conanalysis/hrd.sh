cat *.ddg > ddg
./rosetta_cartddg_filter.py ddg > tmp0
sort -nk1 tmp0 > tmp1
grep -v "#" *.fasta.ca > tmp2
paste tmp2 tmp1 > hybrid_design.sc 
