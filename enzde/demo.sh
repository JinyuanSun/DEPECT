python3 singlepointdesign.py -l RPBE_ref.pdbqt -r 4hs9_p.pdb -bc RPBE_ref.boxcfg -mf zlj.mutfile -m list > log
python3 depect_geo.py depect_enzde.sc zlj.cst RPBE_ref.pdbqt
python3 combine_score.py
echo "Demo succeed!"