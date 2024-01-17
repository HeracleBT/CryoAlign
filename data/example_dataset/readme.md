## dataset architecture
in each directory: example_"source emd id"_"target emd id"

file list:

Initial points and density vectors for density map: emd_\*\*\*\*\_5.txt

Extracted key points: Points_*.xyz

PDB atomic models: *.pdb or *.cif

Ground truth transformation matrix (calculated by MM_align): *_*_mmalign.txt

Transformation matrix of CryoAlign: "source map id"_"target map id".npy
