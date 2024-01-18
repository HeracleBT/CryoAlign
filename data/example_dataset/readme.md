## dataset architecture
In each folder: example_"source emd id"_"target emd id"

File list:

Initial points and density vectors for density map: emd_\*\*\*\*\_5.txt

Extracted key points: Points_*.xyz

PDB atomic models: *.pdb or *.cif

Ground truth transformation matrix (calculated by MM_align): *_*_mmalign.txt

Transformation matrix of CryoAlign: "source map id"_"target map id".npy

The combinations of different keypoint detectors and feature descriptions: "map id"_"description"\_keypoint\_"detector" (Collected descriptions: SHOT, PFH, FPFH, 3DSC, USC, ROPS; Collected detectors: Harris, ISS, SIFT, CryoAlign)
