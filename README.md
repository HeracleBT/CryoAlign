# CryoAlign

Here is an official implementation for CryoAlign. (The current code is a temporary version and will be re-written in a unified language like C++ in the future.)

Advances on cryo-electron imaging technologies have led to a rapidly increasing number of density maps. Alignment and comparison of density maps play a crucial role in interpreting structural information, such as conformational heterogeneity analysis using global alignment and atomic model assembly through local alignment. Here, we propose a fast and accurate global and local cryo-electron microscopy density map alignment method CryoAlign, which leverages local density feature descriptors to capture spatial structure similarities. CryoAlign is the first feature-based EM map alignment tool, in which the employment of feature-based architecture enables the rapid establishment of point pair correspondences and robust estimation of alignment parameters. Extensive experimental evaluations demonstrate the superiority of CryoAlign over the existing methods in both alignment accuracy and speed.

## Operation System

Ubuntu 18.04 or later

## Requirements

The current version needs both python and C++ environment. If you have a Docker environment with GPU, we strongly recommend you to generate image using the dockerfile.

### Python environment:

Python 3.6

open3d 0.13.0

pytorch 1.7

numpy 1.19

scipy 1.5

scikit-learn 0.24.0

TEASER++

### C++ environment:

FFTW

PCL(point cloud libraries)

EIGEN

### Docker

Build image from dockerfile and create a container from the image

```
docker build -f dockerfile -t [image name] .
docker run -it --name [container name] --gpu all [image name]
```

Enter the container, install PCL lib

```
sh /CryoAlign/docker/pcl.sh
```

## Usage

In the "source" directory, there is a main.py file.
```
cd /CryoAlign/source
```

```
Usage: python main.py --data_dir [data dir] --source [SOURCE_MAP.map] --target [TARGET_MAP.map] [(options)]

---options---
--source_contour [float] : Threshold of source density map
--target_contour [float] : Threshold of target density map
--source_pdb [str] : Corresponding pdb of source density map
--source_sup_pdb [str] : Ground truth of transformed source pdb
--voxel [float] : Sampling interval
--seg : Mask version or not
```

## For global and local registration (no mask version)

example: data/emd_3695_3696

run the following code 
```
python main.py --data_dir ../data/emd_3695_3696 --source emd_3695.map --target emd_3696.map --source_contour 0.008 --target_contour 0.002 --source_pdb 5nsr.pdb --source_sup_pdb 5nsr_sup.pdb --voxel 5.0
```

## For local registration (mask version)

example: data/emd_3661_6647

run the following code
```
python main.py --data_dir ../data/emd_3661_6647 --source emd_3661.map --target emd_6647.map --source_contour 0.07 --target_contour 0.017 --source_pdb 5no2.pdb --source_sup_pdb 5no2_5juu_sup.pdb --voxel 5.0 --seg
```

## For atomic modeling

example: data/emd_4775_4776

We extract chains A and B from EMD-4775 (PDB ID: 6rah) and rigidly align them into map EMD-4776 (PDB ID: 6rai).
```
python atomic_model_example.py --stage rigid --visualize
```

The superimposition of point clouds is visualized, in which chain A is colored red, chain B is colored green, and the target map is colored blue. The rigid parameters are also applied to corresponding PDB structures, named "chain_*_rigid.pdb".

### Flexible fitting

Based on above well-assembly of point clouds, CryoAlign can further perform flexible fitting by integrating with Bayesian Coherent Point Drift (O. Hirose, "[Geodesic-Based Bayesian Coherent Point Drift](https://ieeexplore.ieee.org/document/9918058)," [IEEE TPAMI](https://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=34), Oct 2022.).
```
python atomic_model_example.py --stage flexible --visualize
```

For convenience, we provide the result of BCPD, named "output_flexible.txt", and the visualization of transformed point clouds is shown. Furthermore, the point displacement between matched points is performed to the atomic model, named "chain_*_flexible.pdb".