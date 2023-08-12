from argparse import ArgumentParser
from alignment.Registration import direct_alignment, cal_pdb_RMSD, mask_alignment, extract_top_K
from extract_points.Sample_based_VoxEM import Sample_Cluster
import os

"""
shell script
"""
parser = ArgumentParser(description="This script will align the source point cloud to the target one")
parser.add_argument("--data_dir", type=str, default="", help="file dir")
parser.add_argument("--source", type=str, default="", help="source emdb num")
parser.add_argument("--target", type=str, default="", help="target emdb num")
parser.add_argument("--source_contour", type=float, default=5.0, help="author recommend contour level")
parser.add_argument("--target_contour", type=float, default=5.0, help="author recommend contour level")
parser.add_argument("--source_pdb", type=str, default="", help="source pdb name")
parser.add_argument("--source_sup_pdb", type=str, default="", help="transformed source pdb name (ground truth)")
parser.add_argument("--voxel", type=float, default=5.0, help="VOXEL_SIZE")
parser.add_argument("--mask_file", type=str, default="", help='mask file')
parser.add_argument("--seg", action='store_true', default=False, help='segmentation or not')
parser.add_argument("--atom", action='store_true', default=False, help='atomic model or not')

args = parser.parse_args()

data_dir = args.data_dir
source_name = args.source
target_name = args.target
source_contour = args.source_contour
target_contour = args.target_contour
source_pdb_name = args.source_pdb
source_sup_pdb_name = args.source_sup_pdb
VOXEL_SIZE = args.voxel
seg_flag = args.seg
atom_flag = args.atom
mask_file = args.mask_file

"""
sample and clustering
"""
# Sample_Cluster(data_dir, source_name, source_contour, VOXEL_SIZE)
# Sample_Cluster(data_dir, target_name, target_contour, VOXEL_SIZE)

if not seg_flag:

    """
    Gloabl and local registration (No mask version)

    example: python main.py --data_dir ../data/emd_3695_3696 --source emd_3695.map --target emd_3696.map --source_contour 0.008 --target_contour 0.002 --source_pdb 5nsr.pdb --source_sup_pdb 5nsr_sup.pdb --voxel 5.0
    """

    T = direct_alignment(data_dir, source_name, target_name, VOXEL_SIZE)
    print("estiatmed transformation matrix: ", T)

    # source_pdb_dir = "%s/%s" % (data_dir, source_pdb_name)
    # source_sup_dir = "%s/%s" % (data_dir, source_sup_pdb_name)
    # print("RMSD between estiamted transformed PDB and ground truth", cal_pdb_RMSD(source_pdb_dir, source_sup_dir, T))

else:

    """
    Local registration (mask version)

    example: python main.py --data_dir ../data/emd_3661_6647 --source emd_3661.map --target emd_6647.map --source_contour 0.07 --target_contour 0.017 --source_pdb 5no2.pdb --source_sup_pdb 5no2_5juu_sup.pdb --voxel 5.0 --seg
    """
    # mask_alignment(data_dir, source_name, target_name, VOXEL_SIZE)

    record_dir = "%s/record.txt" % data_dir
    record_T_dir = "%s/record_T.npy" % data_dir
    K = 10
    save_dir = "%s/extract_top_%d.txt" % (data_dir, K)
    source_pdb_dir = "%s/%s" % (data_dir, "5no2.pdb")
    source_sup_dir = "%s/%s" % (data_dir, "5no2_5juu_sup.pdb")
    extract_top_K(record_dir, record_T_dir, K, save_dir, source_pdb_dir=source_pdb_dir, source_sup_dir=source_sup_dir)
    
