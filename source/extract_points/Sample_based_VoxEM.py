from .VoxEM import *
from .Supporting import *


def Sample_Cluster(data_dir, map_name, threshold, VOXEL_SIZE):

    # Activate cuda
    cuda = torch.cuda.is_available()
    torch.cuda.set_device(0)
    Tensor = torch.cuda.FloatTensor

    # 0. sampling points and generate vectors

    mrc_file = "%s/%s" % (data_dir, map_name)
    sample_file = "%s/%s_%.2f.txt" % (data_dir, map_name[:-4], VOXEL_SIZE)
    os.system("extract_points/Sample -a %s -t %.4f -s %.2f > %s" % (mrc_file, threshold, VOXEL_SIZE, sample_file))

    # 1. Load Map
    mrcobject = VoxEM()
    mrcobject.IO_ReadMrc(mrc_inputname=mrc_file, voxel_outputname='Default', description_outputname='Default',
                        statistics_outputname='Default')
    # 2. Prune Non-Negative
    mrcobject.Voxel_Prune_RangeZero(lowerbound=threshold, upperbound=None, inputname='Default', outputname='Voxel')
    vox = mrcobject.voxel_workspace['Voxel']
    start = mrcobject.description_workspace['Default']['Start']
    origin = mrcobject.description_workspace['Default']['Origin']
    vox_length = mrcobject.description_workspace['Default']['Angstrom']
    vox_ang = np.array(vox_length).astype(np.float32) / np.array(vox.shape).astype(np.float32)

    sample_points, sample_vector = load_sample_points(sample_file)

    shifting = start * vox_ang + origin
    sample_points = sample_points - shifting

    mrcobject.point_workspace["sample"] = sample_points.T
    print("sample points: ", sample_points.shape)

    # 3. Meanshift

    mrcobject.Point_Create_Meanshift_sample(lower_bound=threshold, window=17, voxel_inputname='Voxel', bandwidth=3.0,
                                            point_outputname='Meanshift', point_inputname='sample', iteration=2000,
                                            convergence=0.000187, step_size=0.05)

    # 4. Dbscan
    mrcobject.Point_Transform_DBSCAN(distance_tolerance=VOXEL_SIZE * 3, clustersize_tolerance=3, inputname='Meanshift', outputname='Meanshift')

    print("meanshift points: ", mrcobject.point_workspace["Meanshift"].shape)

    mrcobject.Point_Transform_DBSCAN(distance_tolerance=VOXEL_SIZE, clustersize_tolerance=1, inputname='Meanshift', outputname='DBSCAN', centroid_only=True)

    print("DBSCAN points: ", mrcobject.point_workspace["DBSCAN"].shape)

    mrcobject.IO_WriteXYZ(point_inputname='DBSCAN', file_outputname='%s/Points_%s_Key.xyz' % (data_dir, map_name[4:8]), atom_name='H', shifting=shifting)

    print("Sampling and clustering success")

