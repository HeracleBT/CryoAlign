import open3d as o3d
import numpy as np
from scipy.spatial import cKDTree
import teaserpp_python
import copy
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture
from scipy.ndimage import affine_transform, map_coordinates
from sklearn.mixture._gaussian_mixture import _compute_precision_cholesky
from sklearn.neighbors import radius_neighbors_graph
import re
import os
from collections import defaultdict
import scipy.special as sp
import ot
import random
import math


def load_data(file_path, file_down_path):
    points = []
    with open(file_path, "r") as f:
        for line in f.readlines()[8:]:
            points.append([float(i) for i in line.strip().split()])
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points)

    points_down = []
    with open(file_down_path, "r") as f:
        for line in f.readlines()[8:]:
            points_down.append([float(i) for i in line.strip().split()])
    pcd_down = o3d.geometry.PointCloud()
    pcd_down.points = o3d.utility.Vector3dVector(points_down)
    points_down = np.array(points).T

    return pcd, pcd_down, points_down


def load_xyz(file_path):
    pcd = o3d.io.read_point_cloud(file_path, format="xyz")
    return pcd


def load_sample(file_path, density=False):
    point_list = []
    vector_list = []
    density_list = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        sample = float(lines[0].strip())
        origin_x, origin_y, origin_z = [float(i) for i in lines[3].strip().split()]
        for i in range(5, len(lines)):
            line = lines[i]
            if i % 2:
                _, x, y, z = line.strip().split()
                point_list.append([float(x) * sample + origin_x, float(y) * sample + origin_y, float(z) * sample + origin_z])
            else:
                v_x, v_y, v_z, d = line.strip().split()
                vector_list.append([float(v_x), float(v_y), float(v_z)])
                density_list.append([float(d)])
    if density:
        return np.array(point_list), np.array(vector_list), np.array(density_list)
    return np.array(point_list), np.array(vector_list)


def writePLY(pcd, file_path, normal=False):
    points = np.array(pcd.points)
    vertex_num = points.shape[0]
    if normal:
        normals = np.array(pcd.normals)
        with open(file_path, "w") as f:
            f.write("ply\n")
            f.write("format ascii 1.0\n")
            f.write("comment Created by Open3D\n")
            f.write("element vertex {}\n".format(vertex_num))
            f.write("property float x\n")
            f.write("property float y\n")
            f.write("property float z\n")
            f.write("property float nx\n")
            f.write("property float ny\n")
            f.write("property float nz\n")
            f.write("end_header\n")
            for i in range(vertex_num):
                f.write("{} {} {} {} {} {}\n".format(points[i, 0], points[i, 1], points[i, 2],
                                                     normals[i, 0], normals[i, 1], normals[i, 2]))
    else:
        with open(file_path, "w") as f:
            f.write("ply\n")
            f.write("format ascii 1.0\n")
            f.write("comment Created by Open3D\n")
            f.write("element vertex {}\n".format(vertex_num))
            f.write("property float x\n")
            f.write("property float y\n")
            f.write("property float z\n")
            f.write("end_header\n")
            for i in range(vertex_num):
                f.write("{} {} {}\n".format(points[i, 0], points[i, 1], points[i, 2]))
    return


def readPLY(file_path):
    point_list = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        for line in lines[8:]:
            # print(line)
            point_list.append(line.strip().split())
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(point_list)
    return pcd


def write_trans(transformation, file):
    with open(file, "w") as f:
        for i in range(4):
            for j in range(4):
                f.write("%.6f" % transformation[i, j])
                f.write("\t")
            f.write("\n")
    return


def read_trans(file):
    transformation = np.zeros((4, 4))
    transformation[3, 3] = 1.0
    with open(file, "r") as f:
        temp = f.readlines()
        for i in range(3):
            for j in range(4):
                transformation[i, j] = float(temp[i].strip().split()[j])
    return transformation


def read_MMAlign_T(file):
    with open(file, "r") as f:
        T = np.zeros((4, 4))
        T[3, 3] = 1.0
        for line in f.readlines():
            if line.startswith("0"):
                num_list = line.strip().split()
                T[0, 3] = num_list[1]
                T[0, 0] = num_list[2]
                T[0, 1] = num_list[3]
                T[0, 2] = num_list[4]
            elif line.startswith("1"):
                num_list = line.strip().split()
                T[1, 3] = num_list[1]
                T[1, 0] = num_list[2]
                T[1, 1] = num_list[3]
                T[1, 2] = num_list[4]
            elif line.startswith("2"):
                num_list = line.strip().split()
                T[2, 3] = num_list[1]
                T[2, 0] = num_list[2]
                T[2, 1] = num_list[3]
                T[2, 2] = num_list[4]
    return T


def getPointsFromPDB(pdbFile='*.pdb', flag=None):
    all_atom_coord = []
    if flag is None:
        flag = pdbFile[-3:]
    if flag == "pdb":
        with open(pdbFile, 'r') as fr:
            for line in fr:
                if line.startswith('ATOM'):
                    temp_atom_coord = line[27:54].strip().split()
                    try:
                        if len(temp_atom_coord) == 3:
                            atom_coord = [float(item) for item in temp_atom_coord if item != '']
                            all_atom_coord.append(atom_coord)
                        elif len(temp_atom_coord) == 2:
                            if len(temp_atom_coord[0]) > len(temp_atom_coord[1]):
                                temp1, temp2 = temp_atom_coord[0][:-8], temp_atom_coord[0][-8:]
                                atom_coord = [float(temp1), float(temp2), float(temp_atom_coord[1])]
                            else:
                                # print(temp_atom_coord)
                                temp1, temp2 = temp_atom_coord[1][:-8], temp_atom_coord[1][-8:]
                                atom_coord = [float(temp_atom_coord[0]), float(temp1), float(temp2)]
                            all_atom_coord.append(atom_coord)
                    # throw away bad coordinate, e.g., '49.680-100.424  42.761'
                    except:
                        print('throw away 3D coordinate %s' %
                              str(line[27:54].strip()))
                        print(pdbFile)
    elif flag == "cif":
        with open(pdbFile, 'r') as fr:
            for line in fr:
                if line.startswith('ATOM'):
                    temp_atom_coord = line.strip().split()
                    # print(temp_atom_coord[10:13])
                    try:
                        atom_coord = [float(item) for item in temp_atom_coord[10:13] if item != '']
                        all_atom_coord.append(atom_coord)
                    # throw away bad coordinate, e.g., '49.680-100.424  42.761'
                    except:
                        print('throw away 3D coordinate %s' % " ".join(temp_atom_coord[10:13]))
                        print(pdbFile)
    if all_atom_coord == []:
        print('Warning: point set is empty!')
    return np.array(all_atom_coord)


def pcd2xyz(pcd):
    return np.asarray(pcd.points).T


def find_knn_cpu(feat0, feat1, knn=1, return_distance=False):
    feat1tree = cKDTree(feat1)
    dists, nn_inds = feat1tree.query(feat0, k=knn, n_jobs=-1)
    if return_distance:
        return nn_inds, dists
    else:
        return nn_inds


def find_correspondences(feats0, feats1, mutual_filter=True):
    nns01 = find_knn_cpu(feats0, feats1, knn=1, return_distance=False)
    corres01_idx0 = np.arange(len(nns01))
    corres01_idx1 = nns01

    # print(min(corres01_idx1))
    # print(max(corres01_idx1))

    if not mutual_filter:
        return corres01_idx0, corres01_idx1

    nns10 = find_knn_cpu(feats1, feats0, knn=1, return_distance=False)
    corres10_idx1 = np.arange(len(nns10))

    if max(corres01_idx1) == feats1.shape[0]:
        corres10_idx0 = np.append(nns10, np.nan)
    else:
        corres10_idx0 = nns10

    mutual_filter = (corres10_idx0[corres01_idx1] == corres01_idx0)
    corres_idx0 = corres01_idx0[mutual_filter]
    corres_idx1 = corres01_idx1[mutual_filter]

    return corres_idx0, corres_idx1


def get_teaser_solver(noise_bound):
    solver_params = teaserpp_python.RobustRegistrationSolver.Params()
    solver_params.cbar2 = 1.0
    # solver_params.cbar2 = 1.0
    solver_params.noise_bound = noise_bound
    solver_params.estimate_scaling = False
    solver_params.inlier_selection_mode = \
        teaserpp_python.RobustRegistrationSolver.INLIER_SELECTION_MODE.PMC_EXACT
    solver_params.rotation_tim_graph = \
        teaserpp_python.RobustRegistrationSolver.INLIER_GRAPH_FORMULATION.CHAIN
    # solver_params.rotation_tim_graph = \
    #     teaserpp_python.RobustRegistrationSolver.INLIER_GRAPH_FORMULATION.COMPLETE
    solver_params.rotation_estimation_algorithm = \
        teaserpp_python.RobustRegistrationSolver.ROTATION_ESTIMATION_ALGORITHM.GNC_TLS
    solver_params.rotation_gnc_factor = 1.4
    # solver_params.rotation_gnc_factor = 1.4
    # solver_params.rotation_max_iterations = 100
    solver_params.rotation_max_iterations = 10000
    solver_params.rotation_cost_threshold = 1e-16
    solver = teaserpp_python.RobustRegistrationSolver(solver_params)
    return solver


def Rt2T(R, t):
    T = np.identity(4)
    T[:3, :3] = R
    T[:3, 3] = t
    return T


def draw_registration_result(source, target, transformation=None):
    source_temp = copy.deepcopy(source)
    target_temp = copy.deepcopy(target)
    source_temp.paint_uniform_color([1, 0.706, 0])
    target_temp.paint_uniform_color([0, 0.651, 0.929])
    # source_temp.paint_uniform_color([0, 0.651, 0.929])
    # target_temp.paint_uniform_color([1, 0.706, 0])
    if transformation is not None:
        source_temp.transform(transformation)
        source_temp.paint_uniform_color([1, 0.706, 0])
        # source_temp.paint_uniform_color([0, 0.651, 0.929])
    o3d.visualization.draw_geometries([source_temp, target_temp],
                                      zoom=0.4559,
                                      front=[0.6452, -0.3036, -0.7011],
                                      lookat=[1.9892, 2.0208, 1.8945],
                                      up=[-0.2779, -0.9482, 0.1556])


def create_GMM(array, sigma, weights=None):
    n_gaussians = len(array)
    covirance = sigma * np.ones(n_gaussians)
    if weights is None:
        weights = (1. / n_gaussians) * np.ones(n_gaussians)
        # print(weights.shape)
    gmm = GaussianMixture(n_components=n_gaussians, covariance_type="spherical")
    gmm.weights_ = weights
    gmm.means_ = array
    gmm.covariances_ = covirance
    gmm.precisions_cholesky_ = _compute_precision_cholesky(covirance, "spherical")
    return gmm


def set_weights(array):
    if len(array) < 1:
        return Nonec
    weights = np.abs(array)
    weights = weights - np.min(weights)
    weights = weights / np.sum(weights)
    return weights


def cal_KL(gmm_p, gmm_q, n_samples):
    X, _ = gmm_p.sample(n_samples)
    log_p_X = gmm_p.score_samples(X)
    log_q_X = gmm_q.score_samples(X)
    return (log_p_X.mean() - log_q_X.mean()) / (-log_q_X.mean())


def cal_JS_KL(gmm_p, gmm_q, n_samples):
    X, _ = gmm_p.sample(n_samples)
    log_p_X = gmm_p.score_samples(X)
    log_q_X = gmm_q.score_samples(X)
    log_mix_X = np.logaddexp(log_p_X, log_q_X)

    Y, _ = gmm_q.sample(n_samples)
    log_p_Y = gmm_p.score_samples(Y)
    log_q_Y = gmm_q.score_samples(Y)
    log_mix_Y = np.logaddexp(log_p_Y, log_q_Y)

    return (log_p_X.mean() - (log_mix_X.mean() - np.log(2))
            + log_q_Y.mean() - (log_mix_Y.mean() - np.log(2))) / 2


def calRMSD(source, target, transformation=None):
    atom_num = source.shape[0]
    if transformation is not None:
        homo_source = np.concatenate([source, np.ones([atom_num, 1])], axis=1)
        transformed_source = np.matmul(transformation, homo_source.T).T
        # print(transformed_source[:2, :])
    else:
        transformed_source = source
    msd = np.sum(np.sqrt(np.reshape(np.sum(np.square(transformed_source[:, :3] - target), axis=1), (atom_num,))))
    return msd / atom_num


def calDist(source, target):
    return np.sqrt(np.sum(np.square(source - target), axis=1))


def getTransformed(array, T):
    if array.shape[1] != T.shape[0]:
        array = np.concatenate([array, np.ones((array.shape[0], 1))], axis=1)
    transformed = np.matmul(T, array.T).T
    return transformed[:, :3]


def eval(A_pcd, B_pcd, threshold, trans):
    return o3d.pipelines.registration.evaluate_registration(A_pcd, B_pcd, threshold, trans)


def save_xyz(pcd, data_dir, PC_FLAG=False):
    XYZTrial = []
    if PC_FLAG:
        Points = np.array(pcd.points)
    else:
        Points = pcd
    for point in Points:
        XYZTrial.append('%.5f        %.5f        %.5f\n' % (point[0], point[1], point[2]))
    with open(data_dir, 'w') as f:
        for point in XYZTrial:
            f.write(point)
    return


def transform_pdb(pdbFile, output_dir, T):

    if pdbFile[-3:] == "pdb":
        with open(output_dir, "w") as of:
            with open(pdbFile, 'r') as fr:
                for line in fr:
                    if line.startswith('ATOM'):
                        temp_atom_coord = line[27:55].strip().split()

                        # throw away bad coordinate, e.g., '49.680-100.424  42.761'
                        if len(temp_atom_coord) == 3:
                            atom_coord = [float(item) for item in temp_atom_coord if item != '']
                        elif len(temp_atom_coord) == 2:
                            if len(temp_atom_coord[0]) > len(temp_atom_coord[1]):
                                temp1, temp2 = temp_atom_coord[0][:-8], temp_atom_coord[0][-8:]
                                atom_coord = [float(temp1), float(temp2), float(temp_atom_coord[1])]
                            else:
                                # print(temp_atom_coord)
                                temp1, temp2 = temp_atom_coord[1][:-8], temp_atom_coord[1][-8:]
                                atom_coord = [float(temp_atom_coord[0]), float(temp1), float(temp2)]
                        transformed_atom_coord = np.matmul(T, np.array(atom_coord + [1.0]))
                        transformed_atom_line = ["%.3f" % i for i in list(transformed_atom_coord[:3])]
                        transformed_atom_line = [" " * (8 - len(i)) + i for i in transformed_atom_line]
                        transformed_line = "".join(transformed_atom_line)
                        transformed_line = " " * (27 - len(transformed_line)) + transformed_line
                        # print(len(transformed_line))
                        transformed_line = line[:27] + transformed_line + line[54:]
                        # print(len(line), len(transformed_line))
                        of.write(transformed_line)
                        # print(len(line), len(transformed_line))
                    else:
                        of.write(line)

    elif pdbFile[-3:] == "cif":
        with open(pdbFile, 'r') as fr:
            for line in fr:
                if line.startswith('ATOM'):
                    temp_atom_coord = line.strip().split()
                    # print(temp_atom_coord[10:13])
                    try:
                        atom_coord = [float(item) for item in temp_atom_coord[10:13] if item != '']
                    # throw away bad coordinate, e.g., '49.680-100.424  42.761'
                    except:
                        print('throw away 3D coordinate %s' % " ".join(temp_atom_coord[10:13]))
                        print(pdbFile)
    return


def mask_sphere_points(pcd, center, radius, return_indices=False):
    raw_points = np.array(pcd.points)
    n_points = raw_points.shape[0]
    tile_center = np.tile(center, (n_points, 1))
    dists = np.sqrt(np.sum(np.square(raw_points - tile_center), axis=1))
    indices = np.where(dists <= radius)[0]
    if return_indices:
        return pcd.select_by_index(indices), indices
    return pcd.select_by_index(indices)


def mask_rectangle_points(pcd, mask, return_indices=False):
    raw_points = np.array(pcd.points)
    min_indices = (raw_points[:, 0] > mask[0, 0]) & (raw_points[:, 1] > mask[1, 0]) & (raw_points[:, 2] > mask[2, 0])
    # print(np.sum(min_indices))
    max_indices = (raw_points[:, 0] < mask[0, 1]) & (raw_points[:, 1] < mask[1, 1]) & (raw_points[:, 2] < mask[2, 1])
    # print(np.sum(max_indices))
    indices = min_indices & max_indices
    indices = np.where(indices == True)[0]
    # print(np.sum(indices))
    if return_indices:
        return pcd.select_by_index(indices), indices
    return pcd.select_by_index(indices)


def txt2pcd(pcd_points, output):

    with open(output, 'w') as f:
        f.writelines("# .PCD v0.7 - Point Cloud Data file format\n")
        f.writelines("VERSION 0.7\n")
        f.writelines("FIELDS x y z\n")
        f.writelines("SIZE 4 4 4\n")
        f.writelines("TYPE F F F\n")
        f.writelines("COUNT 1 1 1\n")
        f.writelines("WIDTH " + str(pcd_points.shape[0]) + "\n")
        f.writelines("HEIGHT 1\n")
        f.writelines("VIEWPOINT 0 0 0 1 0 0 0\n")
        f.writelines("POINTS " + str(pcd_points.shape[0]) + "\n")
        f.writelines("DATA ascii\n")

        '''for i in range(len(xlist)):
            file_to_write.writelines(str(xlist[i]) + " " + str(ylist[i]) + " " + str(zlist[i]) + "\n")'''

        for i in range(pcd_points.shape[0]):
            f.write(" ".join(["%.5f" % j for j in pcd_points[i]]))
            f.write("\n")


def read_features(feature_dir, mode='SHOT', key=False):
    features = []
    with open(feature_dir, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if not key:
                if mode == 'SHOT' or mode == '3DSC' or mode == 'USC':
                    line = lines[i].strip().split('(')[2][:-1]
                else:
                    line = lines[i].strip()[1:-1]
            else:
                if mode == 'SHOT':
                    line = lines[i].strip().split('(')[2][:-1]
                else:
                    line = lines[i].strip()[1:-1]

            features.append([float(k) for k in line.split(',')])
    return np.stack(features, axis=0)


def save_features(features, feature_dir):
    with open(feature_dir, "w") as f:
        for i in range(len(features)):
            histo = []
            for k in features[i]:
                if k == 0.0:
                    histo.append("0")
                else:
                    histo.append(str(k))
            line = "(" + ", ".join(histo) + ")"
            f.write(line)
            f.write("\n")
    return


def cal_SHOT(A_points, A_normals, temp_dir, A_key_points, save_dir=None, radius=25.0):
    points_dir = "%s/points.pcd" % temp_dir
    txt2pcd(A_points, points_dir)
    normal_dir = "%s/normals.txt" % temp_dir
    np.savetxt(normal_dir, A_normals, fmt='%.5f')
    key_points_dir = "%s/key_points.pcd" % temp_dir
    txt2pcd(A_key_points, key_points_dir)
    if save_dir is None:
        feature_dir = "%s/SHOT_features.txt" % temp_dir
    else:
        feature_dir = save_dir
    os.system("alignment/point_cloud_feature %s %s %s %.2f > %s" % (points_dir, normal_dir, key_points_dir, radius, feature_dir))
    A_features = read_features(feature_dir)
    return A_features
