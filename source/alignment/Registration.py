from .Utils import *
import numpy as np
import os
from copy import deepcopy


def Registration_given_feature(data_dir, source_dir, target_dir, source_sample_dir, target_sample_dir, VOXEL_SIZE=5.0, visualize=False, T=None, one_stage=False):

    sample_A_points, sample_A_normals = load_sample(source_sample_dir)
    sample_B_points, sample_B_normals = load_sample(target_sample_dir)
    A_pcd = o3d.geometry.PointCloud()
    A_pcd.points = o3d.utility.Vector3dVector(sample_A_points)
    A_pcd.normals = o3d.utility.Vector3dVector(sample_A_normals)
    B_pcd = o3d.geometry.PointCloud()
    B_pcd.points = o3d.utility.Vector3dVector(sample_B_points)
    B_pcd.normals = o3d.utility.Vector3dVector(sample_B_normals)

    A_key_pcd = load_xyz(source_dir)
    A_keypoint = np.array(A_key_pcd.points)
    # A_key_feats = read_features(source_feature)
    B_key_pcd = load_xyz(target_dir)
    B_keypoint = np.array(B_key_pcd.points)
    # B_key_feats = read_features(target_feature)

    temp_dir = "%s/temp" % data_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    A_key_feats = cal_SHOT(sample_A_points, sample_A_normals, temp_dir, A_keypoint, radius=VOXEL_SIZE * 7)
    B_key_feats = cal_SHOT(sample_B_points, sample_B_normals, temp_dir, B_keypoint, radius=VOXEL_SIZE * 7)

    corrs_A, corrs_B = find_correspondences(A_key_feats, B_key_feats, mutual_filter=True)
    if len(corrs_A) < 3:
        print("mutual failed")
        corrs_A, corrs_B = find_correspondences(A_key_feats, B_key_feats, mutual_filter=False)

    A_corr = A_keypoint[corrs_A, :].T
    B_corr = B_keypoint[corrs_B, :].T

    num_corrs = A_corr.shape[1]
    print(f'generates {num_corrs} putative correspondences.')

    if T is None:
        NOISE_BOUND = VOXEL_SIZE
        teaser_solver = get_teaser_solver(NOISE_BOUND)
        teaser_solver.solve(A_corr, B_corr)
        solution = teaser_solver.getSolution()
        R_teaser = solution.rotation
        t_teaser = solution.translation
        T_teaser = Rt2T(R_teaser, t_teaser)
        init_transformation = T_teaser

        if one_stage:
            if visualize:
                draw_registration_result(A_pcd, B_pcd, init_transformation)

            return init_transformation

        writePLY(A_pcd, "%s/source.ply" % temp_dir)
        writePLY(B_pcd, "%s/target.ply" % temp_dir)
        init_file = "%s/init_trans" % temp_dir
        write_trans(init_transformation, init_file)

        os.system(
            "alignment/ICP %s/target.ply %s/source.ply %s %s/init_trans 3" % (
            temp_dir, temp_dir, temp_dir, temp_dir))

        icp_file = "%s/final_trans" % temp_dir
        T_icp = read_trans(icp_file)

    else:
        T_icp =T

    if visualize:
        draw_registration_result(A_pcd, B_pcd, T_icp)

    os.system("rm -rf %s" % temp_dir)

    return T_icp


def direct_alignment(data_dir, source_name, target_name, VOXEL_SIZE=5):
    source_key_dir = "%s/Points_%s_Key.xyz" % (data_dir, source_name[4:-4])
    target_key_dir = "%s/Points_%s_Key.xyz" % (data_dir, target_name[4:-4])
    source_sample_dir = "%s/%s_%.2f.txt" % (data_dir, source_name[:-4], VOXEL_SIZE)
    target_sample_dir = "%s/%s_%.2f.txt" % (data_dir, target_name[:-4], VOXEL_SIZE)

    T_icp = Registration_given_feature(data_dir, source_key_dir, target_key_dir, source_sample_dir, target_sample_dir, VOXEL_SIZE=VOXEL_SIZE, visualize=False)
    # print("refine registration: ", calRMSD(source_pdb, source_sup, T_icp))
    return T_icp


def cal_pdb_RMSD(source_pdb_dir, source_sup_dir, T):
    source_pdb = getPointsFromPDB(source_pdb_dir)
    source_sup = getPointsFromPDB(source_sup_dir)
    return calRMSD(source_pdb, source_sup, T)


def Registration_mask(data_dir, A_pcd, B_pcd, A_key_pcd, B_key_pcd, A_key_feats, B_key_feats, mask, max_correspondence_dist=10.0, VOXEL_SIZE=5.0, store_partial=False):

    # A_points = np.array(A_pcd.points)
    B_points = np.array(B_pcd.points)
    A_keypoint = np.array(A_key_pcd.points)
    B_keypoint = np.array(B_key_pcd.points)

    """
    mask op
    """
    if mask["name"] == "sphere":
        print("spherical mask")
        center, radius = mask["center"], mask["radius"]
        print("center: ", center)
        mask_B_pcd, mask_indices = mask_sphere_points(B_pcd, center, radius, return_indices=True)
        mask_B_key_pcd, mask_key_indices = mask_sphere_points(B_key_pcd, center, radius, return_indices=True)
    elif mask["name"] == "rectangle":
        B_mask = mask["mask"]
        mask_B_pcd, mask_indices = mask_rectangle_points(B_pcd, B_mask, return_indices=True)
        mask_B_key_pcd, mask_key_indices = mask_rectangle_points(B_key_pcd, B_mask, return_indices=True)
    else:
        mask_B_pcd = B_pcd
        mask_B_key_pcd = B_key_pcd
        mask_indices = np.arange(B_points.shape[0])
        mask_key_indices = np.arange(B_keypoint.shape[0])

    if len(mask_indices) <= 100 or len(mask_key_indices) <= 10:
        print("mask content is small")
        return None, 0.0

    # mask_B_xyz = np.array(mask_B_pcd.points).T
    B_mask_key_feats = B_key_feats[mask_key_indices]
    B_mask_key_points = np.array(mask_B_key_pcd.points)

    # establish correspondences by nearest neighbour search in feature space
    corrs_A, corrs_B = find_correspondences(A_key_feats, B_mask_key_feats, mutual_filter=True)
    if len(corrs_A) < 3:
        corrs_A, corrs_B = find_correspondences(A_key_feats, B_mask_key_feats, mutual_filter=False)

    A_corr = A_keypoint[corrs_A, :].T  # np array of size 3 by num_corrs
    B_corr = B_mask_key_points[corrs_B, :].T  # np array of size 3 by num_corrs

    num_corrs = A_corr.shape[1]
    print(f'generates {num_corrs} putative correspondences.')

    NOISE_BOUND = VOXEL_SIZE
    teaser_solver = get_teaser_solver(NOISE_BOUND)
    teaser_solver.solve(A_corr, B_corr)
    solution = teaser_solver.getSolution()
    R_teaser = solution.rotation
    t_teaser = solution.translation
    T_teaser = Rt2T(R_teaser, t_teaser)

    def cal_score(A_pcd, B_pcd, max_correspondence_dist, final_T):
        A_points = np.array(A_pcd.points)
        eval_metric = eval(A_pcd, B_pcd, max_correspondence_dist, final_T)
        A_transformed = deepcopy(A_pcd)
        A_transformed.transform(final_T)
        correspondence_set = np.array(eval_metric.correspondence_set)

        if len(correspondence_set) < A_points.shape[0] * 0.1:
            print("correspondence_set is small")
            return 0.0

        A_points, B_points = np.array(A_transformed.points), np.array(B_pcd.points)
        A_vector, B_vector = np.array(A_transformed.normals), np.array(B_pcd.normals)

        # sigma = 3.5
        # weights = None
        # A_gmm = create_GMM(A_points, sigma, weights)
        # B_gmm = create_GMM(B_points, sigma, weights)
        # distri_score = cal_JS_KL(A_gmm, B_gmm, 10 ** 4)

        cosin_dist_list = np.diagonal(np.dot(A_vector[correspondence_set[:, 0]], B_vector[correspondence_set[:, 1]].T))
        cosin_dist = cosin_dist_list[cosin_dist_list >= 0.6].shape[0] / cosin_dist_list.shape[0]

        # score = cosin_dist * (1 - distri_score)
        score = cosin_dist

        return score

    init_transformation = T_teaser
    writePLY(A_pcd, "%s/source.ply" % data_dir)
    if store_partial:
        writePLY(mask_B_pcd, "%s/target_partial.ply" % data_dir)

        init_file = "%s/init_trans" % data_dir
        write_trans(init_transformation, init_file)
        os.system(
            "alignment/ICP %s/target_partial.ply %s/source.ply %s %s/init_trans 3" % (
                data_dir, data_dir, data_dir, data_dir))
    else:
        writePLY(B_pcd, "%s/target.ply" % data_dir)

        init_file = "%s/init_trans" % data_dir
        write_trans(init_transformation, init_file)
        os.system(
            "alignment/ICP %s/target.ply %s/source.ply %s %s/init_trans 3" % (
                data_dir, data_dir, data_dir, data_dir))

    icp_file = "%s/final_trans" % data_dir
    T_icp = read_trans(icp_file)
    final_T = T_icp

    score = cal_score(A_pcd, B_pcd, max_correspondence_dist, final_T)
    return final_T, score


def Registration_mask_list(data_dir, source_key_dir, source_sample_dir, target_key_dir, target_sample_dir, VOXEL_SIZE=5.0, store_partial=False):
    sample_A_points, sample_A_normals = load_sample(source_sample_dir)
    sample_B_points, sample_B_normals = load_sample(target_sample_dir)
    A_pcd = o3d.geometry.PointCloud()
    A_pcd.points = o3d.utility.Vector3dVector(sample_A_points)
    A_pcd.normals = o3d.utility.Vector3dVector(sample_A_normals)
    B_pcd = o3d.geometry.PointCloud()
    B_pcd.points = o3d.utility.Vector3dVector(sample_B_points)
    B_pcd.normals = o3d.utility.Vector3dVector(sample_B_normals)
    A_key_pcd = load_xyz(source_key_dir)
    B_key_pcd = load_xyz(target_key_dir)

    temp_dir = "%s/temp" % data_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    A_key_feats = cal_SHOT(sample_A_points, sample_A_normals, temp_dir, np.array(A_key_pcd.points),
                           radius=VOXEL_SIZE * 7.0)
    B_key_feats = cal_SHOT(sample_B_points, sample_B_normals, temp_dir, np.array(B_key_pcd.points),
                           radius=VOXEL_SIZE * 7.0)

    A_min_bound = A_pcd.get_min_bound()
    A_max_bound = A_pcd.get_max_bound()
    A_dist_cor = A_max_bound - A_min_bound
    A_mean_dis = np.mean(A_dist_cor)
    B_min_bound = B_pcd.get_min_bound()
    B_max_bound = B_pcd.get_max_bound()
    B_dist_cor = B_max_bound - B_min_bound
    B_mean_dis = np.mean(B_dist_cor)
    if np.abs(A_mean_dis - B_mean_dis) > 0:
        mask = {}
        radius = np.max(A_dist_cor) * 1.1 / 2
        center = B_min_bound + radius / 10
        terminal = B_max_bound
        step = int(radius // 2) # translation interval of mask

        max_correspondence_dist = 10.0
        store_partial = store_partial

        record_dir = "%s/record.txt" % data_dir
        record_T_dir = "%s/record_T.npy" % data_dir
        record_T = []

        with open(record_dir, 'w') as f:
            f.write("\t".join(["t_x", "t_y", "t_z", "score"]))
            f.write("\n")

            for i in range(0, int(terminal[0]), step):
                for j in range(0, int(terminal[1]), step):
                    for k in range(0, int(terminal[2]), step):
                        temp_center = [center[0] + i, center[1] + j, center[2] + k]

                        mask["name"] = "sphere"
                        mask["center"] = temp_center
                        mask["radius"] = radius

                        final_T, score = Registration_mask(temp_dir, A_pcd, B_pcd, A_key_pcd, B_key_pcd, A_key_feats, B_key_feats, mask, max_correspondence_dist=max_correspondence_dist, VOXEL_SIZE=VOXEL_SIZE, store_partial=store_partial)
                        
                        if final_T is not None:
                            f.write("\t".join(
                                ["%.2f" % i, "%.2f" % j, "%.2f" % k, "%.2f" % score]))
                            f.write("\n")
                            record_T.append(final_T)
                        else:
                            f.write("\t".join(
                                ["%.2f" % i, "%.2f" % j, "%.2f" % k, "0.00"]))
                            f.write("\n")
                            record_T.append(np.identity(4))
        
        np.save(record_T_dir, np.stack(record_T, axis=0))

    else:
        print("please exchange the source map and target one")

    os.system("rm -rf %s" % temp_dir)

    return


def mask_alignment(data_dir, source_name, target_name, VOXEL_SIZE, store_partial=False):
    source_key_dir = "%s/Points_%s_Key.xyz" % (data_dir, source_name[4:-4])
    target_key_dir = "%s/Points_%s_Key.xyz" % (data_dir, target_name[4:-4])
    source_sample_dir = "%s/%s_%.2f.txt" % (data_dir, source_name[:-4], VOXEL_SIZE)
    target_sample_dir = "%s/%s_%.2f.txt" % (data_dir, target_name[:-4], VOXEL_SIZE)
    Registration_mask_list(data_dir, source_key_dir, source_sample_dir, target_key_dir, target_sample_dir, VOXEL_SIZE=5.0, store_partial=store_partial)
    return


def extract_top_K(record_dir, record_T_dir, K, save_dir, source_pdb_dir=None, source_sup_dir=None):
    score_list = []
    with open(record_dir, "r") as f:
        lines = f.readlines()[1:]
        for line in lines:
            score = float(line.strip().split()[-1])
            score_list.append(score)
    lst_with_index = [(idx, val) for idx, val in enumerate(score_list)]
    sorted_list = sorted(lst_with_index, key=lambda x: x[1], reverse=True)
    top_k_index = [x[0] for x in sorted_list[:K]]

    record_T_list = np.load(record_T_dir)
    with open(save_dir, "w") as f:
        f.write("\t".join(["score", "transformation matrix", "RMSD"]))
        f.write("\n")
        for idx in top_k_index:
            score = score_list[idx]
            T = record_T_list[idx]
            f.write("%.2f\t" % score)
            f.write(",".join(["%.4f" % i for i in T.reshape(-1).tolist()]))
            if source_pdb_dir:
                rmsd = cal_pdb_RMSD(source_pdb_dir, source_sup_dir, T)
                f.write("\t%.2f\n" % rmsd)
            else:
                f.write("\n")
    return
