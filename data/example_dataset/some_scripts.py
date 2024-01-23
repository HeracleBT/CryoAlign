import numpy as np
import open3d as o3d
from alignment.Registration import find_correspondences


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


def cal_feature_corres(A_pcd, A_features, B_pcd, B_features, gt_T, accept_thres):
    gt_eval = o3d.pipelines.registration.evaluate_registration(A_pcd, B_pcd, accept_thres, gt_T)
    gt_corres = np.array(gt_eval.correspondence_set)
    corrs_A, corrs_B = find_correspondences(A_features, B_features, mutual_filter=False)
    B_points = np.array(B_pcd.points)
    pre_list = []
    recall_list = []
    match_num = 0
    acceptable_thres = accept_thres
    correct_corres = []
    false_corres = []
    total_match = len(corrs_A)
    if total_match != 0:
        nns01 = np.ones((len(A_features))) * (-1)
        nns01[corrs_A] = corrs_B
        for i in range(len(gt_corres)):
            corres_i, corres_j = gt_corres[i]
            gt_point = B_points[corres_j]
            knn_corres_j = int(nns01[corres_i])
            if knn_corres_j == -1 or knn_corres_j == len(B_features):
                continue
            target_point = B_points[knn_corres_j]
            if np.linalg.norm((gt_point - target_point)) <= acceptable_thres:
                correct_corres.append([corres_i, corres_j])
                match_num += 1
            else:
                false_corres.append([corres_i, corres_j])
        pre_ratio = match_num / total_match
        recall_ratio = match_num / len(A_features)
        pre_list.append(pre_ratio)
        recall_list.append(recall_ratio)
    return pre_list, recall_list
