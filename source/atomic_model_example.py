from argparse import ArgumentParser
from alignment.Registration import Registration_mask, cal_pdb_RMSD
from alignment.Utils import *
import os
from scipy.spatial import cKDTree


def chain_alignment(A_pcd, A_key_pcd, B_pcd, B_key_pcd, temp_dir, VOXEL_SIZE):
    A_key_feats = cal_SHOT(np.array(A_pcd.points), np.array(A_pcd.normals), temp_dir, np.array(A_key_pcd.points),
                           radius=VOXEL_SIZE * 7.0)
    B_key_feats = cal_SHOT(np.array(B_pcd.points), np.array(B_pcd.normals), temp_dir, np.array(B_key_pcd.points),
                           radius=VOXEL_SIZE * 7.0)

    A_min_bound = A_pcd.get_min_bound()
    A_max_bound = A_pcd.get_max_bound()
    A_dist_cor = A_max_bound - A_min_bound
    B_min_bound = B_pcd.get_min_bound()
    B_max_bound = B_pcd.get_max_bound()

    mask = {}
    radius = np.max(A_dist_cor) * 1.1 / 2
    center = B_min_bound + radius / 10
    terminal = B_max_bound
    step = int(radius // 2)  # translation interval of mask

    max_correspondence_dist = 2.0

    best_score = 0.0
    best_T = None

    for i in range(0, int(terminal[0]), step):
        for j in range(0, int(terminal[1]), step):
            for k in range(0, int(terminal[2]), step):
                temp_center = [center[0] + i, center[1] + j, center[2] + k]

                mask["name"] = "sphere"
                mask["center"] = temp_center
                mask["radius"] = radius

                final_T, score = Registration_mask(temp_dir, A_pcd, B_pcd, A_key_pcd, B_key_pcd, A_key_feats,
                                                   B_key_feats, mask,
                                                   max_correspondence_dist=max_correspondence_dist,
                                                   VOXEL_SIZE=VOXEL_SIZE, store_partial=False)
                if score > best_score:
                    best_score = score
                    best_T = final_T

    return best_T


def transform_pdb_chain(pdbFile, output_dir, T):
    if "cif" in pdbFile:
        interval = 9
        with open(output_dir, "w") as of:
            with open(pdbFile, 'r') as fr:
                for line in fr:
                    truncate_p = None
                    concat_p = None
                    if line.startswith('ATOM'):
                        if truncate_p is None:
                            space_count = 0
                            i = 0
                            pre = False  # judge space or not
                            while 1:
                                if line[i] == " ":
                                    pre = True
                                else:
                                    if pre:
                                        # print(i)
                                        space_count += 1
                                    pre = False
                                i += 1
                                if space_count == interval and truncate_p is None:
                                    truncate_p = i
                                if space_count == interval + 3:
                                    concat_p = i
                                    break

                        temp_atom_coord = line.strip().split()
                        atom_coord = [float(item) for item in temp_atom_coord[interval: interval + 3] if item != '']

                        transformed_atom_coord = np.matmul(T, np.array(atom_coord + [1.0]))
                        transformed_atom_line = ["%.3f" % i for i in list(transformed_atom_coord[:3])]
                        transformed_atom_line = [i + " " * (8 - len(i)) for i in transformed_atom_line]
                        transformed_line = "".join(transformed_atom_line)


                        transformed_line = line[:truncate_p - 1] + transformed_line + line[concat_p - 1:]
                        of.write(transformed_line)
                    else:
                        of.write(line)
    elif "pdb" in pdbFile:
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
    return



"""
shell script
"""
parser = ArgumentParser(description="An example of flexible fitting")
parser.add_argument("--stage", type=str, default="rigid", help="initial rigid or flexible fitting")
parser.add_argument("--visualize", action='store_true', default=False, help='point cloud visualization')

args = parser.parse_args()

stage = args.stage
visualize = args.visualize

data_dir = "../data/emd_4775_4776"
chain_A_points, chain_A_normals = load_sample("%s/chain_A.txt" % data_dir)
chain_A_pcd = o3d.geometry.PointCloud()
chain_A_pcd.points = o3d.utility.Vector3dVector(chain_A_points)
chain_A_pcd.normals = o3d.utility.Vector3dVector(chain_A_normals)
chain_A_key_pcd = load_xyz("%s/Points_chain_A_Key.xyz" % data_dir)
chain_B_points, chain_B_normals = load_sample("%s/chain_B.txt" % data_dir)
chain_B_pcd = o3d.geometry.PointCloud()
chain_B_pcd.points = o3d.utility.Vector3dVector(chain_B_points)
chain_B_pcd.normals = o3d.utility.Vector3dVector(chain_B_normals)
chain_B_key_pcd = load_xyz("%s/Points_chain_B_Key.xyz" % data_dir)

target_points, target_normals = load_sample("%s/emd_4776.txt" % data_dir)
target_pcd = o3d.geometry.PointCloud()
target_pcd.points = o3d.utility.Vector3dVector(target_points)
target_pcd.normals = o3d.utility.Vector3dVector(target_normals)
target_key_pcd = load_xyz("%s/Points_4776_Key.xyz" % data_dir)

temp_dir = "%s/temp" % data_dir
if not os.path.exists(temp_dir):
    os.mkdir(temp_dir)

if stage == "rigid":

    chain_A_T = chain_alignment(chain_A_pcd, chain_A_key_pcd, target_pcd, target_key_pcd, temp_dir, 2.0)
    np.save("%s/chain_A.npy" % data_dir, chain_A_T)
    transform_pdb_chain("%s/chain_A.pdb" % data_dir, "%s/chain_A_rigid.pdb" % data_dir, chain_A_T)

    chain_B_T = chain_alignment(chain_B_pcd, chain_B_key_pcd, target_pcd, target_key_pcd, temp_dir, 2.0)
    np.save("%s/chain_B.npy" % data_dir, chain_B_T)
    transform_pdb_chain("%s/chain_B.pdb" % data_dir, "%s/chain_B_rigid.pdb" % data_dir, chain_B_T)

    if visualize:
        chain_A_T = np.load("%s/chain_A.npy" % data_dir)
        chain_B_T = np.load("%s/chain_B.npy" % data_dir)
        chain_A_pcd.paint_uniform_color([1.0, 0, 0])
        chain_B_pcd.paint_uniform_color([0, 1.0, 0])
        target_pcd.paint_uniform_color([0, 0, 1.0])
        chain_A_pcd.transform(chain_A_T)
        chain_B_pcd.transform(chain_B_T)
        o3d.visualization.draw_geometries([chain_A_pcd, chain_B_pcd, target_pcd])


elif stage == "flexible":

    chain_A_T = np.load("%s/chain_A.npy" % data_dir)
    chain_B_T = np.load("%s/chain_B.npy" % data_dir)
    chain_A_pcd.transform(chain_A_T)
    chain_B_pcd.transform(chain_B_T)

    flexible_points_dir = "%s/output_flexible.txt" % data_dir
    flexible_points = np.loadtxt(flexible_points_dir)
    flexible_pcd = o3d.geometry.PointCloud()
    flexible_pcd.points = o3d.utility.Vector3dVector(flexible_points)

    index = 4576
    transformed_A_points = flexible_points[:index]
    transformed_B_points = flexible_points[index:]
    transformed_A_pcd = o3d.geometry.PointCloud()
    transformed_A_pcd.points = o3d.utility.Vector3dVector(transformed_A_points)
    transformed_B_pcd = o3d.geometry.PointCloud()
    transformed_B_pcd.points = o3d.utility.Vector3dVector(transformed_B_points)

    if visualize:
        transformed_A_pcd.paint_uniform_color([1.0, 0, 0])
        transformed_B_pcd.paint_uniform_color([0, 1.0, 0])
        o3d.visualization.draw_geometries([transformed_A_pcd, transformed_B_pcd, target_pcd])

    # A_indices = []
    # for i in range(len(transformed_A_points)):
    #     if transformed_A_points[i, 2] >= 150:
    #         A_indices.append(i)
    # B_indices = []
    # for i in range(len(transformed_B_points)):
    #     if transformed_B_points[i, 2] >= 150:
    #         B_indices.append(i)
    #
    # select_A_points = transformed_A_points[A_indices]
    # select_A_pcd = transformed_A_pcd.select_by_index(A_indices)
    # select_B_points = transformed_B_points[B_indices]
    # select_B_pcd = transformed_B_pcd.select_by_index(B_indices)
    # select_points = np.concatenate([select_A_points, select_B_points], axis=0)
    # # o3d.visualization.draw_geometries([select_A_pcd, select_B_pcd])
    #
    # raw_select_A_pcd = chain_A_pcd.select_by_index(A_indices)
    # raw_select_B_pcd = chain_B_pcd.select_by_index(B_indices)
    # raw_select_A_points = np.array(raw_select_A_pcd.points)
    # raw_select_B_points = np.array(raw_select_B_pcd.points)
    # raw_select_points = np.concatenate([raw_select_A_points, raw_select_B_points], axis=0)
    # # o3d.visualization.draw_geometries([raw_select_A_pcd, raw_select_B_pcd])

    chain_A_points = np.array(chain_A_pcd.points)
    chain_B_points = np.array(chain_B_pcd.points)
    total_points = np.concatenate([chain_A_points, chain_B_points], axis=0)

    point_tree = cKDTree(total_points)
    chain_A_aligned_pdb = "%s/chain_A_rigid.pdb" % data_dir
    chain_A_flexible_pdb = "%s/chain_A_flexible.pdb" % data_dir
    with open(chain_A_flexible_pdb, "w") as of:
        with open(chain_A_aligned_pdb, 'r') as fr:
            for line in fr:
                if line.startswith('ATOM'):
                    temp_atom_coord = line[27:55].strip().split()
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
                    atom_coord = np.array(atom_coord)
                    dists, nn_inds = point_tree.query(atom_coord, k=1, n_jobs=-1)
                    x_disp, y_disp, z_disp = flexible_points[nn_inds] - total_points[nn_inds]
                    transformed_atom_coord = atom_coord + np.array([x_disp, y_disp, z_disp])
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

    chain_B_aligned_pdb = "%s/chain_B_rigid.pdb" % data_dir
    chain_B_flexible_pdb = "%s/chain_B_flexible.pdb" % data_dir
    with open(chain_B_flexible_pdb, "w") as of:
        with open(chain_B_aligned_pdb, 'r') as fr:
            for line in fr:
                if line.startswith('ATOM'):
                    temp_atom_coord = line[27:55].strip().split()
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
                    atom_coord = np.array(atom_coord)
                    dists, nn_inds = point_tree.query(atom_coord, k=1, n_jobs=-1)
                    x_disp, y_disp, z_disp = flexible_points[nn_inds] - total_points[nn_inds]
                    transformed_atom_coord = atom_coord + np.array([x_disp, y_disp, z_disp])
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
