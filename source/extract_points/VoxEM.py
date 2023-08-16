from numpy.core.fromnumeric import shape
import tqdm
import os
from functools import partial
import struct
import numpy as np
import multiprocessing
from sklearn.cluster import DBSCAN
import torch
from .Supporting import *

# =============================
# Voxel Operation Class
# =============================

class VoxEM():

    # =============================
    # Initialise Workspaces
    # =============================

    def __init__(self):

        # Initiate Workspaces
        self.point_workspace = {}
        self.graph_workspace = {}
        self.voxel_workspace = {}
        self.segment_workspace = {}
        # Statistics workspace can be derived entirely from any Voxel workspace
        self.statistics_workspace = {}
        # Description workspace can only be created from reading a proper mrc file.
        self.description_workspace = {}

        return

    # =============================
    # I/O File managementI/O文件管理
    # =============================

    # 1. Read
    def IO_ReadMrc(self, mrc_inputname='XXX.mrc', voxel_outputname='Default', grid_outputname=None,
                   description_outputname='Default', statistics_outputname='Default', nonzerostatistics=True):

        with open(mrc_inputname, 'rb') as fin:

            MRCdata = fin.read()

            # Dimension of column, row and section in unit cell
            nx, ny, nz = struct.unpack_from('<iii', MRCdata, 0)

            # Mode
            mode = struct.unpack_from('<i', MRCdata, 12)

            # Start
            xs, ys, zs = struct.unpack_from('<iii', MRCdata, 16)
            start = (xs, ys, zs)

            # Sampling along axes
            mx, my, mz = struct.unpack_from('<iii', MRCdata, 28)
            unitcelldimension = (mx, my, mz)

            # Cell dimension in angstrom
            X_angstrom, Y_angstrom, Z_angstrom = struct.unpack_from('<fff', MRCdata, 40)
            angstrom = (X_angstrom, Y_angstrom, Z_angstrom)

            # Cell angle in degree
            X_degree, Y_degree, Z_degree = struct.unpack_from('<fff', MRCdata, 52)
            angle = (X_degree, Y_degree, Z_degree)

            # Axis
            MAPC, MAPR, MAPS = struct.unpack_from('<iii', MRCdata, 64)
            axis = (MAPC, MAPR, MAPS)

            # Misc
            ISPG, NSYMBT = struct.unpack_from('<ii', MRCdata, 88)

            # Extra
            EXTRA = struct.unpack_from('<' + 'f' * 12, MRCdata, 96)

            # Origin
            X_origin, Y_origin, Z_origin = struct.unpack_from('<fff', MRCdata, 196)
            origin = (X_origin, Y_origin, Z_origin)

            # character string 'MAP ' to identify file type
            MAP_String = struct.unpack_from('<ssss', MRCdata, 208)

            # Machine Stamp
            MACHST = struct.unpack_from('<BBBB', MRCdata, 212)

            # Number of labels in use
            NLABL = struct.unpack_from('<i', MRCdata, 220)[0]

            # Density
            # The original voxel from the file is stored in self.Voxel
            fin.seek(1024, os.SEEK_SET)
            Voxel = np.fromfile(file=fin, dtype=np.dtype(np.float32)).reshape((nx, ny, nz), order='F')
            fin.close()

            # Workspaces
            self.voxel_workspace[voxel_outputname] = Voxel

            if grid_outputname is not None:
                self.Point_Create_Grid(outputname=grid_outputname)

            if statistics_outputname is not None:
                self.IO_Statistics(voxel_inputname=voxel_outputname, statistics_outputname=statistics_outputname,
                                   nonzerostatistics=nonzerostatistics)

            if description_outputname is not None:
                self.description_workspace[description_outputname] = {"Mode": mode, "Start": start,
                                                                      "UnitCellDim": unitcelldimension,
                                                                      "Angstrom": angstrom, "Angle": angle,
                                                                      "Axis": axis, "Extra": EXTRA, "Ispg": ISPG,
                                                                      "Nsymbt": NSYMBT, "MapString": MAP_String,
                                                                      "Origin": origin, "MachineString": MACHST,
                                                                      "Label": NLABL}

        print("IO_ReadMrc finished")

    def IO_Statistics(self, voxel_inputname='Voxel', statistics_outputname='Default', nonzerostatistics=True):

        vox = self.voxel_workspace[voxel_inputname]

        if nonzerostatistics:
            RMS = np.std(vox[vox > 0.0])
            dmin = np.min(vox[vox > 0.0])
            dmax = np.max(vox[vox > 0.0])
            dmean = np.mean(vox[vox > 0.0])
        else:
            RMS = np.std(vox)
            dmin = np.min(vox)
            dmax = np.max(vox)
            dmean = np.mean(vox)

        # Update Dimension
        dimension = vox.shape
        self.statistics_workspace[statistics_outputname] = {'Rms': RMS, 'Min': dmin, 'Max': dmax, 'Mean': dmean,
                                                            'Dim': dimension}

        print("Update_Statistics finished")
        # print("Rms: %.4f, Min: %.4f, Max: %.4f, Mean: %.4f, Dim: %d*%d*%d" % (RMS, dmin, dmax, dmean, dimension[0], dimension[1], dimension[2]))

        return

    def IO_WriteXYZ(self, point_inputname='Voxel', file_outputname='hello.xyz', atom_name='H', shifting=None):

        Y = self.point_workspace[point_inputname]
        if shifting is not None:
            Y = Y + shifting.reshape(3, 1)
        Y = Y.reshape((3, np.prod(Y.shape[1:])))
        XYZ(Y.T.tolist(), '%s' % (atom_name), '%s' % (file_outputname))
        return

    def Voxel_Prune_RangeZero(self, upperbound=None, lowerbound=None, inputname='Voxel', outputname='Voxel'):

        vox = self.voxel_workspace[inputname]

        if lowerbound != None:
            vox[vox <= lowerbound] = 0.0

        if upperbound != None:
            vox[vox >= upperbound] = 0.0

        self.voxel_workspace[outputname] = vox

        print("Voxel_Prune_RangeZero finished")
        return

    # =========================
    # Point Routine
    # =========================

    def Point_Create_Grid(self, voxel_inputname='Default', point_outputname='Default', description_inputname='Default'):

        # Create grid points from voxel

        # Interval between grid
        # interval = np.array(self.angstrom).astype(np.float32) / np.array(self.Voxel.shape).astype(np.float32)
        interval = np.array(self.description_workspace[description_inputname]['Angstrom']).astype(
            np.float32) / np.array(self.voxel_workspace[voxel_inputname].shape).astype(np.float32)

        point_scatter = np.indices(self.voxel_workspace[voxel_inputname].shape).astype(np.float)
        point_scatter[0, :] = point_scatter[0, :] * interval[0]
        point_scatter[1, :] = point_scatter[1, :] * interval[0]
        point_scatter[2, :] = point_scatter[2, :] * interval[0]

        self.point_workspace[point_outputname] = point_scatter.astype(np.float32)

        return

    def Point_Create_Meanshift_sample(self, lower_bound=0.05, window=7, bandwidth=1.0, iteration=800, step_size=0.01,
                                      convergence=0.000187, voxel_inputname='Voxel', point_outputname='Meanshift',
                                      point_inputname=None, description_inputname='Default'):

        if window % 2 == 0:
            print("User should consider a odd window, which centers the distance window!")

        # Interval between grid
        interval = np.array(self.description_workspace[description_inputname]['Angstrom']).astype(
            np.float32) / np.array(self.voxel_workspace[voxel_inputname].shape).astype(np.float32)
        # interval = np.array([5.0, 5.0, 5.0]).astype(np.float32)

        # Initialisation
        self.Point_Create_Grid()
        grid = self.point_workspace['Default']

        # print(grid.shape)

        vox = self.voxel_workspace[voxel_inputname]

        # print(np.count_nonzero(vox))

        # Background field P
        # Q is the voxel value X is the corresponding coordinates of grid on the voxel
        Q = torch.tensor(vox).cuda()
        Q = Q.repeat(3, 1, 1, 1)
        X = torch.tensor(grid).cuda()
        P = Q * X
        P = P.view(3, 1, *vox.shape)
        Q = Q.view(3, 1, *vox.shape)
        X = X.view(3, 1, *vox.shape)

        # A distance window matrix W to convolve with
        w_x, w_y, w_z = np.mgrid[0: window, 0: window, 0: window]
        w_x = w_x.astype(np.float32) * interval[0]
        w_y = w_y.astype(np.float32) * interval[1]
        w_z = w_z.astype(np.float32) * interval[2]
        centroid = np.array([int(window / 2)] * 3).astype(np.float32) * interval
        w_x = (w_x - centroid[0]) ** 2
        w_y = (w_y - centroid[1]) ** 2
        w_z = (w_z - centroid[2]) ** 2
        W = np.sqrt(w_z + w_y + w_x)

        # Gaussian kernel
        # W = np.exp(-1.5 * ((W / bandwidth)) ** 2) / (bandwidth * math.sqrt(2 * math.pi))
        W = np.exp(-1.5 * ((W / bandwidth)) ** 2)
        W = torch.tensor(W).cuda()
        W = W.repeat(1, 1, 1, 1, 1)

        # Precaclulate Y_i+1 for each grid
        nominator = torch.nn.functional.conv3d(P, W, bias=None, stride=1, padding=int(window / 2), dilation=1, groups=1)
        denominator = torch.nn.functional.conv3d(Q, W, bias=None, stride=1, padding=int(window / 2), dilation=1,
                                                 groups=1)
        denominator += 0.000001
        denominator = 1.000 / denominator

        del P, Q, W, X

        C = nominator * denominator
        C = C[:, 0, :, :, :]
        C = np.array(C.cpu())

        Y = self.point_workspace[point_inputname]

        i = 0
        Y_diff_magnitude = vox.shape[0] ** 3
        while i <= iteration and Y_diff_magnitude >= convergence:
            Y_indice = (Y - interval.reshape(3, 1)) / interval.reshape(3, 1) + 1.0
            Y_indice = np.around(Y_indice).astype(int)
            Y_proposed = C[:, Y_indice[0], Y_indice[1], Y_indice[2]]
            Y_proposed_zero = np.where(np.sum(Y_proposed ** 2, axis=0) < lower_bound ** 2)

            Y_diff = Y - Y_proposed
            Y_diff = Y_diff.T
            Y_diff[Y_proposed_zero[0]] *= 0.0
            Y_diff = Y_diff.T

            Y_diff_magnitude = np.mean(np.sum(Y_diff ** 2, axis=0))

            if i % 50 == 0:
                print("Meanshift: iteration %s convergence %s" % (i, Y_diff_magnitude))

            Y = Y - (Y_diff * step_size)
            # print(Y[0][1227197] / interval[0])
            i += 1

        del C

        self.point_workspace[point_outputname] = Y
        print("Meanshift Completed.")

        return

    def Point_Transform_DBSCAN(self, distance_tolerance=1.9, clustersize_tolerance=13, showtopclusters=5,
                               inputname='Voxel', outputname='DBSCAN', centroid_only=False, n_job=50):

        Y = self.point_workspace[inputname]
        Y = Y.reshape((3, np.prod(Y.shape[1:])))
        Y_size = Y.T.shape[0]
        Y_transpose = Y.T

        clustering = DBSCAN(eps=distance_tolerance, min_samples=clustersize_tolerance).fit(Y.T)
        self.point_workspace["clusters"] = len(set(clustering.labels_))
        print("DBSCan found %s clusters" % self.point_workspace["clusters"])

        if centroid_only:
            cluster_centroid = []

            pool = multiprocessing.Pool(n_job)
            Misc_Centroid_Calculation_partial = partial(Misc_Centroid_Calculation, Y=Y_transpose,
                                                        labels=clustering.labels_)

            cluster_centroid = pool.map(Misc_Centroid_Calculation_partial,
                                        [i for i in tqdm.tqdm(sorted(set(clustering.labels_)))])
            pool.close()

            self.point_workspace[outputname] = np.array(cluster_centroid).T

            return

        # most_common = Counter(clustering.labels_).most_common(showtopclusters)

        self.point_workspace[outputname] = Y

        print("%s grids were reduced to %s grid after this step" % (Y_size, self.point_workspace[outputname].shape[1]))

        return

    def txt2pcd(self, inputname='PCA', file_outputname='%s/FPS_8724.xyz' % ('Storage_Xyz'), atom_name='H'):
        Y = self.point_workspace[inputname]
        Y = Y.reshape((3, np.prod(Y.shape[1:])))
        Y_transpose = Y.T
        Y_size = Y.T.shape[0]
        file_to_write = []

        with open(file_outputname, 'w') as file_to_write:
            file_to_write.writelines("# .PCD v0.7 - Point Cloud Data file format\n")
            file_to_write.writelines("VERSION 0.7\n")
            file_to_write.writelines("FIELDS x y z\n")
            file_to_write.writelines("SIZE 4 4 4\n")
            file_to_write.writelines("TYPE F F F\n")
            file_to_write.writelines("COUNT 1 1 1\n")
            file_to_write.writelines("WIDTH " + str(Y_size) + "\n")
            file_to_write.writelines("HEIGHT 1\n")
            file_to_write.writelines("VIEWPOINT 0 0 0 1 0 0 0\n")
            file_to_write.writelines("POINTS " + str(Y_size) + "\n")
            file_to_write.writelines("DATA ascii\n")

            '''for i in range(len(xlist)):
                file_to_write.writelines(str(xlist[i]) + " " + str(ylist[i]) + " " + str(zlist[i]) + "\n")'''

        XYZ_TXT(Y.T.tolist(), '%s' % (file_outputname))
        return
