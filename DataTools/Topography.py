import array
import time

import numpy as np
import pandas as pd
import os
from DataFactory import CommonUtils, ROOTUtils
from FwdTools import points_columns, centers_columns, normal_columns, plane_equation_columns, intersection_columns, \
    path_length_columns, in_out_columns

from ROOT import TMath, TVector3
import logging

logger = logging.getLogger(__name__)


class Topography:
    def __init__(self, topo_fname: str, itv_x=50., itv_y=50., itv_z=50., topo_name=''):
        self._topo_fname = topo_fname
        self._df_topo = None
        self._itv_x = itv_x
        self._itv_y = itv_y
        self._itv_z = itv_z
        self._check_in_out = False
        self._search_all = False
        self._dic_tri_columns = {}
        self._np_tri = None
        self._topo_name = topo_name

    def get_topo_name(self):
        return self._topo_name

    def get_topo_fname(self):
        return self._topo_fname

    def read_topo_obj(self, check_in_out=False, search_all=False):
        self._search_all = search_all

        tmp_fname = self._topo_fname.replace('.obj', '_tmp.obj')
        os.system('echo "val1 val2 val3" > %s' % tmp_fname)
        command_str = " awk '{if($1 == \"v\") print $1,$2,$3,$4; if($1 == \"vn\" || $1 == \"f\") print $0}' %s >> %s" \
                      % (self._topo_fname, tmp_fname)
        os.system(command_str)

        df_whole = pd.read_csv(tmp_fname, sep=' ', index_col=0, low_memory=False)
        df_vertices = df_whole[df_whole.index == 'v']
        df_normals = df_whole[df_whole.index == 'vn']
        df_faces = df_whole[df_whole.index == 'f']

        df_v_id = self._split_faces(df_faces)
        df_points = self._find_vertices(df_faces_id=df_v_id, df_vertices=df_vertices)
        df_centers = self._calc_center_points(df_points)

        df_plane = self._calc_triangle_plane(df_points)

        df_all = pd.concat([df_points, df_centers, df_plane], axis=1)
        df_all = df_all.reset_index(drop=True)

        # get df_normals
        if check_in_out:
            self._check_in_out = True
            df_vn_id = self._split_faces(df_faces)
            df_normal_points = self._find_vn(df_vn_id=df_vn_id, df_normals=df_normals)
            df_normal_aver = self._calc_center_points(df_normal_points)
            df_normal_aver.columns = normal_columns
            df_normal_aver = df_normal_aver.reset_index(drop=True)

            df_all = pd.concat([df_all, df_normal_aver], axis=1)

        self._df_topo = df_all

        os.remove(tmp_fname)

    def write_topo_csv(self, out_fname=''):
        if not out_fname.__contains__('csv'):
            out_fname = self._topo_fname.replace('.obj', '.csv')
        if self._df_topo is not None:
            self._df_topo.to_csv(out_fname)
        else:
            logger.error('self._df_topo is not None; please read topography first.')

    def read_topo_csv(self):
        if not self._topo_fname.__contains__('.csv'):
            raise Exception('[read_topo_csv] Cannot read %s with pandas.read_csv, '
                            'please use [read_topo_obj] if it is a obj file.' % self._topo_fname)
        self._df_topo = pd.read_csv(self._topo_fname)

    @classmethod
    def _split_faces(cls, df_faces):
        faces1 = df_faces['val1'].str.split('/', expand=True)
        faces1 = faces1.replace('', 0)
        faces1.fillna(value=np.nan, inplace=True)
        faces1 = faces1.replace(np.nan, 0)
        faces1 = faces1[faces1.iloc[:, 0].str.isdigit()].astype('int') - 1

        faces2 = df_faces['val2'].str.split('/', expand=True)
        faces2 = faces2.replace('', 0)
        faces2.fillna(value=np.nan, inplace=True)
        faces2 = faces2.replace(np.nan, 0)
        faces2 = faces2[faces2.iloc[:, 0].str.isdigit()].astype('int') - 1

        faces3 = df_faces['val3'].str.split('/', expand=True)
        faces3 = faces3.replace('', 0)
        faces3.fillna(value=np.nan, inplace=True)
        faces3 = faces3.replace(np.nan, 0)
        faces3 = faces3[faces3.iloc[:, 0].str.isdigit()].astype('int') - 1

        df_split = pd.concat(
            [faces1[0], faces2[0], faces3[0], faces1.iloc[:, -1:], faces2.iloc[:, -1:], faces3.iloc[:, -1:]], axis=1)
        df_split.columns = ['v1', 'v2', 'v3', 'vn1', 'vn2', 'vn3']
        return df_split

    @classmethod
    def _find_vertices(cls, df_faces_id, df_vertices):
        p1 = df_vertices.iloc[df_faces_id['v1']].astype('float')
        p2 = df_vertices.iloc[df_faces_id['v2']].astype('float')
        p3 = df_vertices.iloc[df_faces_id['v3']].astype('float')
        df_points = pd.concat([p1, p2, p3], axis=1)
        df_points.columns = points_columns
        return df_points

    @classmethod
    def _find_vn(cls, df_vn_id, df_normals):
        p1 = df_normals.iloc[df_vn_id['vn1']].astype('float')
        p2 = df_normals.iloc[df_vn_id['vn2']].astype('float')
        p3 = df_normals.iloc[df_vn_id['vn3']].astype('float')
        df_points = pd.concat([p1, p2, p3], axis=1)
        df_points.columns = points_columns
        return df_points

    def _calc_center_points(self, df_points):
        # ctr_x = df_points.apply(
        #     lambda df: self._calc_3average(df[points_columns[0]], df[points_columns[3]], df[points_columns[6]]), axis=1)
        # ctr_y = df_points.apply(
        #     lambda df: self._calc_3average(df[points_columns[1]], df[points_columns[4]], df[points_columns[7]]), axis=1)
        # ctr_z = df_points.apply(
        #     lambda df: self._calc_3average(df[points_columns[2]], df[points_columns[5]], df[points_columns[8]]), axis=1)

        ctr_x = (df_points[points_columns[0]] + df_points[points_columns[3]] + df_points[points_columns[6]]) / 3
        ctr_y = (df_points[points_columns[1]] + df_points[points_columns[4]] + df_points[points_columns[7]]) / 3
        ctr_z = (df_points[points_columns[2]] + df_points[points_columns[5]] + df_points[points_columns[8]]) / 3

        df_center = pd.concat([ctr_x, ctr_y, ctr_z], axis=1)
        df_center.columns = centers_columns
        return df_center

    @classmethod
    def _calc_3average(cls, col1, col2, col3):
        return (col1 + col2 + col3) / 3

    def _calc_triangle_plane(self, df_tri):
        x1 = df_tri[points_columns[0]]
        y1 = df_tri[points_columns[1]]
        z1 = df_tri[points_columns[2]]
        x2 = df_tri[points_columns[3]]
        y2 = df_tri[points_columns[4]]
        z2 = df_tri[points_columns[5]]
        x3 = df_tri[points_columns[6]]
        y3 = df_tri[points_columns[7]]
        z3 = df_tri[points_columns[8]]

        # normal vector(a, b, c)
        a = (y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1)
        b = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1)
        c = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)

        c = c.replace(0., 1e-15)

        d = a * x1 + b * y1 + c * z1
        a0 = d / c
        a1 = -1.0 * a / c
        a2 = -1.0 * b / c

        df_plane = pd.concat([a0, a1, a2], axis=1)
        df_plane.columns = plane_equation_columns
        return df_plane

    @classmethod
    def _select_near_triangles(cls, t_arr, t_dic, px, py, pz=None, itv_x=50., itv_y=50., itv_z=50.):
        if pz is not None:
            _near_index = np.where((abs(t_arr[:, t_dic[centers_columns[0]]] - px) < itv_x)
                                   & (abs(t_arr[:, t_dic[centers_columns[1]]] - py) < itv_y)
                                   & (abs(t_arr[:, t_dic[centers_columns[2]]] - pz) < itv_z))
        else:
            _near_index = np.where((abs(t_arr[:, t_dic[centers_columns[0]]] - px) < itv_x)
                                   & (abs(t_arr[:, t_dic[centers_columns[1]]] - py) < itv_y))
        return _near_index

    @classmethod
    def _check_theta_phi(cls, theta, phi):
        while phi > TMath.Pi():
            phi = phi - TMath.Pi() * 2
        while phi < -1. * TMath.Pi():
            phi = phi + TMath.Pi() * 2
        if not 0 <= theta <= TMath.Pi() / 2:
            raise ValueError('theta < 0 or theta > Pi() : theta = %.3f' % theta)
        return theta, phi

    @classmethod
    def _two_bigger(cls, x1, x2, x3):
        lis = [x1, x2, x3]
        lis.sort()
        return lis[2], lis[1]

    def _intersection_in_triangle(self, tri_index, xi, yi, zi, z0):
        if zi < z0:
            return False

        x1 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[0]]]
        y1 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[1]]]
        z1 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[2]]]
        x2 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[3]]]
        y2 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[4]]]
        z2 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[5]]]
        x3 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[6]]]
        y3 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[7]]]
        z3 = self._np_tri[tri_index, self._dic_tri_columns[points_columns[8]]]

        # size of the triangle
        v_12 = TVector3(x2 - x1, y2 - y1, z2 - z1)
        v_13 = TVector3(x3 - x1, y3 - y1, z3 - z1)
        v_23 = TVector3(x3 - x2, y3 - y2, z3 - z2)
        v_1i = TVector3(xi - x1, yi - y1, zi - z1)
        v_2i = TVector3(xi - x2, yi - y2, zi - z2)
        # v_3i = TVector3(xi - x3, yi - y3, zi - z3)

        tri_size = (v_12.Cross(v_13)).Mag()
        tri_int_1 = (v_12.Cross(v_1i)).Mag()
        tri_int_2 = (v_13.Cross(v_1i)).Mag()
        tri_int_3 = (v_23.Cross(v_2i)).Mag()

        tri_size_1, tri_size_2 = self._two_bigger(tri_int_1, tri_int_2, tri_int_3)

        size_two = tri_size_1 + tri_size_2 - tri_size

        if size_two > 1e-15:
            return False
        return True

    def calc_intersections(self, lis_pos, theta, phi, last_inter=None, lis_xy=None):
        # give the real intersections
        _topo_near_index = None
        if last_inter is not None:
            for i in range(len(last_inter)):
                px = last_inter[i, 0]
                py = last_inter[i, 1]
                pz = last_inter[i, 2]
                if _topo_near_index is None:
                    _topo_near_index = self._select_near_triangles(self._np_tri, self._dic_tri_columns, px, py, pz,
                                                                   self._itv_x, self._itv_y, self._itv_z)
                else:
                    tmp = self._select_near_triangles(self._np_tri, self._dic_tri_columns, px, py, pz,
                                                      self._itv_x, self._itv_y, self._itv_z)
                    _topo_near_index = np.hstack((_topo_near_index, tmp))
            _topo_near_index = np.unique(_topo_near_index)
        elif lis_xy is not None:
            px = lis_xy[0]
            py = lis_xy[1]
            _topo_near_index = self._select_near_triangles(self._np_tri, self._dic_tri_columns, px, py,
                                                           itv_x=self._itv_x, itv_y=self._itv_y, itv_z=self._itv_z)
        else:
            _topo_near_index = self._df_topo.index

        _topo_near = self._df_topo.loc[_topo_near_index]
        _np_topo_near = self._np_tri[_topo_near_index]
        _dic_topo_near = self._dic_tri_columns

        theta, phi = self._check_theta_phi(theta, phi)
        v_ray = TVector3(0, 0, 1)
        v_ray.SetMagThetaPhi(1, theta, phi)
        t1 = v_ray.X()
        t2 = v_ray.Y()
        t3 = v_ray.Z()

        det_x = lis_pos[0]
        det_y = lis_pos[1]
        det_z = lis_pos[2]

        # a0 = _topo_near[plane_equation_columns[0]]
        # a1 = _topo_near[plane_equation_columns[1]]
        # a2 = _topo_near[plane_equation_columns[2]]

        a0 = _np_topo_near[:, _dic_topo_near[plane_equation_columns[0]]]
        a1 = _np_topo_near[:, _dic_topo_near[plane_equation_columns[1]]]
        a2 = _np_topo_near[:, _dic_topo_near[plane_equation_columns[2]]]

        den = a1 * t1 + a2 * t2 - t3
        den = np.where(den == 0., 1e-15, den)

        xi = -1. * (a0 * t1 - a2 * t2 * det_x + t3 * det_x + a2 * t1 * det_y - t1 * det_z) / den
        yi = -1. * (a0 * t2 + a1 * t2 * det_x - a1 * t1 * det_y + t3 * det_y - t2 * det_z) / den
        zi = -1. * (a0 * t3 + a1 * t3 * det_x + a2 * t3 * det_y - a1 * t1 * det_z - a2 * t2 * det_z) / den

        lis_in_tri_index = []
        _inter = None
        for i in range(len(_topo_near)):
            ind = _topo_near.index[i]
            if self._intersection_in_triangle(ind, xi[i], yi[i], zi[i], det_z):
                lis_in_tri_index.append(ind)
                if _inter is None:
                    _inter = np.array([xi[i], yi[i], zi[i]])
                else:
                    _inter = np.vstack((_inter, np.array([xi[i], yi[i], zi[i]])))
        # index_in_tri = self._intersection_in_triangle(tmp_topo_near_inter.index, det_z)

        df_inter_tri = _topo_near.loc[lis_in_tri_index]

        return df_inter_tri, _inter

    def _check_pass_in_out_pairs(self, theta, phi, df_inter_tri, arr_inter, det):
        if len(df_inter_tri) % 2 == 0:
            return True
        if self._check_in_out:
            _in = df_inter_tri[df_inter_tri[in_out_columns] == True]
            _out = df_inter_tri[df_inter_tri[in_out_columns] == False]
            logger.warning(
                'Intersections are not in pairs for %s, %s.\n'
                'Ray direction (%.2f, %.2f)\t pass in: %d, pass out: %d\t\t%s%s'
                % (self._topo_name, det.get_det_name(), theta, phi, len(_in), len(_out), df_inter_tri.__str__(), arr_inter.__str__()))
        else:
            logger.warning('Intersections are not in pairs for %s, %s.\n'
                           'Ray direction (%.2f, %.2f)\n%s%s'
                           % (self._topo_name, det.get_det_name(), theta, phi, df_inter_tri.__str__(), arr_inter.__str__()))
        return False

    @classmethod
    def _check_intersections_in_out(cls, df_inter_tri, theta, phi):
        v_ray = TVector3(0, 0, 1)
        v_ray.SetMagThetaPhi(1, theta, phi)

        t1 = v_ray.X()
        t2 = v_ray.Y()
        t3 = v_ray.Z()

        _normal = df_inter_tri[normal_columns]

        _dot = _normal[normal_columns[0]] * t1 + _normal[normal_columns[1]] * t2 + _normal[normal_columns[2]] * t3

        _in = _dot[0].apply(lambda x: True if x < 0 else False)
        _in.columns = in_out_columns
        return _in

    @classmethod
    def _calc_distance(cls, df_inter, det_x, det_y, det_z):
        del_x = df_inter[intersection_columns[0]] - det_x
        del_y = df_inter[intersection_columns[1]] - det_y
        del_z = df_inter[intersection_columns[2]] - det_z
        path_length = TMath.Sqrt(del_x * del_x + del_y * del_y + del_z * del_z)
        return path_length

    def _calc_path_length_for_one_ray(self, lis_pos, arr_inter, np_in):
        det_x = lis_pos[0]
        det_y = lis_pos[1]
        det_z = lis_pos[2]

        path_length_square = (arr_inter[:, 0] - det_x) * (arr_inter[:, 0] - det_x) \
                             + (arr_inter[:, 1] - det_y) * (arr_inter[:, 1] - det_y) \
                             + (arr_inter[:, 2] - det_z) * (arr_inter[:, 2] - det_z)
        path_length = np.sqrt(path_length_square)

        # df_inter = pd.concat([df_inter, path_length], axis=1)
        # df_inter.columns = intersection_columns + path_length_columns

        lpath = 0.
        if self._check_in_out:
            np_in = np_in.reshape(len(np_in))
            _in = arr_inter[np_in]
            _in = np.sort(_in)
            _out = arr_inter[~np_in]
            _out = np.sort(_out)
            if not len(_in) == len(_out):
                logger.warning('Number of intersections pass in dose not equal to that pass out.\n'
                               'pass in: %d, pass out: %d\n'
                               'The path length would be compute without in_out information' % (len(_in), len(_out)))
            else:
                for i in range(len(_in)):
                    lpath = lpath + (_out[i] - _in[i])
                return lpath

        lpath = 0.
        path_length = np.sort(path_length)
        for i in range(len(path_length)):
            if i % 2 == 0:
                sig = -1.
            else:
                sig = 1.
            lpath = lpath + path_length[i] * sig

        if lpath < 0:
            print("ATTENTION!!!  lpath < 0  for %s  ================\n" % self._topo_name,
                  "intersections are:\n", arr_inter,
                  "\npath_lengths are :\n", path_length)
        return lpath

    def calc_path_length(self, det, h_accpt, opt: str, calc_all=False):
        self._np_tri, self._dic_tri_columns = CommonUtils.cvt_dataframe_to_numpy(self._df_topo)

        if opt.lower().__contains__('phi'):
            mode = 0
        elif opt.lower().__contains__('thetax'):
            mode = 1
        else:
            raise ValueError('opt should be thetaphi pr thetaxy')

        h_ltopo = h_accpt.Clone()
        h_ltopo.Reset("ICES")

        rot_phi = det.get_rotAngle()[2]

        xbins = h_accpt.GetNbinsX()
        ybins = h_accpt.GetNbinsY()

        # count time
        time_start = time.time()
        time_end = time.time()

        for i_xbin in range(1, xbins + 1):

            xx = ROOTUtils.get_bin_position(axis=h_accpt.GetXaxis(), ibin=i_xbin, opt="c")
            t_inter = None
            for i_ybin in range(1, ybins + 1):
                yy = ROOTUtils.get_bin_position(axis=h_accpt.GetYaxis(), ibin=i_ybin, opt="c")

                ibin = h_accpt.GetBin(i_xbin, i_ybin)
                accpt = h_accpt.GetBinContent(ibin)
                if accpt < 1e-15 and not calc_all:
                    continue

                if mode == 0:
                    theta = yy
                    phi = xx + rot_phi
                else:
                    vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
                    theta = vec.Theta()
                    phi = vec.Phi() + rot_phi

                theta, phi = self._check_theta_phi(theta, phi)

                # 射线是否穿过立方体，提前终结

                # calculate intersections
                if t_inter is None:
                    df_inter_tri, _inter = self.calc_intersections(det.get_positions(), theta, phi)
                else:
                    df_inter_tri, _inter = self.calc_intersections(det.get_positions(), theta, phi, t_inter)
                    if self._search_all and _inter is None:
                        df_inter_tri, _inter = self.calc_intersections(det.get_positions(), theta, phi)
                    elif self._search_all and _inter is not None and len(_inter) % 2 != 0:
                        df_inter_tri, _inter = self.calc_intersections(det.get_positions(), theta, phi)

                # df_inter_tri, _inter = self.calc_intersections(det.get_positions(), theta, phi)

                if len(df_inter_tri) == 0:
                    h_ltopo.SetBinContent(ibin, 1e-20)
                    t_inter = None
                    continue

                # add a column of pass in or out if check
                np_in = None
                if self._check_in_out:
                    _in = self._check_intersections_in_out(df_inter_tri, theta, phi)
                    df_inter_tri = pd.concat([df_inter_tri, _in], axis=1)
                    np_in = _in.to_numpy()

                # check if intersections are in pairs (one in for every out), if not discard
                _if_inter_pairs = self._check_pass_in_out_pairs(df_inter_tri=df_inter_tri,
                                                                theta=theta, phi=phi, arr_inter=_inter, det=det)
                if not _if_inter_pairs:
                    # logger.warning('Intersections are not in pair. %s  ----------' %  % self._topo_name)
                    continue

                ltopo = self._calc_path_length_for_one_ray(det.get_positions(), _inter, np_in)

                h_ltopo.SetBinContent(ibin, ltopo)
                h_ltopo.SetBinError(ibin, 0)

                t_inter = _inter

        time_end = time.time()
        print('time cost', time_end - time_start, 's')

        return h_ltopo

    def print_intersection(self, bottom_fname, out_inter_fname, out_obs_fname):
        """
        print_intersection for FillTopoDensity
        @param bottom_fname: 输入底坐标文件名
        @param out_inter_fname: 输出交点文件名
        @param out_obs_fname: 输出obs文件名
        @return:
        """
        self._np_tri, self._dic_tri_columns = CommonUtils.cvt_dataframe_to_numpy(self._df_topo)

        theta = phi = 0.

        fout_inter = open(out_inter_fname, "w")
        fout_obs = open(out_obs_fname, "w")

        fout_obs.write('directions:\n')
        i = 1
        for line in open(bottom_fname):
            det_x, det_y, det_z = self._get_pos_from_txt(line)
            df_inter_tri, _inter = self.calc_intersections([det_x, det_y, det_z], theta, phi, lis_xy=[det_x, det_y])

            obs_print = '%d %.2f %.2f %.2f %.2f %.2f -10 -10 -10' % (i, det_x, det_y, det_z, theta, phi)
            str_print = ''

            if _inter is None:
                continue

            for i_arr in _inter.flatten():
                str_print = str_print + '%f ' % i_arr

            fout_inter.write(str_print)
            fout_inter.write('\n')
            fout_obs.write(obs_print)
            fout_obs.write('\n')

            i += 1

        return 0

    @classmethod
    def _get_pos_from_txt(cls, line):
        lis_line = line.split()
        det_x, det_y, det_z = float(lis_line[0]), float(lis_line[1]), float(lis_line[2])
        return det_x, det_y, det_z
