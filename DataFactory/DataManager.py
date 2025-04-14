import logging
import os

import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count, Process
import time, psutil

from scipy.stats import expon

from DataFactory import Setting, PhysModel, FluxModel, CommonUtils, ROOTUtils
from DataTools import Detector, AngleMesh, Topography

from FwdTools import Topo2Ltopo, Ltopo2Leff, Leff2Res, ResolutionTools
from DataAnalysisTools import MuonStatistics, GetRatioLeff, PhysFitting, \
    DatabaseInterface, ComparingHists, RayAnalyzer, TimeSeriesTools

from DataAnalysisTools.DataframeProcessor import DataframeProcessor
from DataAnalysisTools.DataframeAnalyzer import DataframeAnalyzer
from DataAnalysisTools.DataFormatTransformer import DataFormatTransformer
from DataAnalysisTools import CalculatePositionsDF
from DataTools.PrintTObject import PrintTObject
from DataTools import PaperDraw

from ROOT import TString, TFile

logger = logging.getLogger(__name__)


class DataManager:
    def __init__(self, setting: Setting):
        self._setting = setting

        self._bins_dic = self._setting.get_bins()
        """
        classes and tools
        """

    def _make_dir_if_not_exists(self, _dir):
        if not CommonUtils.path_isabs(_dir):
            _dir = os.path.join(self._setting.get_global_dir(), _dir)
        if not os.path.exists(_dir):
            os.mkdir(_dir)

    @classmethod
    def _print_root_in_pdf(cls, out_fname, lis_draw_option):
        pdf_name = out_fname.replace('.root', '.pdf')
        _print_Tobject = PrintTObject()

        dic_draw_option = {}
        if not lis_draw_option == '':
            for i_draw_opt in range(len(lis_draw_option)):
                [class_name, draw_option] = lis_draw_option[i_draw_opt].split(':')
                dic_draw_option[class_name] = draw_option
            _print_Tobject.PrintTFile2PDF(out_fname, pdf_name, dic_draw_option)
        else:
            _print_Tobject.PrintTFile2PDF(out_fname, pdf_name)

    def Topo2Ltopo(self, kwargs):
        """
        Forward modeling step 1
        @param kwargs: OutDir; Input.Pool; Input.DetNames; Input.TopoNames (topography (obj file))
        @return:
        """
        logger.info('Run Topo2Ltopo')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs in the config file
        n_pools = int(kwargs.get('Input.Pool', 4))
        if n_pools == '':
            n_pools = 4
        det_names = kwargs.get('Input.DetNames')
        if det_names is None or not isinstance(det_names, list):
            raise Exception('Topo2Ltopo Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(det_names)

        topo_keys = kwargs.get('Input.TopoNames')
        if topo_keys is None:
            raise Exception('Topo2Ltopo Input.TopoNames is None.')

        calc_all_str = kwargs.get('Input.CalcAll', '')
        if calc_all_str.lower() != 'true':
            calc_all = False
        else:
            calc_all = True

        # prepare topography info and object
        dic_topo_info = self._setting.get_topos_info()
        fin_dir = os.path.join(self._setting.get_global_dir(), dic_topo_info.get('Directory', ''))
        _topo2ltopo = Topo2Ltopo.Topo2Ltopo()

        # create a Process Pool
        print("CPU内核数:{}".format(cpu_count()))
        print('当前母进程: {}'.format(os.getpid()))

        start = time.time()
        p = Pool(n_pools)

        for i_topo in topo_keys:
            # get topography info
            i_topo_info = dic_topo_info.get(i_topo)
            if i_topo_info is None:
                Exception('Cannot find %s in given topography info file' % i_topo)

            fin_name = os.path.join(fin_dir, i_topo_info.get('ObjFile'))

            search_all = eval(i_topo_info.get('SearchAll', 'True'))
            if search_all == '':
                search_all = True
            search_dimension = i_topo_info.get('SearchDimension')
            if search_dimension is not None and search_dimension != '':
                search_dimension = [float(i) for i in search_dimension]

            if fin_name is None:
                Exception('Cannot get ObjFile in given topography info file')

            # begin to do the calculation
            for i_det in _v_dets:
                logger.info('Run _multi_progressing_Topo2Ltopo for %s detector %s' % (i_topo, i_det.get_det_name()))

                out_fname = os.path.basename(fin_name).replace('.obj', '_%s.root' % i_det.get_det_name()).replace(
                    '.csv', '_%s.root' % i_det.get_det_name())
                out_fname = os.path.join(_directory, out_fname)

                if search_dimension is not None and search_dimension != '':
                    _topo = Topography.Topography(fin_name, search_dimension[0], search_dimension[1],
                                                  search_dimension[2], topo_name=i_topo)
                else:
                    _topo = Topography.Topography(fin_name, topo_name=i_topo)
                _topo.read_topo_obj(search_all=search_all)

                p.apply_async(func=_topo2ltopo.write_ltopo, args=(_topo, i_det, out_fname, calc_all))

        p.close()
        p.join()

        end = time.time()
        print("总共用时{}秒".format((end - start)))

    def Ltopo2Leff(self, kwargs):
        """
        Forward modeling step 2
        @param kwargs: OutDir; Input.LtopoDir; Input.Expressions; Output.Prefixes
        @return:
        """
        logger.info('Run Ltopo2Leff')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs in the config file
        det_names = kwargs.get('Input.DetNames')
        if det_names is None or not isinstance(det_names, list):
            raise Exception('Ltopo2Leff Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(det_names)

        input_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.LtopoDir', ''))
        _v_expr = kwargs.get('Input.Expressions')
        _v_outfile_prefix = kwargs.get('Output.Prefixes')
        if not len(_v_expr) == len(_v_outfile_prefix):
            raise Exception('The number of Expressions dose not equal to that of output Prefixes.')

        # begin to calculate Ltopo2Leff
        for i in range(len(_v_expr)):
            for i_det in _v_dets:
                expr = _v_expr[i]
                leff_fname = os.path.join(_directory,
                                          _v_outfile_prefix[i] + '_Leff_' + i_det.get_det_name() + '.root')
                _ltopo2leff = Ltopo2Leff.Ltopo2Leff(self._setting.get_topos_info(), input_dir)
                _ltopo2leff.write_Leff_from_expression(expr, leff_fname, i_det.get_det_name())

    def Leff2Res(self, kwargs):
        """
        Forward modeling step 3
        @param kwargs: OutDir; Input.DetNames; Input.PhysModelNames; Input.Epsilon; Input.LeffDir; Input.Prefixes
        @return:
        """
        logger.info('Run Leff2Res')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get infos from its config file
        dics_phys_model = self._setting.get_phys_model_info()

        # get kwargs in the config file
        det_names = kwargs.get('Input.DetNames')
        if det_names is None or not isinstance(det_names, list):
            raise Exception('Leff2Res Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(det_names)

        _v_phys_model_name = kwargs.get('Input.PhysModelNames')
        if _v_phys_model_name is None:
            raise Exception('Leff2Res Input.PhysModelNames is None.')
        eps = kwargs.get('Input.Epsilon', '')
        if eps == '':
            eps = 0.0002

        input_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.LeffDir', ''))

        _v_infile_prefix = kwargs.get('Input.Prefixes')

        for i_phys_model_name in _v_phys_model_name:
            _phys_model = self._get_phys_model(i_phys_model_name)

            for in_prefix in _v_infile_prefix:
                for i_det in _v_dets:
                    _size_plane_det = i_det.get_size_plane_det()
                    area = _size_plane_det[0] * _size_plane_det[1] / 100  # /100 means conversion mm into cm

                    _leff2res = Leff2Res.Leff2Res()

                    leff_fname = in_prefix + '_Leff_' + i_det.get_det_name() + '.root'
                    leff_fname = os.path.join(input_dir, leff_fname)

                    # calculate ratio and flux
                    # calculate ratio
                    _leff2res.calc_ratio(_phys_model, leff_fname, eps=float(eps))
                    # calculate flux
                    _leff2res.calc_flux(_phys_model, leff_fname, eps=float(eps))
                    # calculate flux times solid angle
                    solid_angle_info = {'Output.FileName': os.path.join(self._setting.get_global_dir(),
                                                                        kwargs.get('Input.SolidAngleFile'))}
                    _h_thetaphi_solid_angle, _h_thetaxy_solid_angle = self.GetSolidAngle(solid_angle_info)
                    _leff2res.calc_flux_x_SldAgl(_h_thetaphi_solid_angle, 'thetaphi')
                    _leff2res.calc_flux_x_SldAgl(_h_thetaxy_solid_angle, 'thetaxy')

                    # calculate counting rate and counts
                    _h_accept_thetaphi, _h_accept_thetaxy = i_det.get_acceptance_plane_det()
                    _leff2res.calc_counting_rate(_h_accept_thetaphi, area, 'thetaphi')
                    _leff2res.calc_counts(i_det.get_meas_time(), 'thetaphi')
                    _leff2res.calc_counting_rate(_h_accept_thetaxy, area, 'thetaxy')
                    _leff2res.calc_counts(i_det.get_meas_time(), 'thetaxy')

                    res_fname = in_prefix + '_' + i_phys_model_name + '_Res_' + i_det.get_det_name() + '.root'
                    res_fname = os.path.join(_directory, res_fname)

                    _leff2res.write_pred_res(res_fname)

                    # background calculation
                    _leff2res_bg = Leff2Res.Leff2Res()
                    bg_meas_time = i_det.get_background_meas_time()
                    # calculate flux
                    _h_accept_thetaphi, _h_accept_thetaxy = i_det.get_acceptance_plane_det()
                    _leff2res_bg.calc_flux_free_sky(_phys_model, _h_accept_thetaphi, _h_accept_thetaxy, eps=float(eps))
                    # calculate flux times solid angle
                    solid_angle_info = {'Output.FileName': os.path.join(self._setting.get_global_dir(),
                                                                        kwargs.get('Input.SolidAngleFile'))}
                    _h_thetaphi_solid_angle, _h_thetaxy_solid_angle = self.GetSolidAngle(solid_angle_info)
                    _leff2res_bg.calc_flux_x_SldAgl(_h_thetaphi_solid_angle, 'thetaphi')
                    _leff2res_bg.calc_flux_x_SldAgl(_h_thetaxy_solid_angle, 'thetaxy')

                    # calculate counting rate and counts
                    _leff2res_bg.calc_counting_rate(_h_accept_thetaphi, area, 'thetaphi')
                    _leff2res_bg.calc_counts(bg_meas_time, 'thetaphi')
                    _leff2res_bg.calc_counting_rate(_h_accept_thetaxy, area, 'thetaxy')
                    _leff2res_bg.calc_counts(bg_meas_time, 'thetaxy')

                    res_fname = in_prefix + '_' + i_phys_model_name + '_BackgroundRes_' + i_det.get_det_name() + '.root'
                    res_fname = os.path.join(_directory, res_fname)

                    _leff2res_bg.write_pred_res(res_fname)

    def ExpectMeasTime(self, kwargs):
        """
        print simLeff as an obs.dat
        @param kwargs: OutDir; Input.DetNames; Input.LeffDir; Input.LeffPrefixes; Input.AngleCoors; Input.LeffUncertaintyLevel; Output.FileNames
        @return:
        """
        logger.info('Run ExpectMeasTime')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get infos from its config file
        input_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.ResDir', ''))
        sig = eval(kwargs.get('Input.Significance', "2."))
        print("Set the significance level as ", sig)

        det_names = kwargs.get('Input.DetNames')
        if det_names is None or not isinstance(det_names, list):
            raise Exception('ExpectMeasTime Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(det_names)

        lis_in_prefix = kwargs.get('Input.PrefixesExpressions', [])

        lis_fout_name = kwargs.get('Output.FileNames', '')
        if lis_fout_name == '':
            lis_fout_name = []
            for i in range(len(lis_in_prefix)):
                lis_fout_name.append('%s_measTime.root' % lis_in_prefix[i])
        lis_fout_name = [os.path.join(_directory, i) for i in lis_fout_name]

        if len(lis_in_prefix) != len(lis_fout_name):
            raise Exception('ExpectMeasTime len(lis_in_prefix) != len(lis_fout_name)')

        # begin to calc
        _resolution_tools = ResolutionTools.ResolutionTools()

        for i in range(len(lis_in_prefix)):
            _prefix_exp = lis_in_prefix[i]
            _pre_with = _prefix_exp.split('-')[0]
            _pre_without = _prefix_exp.split('-')[1]

            _fout_name = lis_fout_name[i]

            if os.path.exists(_fout_name):
                os.remove(_fout_name)

            lis_f_with_name = []
            lis_f_without_name = []
            lis_suffix = []
            for i_det in _v_dets:
                leff_fname_with = _pre_with + '_Res_' + i_det.get_det_name() + '.root'
                leff_fname_with = os.path.join(input_dir, leff_fname_with)
                lis_f_with_name.append(leff_fname_with)

                leff_fname_without = _pre_without + '_Res_' + i_det.get_det_name() + '.root'
                leff_fname_without = os.path.join(input_dir, leff_fname_without)
                lis_f_without_name.append(leff_fname_without)

                suffix = _pre_with + "_" + i_det.get_det_name() + "-" + _pre_without + "_" + i_det.get_det_name()
                lis_suffix.append(suffix)

            _resolution_tools.ExpectMeasTime(lis_f_without_name, lis_f_with_name, lis_suffix, _fout_name, sig)

    def PrintSimObsData(self, kwargs):
        """
        print simLeff as an obs.dat
        @param kwargs: OutDir; Input.DetNames; Input.LeffDir; Input.LeffPrefixes; Input.AngleCoors; Input.LeffUncertaintyLevel; Output.FileNames
        @return:
        """
        logger.info('Run PrintSimObsData')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get infos from its config file
        det_names = kwargs.get('Input.DetNames')
        if det_names is None or not isinstance(det_names, list):
            raise Exception('PrintSimObsData Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(det_names)

        _in_dir = kwargs.get('Input.LeffDir', '')
        _in_dir = os.path.join(self._setting.get_global_dir(), _in_dir)
        lis_in_prefix = kwargs.get('Input.LeffPrefixes', [])
        lis_angle_coor = kwargs.get('Input.AngleCoors', ['thetaphi'])
        leff_uncertainty_level = kwargs.get('Input.LeffUncertaintyLevel', 0.05)
        if leff_uncertainty_level == '':
            leff_uncertainty_level = 0.05
        leff_uncertainty_level = float(leff_uncertainty_level)

        if_print_zero = eval(kwargs.get('Input.PrintZero', 'True'))

        lis_fout_name = kwargs.get('Output.FileNames', '')
        if lis_fout_name == '':
            lis_fout_name = []
            for i in range(len(lis_in_prefix)):
                lis_fout_name.append('%s.dat' % lis_in_prefix[i])
        lis_fout_name = [os.path.join(_directory, i) for i in lis_fout_name]

        if len(lis_in_prefix) != len(lis_fout_name):
            raise Exception('PrintSimObsData len(lis_in_prefix) != len(lis_fout_name)')

        # begin to print
        _print_obj = PrintTObject()

        for i in range(len(lis_in_prefix)):
            _in_prefix = lis_in_prefix[i]
            _fout_name = lis_fout_name[i]

            if os.path.exists(_fout_name):
                os.remove(_fout_name)

            print('PrintLeffHist2Obs to ', _fout_name)
            fout = open(_fout_name, 'w')
            print('PrintLeffHist2Obs', file=fout)
            for i_det in _v_dets:
                fin_name = _in_prefix + '_Leff_' + i_det.get_det_name() + '.root'
                fin_name = os.path.join(_in_dir, fin_name)

                for _angle_coor in lis_angle_coor:
                    _print_obj.PrintLeffHist2Obs(fin_name, i_det, fout, _angle_coor, leff_uncertainty_level,
                                                 if_print_zero=if_print_zero)

    def _get_detector_list(self, lis_det_name):
        """
        create a list of detector objects
        @param lis_det_name: list of DetNames
        @return: list of Detector.Detector objects
        """
        _v_dets = []
        for i_det_name in lis_det_name:
            logger.info('Get detector information: ' + i_det_name)

            # get dictionary of the detector information
            dic_det_info = self._setting.get_dets_info().get(i_det_name)
            if dic_det_info is None:
                raise Exception('Cannot get detector information of %s from %s'
                                % (i_det_name, self._setting.get_det_info_fname()))

            # get the directory for the detector information and add it into det_info dictionary
            t_dir = self._setting.get_dets_info().get('Directory', '')
            if not CommonUtils.path_isabs(t_dir):
                t_dir = os.path.join(self._setting.get_global_dir(), t_dir)
            dic_det_info['GlobalDir'] = t_dir

            # get bins from global for creating histograms in detector functions (acceptance hists)
            bins_dic = self._setting.get_bins()

            # create a detector object
            det = Detector.Detector(det_name=i_det_name, dic_det_info=dic_det_info, bins_dic=bins_dic)
            _v_dets.append(det)
        return _v_dets

    def _get_phys_model(self, phys_name):
        all_phys_model_info = self._setting.get_phys_model_info()
        dic_phys_model = all_phys_model_info.get(phys_name)
        if dic_phys_model is None:
            raise Exception('Cannot find PhysModel %s' % phys_name)

        if_formula = True
        for i in dic_phys_model.keys():
            if i.lower().__contains__('hist'):
                if_formula = False
        if if_formula:
            phys_model = PhysModel.PhysModelFromFormula(dic_phys_model)
        else:
            if not CommonUtils.path_isabs(dic_phys_model['FluxModelRootName']):
                dic_phys_model['FluxModelRootName'] = os.path.join(self._setting.get_global_dir(),
                                                                   dic_phys_model['FluxModelRootName'])
            phys_model = PhysModel.PhysModelFromHist(dic_phys_model)

        phys_model.print_physModel()

        return phys_model

    def GetDetAcceptance(self, kwargs):
        """
        create a list of detector objects
        @param kwargs: OutDir; Input.DetNames; Output.Prefix
        @return:
        """
        logger.info('Run GetDetAcceptance')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs in the config file
        det_names = kwargs.get('Input.DetNames')
        if det_names is None or not isinstance(det_names, list):
            raise Exception('GetDetAcceptance Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(det_names)
        out_prefix = os.path.join(_directory, kwargs.get('Output.Prefix', ''))

        for i_det in _v_dets:
            det_type_id, det_type_str = i_det.get_det_type()
            if det_type_str.__contains__('plane'):
                i_det.get_acceptance_plane_det(out_prefix)
            else:
                raise Exception('This function only support planar detector.'
                                'Please expand this function by yourself.')
            # if det_type_str.__contains__('borehole'):

    def GetSolidAngle(self, kwargs):
        """
        calculate solid angle and output with root file
        @param kwargs: OutDir; Input.NRand; Output.FileName
        @return:
        """
        logger.info('Run GetSolidAngle.')

        if_enforce = False
        if kwargs.get('OutDir') is not None:
            if_enforce = True

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        nrand = int(float(kwargs.get('Input.NRand', 1e7)))
        out_file = kwargs.get('Output.FileName', 'size_of_angle.root')
        if out_file == '':
            out_file = 'size_of_angle.root'
        out_file = os.path.join(_directory, out_file)

        _angleMesh = AngleMesh.AngleMesh()

        _angleMesh.write_solid_angle(bins_dic=self._bins_dic, solid_angle_fname=out_file,
                                     nrand=nrand, enforce=if_enforce)
        _h_thetaphi_solid_angle = _angleMesh.get_thetaphi_solid_angle().Clone()
        _h_thetaxy_solid_angle = _angleMesh.get_thetaxy_solid_angle().Clone()

        return _h_thetaphi_solid_angle, _h_thetaxy_solid_angle

    def PrintDevicesInfo(self, kwargs):
        """
        print the information of all devices and output as csv
        @param kwargs: OutDir; Output.DeciveInfoFile
        @return:
        """
        logger.info('Run PrintDevicesInfo.')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        outf_name = kwargs.get('Output.DeciveInfoFile', 'device_info.csv')
        outf_name = os.path.join(_directory, outf_name)

        _database_interface = DatabaseInterface.DatabaseInterface()
        _database_interface.PrintDevicesInfo(outf_name=outf_name)

    def PrintOnlineTablesInfo(self, kwargs):
        """
        print the information of all online tables for certain serials and output as csv
        @param kwargs: OutDir; Input.serialNumbers; Input.serialNames; Output.Prefix
        @return:
        """
        logger.info('Run PrintOnlineTablesInfo.')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        lis_serial_number = kwargs.get('Input.serialNumbers')
        if lis_serial_number is None:
            raise Exception('Cannot get serialNumbers of PrintTablesInfo.')
        lis_serial_name = kwargs.get('Input.serialNames')
        if not len(lis_serial_name) == len(lis_serial_name):
            raise Exception('The number of serialNumbers must equal to theta of serialNames.')

        outf_prefix = os.path.join(_directory, kwargs.get('Output.Prefix'))

        # create a PrintDatabaseInfo object
        _database_interface = DatabaseInterface.DatabaseInterface()

        # print tables
        for i in range(len(lis_serial_number)):
            serial_number = lis_serial_number[i]
            serial_name = lis_serial_name[i]

            outf_name = '%s_table_%s.csv' % (outf_prefix, serial_name)
            _database_interface.PrintOnlineTablesInfo(serial_number=serial_number, outf_name=outf_name,
                                                      serial_name=serial_name)

    def PrintOfflineTablesInfo(self, kwargs):
        """
        print the information of all offline tables for certain serials and output as csv
        @param kwargs: OutDir; Input.Hdf5Dir; Output.FileName
        @return:
        """
        logger.info('Run PrintOfflineTablesInfo.')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = kwargs.get('Input.Hdf5Dir')
        h5_dir = os.path.join(self._setting.get_global_dir(), h5_dir)

        outf_name = os.path.join(_directory, kwargs.get('Output.FileName'))

        # create a PrintDatabaseInfo object
        _database_interface = DatabaseInterface.DatabaseInterface()

        _database_interface.PrintOfflineTablesInfo(h5_dir, outf_name)

    def DownloadHdf5(self, kwargs):
        """
        download data from database and save as hdf5 files
        @param kwargs: OutDir; Input.TablesInfoFile; Input.PageCapacity
        @return:
        """
        logger.info('Run DownloadTables.')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        table_file = kwargs.get('Input.TablesInfoFile')
        table_file = os.path.join(self._setting.get_global_dir(), table_file)
        page_capacity = int(kwargs.get('Input.PageCapacity'))

        # read table names for converting
        df_table_infos = pd.read_csv(table_file)

        _database_interface = DatabaseInterface.DatabaseInterface()
        _database_interface.DownloadHdf5(df_table_infos, _directory, page_capacity)

    def FindBarChannelRelation(self, kwargs):
        """
        find the bar-channel relation for hdf5 files
        :param kwargs: hdf5 file names
                OutDir; Input.Hdf5Dir; Input.Hdf5Names; Output.CsvNames; Input.NumOfBars (option)
        :return: csv files
                OutDir/CsvNames
        """
        logger.info('Run FindBarChannelRelation.')

        # create a directory if not exist
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fnames = kwargs.get('Input.Hdf5Names')

        lis_out_fname = kwargs.get('Output.CsvNames', '')
        if lis_out_fname == '':
            lis_out_fname = [os.path.splitext(i)[0] + '.csv' for i in lis_h5_fnames]

        lis_num_bars = kwargs.get('Input.NumOfBars', '')
        if lis_num_bars == '':
            lis_num_bars = None
        elif isinstance(lis_num_bars, str):
            n_bars = int(lis_num_bars)
            lis_num_bars = [n_bars for i in range(len(lis_h5_fnames))]

        if len(lis_h5_fnames) != len(lis_out_fname):
            raise Exception(
                'FindBarChannelRelation: number of Input.Hdf5Names is not equal to that of Output.CsvNames.')

        for i in range(len(lis_h5_fnames)):
            h5_name = lis_h5_fnames[i]
            (t_dir, t_h5_name) = os.path.split(h5_name)
            h5_dir = os.path.join(h5_dir, t_dir)
            t_h5_name = CommonUtils.splitext_h5_fname(t_h5_name)

            t_csv_fname = lis_out_fname[i]
            t_csv_fname = os.path.join(_directory, t_csv_fname)

            if lis_num_bars is None:
                DataframeAnalyzer.FindBarChannelRelation(t_h5_name, t_csv_fname, h5_dir)
            else:
                n_bars = lis_num_bars[i]
                DataframeAnalyzer.FindBarChannelRelation(t_h5_name, t_csv_fname, h5_dir, n_bars)

    def SelectWithCondition(self, kwargs):
        """
        select data with following conditions and save in a new h5 file
                    start_time, end_time, start_trigger, end_trigger
        @param kwargs: OutDir; Input.Hdf5Dir; Input.Hdf5Names; Output.Hdf5Name
                    Input.StartTime, Input.EndTime, Input.StartTrigger, Input.EndTrigger
        @return: csv files
                    OutDir/Hdf5Name
        """
        logger.info('Run SelectWithCondition.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        h5_name = kwargs.get('Input.Hdf5Name')
        h5_name = CommonUtils.splitext_h5_fname(h5_name)

        key_val = {}

        start_time = kwargs.get('Input.StartTime', '')
        if not start_time == '':
            key_val['start_time'] = start_time
        end_time = kwargs.get('Input.EndTime', '')
        if not end_time == '':
            key_val['end_time'] = end_time
        start_trigger = kwargs.get('Input.StartTrigger', '')
        if not start_trigger == '':
            key_val['start_trigger'] = start_trigger
        end_trigger = kwargs.get('Input.EndTrigger', '')
        if not end_trigger == '':
            key_val['end_trigger'] = end_trigger

        out_fname = os.path.join(_directory, kwargs.get('Output.Hdf5Name'))
        if not out_fname.__contains__('.h5'):
            out_fname += '.h5'
        DataframeProcessor.SelectWithCondition(h5_name, out_fname, h5_dir, key_val)

    def PrintHdf5ToFile(self, kwargs):
        """
        print hdf5 to file
        @param kwargs: OutDir; Input.Hdf5Dir; Input.Hdf5Name; Output.FileName
        @return: OutDir/csvFile
        """
        logger.info('Run PrintHdf5ToFile.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        h5_name = kwargs.get('Input.Hdf5Name')

        out_fname = kwargs.get('Output.FileName', 'printHdf5')
        if not out_fname.__contains__('.'):
            out_fname += '.txt'
        out_fname = os.path.join(_directory, out_fname)

        DataframeProcessor.PrintHdf5ToFile(h5_name, out_fname, h5_dir)

    def ReplaceChannel(self, kwargs):
        """
        replace channel and save in another hdf5 file
        @param kwargs: OutDir; Input.Hdf5Dir; Input.Hdf5Names;
            Input.BarChannelRelationDir; Input.BarChannelRelationNames; Output.Hdf5Names
        @return: OutDir/Hdf5Names
        """
        logger.info('Run ReplaceChannel.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_name = [os.path.join(h5_dir, i) for i in kwargs.get('Input.Hdf5Names')]

        relation_info_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.BarChannelRelationDir', ''))
        lis_relation_info_fname = kwargs.get('Input.BarChannelRelationNames', '')
        lis_relation_info_fname = [os.path.join(relation_info_dir, i) for i in lis_relation_info_fname]

        lis_out_fname = kwargs.get('Output.Hdf5Names', '')
        if lis_out_fname == '':
            lis_out_fname = None
        else:
            lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        if len(lis_h5_name) != len(lis_relation_info_fname):
            raise Exception('ReplaceChannel: number of Input.Hdf5Names '
                            'and Input.BarChannelRelationFiles must be the same.')

        for i in range(len(lis_h5_name)):
            h5_fname = lis_h5_name[i]
            (h5_dir, h5_name) = os.path.split(h5_fname)
            relation_info_fname = lis_relation_info_fname[i]

            if lis_out_fname is not None:
                out_fname = lis_out_fname[i]
            else:
                out_fname = os.path.join(_directory, 'reChn_%s' % h5_name)
            if not out_fname.__contains__('.h5'):
                out_fname += '.h5'

            DataframeProcessor.ReplaceChannel(h5_name, out_fname, h5_dir, relation_info_fname=relation_info_fname)

    def FiltrateAllHit(self, kwargs):
        """
        filtrate all boards hit events
        @param kwargs: OutDir; Input.CheckContinueTrigger; Input.Hdf5Dir; Input.Hdf5Names;
            Output.Hdf5Name; Output.LoggerName (optional)
        @return: OutDir/Hdf5Name
        """
        logger.info('Run FiltrateAllHit.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        check_continue_trigger = kwargs.get('Input.CheckContinueTrigger', 'False')
        check_continue_trigger = eval(check_continue_trigger)

        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        if isinstance(lis_h5_fname, str) and lis_h5_fname.__contains__("re:"):
            pattern = lis_h5_fname.replace("re:", "")
            lis_h5_fname = CommonUtils.search_files_with_re(h5_dir, pattern)
        else:
            lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        lis_out_fname = kwargs.get('Output.Hdf5Name', '')
        if lis_out_fname == '':
            lis_out_fname = None

        lis_log_fname = kwargs.get('Output.LoggerName', '')
        if lis_log_fname == '':
            lis_log_fname = None

        for i in range(len(lis_h5_fname)):
            h5_fname = lis_h5_fname[i]
            (h5_dir, h5_name) = os.path.split(h5_fname)
            (h5_name, file_type) = os.path.splitext(h5_name)

            if lis_out_fname is None:
                out_fname = 'allHit_%s.h5' % h5_name
            else:
                out_fname = lis_out_fname[i]
            out_fname = os.path.join(_directory, out_fname)

            if lis_log_fname is None:
                out_log_fname = './log/allHit_%s.log' % h5_name
            else:
                out_log_fname = lis_log_fname[i]
            out_log_fname = os.path.join(_directory, out_log_fname)
            (log_dir, log_file) = os.path.split(out_log_fname)
            if not os.path.exists(log_dir):
                os.mkdir(log_dir)

            DataframeProcessor.FiltrateAllHit(h5_name, out_fname, h5_dir, out_log_fname, check_continue_trigger)

    def DeleteGivenChannels(self, kwargs):
        """
        delete given channels in dataframe and save in another h5 file
        @param kwargs: OutDir; Input.Hdf5Dir; Input.Hdf5Name; Input.DeleteChannelNumbers;
            Output.Hdf5Name; Output.LoggerName (optional)
        @return: OutDir/Hdf5Name
        """
        logger.info('Run DeleteGivenChannel.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        h5_name = kwargs.get('Input.Hdf5Name')

        lis_del_number = kwargs.get('Input.DeleteChannelNumbers')

        out_fname = kwargs.get('Output.Hdf5Name', 'delChn_%s' % h5_name)
        if out_fname == '':
            out_fname = 'delChn_%s' % h5_name
        if not out_fname.__contains__('.h5'):
            out_fname += '.h5'
        out_fname = os.path.join(_directory, out_fname)

        out_log_fname = kwargs.get('Output.LoggerName')
        if out_log_fname is None or out_log_fname == '':
            DataframeProcessor.DeleteGivenChannels(h5_name, out_fname, h5_dir, lis_del_number)
        else:
            DataframeProcessor.DeleteGivenChannels(h5_name, out_fname, h5_dir, lis_del_number, out_log_fname)

    def Set0Randomly(self, kwargs):
        """
        set adc as 0 randomly for certain channels and save in another h5 file
        @param kwargs: OutDir; Input.Hdf5Dir; Input.Hdf5Name; Input.ChannelNumbers; Input.Opt; Input.Probability
            Output.Hdf5Name;
        @return: OutDir/Hdf5Name
        """
        logger.info('Run Set0Randomly.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        h5_name = kwargs.get('Input.Hdf5Name')

        lis_del_number = kwargs.get('Input.ChannelNumbers')

        out_fname = kwargs.get('Output.Hdf5Name', 'delChn_%s' % h5_name)
        if out_fname == '':
            out_fname = 'delChn_%s' % h5_name
        if not out_fname.__contains__('.h5'):
            out_fname += '.h5'
        out_fname = os.path.join(_directory, out_fname)

        opt = eval(kwargs.get('Input.Opt', '1'))
        prob = eval(kwargs.get('Input.Probability', '0.5'))

        DataframeProcessor.Set0Randomly(h5_name, out_fname, h5_dir, lis_del_number, opt, prob)

    def FindBaseline(self, kwargs):
        """
        find the baseline with the given baseline csv file
        @param kwargs: OutDir; Input.Hdf5Dir; Input.Hdf5Names; Output.CsvNames; Output.RootNames
        @return: OutDir/RootNames
        """
        logger.info('Run FindBaseline.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fnames = kwargs.get('Input.Hdf5Names')

        lis_csv_fname = kwargs.get('Output.CsvNames', '')
        if lis_csv_fname == '':
            lis_csv_fname = [os.path.splitext(i)[0] + '_baseline.csv' for i in lis_h5_fnames]

        lis_root_fname = kwargs.get('Output.RootNames', '')
        if lis_root_fname == '':
            lis_root_fname = [os.path.splitext(i)[0] + '_baseline.root' for i in lis_h5_fnames]

        if len(lis_h5_fnames) != len(lis_csv_fname) or len(lis_h5_fnames) != len(lis_root_fname):
            raise Exception('FindBaseline: number of Input.Hdf5Names, '
                            'Output.CsvNames and Output.RootNames must be the same.')

        for i in range(len(lis_h5_fnames)):
            h5_name = lis_h5_fnames[i]
            (t_dir, t_h5_name) = os.path.split(h5_name)
            h5_dir = os.path.join(h5_dir, t_dir)
            t_h5_name = CommonUtils.splitext_h5_fname(t_h5_name)

            t_csv_fname = lis_csv_fname[i]
            t_csv_fname = os.path.join(_directory, t_csv_fname)

            t_root_fname = lis_root_fname[i]
            if not t_root_fname.__contains__('.root'):
                t_root_fname = t_root_fname + '.root'
            t_root_fname = os.path.join(_directory, t_root_fname)

            DataframeAnalyzer.FindBaseline(t_h5_name, t_csv_fname, h5_dir, t_root_fname)

    def DrawBaselineList(self, kwargs):
        """
        compare the baseline from different files
        @param kwargs: OutDir; Input.BaselineDir; Input.BaselineCsvNames; Output.PdfName
        @return: OutDir/PdfName, RootName
        """
        logger.info('Run DrawBaselineList.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        csv_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.BaselineDir', ''))
        lis_csv_fnames = kwargs.get('Input.BaselineCsvNames')
        lis_csv_fnames = [os.path.join(csv_dir, i) for i in lis_csv_fnames]

        pdf_name = kwargs.get('Output.PdfName', '')
        if pdf_name == '':
            raise Exception('DrawBaselineList must give a Output.PdfName')
        pdf_name = os.path.join(_directory, pdf_name)

        DataframeAnalyzer.DrawBaselineList(lis_csv_fnames, pdf_name)

        _print_Tobject = PrintTObject()
        root_fname = os.path.splitext(pdf_name)[0] + '.root'
        _print_Tobject.PrintTFile2PDF(root_fname, pdf_name)

    def FindThreshold(self, kwargs):
        """
        find the threshold with the given threshold file
        @param kwargs: OutDir; Input.Hdf5Dir; Input.Hdf5Names; Input.BaselineDir; Input.BaselineNames;
            Input.SmoothLevels (default 2); Input.ADCOffset (default 3); Output.CsvNames; Output.RootNames
        @return: OutDir/CsvNames, RootNames
        """
        logger.info('Run FindThreshold.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fnames = kwargs.get('Input.Hdf5Names')
        lis_h5_fnames = [os.path.join(h5_dir, i) for i in lis_h5_fnames]

        baseline_dir = kwargs.get('Input.BaselineDir')
        baseline_dir = os.path.join(self._setting.get_global_dir(), baseline_dir)
        lis_baseline_fnames = kwargs.get('Input.BaselineNames')
        lis_baseline_fnames = [os.path.join(baseline_dir, i) for i in lis_baseline_fnames]

        lis_smooth_level = kwargs.get('Input.SmoothLevels', '')
        if lis_smooth_level == '':
            lis_smooth_level = []
            for i in range(len(lis_h5_fnames)):
                lis_smooth_level.append(2)
        if isinstance(lis_smooth_level, str):
            smooth_level = int(lis_smooth_level)
            lis_smooth_level = []
            for i in range(len(lis_h5_fnames)):
                lis_smooth_level.append(smooth_level)
        lis_smooth_level = [int(i) for i in lis_smooth_level]

        offset = kwargs.get('Input.ADCOffset', '')
        if offset == '':
            offset = 3
        offset = int(offset)

        lis_csv_fname = kwargs.get('Output.CsvNames', '')
        if lis_csv_fname == '':
            lis_csv_fname = None
        else:
            lis_csv_fname = [os.path.join(_directory, i) for i in lis_csv_fname]

        lis_root_fname = kwargs.get('Output.RootNames', '')
        if lis_root_fname == '':
            lis_root_fname = None
        else:
            lis_root_fname = [os.path.join(_directory, i) for i in lis_root_fname]

        if len(lis_h5_fnames) != len(lis_baseline_fnames):
            raise Exception('FindThreshold: number of Input.Hdf5Names and Input.BaselineNames'
                            ' must be the same.\nHdf5Names: %d, BaselineNames: %d'
                            % (len(lis_h5_fnames), len(lis_baseline_fnames)))

        for i in range(len(lis_h5_fnames)):
            h5_name = lis_h5_fnames[i]
            (t_dir, t_h5_name) = os.path.split(h5_name)
            h5_dir = os.path.join(h5_dir, t_dir)
            t_h5_name = CommonUtils.splitext_h5_fname(t_h5_name)

            t_baseline_fname = lis_baseline_fnames[i]

            if lis_csv_fname is not None:
                t_csv_fname = lis_csv_fname[i]
                t_csv_fname = os.path.join(_directory, t_csv_fname)
            else:
                t_csv_fname = t_h5_name + '_threshold.csv'
                t_csv_fname = os.path.join(_directory, t_csv_fname)

            if lis_csv_fname is not None:
                t_root_fname = lis_root_fname[i]
                t_root_fname = os.path.join(_directory, t_root_fname)
            else:
                t_root_fname = t_h5_name + '_threshold.root'
                t_root_fname = os.path.join(_directory, t_root_fname)

            smooth_level = lis_smooth_level[i]

            DataframeAnalyzer.FindThreshold(h5_fname=t_h5_name, directory=h5_dir, baseline_csv=t_baseline_fname,
                                            csv_out_fname=t_csv_fname, root_fname=t_root_fname,
                                            smooth_level=smooth_level, offset=offset)

    def DrawEnergySpectrum(self, kwargs):
        """draw energy spectrum with given hdf5 file and key name"""
        logger.info('Run DrawEnergySpectrum.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fnames = kwargs.get('Input.Hdf5Names')
        lis_h5_fnames = [os.path.join(h5_dir, i) for i in lis_h5_fnames]

        lis_root_fname = kwargs.get('Output.RootNames', '')
        if lis_root_fname == '':
            lis_root_fname = None
        else:
            lis_root_fname = [os.path.join(_directory, i) for i in lis_root_fname]

        key_name = kwargs.get('Input.DrawKeyName')

        for i in range(len(lis_h5_fnames)):
            h5_name = lis_h5_fnames[i]
            (t_dir, t_h5_name) = os.path.split(h5_name)
            h5_dir = os.path.join(h5_dir, t_dir)
            t_h5_name = CommonUtils.splitext_h5_fname(t_h5_name)

            if lis_root_fname is not None:
                t_root_fname = lis_root_fname[i]
            else:
                t_root_fname = t_h5_name + '_energySpectrum.root'
                t_root_fname = os.path.join(_directory, t_root_fname)

            DataframeAnalyzer.DrawEnergySpectrum(h5_fname=t_h5_name, directory=h5_dir,
                                                 root_fname=t_root_fname, key_name=key_name)

    def SpectrumTimeSlice(self, kwargs):
        logger.info('Run SpectrumTimeSlice.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fnames = kwargs.get('Input.Hdf5Names')

        lis_root_fout_name = kwargs.get('Output.RootNames', '')
        if lis_root_fout_name == '':
            lis_root_fout_name = ['specSlice_' + CommonUtils.splitext_h5_fname(i) + '.root' for i in lis_h5_fnames]
            lis_root_fout_name = [os.path.join(_directory, i) for i in lis_root_fout_name]
        if len(lis_h5_fnames) != len(lis_root_fout_name):
            raise Exception('Number of hdf5 names and root names must be the same.')

        lis_h5_fnames = [os.path.join(h5_dir, i) for i in lis_h5_fnames]

        lis_time_interval = kwargs.get('Input.TimeIntervals', '')
        if lis_time_interval == '':
            raise Exception('Input.TimeIntervals = ??')
        if isinstance(lis_time_interval, str):
            time_interval = CommonUtils.time_convert(lis_time_interval)
            lis_time_interval = []
            for i in range(len(lis_h5_fnames)):
                lis_time_interval.append(time_interval)
        else:
            lis_time_interval = [CommonUtils.time_convert(i) for i in lis_time_interval]

        if_print_pdf = kwargs.get('Input.PrintPdf', 'True')
        if_print_pdf = eval(if_print_pdf)
        _print_Tobject = PrintTObject()

        for i in range(len(lis_h5_fnames)):
            h5_fname = os.path.splitext(lis_h5_fnames[i])[0] + '.h5'
            time_interval = lis_time_interval[i]
            root_fout_name = lis_root_fout_name[i]

            _dataframeAnalyzer = DataframeAnalyzer()
            df_meas_time = _dataframeAnalyzer.get_df_meas_time(h5_fname)
            lis_start_time, lis_end_time = _dataframeAnalyzer.get_time_list(df_meas_time, time_interval)

            _dataframeProcessor = DataframeProcessor()
            str_hadd = 'hadd -f'
            str_lis_root_fname = ''
            str_lis_h5_fname = ''
            for i_time in range(len(lis_start_time)):
                start_time = lis_start_time[i_time]
                end_time = lis_end_time[i_time]
                dic_time = {'i_time': i_time,
                            'start_time': str(start_time),
                            'end_time': str(end_time)}

                h5_out_name = h5_fname.replace('.h5', '_split%d.h5' % i_time)
                (h5_dir, h5_name) = os.path.split(h5_fname)
                dict_condition = _dataframeProcessor.SelectWithCondition(h5_fname, h5_out_name, h5_dir, dic_time)

                # update end_time
                dic_time['end_time'] = dict_condition['end_time']

                (h5_dir, h5_name) = os.path.split(h5_out_name)
                root_split_name = h5_out_name.replace('.h5', '.root')
                _dataframeAnalyzer.DrawEnergySpectrum(h5_name, h5_dir, root_split_name, 'data',
                                                      dic_split=dic_time)

                str_lis_root_fname = str_lis_root_fname + ' ' + root_split_name
                str_lis_h5_fname = str_lis_h5_fname + ' ' + h5_out_name

            str_hadd += ' ' + root_fout_name + str_lis_root_fname
            os.system(str_hadd)
            os.system('rm %s %s' % (str_lis_root_fname, str_lis_h5_fname))

            # print pdf
            if if_print_pdf:
                lis_draw_option = kwargs.get('Input.DrawOptions', '')
                self._print_root_in_pdf(root_fout_name, lis_draw_option)

    def ReplaceCsvChannel(self, kwargs):
        """replace channel for csv file and save in another csv file"""
        logger.info('Run ReplaceCsvChannel.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        bar_chn_rel_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.BarChannelRelation', ''))
        lis_csv_fname = kwargs.get('Input.CsvNames')
        lis_csv_fname = [os.path.join(self._setting.get_global_dir(), i) for i in lis_csv_fname]

        lis_out_fname = kwargs.get('Output.CsvNames')
        lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        if len(lis_csv_fname) != len(lis_out_fname):
            raise Exception('ReplaceCsvChannel: number of Input.CsvNames dose not equal to Output.CsvNames.')

        for i in range(len(lis_csv_fname)):
            csv_fname = lis_csv_fname[i]
            out_fname = lis_out_fname[i]
            CommonUtils.replace_csv_channel(csv_fname, out_fname, bar_chn_rel_fname)

    def MergeHdf5(self, kwargs):
        """merge hdf5 and save as another hdf5 file"""
        logger.info('Run MergeHdf5.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        out_fname = kwargs.get('Output.Hdf5Name')
        if not out_fname.__contains__('.h5'):
            out_fname += '.h5'
        out_fname = os.path.join(_directory, out_fname)

        log_fname = kwargs.get('Output.LoggerName', '')
        if not log_fname == '':
            log_fname = os.path.join(_directory, log_fname)

        if log_fname == '':
            DataframeProcessor.MergeHdf5(lis_h5_fname, out_fname)
        else:
            DataframeProcessor.MergeHdf5(lis_h5_fname, out_fname, log_fname)

    def CutWithThreshold(self, kwargs):
        """cut dataframe with the given threshold"""
        logger.info('Run CutWithThreshold.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        threshold_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.ThresholdDir'))
        lis_threshold_fname = kwargs.get('Input.ThresholdNames')
        lis_threshold_fname = [os.path.join(threshold_dir, i) for i in lis_threshold_fname]

        lis_out_fname = kwargs.get('Output.Hdf5Names', '')
        if lis_out_fname == '':
            lis_out_fname = None
        else:
            lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        if len(lis_h5_fname) != len(lis_threshold_fname):
            raise Exception('CutWithThreshold: the number of Input.Hdf5Names and '
                            'Input.ThresholdNames are not the same.')

        # cut with certanin threshold
        for i in range(len(lis_h5_fname)):
            (h5_dir, h5_fname) = os.path.split(lis_h5_fname[i])
            thre_fname = lis_threshold_fname[i]
            if lis_out_fname is not None:
                out_fname = lis_out_fname[i]
            else:
                out_fname = 'cutThre_' + h5_fname
                out_fname = os.path.join(_directory, out_fname)
            if not out_fname.__contains__('.h5'):
                out_fname += '.h5'

            DataframeProcessor.CutWithThreshold(h5_fname, out_fname, h5_dir, thre_fname)

    def SubtractBaseline(self, kwargs):
        """subtract dataframe to the given baseline"""
        logger.info('Run SubtractBaseline.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        baseline_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.BaselineDir', ''))
        lis_baseline_fname = kwargs.get('Input.BaselineNames')
        lis_baseline_fname = [os.path.join(baseline_dir, i) for i in lis_baseline_fname]

        lis_out_fname = kwargs.get('Output.Hdf5Names', '')
        if lis_out_fname == '':
            lis_out_fname = None
        else:
            lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        if len(lis_h5_fname) != len(lis_baseline_fname):
            raise Exception('SubtractBaseline: the number of Input.Hdf5Names and '
                            'Input.BaselineNames are not the same.')

        # cut with certain threshold
        for i in range(len(lis_h5_fname)):
            (h5_dir, h5_fname) = os.path.split(lis_h5_fname[i])
            baseline_fname = lis_baseline_fname[i]
            if lis_out_fname is not None:
                out_fname = lis_out_fname[i]
            else:
                out_fname = os.path.join(_directory, 'subBl_' + h5_fname)
            if not out_fname.__contains__('.h5'):
                out_fname += '.h5'

            DataframeProcessor.SubtractBaseline(h5_fname, out_fname, h5_dir, baseline_fname)

    def TimeSlice(self, kwargs):
        """time slice"""
        logger.info('Run TimeSlice.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        if isinstance(lis_h5_fname, str) and lis_h5_fname.__contains__("re:"):
            pattern = lis_h5_fname.replace("re:", "")
            lis_h5_fname = CommonUtils.search_files_with_re(h5_dir, pattern)
        else:
            lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        lis_time_interval = kwargs.get('Input.TimeIntervals', '')
        if lis_time_interval == '':
            raise Exception('Input.TimeIntervals = ??')
        if isinstance(lis_time_interval, str):
            time_interval = CommonUtils.time_convert(lis_time_interval)
            lis_time_interval = []
            for i in range(len(lis_h5_fname)):
                lis_time_interval.append(time_interval)
        else:
            lis_time_interval = [CommonUtils.time_convert(i) for i in lis_time_interval]

        lis_out_fname = kwargs.get('Output.RootNames', '')
        if lis_out_fname == '':
            lis_out_fname = None
        else:
            lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        if lis_out_fname is not None and len(lis_out_fname) != len(lis_h5_fname):
            raise Exception('TimeSlice: the number of Input.Hdf5Names '
                            'and Output.RootNames are not the same.')

        pageCapacity = kwargs.get('Input.PageCapacity', '')
        if pageCapacity == '':
            pageCapacity = 100000
        else:
            pageCapacity = int(pageCapacity)

        checkFactor = kwargs.get('Input.CheckFactor', '')
        if checkFactor == '':
            checkFactor = 3
        else:
            checkFactor = int(checkFactor)

        if_temp_str = kwargs.get('Input.IfTemperature', '')
        if if_temp_str.lower() == 'true':
            if_temperature = True
        else:
            if_temperature = False

        if_print_pdf = kwargs.get('Input.PrintPdf', 'True')
        if_print_pdf = eval(if_print_pdf)

        # cut with certain threshold
        _df_analyzer = DataframeAnalyzer()
        _print_Tobject = PrintTObject()
        for i in range(len(lis_h5_fname)):
            (h5_dir, h5_fname) = os.path.split(lis_h5_fname[i])
            time_interval = lis_time_interval[i]

            if lis_out_fname is None:
                (h5_name, file_type) = os.path.splitext(h5_fname)
                out_fname = 'timeSlice_' + h5_name
                out_fname = os.path.join(_directory, out_fname)
            else:
                out_fname = lis_out_fname[i]
            if not out_fname.__contains__('.root'):
                out_fname += '.root'

            _df_analyzer.TimeSlice(h5_fname, h5_dir, out_fname, time_interval, pageCapacity,
                                   checkFactor=checkFactor, if_temperature=if_temperature)

            if if_print_pdf:
                lis_draw_option = kwargs.get('Input.DrawOptions', '')
                self._print_root_in_pdf(out_fname, lis_draw_option)

    def PrintTFile2PDF(self, kwargs):
        """print hists, graphs, canvas to pdf file"""
        logger.info('Run PrintTFile2PDF.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = kwargs.get('Input.RootDir', '')
        root_dir = os.path.join(self._setting.get_global_dir(), root_dir)
        lis_root_fname = kwargs.get('Input.RootNames')
        lis_root_fname = [os.path.join(root_dir, i) for i in lis_root_fname]

        lis_out_fname = kwargs.get('Output.PDFNames', '')
        if lis_out_fname == '':
            lis_out_fname = []
            for i_root_name in lis_root_fname:
                (i_dir, i_name) = os.path.split(i_root_name)
                (i_name, file_type) = os.path.splitext(i_name)
                lis_out_fname.append('%s.pdf' % i_name)
        lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        if len(lis_out_fname) != len(lis_root_fname):
            raise Exception('PrintTFile2PDF: the number of Input.RootNames and Output.PDFNames are not the same.')

        lis_draw_option = kwargs.get('Input.DrawOptions', '')
        dic_draw_option = {}
        if not lis_draw_option == '':
            for i in range(len(lis_draw_option)):
                [class_name, draw_option] = lis_draw_option[i].split(':')
                dic_draw_option[class_name] = draw_option
        else:
            dic_draw_option = None

        if_only_canvas = kwargs.get('Input.OnlyCanvas', 'False')
        if_only_canvas = eval(if_only_canvas)

        # begin to print pdf
        _print_TObject = PrintTObject()
        for i in range(len(lis_root_fname)):
            root_fname = lis_root_fname[i]
            pdf_fname = lis_out_fname[i]

            if dic_draw_option is None:
                _print_TObject.PrintTFile2PDF(root_fname, pdf_fname, if_only_canvas=if_only_canvas)
            else:
                _print_TObject.PrintTFile2PDF(root_fname, pdf_fname, dic_draw_option, if_only_canvas=if_only_canvas)

    def PrintTH2WithError(self, kwargs):
        """print hists with error to pdf file"""
        logger.info('Run PrintTH2WithError.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = kwargs.get('Input.RootDir', '')
        root_dir = os.path.join(self._setting.get_global_dir(), root_dir)
        lis_root_fname = kwargs.get('Input.RootNames')
        lis_root_fname = [os.path.join(root_dir, i) for i in lis_root_fname]

        lis_hname = kwargs.get('Input.HistNames', '')

        if len(lis_hname) != len(lis_root_fname):
            raise Exception('PrintTH2WithError: the number of Input.RootNames and Input.HistNames are not the same.')

        pdf_fname = kwargs.get('Output.PDFName', '')
        if pdf_fname == '':
            pdf_fname = 'hists_with_err.pdf'
        pdf_fname = os.path.join(_directory, pdf_fname)

        draw_option = kwargs.get('Input.DrawOption', '1')
        if_setoptstat = kwargs.get('Input.SetOptStat', '0')

        # begin to print pdf
        _print_TObject = PrintTObject()
        _print_TObject.PrintTH2WithError(lis_root_fname, lis_hname, pdf_fname, draw_option, if_setoptstat=if_setoptstat)

    def DrawBaselineThreshold(self, kwargs):
        """draw baseline and threshold on the energy spectra"""
        logger.info('Run DrawBaselineThreshold.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = kwargs.get('Input.RootDir', '')
        root_dir = os.path.join(self._setting.get_global_dir(), root_dir)
        lis_root_fname = kwargs.get('Input.RootNames')
        lis_root_fname = [os.path.join(root_dir, i) for i in lis_root_fname]

        baseline_dir = kwargs.get('Input.BaselineDir', '')
        baseline_dir = os.path.join(self._setting.get_global_dir(), baseline_dir)
        lis_baseline_fname = kwargs.get('Input.BaselineNames')
        lis_baseline_fname = [os.path.join(baseline_dir, i) for i in lis_baseline_fname]

        threshold_dir = kwargs.get('Input.ThresholdDir', '')
        threshold_dir = os.path.join(self._setting.get_global_dir(), threshold_dir)
        lis_threshold_fname = kwargs.get('Input.ThresholdNames')
        lis_threshold_fname = [os.path.join(threshold_dir, i) for i in lis_threshold_fname]

        lis_out_fname = kwargs.get('Output.PDFNames', '')
        if lis_out_fname == '':
            lis_out_fname = None
        else:
            lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        if len(lis_baseline_fname) != len(lis_root_fname) \
                or len(lis_threshold_fname) != len(lis_baseline_fname):
            raise Exception('DrawBaselineThreshold: the number of Input.RootNames, Input.BaselineNames '
                            'Input.ThresholdNames and Output.PDFNames are not the same.\n'
                            'Input.RootNames %d, Input.BaselineNames %d, Input.ThresholdNames %d' %
                            (len(lis_root_fname), len(lis_baseline_fname), len(lis_threshold_fname)))

        # begin to print pdf
        for i in range(len(lis_root_fname)):
            root_fname = lis_root_fname[i]
            if lis_out_fname is not None:
                pdf_fname = lis_out_fname[i]
            else:
                (root_dir, root_name) = os.path.split(root_fname)
                (root_name, file_type) = os.path.splitext(root_name)
                pdf_fname = root_name + '.pdf'
                pdf_fname = os.path.join(_directory, pdf_fname)
            baseline_fname = lis_baseline_fname[i]
            threshold_fname = lis_threshold_fname[i]

            PrintTObject.DrawBaselineThreshold(root_fname, pdf_fname, baseline_fname, threshold_fname)

    def SeparateFalseTrigger_Plane(self, kwargs):
        """divide false trigger event to another key in the h5 file"""
        logger.info('Run SeparateFalseTrigger_Plane.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        if isinstance(lis_h5_fname, str) and lis_h5_fname.__contains__("re:"):
            pattern = lis_h5_fname.replace("re:", "")
            lis_h5_fname = CommonUtils.search_files_with_re(h5_dir, pattern)
        else:
            lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        lis_h5_out_fname = kwargs.get('Output.Hdf5Names', '')
        if lis_h5_out_fname == '':
            lis_h5_out_fname = None
        else:
            lis_h5_out_fname = [os.path.join(_directory, i) for i in lis_h5_out_fname]

        lis_log_fname = kwargs.get('Output.LoggerName', '')
        if lis_log_fname == '':
            lis_log_fname = None

        for i in range(len(lis_h5_fname)):
            h5_fname = lis_h5_fname[i]
            (h5_dir, h5_name) = os.path.split(h5_fname)
            (h5_name, file_type) = os.path.splitext(h5_name)

            if lis_h5_out_fname is not None:
                out_fname = lis_h5_out_fname[i]
                if not out_fname[-3:] == ".h5":
                    out_fname = out_fname + ".h5"
            else:
                out_fname = 'sep_' + h5_name + '.h5'
                out_fname = os.path.join(_directory, out_fname)

            if lis_log_fname is None:
                out_log_fname = './log/sep_%s.log' % h5_name
            else:
                out_log_fname = lis_log_fname[i]
            out_log_fname = os.path.join(_directory, out_log_fname)
            (log_dir, log_file) = os.path.split(out_log_fname)
            if not os.path.exists(log_dir):
                os.mkdir(log_dir)

            DataframeProcessor.SeparateFalseTrigger_Plane(h5_name, out_fname, h5_dir, out_log_fname)

    def AddUniformTimeTag(self, kwargs):
        """reconstruct uniform timeTag into the h5 file"""
        logger.info('Run AddUniformTimeTag.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        if isinstance(lis_h5_fname, str) and lis_h5_fname.__contains__("re:"):
            pattern = lis_h5_fname.replace("re:", "")
            lis_h5_fname = CommonUtils.search_files_with_re(h5_dir, pattern)
        else:
            lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        lis_h5_out_fname = kwargs.get('Output.Hdf5Names', '')
        if lis_h5_out_fname == '':
            lis_h5_out_fname = None
        else:
            lis_h5_out_fname = [os.path.join(_directory, i) for i in lis_h5_out_fname]

        lis_time_interval = kwargs.get('Input.TimeDeltas', '')
        if isinstance(lis_time_interval, str):
            time_interval = CommonUtils.time_convert(lis_time_interval)
            if time_interval < 0:
                raise Exception('AddUniformTimeTag: TimeDelta < 0 or cannot identify')
            lis_time_interval = []
            for i in range(len(lis_h5_fname)):
                lis_time_interval.append(time_interval)
        else:
            lis_time_interval = [CommonUtils.time_convert(i) for i in lis_time_interval]

        lis_bgtime = kwargs.get('Input.BeginTimes', '')
        if isinstance(lis_bgtime, str) and lis_bgtime == "0":
            lis_bgtime = [0 for i in lis_time_interval]
        else:
            lis_bgtime = [CommonUtils.time_convert(i) for i in lis_bgtime]

        for i in range(len(lis_h5_fname)):
            h5_fname = lis_h5_fname[i]
            (h5_dir, h5_name) = os.path.split(h5_fname)
            (h5_name, file_type) = os.path.splitext(h5_name)

            if lis_h5_out_fname is not None:
                out_fname = lis_h5_out_fname[i]
            else:
                out_fname = 'uniTimeTag_' + h5_name + '.h5'
                out_fname = os.path.join(_directory, out_fname)

            time_delta = lis_time_interval[i]
            begin_time = lis_bgtime[i]

            DataframeProcessor.AddUniformTimeTag(h5_name, out_fname, h5_dir, time_delta, begin_time)

    def CalcHitPoints(self, kwargs):
        """calculate every muon hit point for given h5 file
        and save into tree"""
        logger.info('Run CalcHitPoints.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        lis_det_name = kwargs.get('Input.DetNames')
        if lis_det_name is None or not isinstance(lis_det_name, list):
            raise Exception('CalcHitPoints Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(lis_det_name)

        lis_calc_multi_trigger = kwargs.get('Input.CalcMultiTriggers', '')
        if lis_calc_multi_trigger == '':
            lis_time_interval = 'True'
        if isinstance(lis_calc_multi_trigger, str):
            calc_multi_trigger = eval(lis_calc_multi_trigger)
            lis_calc_multi_trigger = []
            for i in range(len(lis_h5_fname)):
                lis_calc_multi_trigger.append(calc_multi_trigger)
        else:
            lis_calc_multi_trigger = [eval(i) for i in lis_calc_multi_trigger]

        lis_use_det_meas_time = kwargs.get('Input.UseDetMeasTime', '')
        if lis_use_det_meas_time == '':
            lis_use_det_meas_time = 'False'
        if isinstance(lis_use_det_meas_time, str):
            use_det_meas_time = eval(lis_use_det_meas_time)
            lis_use_det_meas_time = []
            for i in range(len(lis_h5_fname)):
                lis_use_det_meas_time.append(use_det_meas_time)
        else:
            lis_use_det_meas_time = [eval(i) for i in lis_use_det_meas_time]

        lis_root_fname = kwargs.get('Output.RootNames', '')
        if lis_root_fname == '':
            lis_root_fname = None
        else:
            lis_root_fname = [os.path.join(_directory, i) for i in lis_root_fname]

        _dataframe_analyzer = DataframeAnalyzer()
        for i in range(len(lis_h5_fname)):
            h5_fname = lis_h5_fname[i]
            (h5_dir, h5_name) = os.path.split(h5_fname)
            (h5_name, file_type) = os.path.splitext(h5_name)

            det = _v_dets[i]
            calc_multi_trigger = lis_calc_multi_trigger[i]
            use_det_meas_time = lis_use_det_meas_time[i]

            if lis_root_fname is not None:
                root_fname = lis_root_fname[i]
            else:
                root_fname = h5_name + '.root'
                root_fname = os.path.join(_directory, root_fname)

            _dataframe_analyzer.CalcHitPoints_4Plane(h5_name, h5_dir, root_fname, det,
                                                     calc_multi_trigger, use_det_meas_time)

    def FixHdf5Data(self, kwargs):
        """fix hdf5 data, the given file will be backup and modified"""
        logger.info('Run FixHdf5Data.')

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        _data_processor = DataframeProcessor()
        for h5_fname in lis_h5_fname:
            (h5_dir, h5_name) = os.path.split(h5_fname)
            _data_processor.FixHdf5Data(h5_name, h5_dir)

    def AddKeyMeasTime(self, kwargs):
        """add a key meas_time in the hdf5 file"""
        logger.info('Run AddKeyMeasTime.')

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        _data_processor = DataframeProcessor()
        for h5_fname in lis_h5_fname:
            (h5_dir, h5_name) = os.path.split(h5_fname)
            _data_processor.AddKeyMeasTime(h5_name, h5_dir)

    def PrintAllCountsMeasTime(self, kwargs):
        """print counts, count rate and measurement time for all files in given directory
                and save as a csv file"""
        logger.info('Run PrintAllCountsMeasTime.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        tree_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootDir', ''))

        fout_name = kwargs.get('Output.CsvName', '')
        if fout_name == '':
            fout_name = kwargs.get('Input.TreeRootDir') + '.csv'
        fout_name = os.path.join(_directory, fout_name)

        PrintTObject.PrintAllCountsMeasTime(tree_dir, fout_name)

    def DoAllStatistics(self, kwargs):
        """count muons in each bin and save in a root file, procedures includes:
            DoStatistics1D_1Layer, DoStatistics1D_2Layer, DoStatistics1D_4Layer, DoStatistics2D"""
        logger.info('Run DoAllStatistics.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        tree_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootDir', ''))
        lis_tr_fname = kwargs.get('Input.TreeRootNames')
        lis_tr_fname = [os.path.join(tree_dir, i) for i in lis_tr_fname]

        lis_det_name = kwargs.get('Input.DetNames')
        if lis_det_name is None or not isinstance(lis_det_name, list):
            raise Exception('DoAllStatistics Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(lis_det_name)

        lis_stats_fname = kwargs.get('Output.StatsNames', '')
        if lis_stats_fname == '':
            lis_stats_fname = None
        else:
            lis_stats_fname = [os.path.join(_directory, i) for i in lis_stats_fname]

        solid_angle_fname = kwargs.get('Input.SolidAngleName', '')
        solid_angle_fname = os.path.join(self._setting.get_global_dir(), solid_angle_fname)

        if_print_pdf = kwargs.get('Input.PrintPdf', 'True')
        if_print_pdf = eval(if_print_pdf)

        if len(lis_tr_fname) != len(lis_det_name):
            raise Exception('MuonStatistics: number of Input.TreeRootNames and Input.DetNames '
                            'must be the same\nTreeRootNames: %d, DetNames: %d.'
                            % (len(lis_tr_fname), len(lis_det_name)))

        _statistic = MuonStatistics.MuonStatistics()
        for i in range(len(lis_tr_fname)):
            tree_fname = lis_tr_fname[i]
            det = _v_dets[i]

            if lis_stats_fname is not None:
                stats_fname = lis_stats_fname[i]
            else:
                stats_fname = 'stats_' + tree_fname
                stats_fname = os.path.join(stats_fname)

            if os.path.exists(stats_fname):
                os.remove(stats_fname)

            # DoStatistics1D
            _statistic.DoStatistics1D_1Layer(self._setting.get_bins(), tree_fname, stats_fname)
            _statistic.DoStatistics1D_2Layer(self._setting.get_bins(), tree_fname, stats_fname)
            _statistic.DoStatistics1D_4Layer(self._setting.get_bins(), tree_fname, stats_fname)

            # DoStatistics2D
            dic_kwargs = {'Output.FileName': solid_angle_fname}
            h_thetaphi_solid_angle, h_thetaxy_solid_angle = self.GetSolidAngle(dic_kwargs)
            _statistic.DoStatistics2D(self._setting.get_bins(), tree_fname, stats_fname,
                                      det, h_thetaphi_solid_angle, h_thetaxy_solid_angle)

            if if_print_pdf:
                lis_draw_option = kwargs.get('Input.DrawOptions', '')
                self._print_root_in_pdf(stats_fname, lis_draw_option)

    def GetMeasRatio(self, kwargs):
        """Calculate measured ratio and write into root file"""
        logger.info('Run GetMeasRatio.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        target_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TargetRootDir', ''))
        lis_tg_fname = kwargs.get('Input.TargetRootNames')
        lis_tg_fname = [os.path.join(target_dir, i) for i in lis_tg_fname]
        background_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.BackgroundRootDir', ''))
        lis_bg_fname = kwargs.get('Input.BackgroundRootNames')
        lis_bg_fname = [os.path.join(background_dir, i) for i in lis_bg_fname]

        lis_ratio_cut = kwargs.get('Input.StatsCutNames', '')
        if lis_ratio_cut == '':
            lis_ratio_cut = ['TOT']

        lis_scale = kwargs.get('Input.RatioScale', '')
        if lis_scale == '':
            lis_scale = None
        elif isinstance(lis_scale, str):
            scale = eval(lis_scale)
            lis_scale = []
            for i in lis_tg_fname:
                lis_scale.append(scale)
        else:
            lis_scale = [eval(i) for i in lis_scale]

        lis_ratio_fname = kwargs.get('Output.RatioRootNames')
        lis_ratio_fname = [os.path.join(_directory, i) for i in lis_ratio_fname]

        lis_log_fname = kwargs.get('Output.LogNames', '')
        if lis_log_fname == '':
            lis_log_fname = None
        else:
            lis_log_fname = [os.path.join(_directory, './log' + i) for i in lis_log_fname]

        # check length of lists
        if len(lis_tg_fname) != len(lis_bg_fname) or len(lis_ratio_fname) != len(lis_tg_fname):
            raise Exception('GetMeasRatio: number of Input.TargetRootNames, Input.BackgroundRootNames '
                            'and Input.RatioNames must be the same\n'
                            'TargetRootNames: %d, BackgroundRootNames: %d, RatioNames: %d'
                            % (len(lis_tg_fname), len(lis_bg_fname), len(lis_ratio_fname)))

        # begin to calculate ratio
        _get_ratio_Leff = GetRatioLeff.GetRatioLeff()
        for i in range(len(lis_tg_fname)):
            tg_fname = lis_tg_fname[i]
            bg_fname = lis_bg_fname[i]
            ratio_fname = lis_ratio_fname[i]
            if lis_log_fname is not None:
                log_fname = lis_log_fname[i]
            else:
                (ratio_out_dir, ratio_fout_fname) = os.path.split(ratio_fname)
                log_dir = os.path.join(ratio_out_dir, 'log')
                if not os.path.exists(log_dir):
                    os.mkdir(log_dir)
                log_fname = os.path.join(log_dir, ratio_fout_fname.replace('.root', '.log'))

            if os.path.exists(ratio_fname):
                os.remove(ratio_fname)

            if lis_scale is None:
                _get_ratio_Leff.GetMeasRatio(tg_stats_fname=tg_fname, bg_stats_fname=bg_fname, ratio_fname=ratio_fname,
                                             ratio_mode=lis_ratio_cut, log_fname=log_fname)
            else:
                scale = lis_scale[i]
                _get_ratio_Leff.GetMeasRatio(tg_stats_fname=tg_fname, bg_stats_fname=bg_fname, ratio_fname=ratio_fname,
                                             ratio_mode=lis_ratio_cut, log_fname=log_fname, ratio_scale=1. / scale)

    def GetMeasLeff(self, kwargs):
        """Calculate measured effective length with given physical model and write into root file"""
        logger.info('Run GetMeasLeff.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        ratio_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.RatioRootDir', ''))
        lis_ratio_fname = kwargs.get('Input.RatioRootNames')
        lis_ratio_fname = [os.path.join(ratio_dir, i) for i in lis_ratio_fname]

        phys_model_name = kwargs.get('Input.PhysModelName', '')
        if phys_model_name == '' or not isinstance(phys_model_name, str):
            raise Exception('GetMeasLeff need one Input.PhysModelName.')

        eps = kwargs.get('Input.Epsilon', '')
        if eps == '':
            eps = 0.0002

        lis_leff_fname = kwargs.get('Output.LeffRootNames')
        lis_leff_fname = [os.path.join(_directory, i) for i in lis_leff_fname]

        lis_log_fname = kwargs.get('Output.LogNames', '')
        if lis_log_fname == '':
            lis_log_fname = None
        else:
            lis_log_fname = [os.path.join(_directory, './log' + i) for i in lis_log_fname]

        # check length of lists
        if len(lis_ratio_fname) != len(lis_leff_fname):
            raise Exception('GetMeasRatio: number of Input.RatioRootNames and Input.LeffRootNames must be the same\n'
                            'RatioRootNames: %d, LeffRootNames: %d' % (len(lis_ratio_fname), len(lis_leff_fname)))

        # begin to calculate ratio
        phys_model = self._get_phys_model(phys_model_name)

        _get_ratio_Leff = GetRatioLeff.GetRatioLeff()
        for i in range(len(lis_ratio_fname)):
            ratio_fname = lis_ratio_fname[i]
            leff_fname = lis_leff_fname[i]
            if lis_log_fname is not None:
                log_fname = lis_log_fname[i]
            else:
                (leff_out_dir, leff_fout_fname) = os.path.split(leff_fname)
                log_dir = os.path.join(leff_out_dir, 'log')
                if not os.path.exists(log_dir):
                    os.mkdir(log_dir)
                log_fname = os.path.join(log_dir, leff_fout_fname.replace('.root', '.log'))

            if os.path.exists(leff_fname):
                os.remove(leff_fname)

            _get_ratio_Leff.GetMeasLeff(ratio_fname, leff_fname, phys_model, log_fname=log_fname, eps=float(eps))

    def CompareHistList(self, kwargs):
        """compare given hists, and draw their significance, KSTest, QQplot in canvas"""
        logger.info('Run CompareHistList.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.HistRootDir'))
        lis_fin_names = [os.path.join(fin_dir, i) for i in kwargs.get('Input.HistRootNames')]
        lis_hist_names = kwargs.get('Input.HistNames')
        lis_scale = kwargs.get('Input.Scales')
        if lis_scale == '1' or lis_scale == '1.' or lis_scale == '1.0':
            lis_scale = [1. for i in lis_fin_names]

        if len(lis_fin_names) != len(lis_hist_names) or len(lis_fin_names) != len(lis_scale):
            raise Exception('number of Input.HistRootNames, Input.HistNames and Input.Scales'
                            'must be the same.\nHistRootNames: %d; HistNames: %d; Scales: %d'
                            % (len(lis_fin_names), len(lis_hist_names), len(lis_scale)))

        if_print_pdf = kwargs.get('Input.PrintPdf', 'True')
        if_print_pdf = eval(if_print_pdf)

        fout_name = kwargs.get('Output.CompareRootName', '')
        if (not isinstance(fout_name, str)) or fout_name == '':
            raise Exception('CompareHistList: need one Output.CompareRootName.')
        fout_name = os.path.join(_directory, fout_name)

        if os.path.exists(fout_name):
            os.remove(fout_name)

        _comparing_tool = ComparingHists.ComparingHists()
        _comparing_tool.CompareHistList(lis_fin_names, lis_hist_names, lis_scale, fout_name)

        if if_print_pdf:
            lis_draw_option = kwargs.get('Input.DrawOptions', '')
            self._print_root_in_pdf(fout_name, lis_draw_option)

    def CalcSignificance(self, kwargs):
        """compare given hists, and draw their significance, KSTest, QQplot in canvas"""
        logger.info('Run CalcSignificance.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        lis_expr = kwargs.get('Input.Expressions')
        dic_input = {key.replace('Input.', ''): kwargs[key] for key in kwargs.keys()
                     if key.__contains__('Input.') and not key.__contains__('Expressions')}
        for key in dic_input.keys():
            val = dic_input[key]
            if key.__contains__('RootName'):
                dic_input[key] = os.path.join(self._setting.get_global_dir(), val)
        lis_fout_name = kwargs.get('Output.SigRootNames')
        lis_fout_name = [os.path.join(_directory, i) for i in lis_fout_name]

        if_smooth = eval(kwargs.get('Input.IfSmooth', 'True'))

        if len(lis_expr) != len(lis_fout_name):
            raise Exception('CalcSignificance: number of Input.Expressions and Output.SigRootNames must'
                            'be the same.\nExpressions: %d; SigRootNames: %d.' % (len(lis_expr), len(lis_fout_name)))

        _comparing_tool = ComparingHists.ComparingHists()
        for i in range(len(lis_expr)):
            i_expr = lis_expr[i]
            fout_name = lis_fout_name[i]

            if os.path.exists(fout_name):
                os.remove(fout_name)

            _comparing_tool.CalcSignificance(i_expr, dic_input, fout_name, if_smooth)

    def CutSignificance(self, kwargs):
        """cut significance with given threshold"""
        logger.info('Run CutSignificance.')

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.RootDir'))
        lis_fname = [os.path.join(root_dir, i) for i in kwargs.get('Input.RootNames')]
        h_name = kwargs.get('Input.HistName')
        high_thr = eval(kwargs.get('Input.SigHighThreshold'))
        low_thr = eval(kwargs.get('Input.SigLowThreshold'))

        _comparingHist = ComparingHists.ComparingHists()

        for i_fname in lis_fname:
            _comparingHist.CutSignificance(i_fname, h_name, low_thr, high_thr)

    def Slice2D(self, kwargs):
        """compare experimental and predicted results, and show in canvas"""
        logger.info('Run Slice2D.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        dic_input = {key.replace('Input.', ''): kwargs[key] for key in kwargs.keys()
                     if key.__contains__('Input.')}
        for key in dic_input.keys():
            val = dic_input[key]
            if not (key.__contains__('RootName') or key.__contains__('StatsName')):
                continue
            if isinstance(val, str):
                dic_input[key] = os.path.join(self._setting.get_global_dir(), val)
            else:
                dic_input[key] = [os.path.join(self._setting.get_global_dir(), i) for i in val]

        lis_coor_name = kwargs.get('Output.CoordNames', '')
        if lis_coor_name == '':
            lis_coor_name = ['thetaphi', 'thetaxy']

        if_print_pdf = kwargs.get('Input.PrintPdf', 'True')
        if_print_pdf = eval(if_print_pdf)

        fout_name = kwargs.get('Output.Slice2DName', '')
        if (not isinstance(fout_name, str)) or fout_name == '':
            raise Exception('Slice2D: need one Output.Slice2DName.')
        fout_name = os.path.join(_directory, fout_name)

        if os.path.exists(fout_name):
            os.remove(fout_name)

        _comparing_tool = ComparingHists.ComparingHists()
        for i in lis_coor_name:
            if (not i.__contains__('phi')) and (not i.__contains__('x')):
                raise Exception('Cannot identify CoordName %s' % i)
            _comparing_tool.Slice2D(dic_input, fout_name, i)

            if if_print_pdf:
                lis_draw_option = kwargs.get('Input.DrawOptions', '')
                self._print_root_in_pdf(fout_name, lis_draw_option)

    def TargetBackgroundSlice(self, kwargs):
        """create hstack of target and background measurement"""
        logger.info('Run TargetBackgroundSlice.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        dic_input = {key.replace('Input.', ''): kwargs[key] for key in kwargs.keys()
                     if key.__contains__('Input.')}
        for key in dic_input.keys():
            val = dic_input[key]
            if (key.__contains__('RootName') or key.__contains__('StatsName')) and val != '':
                dic_input[key] = os.path.join(self._setting.get_global_dir(), val)

        lis_coor_name = kwargs.get('Output.CoordNames', '')
        if lis_coor_name == '':
            lis_coor_name = ['thetaphi', 'thetaxy']

        fout_name = kwargs.get('Output.SliceName', '')
        if (not isinstance(fout_name, str)) or fout_name == '':
            raise Exception('TargetBackgroundSlice: need one Output.SliceName.')
        fout_name = os.path.join(_directory, fout_name)

        if os.path.exists(fout_name):
            os.remove(fout_name)

        _comparing_tool = ComparingHists.ComparingHists()
        for i in lis_coor_name:
            if (not i.__contains__('phi')) and (not i.__contains__('x')):
                raise Exception('Cannot identify CoordName %s' % i)
            _comparing_tool.TargetBackgroundSlice(dic_input, fout_name, i)

    def MeasSimLeffSlice(self, kwargs):
        """create hstack of Leff of prediction and measurement"""
        logger.info('Run MeasSimLeffSlice.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        dic_input = {key.replace('Input.', ''): kwargs[key] for key in kwargs.keys()
                     if key.__contains__('Input.')}
        for key in dic_input.keys():
            val = dic_input[key]
            if (key.__contains__('RootName') or key.__contains__('StatsName')) and val != '':
                if isinstance(val, str):
                    dic_input[key] = os.path.join(self._setting.get_global_dir(), val)
                else:
                    dic_input[key] = [os.path.join(self._setting.get_global_dir(), i) for i in val]

        lis_coor_name = kwargs.get('Output.CoordNames', '')
        if lis_coor_name == '':
            lis_coor_name = ['thetaphi', 'thetaxy']

        fout_name = kwargs.get('Output.SliceName', '')
        if (not isinstance(fout_name, str)) or fout_name == '':
            raise Exception('MeasSimLeffSlice: need one Output.SliceName.')
        fout_name = os.path.join(_directory, fout_name)

        if os.path.exists(fout_name):
            os.remove(fout_name)

        _comparing_tool = ComparingHists.ComparingHists()
        for i in lis_coor_name:
            if (not i.__contains__('phi')) and (not i.__contains__('x')):
                raise Exception('Cannot identify CoordName %s' % i)
            _comparing_tool.MeasSimLeffSlice(dic_input, fout_name, i)

    def PrintSliceTH2(self, kwargs):
        """slice a th2 and print in a pdf file"""
        logger.info('Run PrintSliceTH2.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        hist_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TH2RootName'))
        hist_name = kwargs.get('Input.HistName', '')
        pdf_name = os.path.join(_directory, kwargs.get('Output.PdfName'))
        fout_name = pdf_name.replace('.pdf', '.root')

        if os.path.exists(fout_name):
            os.remove(fout_name)

        _comparing_tool = ComparingHists.ComparingHists()
        _comparing_tool.SliceOneTH2(hist_fname, hist_name, fout_name)

        # noinspection DuplicatedCode
        lis_draw_option = kwargs.get('Input.DrawOptions', '')
        self._print_root_in_pdf(fout_name, lis_draw_option)

    def FillHistsIntoTree(self, kwargs):
        """fill hists into tree for ray processing and save it in a root file"""
        logger.info('Run FillHistsIntoTree.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        dic_input = {key.replace('Input.', ''): kwargs[key] for key in kwargs.keys()
                     if key.__contains__('Input.')}
        for key in dic_input.keys():
            val = dic_input[key]
            if key.__contains__('RootName') and val != '':
                if isinstance(val, str):
                    dic_input[key] = os.path.join(self._setting.get_global_dir(), val)
                else:
                    dic_input[key] = [os.path.join(self._setting.get_global_dir(), i) for i in val]

        lis_det_name = kwargs.get('Input.DetNames')
        if lis_det_name is None or not isinstance(lis_det_name, list):
            raise Exception('FillHistsIntoTree Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(lis_det_name)

        dic_input['Detectors'] = _v_dets

        fout_name = os.path.join(_directory, kwargs.get('Output.TreeRootName'))

        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.FillHistsIntoTree(dic_input, fout_name)

    def MergeTreeBranch(self, kwargs):
        """fill hists into tree for ray processing and save it in a root file"""
        logger.info('Run MergeTreeBranch.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootDir', ''))
        lis_tree_fname = [os.path.join(root_dir, i) for i in kwargs.get('Input.TreeRootNames')]

        fout_name = kwargs.get('Output.TreeRootName')
        fout_name = os.path.join(_directory, fout_name)

        log_dir = os.path.join(_directory, 'log')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        log_fname = os.path.join(log_dir, os.path.splitext(os.path.split(fout_name)[1])[0] + '.log')

        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.MergeTreeBranch(lis_tree_fname, fout_name, log_fname)

    def ContinueTree(self, kwargs):
        """continue trees in files"""
        logger.info('Run ContinueTree.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootDir', ''))
        lis_tree_fname = [os.path.join(root_dir, i) for i in kwargs.get('Input.TreeRootNames')]
        lis_tree_name = kwargs.get('Input.TreeName', '')

        fout_name = kwargs.get('Output.TreeRootName')
        fout_name = os.path.join(_directory, fout_name)

        log_dir = os.path.join(_directory, 'log')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        log_fname = os.path.join(log_dir, os.path.splitext(os.path.split(fout_name)[1])[0] + '.log')

        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.ContinueTree(lis_tree_fname, lis_tree_name, fout_name, log_fname)

    def CopyTreeWithSelection(self, kwargs):
        """merge trees in given files and do selection"""
        logger.info('Run CopyTreeWithSelection.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        selection = kwargs.get('Input.Selection')
        if selection is None or not isinstance(selection, str):
            raise Exception('CopyTreeWithSelection Input.Selection is must be a string.')

        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootDir', ''))
        lis_tree_fname = [os.path.join(root_dir, i) for i in kwargs.get('Input.TreeRootNames')]

        fout_name = kwargs.get('Output.TreeRootName')
        fout_name = os.path.join(_directory, fout_name)

        log_dir = os.path.join(_directory, 'log')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        log_fname = os.path.join(log_dir, os.path.splitext(os.path.split(fout_name)[1])[0] + '.log')

        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.CopyTreeWithSelection(selection, lis_tree_fname, fout_name, log_fname)

    def SampleFromMeasCount(self, kwargs):
        """sample from measured fluxes with poisson distribution"""
        logger.info('Run SampleFromMeasFlux.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.MeasTreeRootDir', ''))
        lis_tree_fname = [os.path.join(root_dir, i) for i in kwargs.get('Input.MeasTreeRootNames')]

        phys_model_name = kwargs.get('Input.PhysModelName', '')
        if phys_model_name == '' or not isinstance(phys_model_name, str):
            raise Exception('SampleFromMeasFlux need one Input.PhysModelName.')

        phys_model = self._get_phys_model(phys_model_name)

        eps = kwargs.get('Input.Epsilon', '')
        if eps == '':
            eps = 0.0002

        dic_input = {'TreeRootNames': lis_tree_fname,
                     'PhysModel': phys_model,
                     'Epsilon': eps}

        fout_name = kwargs.get('Output.TreeRootName')
        fout_name = os.path.join(_directory, fout_name)

        log_dir = os.path.join(_directory, 'log')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        log_fname = os.path.join(log_dir, os.path.splitext(os.path.split(fout_name)[1])[0] + '.log')

        # begin to sample from measurement
        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.SampleFromMeasCount(dic_input, fout_name, log_fname)

    def SampleFromSimCount(self, kwargs):
        """sample from predicted fluxes with poisson distribution"""
        logger.info('Run SampleFromSimFlux.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.SimTreeRootDir', ''))
        tree_fname = os.path.join(root_dir, kwargs.get('Input.SimTreeRootName'))

        phys_model_name = kwargs.get('Input.PhysModelName', '')
        if phys_model_name == '' or not isinstance(phys_model_name, str):
            raise Exception('SampleFromSimFlux need one Input.PhysModelName.')

        phys_model = self._get_phys_model(phys_model_name)

        eps = kwargs.get('Input.Epsilon', '')
        if eps == '':
            eps = 0.0002

        dic_input = {'TreeRootName': tree_fname,
                     'PhysModel': phys_model,
                     'Epsilon': eps}

        fout_name = kwargs.get('Output.TreeRootName')
        fout_name = os.path.join(_directory, fout_name)

        log_dir = os.path.join(_directory, 'log')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        log_fname = os.path.join(log_dir, os.path.splitext(os.path.split(fout_name)[1])[0] + '.log')

        # begin to sample from prediction
        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.SampleFromSimCount(dic_input, fout_name, log_fname)

    def PrintTree2Csv(self, kwargs):
        """print all rays in given tree"""
        logger.info('Run PrintTree2Csv.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        tree_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootName'))
        tree_name = kwargs.get('Input.TreeName', '')
        lis_var_name = kwargs.get('Input.VariableNames')
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName'))

        PrintTObject.PrintTree2Csv(tree_fname, tree_name, lis_var_name, fout_name)

    def TreeThetaxy2Thetaphi(self, kwargs):
        """convert tree branch from thetaxy to thetaphi"""
        logger.info('Run TreeThetaxy2Thetaphi.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        tree_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootName', ''))
        tree_name = kwargs.get('Input.TreeName', '')
        keep_thetaxy = eval(kwargs.get('Input.KeepThetaxy', 'True'))
        fout_name = os.path.join(_directory, kwargs.get('Output.TreeRootName'))

        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.TreeThetaxy2Thetaphi(tree_fname, tree_name, fout_name, if_keep_thetaxy=keep_thetaxy)

    def PhysicalModelFitting(self, kwargs):
        """fitting physical model parameters with given tree
            (the given parameters are used as initialization)"""
        logger.info('Run PhysicalModelFitting.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        tree_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootName'))
        tree_name = kwargs.get('Input.TreeName', '')
        if tree_name == '':
            tree_name = 'tr_selected'
        eps = kwargs.get('Input.Epsilon', '')
        if eps == '':
            eps = 0.0002
        rob = kwargs.get('Input.Robust', '')
        if rob == '':
            rob = 1.
        lis_fix_number = kwargs.get('Input.FixNumbers', '')
        try:
            lis_fix_number = [int(i) for i in lis_fix_number]
            print('Fix Numbers: ', lis_fix_number)
        except:
            lis_fix_number = None

        fout_name = os.path.join(_directory, kwargs.get('Output.LogName'))
        out_phys_model_name = kwargs.get('Output.PhysModelName', 'PF')

        _phys_fitting = PhysFitting.PhysFitting()

        phys_model_name = kwargs.get('Input.PhysModelName', '')
        if phys_model_name == '' or not isinstance(phys_model_name, str):
            raise Exception('SampleFromSimFlux  need one Input.PhysModelName.')

        phys_model = self._get_phys_model(phys_model_name)

        if lis_fix_number is None:
            _phys_fitting.fit_parameters(tree_fname, phys_model, tree_name,
                                         eps=float(eps), fout_name=fout_name, rob=float(rob))
        else:
            _phys_fitting.fit_parameters(tree_fname, phys_model, tree_name,
                                         eps=float(eps), fout_name=fout_name,
                                         rob=float(rob), lis_fix_number=lis_fix_number)
        _phys_fitting.write_parameters(fout_name, out_phys_model_name)

    def CalcEminFromRatioTheta(self, kwargs):
        """calculate Emin from the given ratio and theta"""
        logger.info('Run CalcEminFromRatioTheta.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CsvName', ''))
        df = pd.read_csv(fin_name, index_col=False)
        sep = ','
        if len(df.columns) == 1:
            df = pd.read_csv(fin_name, index_col=False, sep=' ')
            sep = ' '
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName', ''))

        phys_model_name = kwargs.get('Input.PhysModelName', '')

        phys_model = self._get_phys_model(phys_model_name)

        for col in df.columns:
            if col.lower().__contains__('ratio'):
                ratio_name = col

                arr_emin = np.zeros(df.shape[0])
                for i in range(df.shape[0]):
                    theta = df['Theta'].iloc[i]
                    ratio = df[ratio_name].iloc[i]
                    arr_emin[i] = phys_model.get_Emin_from_ratio(theta=theta, ratio=ratio)

                df['%s_Emin' % ratio_name] = arr_emin

        df.to_csv(fout_name, index=False, sep=sep)

    def CalcLeffFromRatioTheta(self, kwargs):
        """calculate Emin from the given ratio and theta"""
        logger.info('Run CalcEminFromRatioTheta.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CsvName', ''))
        df = pd.read_csv(fin_name, index_col=False)
        sep = ','
        if len(df.columns) == 1:
            df = pd.read_csv(fin_name, index_col=False, sep=' ')
            sep = ' '
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName', ''))

        phys_model_name = kwargs.get('Input.PhysModelName', '')

        phys_model = self._get_phys_model(phys_model_name)

        for col in df.columns:
            if col.lower().__contains__('ratio'):
                ratio_name = col

                arr_emin = np.zeros(df.shape[0])
                arr_leff = np.zeros(df.shape[0])
                for i in range(df.shape[0]):
                    theta = df['Theta'].iloc[i]
                    ratio = df[ratio_name].iloc[i]
                    arr_emin[i] = phys_model.get_Emin_from_ratio(theta=theta, ratio=ratio)
                    arr_leff[i] = phys_model.get_Leff_from_Emin(arr_emin[i], phys_model.get_parameters())

                df['%s_Emin' % ratio_name] = arr_emin
                df['%s_Leff' % ratio_name] = arr_leff

        df.to_csv(fout_name, index=False, sep=sep)

    def CalcRatioFromThetaLeff(self, kwargs):
        """calculate ratio from the given leff and theta"""
        logger.info('Run CalcRatioFromThetaLeff.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CsvName', ''))
        df = pd.read_csv(fin_name, index_col=False)
        sep = ','
        if len(df.columns) == 1:
            df = pd.read_csv(fin_name, index_col=False, sep=' ')
            sep = ' '
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName', ''))

        phys_model_name = kwargs.get('Input.PhysModelName', '')

        phys_model = self._get_phys_model(phys_model_name)

        lis_leff_fit = kwargs.get('Input.LeffColumnNames', '')

        for col in lis_leff_fit:
            leff_name = col
            arr_ratio = np.zeros(df.shape[0])
            arr_flux = np.zeros(df.shape[0])

            for i in range(df.shape[0]):
                theta = df['Theta'].iloc[i]
                leff = df[leff_name].iloc[i]
                arr_ratio[i] = phys_model.get_ratio_from_Leff(theta=theta, leff=leff)
                arr_flux[i] = phys_model.get_flux_from_Leff(theta=theta, leff=leff)

            df['%s_Ratio' % leff_name] = arr_ratio
            df['%s_Flux' % leff_name] = arr_flux

        df.to_csv(fout_name, index=False, sep=sep)

    def CalcRatioFromThetaEmin(self, kwargs):
        """calculate ratio from the given emin and theta"""
        logger.info('Run CalcRatioFromThetaEmin.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CsvName', ''))
        df = pd.read_csv(fin_name, index_col=False)
        sep = ','
        if len(df.columns) == 1:
            df = pd.read_csv(fin_name, index_col=False, sep=' ')
            sep = ' '
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName', ''))

        phys_model_name = kwargs.get('Input.PhysModelName', '')

        phys_model = self._get_phys_model(phys_model_name)

        lis_emin_fit = kwargs.get('Input.EminColumnNames', '')

        for col in lis_emin_fit:
            emin_name = col
            arr_ratio = np.zeros(df.shape[0])

            for i in range(df.shape[0]):
                theta = df['Theta'].iloc[i]
                emin = df[emin_name].iloc[i]
                arr_ratio[i] = phys_model.get_ratio_from_Emin(theta=theta, emin=emin)

            df['%s_Ratio' % emin_name] = arr_ratio

        df.to_csv(fout_name, index=False, sep=sep)

    def GetFluxModelRatioHists(self, kwargs):
        """calculate ratio on every theta and energy of flux models"""
        logger.info('Run GetFluxModelRatioHists.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fout_name = os.path.join(_directory, kwargs.get('Output.RootName', ''))

        lis_phys_model_name = kwargs.get('Input.PhysModelNames', '')

        lis_phys_model = []
        for i in lis_phys_model_name:
            phys_model = self._get_phys_model(i)
            lis_phys_model.append(phys_model)

            print('Get %s' % i)

        fout = TFile(fout_name, 'recreate')
        for i, i_name in zip(lis_phys_model, lis_phys_model_name):
            h = i.get_ratio_hist()
            h.SetName(h.GetName() + '_%s' % i_name)
            h.Write()

        fout.Save()
        fout.Close()

    def DrawFluxModel(self, kwargs):
        """draw flux model in TF2 and save in a root file"""
        logger.info('Run DrawFluxModel.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        lis_flux_model = kwargs.get('Input.FluxModelIDs')
        fout_name = os.path.join(_directory, kwargs.get('Output.RootName', ''))

        dic_xy_range = {'e_low': kwargs.get('Input.EnergyLow', ''),
                        'e_up': kwargs.get('Input.EnergyUp', ''),
                        'theta_low': kwargs.get('Input.ThetaLow', ''),
                        'theta_up': kwargs.get('Input.ThetaUp', '')
                        }
        for i_key in dic_xy_range.keys():
            if dic_xy_range.get(i_key) == '':
                dic_xy_range.__delitem__(i_key)

        if os.path.exists(fout_name):
            os.remove(fout_name)

        for flux_id in lis_flux_model:
            dic_phys_model = {'FluxModelFromFormula': flux_id}
            _flux_model = FluxModel.FluxModelFromFormula(dic_phys_model)

            _flux_model.DrawFluxModel(out_fname=fout_name, dic_xy_range=dic_xy_range)

    def SampleFromFluxModel(self, kwargs):
        """sample events from flux model in TH2F and save in a root file"""
        logger.info('Run SampleFromFluxModel.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        lis_flux_model = kwargs.get('Input.FluxModelIDs')
        fout_name = os.path.join(_directory, kwargs.get('Output.RootName', ''))
        n_entries = int(eval(kwargs.get('Input.NEntry', '10000')))
        lis_scale = kwargs.get('Input.Scales', '')
        if lis_scale == '':
            lis_scale = ['1' for i in lis_flux_model]
        elif isinstance(lis_scale, str):
            scale = lis_scale
            lis_scale = [scale for i in lis_flux_model]
        if isinstance(lis_flux_model, list):
            lis_scale = [eval(i) for i in lis_scale]

        dic_xy_range = {'e_low': kwargs.get('Input.EnergyLow', ''),
                        'e_up': kwargs.get('Input.EnergyUp', ''),
                        'theta_low': kwargs.get('Input.ThetaLow', ''),
                        'theta_up': kwargs.get('Input.ThetaUp', ''),
                        'e_bins': kwargs.get('Input.EnergyBin', ''),
                        'theta_bins': kwargs.get('Input.ThetaBin', '')
                        }
        for i_key in dic_xy_range.keys():
            if dic_xy_range.get(i_key) == '':
                dic_xy_range.__delitem__(i_key)

        if os.path.exists(fout_name):
            os.remove(fout_name)

        for flux_id, scale in zip(lis_flux_model, lis_scale):
            dic_phys_model = {'FluxModelFromFormula': flux_id}
            _flux_model = FluxModel.FluxModelFromFormula(dic_phys_model)

            _flux_model.SampleFromFluxModel(out_fname=fout_name, dic_xy_range=dic_xy_range, n_entries=n_entries,
                                            scale=scale)

    def DataframeCalculation(self, kwargs):
        """calculate ratio from the given leff and theta"""
        logger.info('Run DataframeCalculation.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CsvName', ''))
        df = pd.read_csv(fin_name, index_col=False)
        if len(df.columns) == 1:
            df = pd.read_csv(fin_name, index_col=False, sep=' ')
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName', ''))

        lis_expr = kwargs.get('Input.Expressions', '')

        for expr in lis_expr:
            [y, fx] = expr.split('=')
            df[y] = df.eval(fx)

        df.to_csv(fout_name, index=False, sep=' ')

    def RunInversion(self, kwargs):
        """func(xxx.yalm)"""
        logger.info('Run RunInversion.')

        cfg_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CfgName'))

        from Muon_Imaging_Algorithm.inversion import start_inversion_by_setting

        start_inversion_by_setting(cfg_name)

    @classmethod
    def _generate_topo_for_ubc(cls, origin_of_coordinates, x_directions, y_directions, z_directions, out_fname):
        def _get_random_points_one_face(dic_points_setting: dict, n_points: int):
            def _get_one_range_points(range, n_points):
                sc = 10
                p1 = expon.rvs(scale=(range[1] - range[0]) / sc, loc=range[0], size=n_points)
                p1 = np.where(p1 < range[1], p1, range[1])
                p2 = (range[0] + range[1]) - p1
                return p1, p2

            dic_one_dim = {'x': [], 'y': [], 'z': []}

            if not isinstance(dic_points_setting['x'], list):
                dic_one_dim['x'].append(np.array([dic_points_setting['x']] * n_points))
                dic_one_dim['x'].append(np.array([dic_points_setting['x']] * n_points))
            else:
                x1, x2 = _get_one_range_points(dic_points_setting['x'], n_points)
                dic_one_dim['x'].append(x1)
                dic_one_dim['x'].append(x2)
            if not isinstance(dic_points_setting['y'], list):
                dic_one_dim['y'].append(np.array([dic_points_setting['y']] * n_points))
                dic_one_dim['y'].append(np.array([dic_points_setting['y']] * n_points))
            else:
                x1, x2 = _get_one_range_points(dic_points_setting['y'], n_points)
                dic_one_dim['y'].append(x1)
                dic_one_dim['y'].append(x2)
            if not isinstance(dic_points_setting['z'], list):
                dic_one_dim['z'].append(np.array([dic_points_setting['z']] * n_points))
                dic_one_dim['z'].append(np.array([dic_points_setting['z']] * n_points))
            else:
                x1, x2 = _get_one_range_points(dic_points_setting['z'], n_points)
                dic_one_dim['z'].append(x1)
                dic_one_dim['z'].append(x2)

            xyz = None
            for ix in dic_one_dim['x']:
                for iy in dic_one_dim['y']:
                    for iz in dic_one_dim['z']:
                        if xyz is None:
                            xyz = np.concatenate((np.array([ix]).T, np.array([iy]).T, np.array([iz]).T), 1)
                        else:
                            new = np.concatenate((np.array([ix]).T, np.array([iy]).T, np.array([iz]).T), 1)
                            xyz = np.concatenate((xyz, new), 0)

            return xyz

        xmin = origin_of_coordinates[0]
        ymin = origin_of_coordinates[1]
        zmin = origin_of_coordinates[2]

        xmax = xmin
        for i_lis in x_directions:
            xmax = xmax + i_lis[0] * i_lis[1]

        ymax = ymin
        for i_lis in y_directions:
            ymax = ymax + i_lis[0] * i_lis[1]

        zmax = zmin
        for i_lis in z_directions:
            zmax = zmax + i_lis[0] * i_lis[1]

        dic_points = {'top': _get_random_points_one_face(
            {'x': [xmin, xmax], 'y': [ymin, ymax], 'z': zmax}, 1000),
            'left': _get_random_points_one_face({'x': xmin, 'y': [ymin, ymax], 'z': [zmin, zmax]},
                                                400),
            'right': _get_random_points_one_face({'x': xmax, 'y': [ymin, ymax], 'z': [zmin, zmax]},
                                                 400),
            'front': _get_random_points_one_face({'x': [xmin, xmax], 'y': ymin, 'z': [zmin, zmax]},
                                                 400),
            'back': _get_random_points_one_face({'x': [xmin, xmax], 'y': ymax, 'z': [zmin, zmax]},
                                                400)}

        xyz = np.concatenate((dic_points['top'],
                              dic_points['left'],
                              dic_points['right'],
                              dic_points['front'],
                              dic_points['back']), 0)

        np.savetxt(out_fname, xyz, delimiter=" ", fmt='%.6f', header='%d' % xyz.shape[0])

    @classmethod
    def _extract_mesh_expression(cls, expr):
        directions = []
        lis_term = expr.split('+')

        for it in lis_term:
            directions.append([eval(i) for i in it.split('*')])

        return directions

    def GenerateMesh(self, kwargs):
        """generate mesh tools"""
        logger.info('Run GenerateMesh.')

        from Muon_Imaging_Algorithm.InvDataTools.MeshTools import MeshTools

        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        origin_of_coordinates = [eval(i) for i in kwargs.get('Input.Origin')]

        x_expr = kwargs.get('Input.XExpression')
        x_directions = self._extract_mesh_expression(x_expr)
        y_expr = kwargs.get('Input.YExpression')
        y_directions = self._extract_mesh_expression(y_expr)
        z_expr = kwargs.get('Input.ZExpression')
        z_directions = self._extract_mesh_expression(z_expr)

        msh_name = os.path.join(_directory, kwargs.get('Output.MeshName'))

        res = MeshTools.generate_mesh_file(origin_of_coordinates, x_directions, y_directions, z_directions, msh_name)
        print(res)

        out_fname = os.path.splitext(msh_name)[0] + '_topo.txt'
        self._generate_topo_for_ubc(origin_of_coordinates, x_directions, y_directions, z_directions, out_fname)

    def FillTopoDensity(self, kwargs):
        """1.  func1(msh_name, f_mesh_pos)
        for topo in lis_topo:
        2. topo.print_intersection(f_mesh_pos, f_inter, f_obs)
        3. func2(f_inter, f_obs, msh_name, ref: float, bnd: list[low, up], ref_name, bnd_name)"""
        logger.info('Run FillTopoDensity')

        from Muon_Imaging_Algorithm.inversion import get_bottom_cells_xyz_coordinate, \
            build_refs_bounds_file_from_intersection_information

        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        msh_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.MeshName'))
        strategy: str = kwargs.get('Input.Strategy', '')

        air_ref = kwargs.get('Input.AirRef', '')
        if air_ref == '':
            air_ref = 0.
        air_ref = float(air_ref)
        air_bnd = kwargs.get('Input.AirBound', '')
        if not isinstance(air_bnd, list) or (isinstance(air_bnd, list) and len(air_bnd) != 2):
            air_bnd = [-0.00001, 0.00001]
            logger.warning('Use default air density bound:-0.00001, 0.00001')
        air_bnd = (float(air_bnd[0]), float(air_bnd[1]))

        bottom_name = os.path.join(_directory, kwargs.get('Output.BottomCoordName'))
        out_inter_fname = os.path.join(_directory, kwargs.get('Output.IntersectionName'))
        out_obs_fname = os.path.join(_directory, kwargs.get('Output.DirectionName'))
        ref_fname = os.path.join(_directory, kwargs.get('Output.DensityRefName'))
        bnd_fname = os.path.join(_directory, kwargs.get('Output.DensityBoundName'))

        try:
            strategy = tuple(int(i) for i in strategy)
            # if not (isinstance(strategy, tuple) and len(strategy) == 3):
            #     strategy = ''
            if not (len(strategy) == 3):
                raise Exception("strategy长度不是3")
        except Exception as e:
            strategy = ''
            print(e)

        if strategy == '':
            get_bottom_cells_xyz_coordinate(msh_name, bottom_name)
        else:
            get_bottom_cells_xyz_coordinate(msh_name, bottom_name, strategy)

        lis_topo_names = kwargs.get('Input.TopoNames')
        dic_topo_info = self._setting.get_topos_info()
        topo_dir = os.path.join(self._setting.get_global_dir(), dic_topo_info.get('Directory', ''))

        if os.path.exists(ref_fname):
            os.remove(ref_fname)
        if os.path.exists(bnd_fname):
            os.remove(bnd_fname)

        for i_topo in lis_topo_names:
            den_ref = kwargs.get('%s.DensityRef' % i_topo, '')
            if den_ref == '':
                den_ref = 2.65
                logger.warning('Use default density reference for %s: 2.65.' % i_topo)
            den_ref = float(den_ref)
            den_bnd = kwargs.get('%s.DensityBound' % i_topo, '')
            if not isinstance(den_bnd, list) or (isinstance(den_bnd, list) and len(den_bnd) != 2):
                raise Exception('Use default density bound for %s: 0,3.' % i_topo)
            den_bnd = (float(den_bnd[0]), float(den_bnd[1]))

            # get topography info
            i_topo_info = dic_topo_info.get(i_topo)
            if i_topo_info is None:
                Exception('Cannot find %s in given topography info file' % i_topo)

            fin_name = os.path.join(topo_dir, i_topo_info.get('ObjFile'))

            search_all = eval(i_topo_info.get('SearchAll', 'True'))
            if search_all == '':
                search_all = True
            search_dimension = i_topo_info.get('SearchDimension', '')
            if search_dimension != '':
                search_dimension = [float(i) for i in search_dimension]

            if search_dimension != '':
                _topo = Topography.Topography(fin_name, search_dimension[0], search_dimension[1],
                                              search_dimension[2], topo_name=i_topo)
            else:
                _topo = Topography.Topography(fin_name, topo_name=i_topo)
            _topo.read_topo_obj(search_all=search_all)

            _topo.print_intersection(bottom_name, out_inter_fname=out_inter_fname, out_obs_fname=out_obs_fname)

            build_refs_bounds_file_from_intersection_information(intersection_points_file=out_inter_fname,
                                                                 obs_file=out_obs_fname,
                                                                 mesh_file=msh_name,
                                                                 refs_file=ref_fname,
                                                                 bounds_file=bnd_fname,
                                                                 refs_value=den_ref,
                                                                 bounds_value=den_bnd,
                                                                 default_refs_value=air_ref,
                                                                 default_bounds_value=air_bnd)

    """utils"""

    def Hist2Csv(self, kwargs):
        """print content and error of a hist into given file"""
        logger.info('Run Hist2Csv.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        hist_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.HistRootName', ''))
        hist_name = kwargs.get('Input.HistName', '')
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName', ''))

        if os.path.exists(fout_name):
            os.remove(fout_name)

        lis_h_name = ROOTUtils.find_all_keys_in_file(hist_fname, hist_name)
        if len(lis_h_name) < 1:
            raise Exception('Cannot find %s in %s' % (hist_name, hist_fname))

        for h_name in lis_h_name:
            (t_fout_name, fout_type) = os.path.splitext(fout_name)
            t_fout_name = t_fout_name + '_%s' % h_name + fout_type
            ROOTUtils.hist2csv(hist_fname, h_name, t_fout_name)

    def CutCsvWithSelection(self, kwargs):
        """select events in the given dataframe"""
        logger.info('Run CutCsvWithSelection.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        csv_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CsvName', ''))
        selection = kwargs.get('Input.HistNames', '')
        fout_name = os.path.join(_directory, kwargs.get('Output.CsvName', ''))

        df = pd.read_csv(csv_name, index_col=False)
        df = CommonUtils.cut_dataframe_with_selection(df, selection)

        df.to_csv(fout_name, index_label=False)

    def DrawGraphWithCsv(self, kwargs):
        """draw a graph with the given csv file and the column names"""
        logger.info('Run DrawGraphWithCsv.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        csv_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.CsvName', ''))
        lis_col_name = kwargs.get('Input.ColumnNames')
        if not isinstance(lis_col_name, list):
            raise Exception('Input ColumnNames as a list.')
        gr_name = kwargs.get('Output.GraphName', '')
        fout_name = os.path.join(_directory, kwargs.get('Output.RootName', ''))

        fit_func = kwargs.get('Input.FitFunction', '')

        # for itheta draw
        # for i in range(1, 31):
        #     i_csv_name = csv_name.replace('itheta1', 'itheta%d' % i)
        #     i_gr_name = gr_name.replace('itheta1', 'itheta%d' % i)
        #     Utils.draw_graph_with_csv(i_csv_name, lis_col_name, fout_name, i_gr_name)

        if fit_func == '':
            ROOTUtils.draw_graph_with_csv(csv_name, lis_col_name, fout_name, gr_name)
        else:
            ROOTUtils.draw_graph_with_csv(csv_name, lis_col_name, fout_name, gr_name, fit_func)

    """hist operation"""

    def HistMultiple(self, kwargs):
        """hists multiple a hist"""
        logger.info('Run HistMultiple.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        lis_fname = [os.path.join(self._setting.get_global_dir(), i) for i in kwargs.get('Input.RootNames')]
        lis_h_name = kwargs.get('Input.HistNames')

        multi_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.MultiRootName', ''))
        multi_h_name = kwargs.get('Input.MultiHistName')

        multi_h_name = ROOTUtils.find_key_in_file(multi_fname, multi_h_name)
        fin1 = TFile(multi_fname, 'read')
        h = fin1.Get(multi_h_name)

        for i_fname, i_hname in zip(lis_fname, lis_h_name):
            i_hname = ROOTUtils.find_key_in_file(i_fname, i_hname)
            fin = TFile(i_fname, 'update')
            ih = fin.Get(i_hname)

            h_multi = ih.Clone()
            h_multi.SetName(ih.GetName() + '_multi_' + multi_h_name)
            h_multi.Multiply(h)

            fin.cd()
            h_multi.Write()

            fin.Close()

        fin1.Close()

    def HistAdd(self, kwargs):
        """hists add a hist"""
        logger.info('Run HistAdd.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        lis_fname = [os.path.join(self._setting.get_global_dir(), i) for i in kwargs.get('Input.RootNames')]
        lis_h_name = kwargs.get('Input.HistNames')

        add_fname = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.AddRootName', ''))
        add_h_name = kwargs.get('Input.AddHistName')
        scale = eval(kwargs.get('Input.Scale'))

        fout_name = kwargs.get('Output.RootName', '')
        if isinstance(fout_name, str) and fout_name != '':
            fout_name = os.path.join(_directory, fout_name)
        else:
            fout_name = None

        add_h_name = ROOTUtils.find_key_in_file(add_fname, add_h_name)
        fin1 = TFile(add_fname, 'read')
        h = fin1.Get(add_h_name)

        for i_fname, i_hname in zip(lis_fname, lis_h_name):
            i_hname = ROOTUtils.find_key_in_file(i_fname, i_hname)
            fin = TFile(i_fname, 'update')
            ih = fin.Get(i_hname)

            h_add = ih.Clone()
            h_add.SetName(ih.GetName() + '_add_' + add_h_name)
            h_add.Add(h, scale)

            if fout_name is not None:
                fout = TFile(fout_name, 'recreate')
                fout.cd()
            else:
                fin.cd()
            h_add.Write()

            fin.Close()

        fin1.Close()

    def DrawSigPNG4Paper(self, kwargs):
        """output png of significance for paper publishing"""
        logger.info('Run DrawSigPNG4Paper.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        f_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.RootDir', ''))
        lis_fname = kwargs.get('Input.RootNames')
        lis_fout_name = kwargs.get('Output.PNGNames', '')
        if not isinstance(lis_fout_name, list):
            lis_fout_name = [i.replace('.root', '.png') for i in lis_fname]
        lis_fout_name = [os.path.join(_directory, i) for i in lis_fout_name]
        lis_fname = [os.path.join(f_dir, i) for i in lis_fname]
        lis_h_name = kwargs.get('Input.HistNames')
        lis_det_name = kwargs.get('Input.DetNames')
        lis_det = self._get_detector_list(lis_det_name)

        high_sig = eval(kwargs.get('Input.HighSig', '2'))
        low_sig = eval(kwargs.get('Input.LowSig', '-2'))
        xrange = kwargs.get('Input.XRangeUser', '')
        if isinstance(xrange, list):
            xrange = [eval(i) for i in xrange]
        else:
            xrange = None
        yrange = kwargs.get('Input.YRangeUser', '')
        if isinstance(yrange, list):
            yrange = [eval(i) for i in yrange]
        else:
            yrange = None

        if len(lis_fname) != len(lis_fout_name) \
                or len(lis_h_name) != len(lis_fname) \
                or len(lis_det) != len(lis_fname):
            raise Exception('Please keep the numbers of '
                            'RootNames, HistNames, DetNames, PNGNames are the same.')

        _paper_draw = PaperDraw.PaperDraw()
        for i in range(len(lis_fname)):
            fin_name = lis_fname[i]
            h_name = lis_h_name[i]
            png_name = lis_fout_name[i]
            det = lis_det[i]
            rot_phi = det.get_rotAngle()[2]

            _paper_draw.DrawSigPNG4Paper(fin_name, h_name, png_name,
                                         rot_phi, xrange, yrange, low_sig, high_sig)

    def DrawRatioLeff4Paper(self, kwargs):
        """output png of ratio/leff for paper publishing"""
        logger.info('Run DrawRatioLeff4Paper.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        meas_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.MeasRootDir', ''))
        lis_meas_fname = kwargs.get('Input.MeasRootNames')
        lis_fout_name = kwargs.get('Output.PNGNames', '')
        if not isinstance(lis_fout_name, list):
            lis_fout_name = [i.replace('.root', '.png') for i in lis_meas_fname]
        lis_fout_name = [os.path.join(_directory, i) for i in lis_fout_name]
        lis_meas_fname = [os.path.join(meas_dir, i) for i in lis_meas_fname]
        sim_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.SimRootDir', ''))
        lis_sim_fname = [os.path.join(sim_dir, i) for i in kwargs.get('Input.SimRootNames')]
        lis_det_name = kwargs.get('Input.DetNames')
        lis_det = self._get_detector_list(lis_det_name)

        lis_cont = kwargs.get('Input.ContLists', '')
        if isinstance(lis_cont, list):
            lis_cont = [eval(i) for i in lis_cont]
        else:
            lis_cont = None
        xrange = kwargs.get('Input.XRangeUser', '')
        if isinstance(xrange, list):
            xrange = [eval(i) for i in xrange]
        else:
            xrange = None
        yrange = kwargs.get('Input.YRangeUser', '')
        if isinstance(yrange, list):
            yrange = [eval(i) for i in yrange]
        else:
            yrange = None

        if len(lis_meas_fname) != len(lis_fout_name) \
                or len(lis_sim_fname) != len(lis_meas_fname) \
                or len(lis_det) != len(lis_meas_fname):
            raise Exception('Please keep the numbers of '
                            'MeasRootNames, SimRootNames, DetNames, PNGNames are the same.')

        _paper_draw = PaperDraw.PaperDraw()
        for i in range(len(lis_meas_fname)):
            fin_meas_name = lis_meas_fname[i]
            h_meas_name = 'htheta_phi'
            fin_sim_name = lis_sim_fname[i]
            h_sim_name = 'htheta_phi'
            png_name = lis_fout_name[i]
            det = lis_det[i]
            rot_phi = det.get_rotAngle()[2]

            _paper_draw.DrawRatioLeff4Paper(fin_meas_name, h_meas_name, fin_sim_name, h_sim_name, png_name, rot_phi,
                                            lis_x_range=xrange, lis_y_range=yrange, lis_contours=lis_cont)

    def TH1sComparison(self, kwargs):
        """output png of one-dimensional hists comparison"""
        logger.info('Run TH1sComparison.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        f_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.HistRootDir', ''))
        lis_fname = kwargs.get('Input.HistRootNames')
        lis_fname = [os.path.join(f_dir, i) for i in lis_fname]
        lis_hname = kwargs.get('Input.HistNames')
        png_name = os.path.join(_directory, kwargs.get('Output.PNGName', ''))

        lis_scale = kwargs.get('Input.Scales', '')
        if isinstance(lis_scale, list):
            lis_scale = [eval(i) for i in lis_scale]
            if len(lis_scale) != len(lis_fname):
                raise Exception('Please keep the numbers of RootNames and Scales the same.')
        else:
            lis_scale = None
        lis_label = kwargs.get('Input.Labels', '')
        if not isinstance(lis_label, list):
            lis_label = None
        else:
            if len(lis_label) != len(lis_fname):
                raise Exception('Please keep the numbers of RootNames and Labels the same.')
        legend_columns = eval(kwargs.get('Input.LegendCol', '3'))
        xtitle = kwargs.get('Input.XaxisTitle', '')
        ytitle = kwargs.get('Input.YaxisTitle', '')
        horg = eval(kwargs.get('Input.Hist2Graph', 'True'))

        if len(lis_fname) != len(lis_hname):
            raise Exception('Please keep the numbers of '
                            'RootNames, HistNames are the same.')

        _paper_draw = PaperDraw.PaperDraw()
        _paper_draw.TH1sComparison(lis_fname, lis_hname, png_name, lis_scale, lis_label, legend_columns,
                                   xtitle, ytitle, h_or_g=horg)

    def DrawMultiGraphs(self, kwargs):
        """output png of one-dimensional hists comparison"""
        logger.info('Run DrawMultiGraphs.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        f_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.GraphRootDir', ''))
        lis_fname = kwargs.get('Input.GraphRootNames')
        lis_fname = [os.path.join(f_dir, i) for i in lis_fname]
        lis_hname = kwargs.get('Input.GraphNames')
        png_name = os.path.join(_directory, kwargs.get('Output.PNGName', ''))

        lis_label = kwargs.get('Input.Labels', '')
        if not isinstance(lis_label, list):
            lis_label = None
        else:
            if len(lis_label) != len(lis_fname):
                raise Exception('Please keep the numbers of RootNames and Labels the same.')
        legend_columns = eval(kwargs.get('Input.LegendCol', '3'))
        xtitle = kwargs.get('Input.XaxisTitle', '')
        ytitle = kwargs.get('Input.YaxisTitle', '')

        if len(lis_fname) != len(lis_hname):
            raise Exception('Please keep the numbers of '
                            'RootNames, HistNames are the same.')

        _paper_draw = PaperDraw.PaperDraw()
        _paper_draw.DrawMultiGraphs(lis_fname, lis_hname, png_name, lis_label, legend_columns,
                                    xtitle, ytitle)

    # test
    def ComparePhysicalModelAndMeasuredData(self, kwargs):
        logger.info('Run ComparePhysicalModelAndMeasuredData.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_name = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootName', ''))
        tr_name = kwargs.get('Input.TreeName')
        lis_det_names = kwargs.get('Input.DetNames', '')
        lis_det = self._get_detector_list(lis_det_names)
        lis_det_id = [i.get_det_id() for i in lis_det]
        lis_phys_names = kwargs.get('Input.PhysModelNames', '')
        lis_phys_model = []
        for i in lis_phys_names:
            phys_model = self._get_phys_model(i)
            lis_phys_model.append(phys_model)

        fout_name = os.path.join(_directory, kwargs.get('Output.RootName', ''))
        _ray_analyzer = RayAnalyzer.RayAnalyzer()
        _ray_analyzer.ComparePhysicalModelAndMeasuredData(fin_name, tr_name, lis_det_id, lis_phys_model, fout_name,
                                                          self._bins_dic)

        lis_draw_option = kwargs.get('Input.DrawOptions', '')
        self._print_root_in_pdf(fout_name, lis_draw_option)

    def PrintSpectrumSlides(self, kwargs):
        logger.info('Run PrintSpectrumSlides.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.SpectraSlideRootDir', ''))
        fin_name = os.path.join(fin_dir, kwargs.get('Input.SpectraSlideRootName', ''))
        lis_slide_number = kwargs.get('Input.SlideNumbers')
        if_smooth = eval(kwargs.get('Input.IfSmooth'))

        fout_name = os.path.join(_directory, kwargs.get('Output.RootName', ''))

        _comparing_tool = ComparingHists.ComparingHists()
        _comparing_tool.PrintSpectrumSlides(fin_name, lis_slide_number, fout_name, if_smooth)

        lis_draw_option = kwargs.get('Input.DrawOptions', '')
        self._print_root_in_pdf(fout_name, lis_draw_option)

    def Convert2OneEventPerRow(self, kwargs):
        logger.info('Run Convert2OneEventPerRow.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        fin_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_fin_name = kwargs.get('Input.Hdf5Names')
        lis_fout_name = kwargs.get('Output.RootNames', '')
        if not isinstance(lis_fout_name, list):
            lis_fout_name = [i.replace('.h5', '.root') for i in lis_fin_name]
            lis_fout_name = [i + '.h5' for i in lis_fout_name if not i.__contains__('.h5')]
        lis_fout_name = [os.path.join(_directory, i) for i in lis_fout_name]
        lis_fin_name = [os.path.join(fin_dir, i) for i in lis_fin_name]

        tr_name = kwargs.get('Input.TreeName', 'tr_data')

        if len(lis_fin_name) != len(lis_fout_name):
            raise Exception('CalcSignificance: number of Input.Hdf5Names and Output.RootNames must'
                            'be the same.\nHdf5Names: %d; RootNames: %d.' % (len(lis_fin_name), len(lis_fout_name)))

        _data_preprocessing_root = DataFormatTransformer()
        for i in range(len(lis_fin_name)):
            in_fname = lis_fin_name[i]
            out_fname = lis_fout_name[i]
            _data_preprocessing_root.Convert2OneEventPerRow(in_fname, out_fname, tr_name)

    def CalcHitPoints_NPlanes(self, kwargs):
        """calculate every muon hit point for given h5 file
        and save into tree"""
        logger.info('Run CalcHitPoints.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_h5_fname = kwargs.get('Input.Hdf5Names')
        lis_h5_fname = [os.path.join(h5_dir, i) for i in lis_h5_fname]

        lis_det_name = kwargs.get('Input.DetNames')
        if lis_det_name is None or not isinstance(lis_det_name, list):
            raise Exception('CalcHitPoints Input.DetNames is must be a list.')
        _v_dets = self._get_detector_list(lis_det_name)

        lis_calc_multi_trigger = kwargs.get('Input.CalcMultiTriggers', '')
        if lis_calc_multi_trigger == '':
            lis_time_interval = 'True'
        if isinstance(lis_calc_multi_trigger, str):
            calc_multi_trigger = eval(lis_calc_multi_trigger)
            lis_calc_multi_trigger = []
            for i in range(len(lis_h5_fname)):
                lis_calc_multi_trigger.append(calc_multi_trigger)
        else:
            lis_calc_multi_trigger = [eval(i) for i in lis_calc_multi_trigger]

        lis_use_det_meas_time = kwargs.get('Input.UseDetMeasTime', '')
        if lis_use_det_meas_time == '':
            lis_use_det_meas_time = 'False'
        if isinstance(lis_use_det_meas_time, str):
            use_det_meas_time = eval(lis_use_det_meas_time)
            lis_use_det_meas_time = []
            for i in range(len(lis_h5_fname)):
                lis_use_det_meas_time.append(use_det_meas_time)
        else:
            lis_use_det_meas_time = [eval(i) for i in lis_use_det_meas_time]

        lis_root_fname = kwargs.get('Output.RootNames', '')
        if lis_root_fname == '':
            lis_root_fname = None
        else:
            lis_root_fname = [os.path.join(_directory, i) for i in lis_root_fname]

        for i in range(len(lis_h5_fname)):
            h5_fname = lis_h5_fname[i]

            det = _v_dets[i]
            calc_multi_trigger = lis_calc_multi_trigger[i]
            use_det_meas_time = lis_use_det_meas_time[i]

            _calc = CalculatePositionsDF.CalculatePositionsDF_PlaneA(det)

            if lis_root_fname is not None:
                root_fname = lis_root_fname[i]
            else:
                root_fname = os.path.splitext(os.path.split(h5_fname)[1])[0] + '.root'
                root_fname = os.path.join(_directory, root_fname)

            _calc.CalcHitPoints_NPlane(h5_fname, root_fname, det,
                                       calc_multi_trigger, use_det_meas_time)

    def PrintMeasTime2Csv(self, kwargs):
        """print df_meas_time to csv file"""
        logger.info('Run PrintMeasTime2Csv.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        h5_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.Hdf5Dir', ''))
        lis_in_fname = kwargs.get('Input.Hdf5Names')
        lis_in_fname = [os.path.join(h5_dir, i) for i in lis_in_fname]

        for fin_name in lis_in_fname:
            if not fin_name.__contains__(".h5"):
                fin_name += ".h5"
            fout_name = os.path.split(fin_name)[1].replace(".h5", "_meas_time.csv")
            fout_name = os.path.join(_directory, fout_name)

            df_meas_time = DataframeAnalyzer.get_df_meas_time(fin_name)
            df_meas_time.to_csv(fout_name)

    def AddSelfDefinedStatistics(self, kwargs):
        """calculate every muon hit point for given h5 file
        and save into tree"""
        logger.info('Run AddSelfDefinedStatistics.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.RootDir', ''))
        lis_in_fname = kwargs.get('Input.RootNames')
        lis_in_fname = [os.path.join(root_dir, i) for i in lis_in_fname]

        csv_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.MeasTimeCsvDir', ''))
        lis_csv_fname = kwargs.get('Input.MeasTimeCsvNames')
        lis_csv_fname = [os.path.join(csv_dir, i) for i in lis_csv_fname]

        lis_var_exp = kwargs.get('Input.VariableExpressions')
        lis_bin_exp = kwargs.get('Input.BinNameExpressions')
        lis_cut_exp = kwargs.get('Input.CutExpressions')
        lis_cut_name = kwargs.get('Input.CutNames')

        # check the length of each list
        if len(lis_var_exp) != len(lis_bin_exp):
            raise Exception("The length of VariableExpressions and BinNameExpressions should be the same.")
        if len(lis_cut_exp) != len(lis_cut_name):
            raise Exception("The length of CutExpressions and CutExpressions should be the same.")

        time_interval = kwargs.get('Input.TimeInterval', '')
        time_interval_int = CommonUtils.time_convert(time_interval)
        if_print_pdf = kwargs.get('Input.PrintPdf', 'True')
        if_print_pdf = eval(if_print_pdf)

        lis_out_fname = kwargs.get('Output.RootNames', '')
        if lis_out_fname == '':
            lis_out_fname = [i.replace(".root", "_stats.root") for i in lis_in_fname]
        else:
            lis_out_fname = [os.path.join(_directory, i) for i in lis_out_fname]

        _statistic = MuonStatistics.MuonStatistics()

        for i in range(len(lis_in_fname)):
            in_fname = lis_in_fname[i]
            out_fname = lis_out_fname[i]
            try:
                csv_fname = lis_csv_fname[i]
            except:
                csv_fname = ''
            _statistic.AddSelfDefinedStatistics(self._setting.get_bins(), in_fname, out_fname,
                                                lis_var_exp, lis_bin_exp,
                                                lis_cut_exp, lis_cut_name, csv_fname, time_interval_int)

            if if_print_pdf:
                lis_draw_option = kwargs.get('Input.DrawOptions', '')
                self._print_root_in_pdf(out_fname, lis_draw_option)

    def SelectBinsInHist(self, kwargs):
        """select bins in hists"""
        logger.info('Run SelectBinsInHist.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.RootDir', ''))
        lis_in_fname = kwargs.get('Input.RootNames')
        lis_in_fname = [os.path.join(root_dir, i) for i in lis_in_fname]
        lis_hist_name = kwargs.get('Input.HistNames')

        sel_exp = kwargs.get('Input.SelectionExpression')

        for i_fname in lis_in_fname:
            fin = TFile(i_fname, "read")
            fout = TFile(i_fname.repalce(".root", "_bool.root"), "recreate")
            for hist_name in lis_hist_name:
                h = fin.Get(hist_name)
                h_bool = ROOTUtils.convert_hist_to_bool_with_selection(h, sel_exp)
                h_bool.SetName(h_bool.GetName() + "_bool")

                fout.cd()
                h_bool.Write()
            fout.Save()
            fout.Close()
            fin.Close()

    def SelectTreeEventWithHist(self, kwargs):
        """select tree events in hists"""
        logger.info('Run SelectTreeEventWithHist.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        tree_file = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TreeRootName', ''))
        hist_file = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.HistRootName', ''))
        tree_name = kwargs.get('Input.TreeName')
        hist_name = kwargs.get('Input.HistName')

        sel_exp = kwargs.get('Input.ValueSelectionExpression')
        lis_xyname = kwargs.get('Input.XYNames')

        fin1 = TFile(hist_file, "read")
        h = fin1.Get(hist_name)

        fin2 = TFile(tree_file, "read")
        tree = fin2.Get(tree_name)

        tr = ROOTUtils.select_events_with_hist2d_selection(tree, lis_xyname[0], lis_xyname[1], h, sel_exp)

        fout_name = kwargs.get('Output.RootName')
        fout_name = os.path.join(_directory, fout_name)
        fout = TFile(fout_name, "recreate")

        tr.Write()

    def MergeRootMeasTime(self, kwargs):
        """if the tfile was merged by hadd, and there are more than one meas_time in it,
        use this function to merge the meas_time"""
        logger.info('Run MergeRootMeasTime.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.RootDir', ''))
        lis_in_fname = kwargs.get('Input.RootNames')
        lis_in_fname = [os.path.join(root_dir, i) for i in lis_in_fname]

        for fname in lis_in_fname:
            ROOTUtils.merge_meas_time_in_tfile(fname)

    def MergeEventsTreeTFile(self, kwargs):
        """merge tree of muon events and the meas_time"""
        logger.info('Run MergeEventsTreeTFile.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        root_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.RootDir', ''))
        lis_in_fname = kwargs.get('Input.RootNames')
        lis_in_fname = [os.path.join(root_dir, i) for i in lis_in_fname]

        fout_name = kwargs.get('Output.RootName', "")
        if fout_name == "":
            raise Exception("Please set an output root file name: Output.RootName = ?")
        fout_name = os.path.join(_directory, fout_name)

        ROOTUtils.merge_root_files(lis_in_fname, fout_name, if_delete=False)
        ROOTUtils.merge_meas_time_in_tfile(fout_name)

    def GetCountsTimeSeries(self, kwargs):
        """get counts time series with selection of leff"""
        logger.info('Run GetCountsTimeSeries.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        hts_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TimeSeriesRootDir', ''))
        lis_hts_fname = kwargs.get('Input.TimeSeriesRootNames')
        lis_hts_fname = [os.path.join(hts_dir, i) for i in lis_hts_fname]

        leff_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.LeffRootDir', ''))
        lis_leff_fname = kwargs.get('Input.LeffRootNames')
        lis_leff_fname = [os.path.join(leff_dir, i) for i in lis_leff_fname]

        if len(lis_hts_fname) != len(lis_leff_fname):
            raise Exception("The length of Input.TimeSeriesRootNames and Input.LeffRootNames"
                            "should be the same.")

        hts_reg = kwargs.get('Input.HistRegularExpression')
        sel_exp = kwargs.get('Input.LeffSelectionExpression')

        lis_fout_name = kwargs.get('Output.RootNames', "")
        if not isinstance(lis_fout_name, list):
            lis_fout_name = [os.path.split(i)[1].replace(".root", "_sel.root") for i in lis_hts_fname]
        lis_fout_name = [os.path.join(_directory, i) for i in lis_fout_name]
        if len(lis_hts_fname) != len(lis_fout_name):
            raise Exception("The length of Input.TimeSeriesRootNames and Output.RootNames"
                            "should be the same.")

        _ts_tool = TimeSeriesTools.TimeSeriesTools()

        for i in range(len(lis_hts_fname)):
            hts_fname = lis_hts_fname[i]
            leff_file = lis_leff_fname[i]
            fout_name = lis_fout_name[i]

            _ts_tool.GetCountsTimeSeries(hts_fname, hts_reg, leff_file, sel_exp, fout_name)

        lis_draw_option = kwargs.get('Input.DrawOptions', '')
        self._print_root_in_pdf(fout_name, lis_draw_option)

    def ScaleHCountsByTotalAirCounts(self, kwargs):
        """scale the hist of counts by total counts from air"""
        logger.info('Run ScaleHCountsByTotalAirCounts.')

        # create a directory if not exist
        # noinspection DuplicatedCode
        _directory = os.path.join(self._setting.get_global_dir(), kwargs.get('OutDir', ''))
        self._make_dir_if_not_exists(_directory)

        # get kwargs from config file
        hts_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.TimeSeriesRootDir', ''))
        lis_hts_fname = kwargs.get('Input.TimeSeriesRootNames')
        lis_hts_fname = [os.path.join(hts_dir, i) for i in lis_hts_fname]

        leff_dir = os.path.join(self._setting.get_global_dir(), kwargs.get('Input.LeffRootDir', ''))
        lis_leff_fname = kwargs.get('Input.LeffRootNames')
        lis_leff_fname = [os.path.join(leff_dir, i) for i in lis_leff_fname]

        if len(lis_hts_fname) != len(lis_leff_fname):
            raise Exception("The length of Input.TimeSeriesRootNames and Input.LeffRootNames"
                            "should be the same.")
