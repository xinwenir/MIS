import json
from ROOT import TMath, TString
from configobj import ConfigObj
from DataFactory import CommonUtils, ROOTUtils
import os
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

my_tz = "Asia/Shanghai"


class Setting:
    """read the global config file"""

    def __init__(self, cfg_fname):

        self.cfg = ConfigObj(cfg_fname)
        self.dic_global = self.cfg.get('global', {})

        tz = self.dic_global.get('TimeZone', '')
        if tz != '':
            global my_tz
            my_tz = tz

    def get_date(self):
        _date = self.dic_global.get('Date', 'today')
        if _date.__contains__('today'):
            _date = datetime.today().strftime("%b-%d-%Y")
        return _date

    def get_global_dir(self):
        return self.dic_global.get('GlobalDir', os.getcwd())

    def get_Tag(self):
        return self.dic_global.get('Tag', '')

    def get_description(self):
        return self.dic_global.get('Description', '')

    def get_debug_level(self):
        return self.dic_global.get('DebugLevel', 'DEBUG')

    def get_log_fname(self):
        log_fname = self.dic_global.getfloat('LogFile', '')
        if not CommonUtils.path_isabs(log_fname):
            log_fname = os.path.join(self.get_global_dir(), log_fname)
        return log_fname

    @classmethod
    def _get_bin_strings(cls, _str_bins):
        _min = _max = 0
        try:
            _str_pi = str(TMath.Pi())
            _min = eval(_str_bins[1].replace('Pi()', _str_pi))
            _max = eval(_str_bins[2].replace('Pi()', _str_pi))
        except Exception as e:
            print(_str_bins)

        return [int(_str_bins[0]), _min, _max]

    def get_bins(self):
        # get bins information or give default values
        str_xbins = self.dic_global.get('XBins', '200,-250,250'.split(','))
        x_bins = self._get_bin_strings(str_xbins)

        str_ybins = self.dic_global.get('YBins', '200,-250,250'.split(','))
        y_bins = self._get_bin_strings(str_ybins)

        str_zbins = self.dic_global.get('ZBins', '200,1,-1'.split(','))
        z_bins = self._get_bin_strings(str_zbins)

        str_theta_bins = self.dic_global.get('ThetaBins', '60,0,Pi()/2.'.split(','))
        theta_bins = self._get_bin_strings(str_theta_bins)

        str_phi_bins = self.dic_global.get('PhiBins', '200,-Pi(),Pi()'.split(','))
        phi_bins = self._get_bin_strings(str_phi_bins)

        str_thetax_bins = self.dic_global.get('ThetaXBins', '200,-Pi()/2.,Pi()/2.'.split(','))
        thetax_bins = self._get_bin_strings(str_thetax_bins)

        str_thetay_bins = self.dic_global.get('ThetaYBins', '200,-Pi()/2.,Pi()/2.'.split(','))
        thetay_bins = self._get_bin_strings(str_thetay_bins)

        bins_dic = {'x_bins': x_bins,
                    'y_bins': y_bins,
                    'z_bins': z_bins,
                    'theta_bins': theta_bins,
                    'phi_bins': phi_bins,
                    'thetax_bins': thetax_bins,
                    'thetay_bins': thetay_bins,
                    'x_title': "x [mm]",
                    'y_title': "y [mm]",
                    'z_title': "z [mm]",
                    'theta_title': "#theta [rad]",
                    'phi_title': "#phi [rad]",
                    'thetax_title': "#theta_{x} [rad]",
                    'thetay_title': "#theta_{y} [rad]"
                    }
        return bins_dic

    def get_det_info_fname(self):
        det_info_fname = self.dic_global.get('DetInfoFile', '')
        if not CommonUtils.path_isabs(det_info_fname):
            det_info_fname = os.path.join(self.get_global_dir(), det_info_fname)
        return det_info_fname

    def get_dets_info(self):
        det_info_fname = self.dic_global.get('DetInfoFile', '')
        if not CommonUtils.path_isabs(det_info_fname):
            det_info_fname = os.path.join(self.get_global_dir(), det_info_fname)
        return ConfigObj(det_info_fname)

    def get_topos_info(self):
        topo_info_fname = self.dic_global.get('TopoInfoFile', '')
        if not CommonUtils.path_isabs(topo_info_fname):
            topo_info_fname = os.path.join(self.get_global_dir(), topo_info_fname)
        return ConfigObj(topo_info_fname)

    def get_phys_model_info(self):
        phys_model_fname = self.dic_global.get('PhysModelFile', '')
        if not CommonUtils.path_isabs(phys_model_fname):
            phys_model_fname = os.path.join(self.get_global_dir(), phys_model_fname)
        return ConfigObj(phys_model_fname)

    def get_step_lists(self):
        lis_steps = self.dic_global.get('StepNames')
        if lis_steps is None:
            raise Exception('No steps are added')
        return lis_steps

    def get_step_info(self, step_name: str):
        dic_step = self.dic_global.get(step_name)
        if dic_step is None:
            raise Exception('Step %s dose not exist.' % step_name)
        return dic_step
