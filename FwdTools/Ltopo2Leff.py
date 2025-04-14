import os.path

from ROOT import TH2F, TFile, gROOT
import re

from FwdTools import dic_hist_name
from DataFactory import ROOTUtils
import logging

logger = logging.getLogger(__name__)


class Ltopo2Leff:
    def __init__(self, dic_ltopo_fname: dict, input_dir):
        self._dic_ltopo_fname = dic_ltopo_fname
        self._input_dir = input_dir

    @classmethod
    def _split_expression(cls, expr: str):
        # lis_topo = re.findall(r'[A-Za-z_]+[\d]*', expr)
        lis_topo = re.findall(r'\*[\w]+', expr)
        lis_topo = [i[1:] for i in lis_topo]

        lis_scale = []
        for i in range(len(lis_topo)):
            if i == 0:
                begin_index = 0
            else:
                begin_index = expr.find(lis_topo[i - 1]) + len(lis_topo[i - 1])
            _scale = expr[begin_index:expr.find(lis_topo[i])].removesuffix('*')
            lis_scale.append(float(_scale))
        return lis_topo, lis_scale

    def write_Leff_from_expression(self, expr: str, leff_fname: str, det_name: str):
        logger.info('Writing Leff from expression "%s" for %s' % (expr, det_name))
        lis_topo, lis_scale = self._split_expression(expr)

        htheta_phi_leff = hthetax_thetay_leff = None
        for i in range(len(lis_topo)):
            _topo_info = self._dic_ltopo_fname.get(lis_topo[i])
            if _topo_info is None:
                raise Exception('Cannot find %s in topoInfo' % lis_topo[i])
            _ltopo_fname = _topo_info.get('ObjFile'). \
                replace('.obj', '_%s.root' % det_name).replace('.csv', '_%s.root' % det_name)
            _ltopo_fname = os.path.join(self._input_dir, _ltopo_fname)
            hthetaphi_name = ROOTUtils.find_key_in_file(_ltopo_fname, dic_hist_name['hthetaphi_ltopo'])
            hthetaxy_name = ROOTUtils.find_key_in_file(_ltopo_fname, dic_hist_name['hthetaxy_ltopo'])

            if hthetaphi_name == '' or hthetaphi_name == '':
                raise Exception('Cannot find %s or %s in file %s.' %
                                (dic_hist_name['hthetaphi_ltopo'], dic_hist_name['hthetaxy_ltopo'], _ltopo_fname))

            _f = TFile(_ltopo_fname, 'read')
            _hthetaphi = _f.Get(hthetaphi_name)
            _hthetaxy = _f.Get(hthetaxy_name)

            if _hthetaphi is None or _hthetaxy is None:
                raise Exception('Did not get %s and %s from %s' % (hthetaphi_name, hthetaxy_name, _ltopo_fname))

            if i == 0:
                htheta_phi_leff = _hthetaphi.Clone()
                htheta_phi_leff.Scale(lis_scale[i])
                hthetax_thetay_leff = _hthetaxy.Clone()
                hthetax_thetay_leff.Scale(lis_scale[i])

                htheta_phi_leff.SetDirectory(gROOT)
                hthetax_thetay_leff.SetDirectory(gROOT)
            else:
                htheta_phi_leff.Add(_hthetaphi, lis_scale[i])
                hthetax_thetay_leff.Add(_hthetaxy, lis_scale[i])

        htheta_phi_leff.SetNameTitle(dic_hist_name['hthetaphi_leff'],
                                     '%s; phi [rad]; theta [rad]; Leff [m g/cm^{3}]' % dic_hist_name['hthetaphi_leff'])
        hthetax_thetay_leff.SetNameTitle(dic_hist_name['hthetaxy_leff'],
                                         '%s; phi [rad]; theta [rad]; Leff [m g/cm^{3}]' % dic_hist_name['hthetaxy_leff'])

        # htheta_phi_leff.SetMinimum(1e-15)
        # hthetax_thetay_leff.SetMinimum(1e-15)

        fout = TFile(leff_fname, 'recreate')
        htheta_phi_leff.Write()
        hthetax_thetay_leff.Write()
        fout.Save()
        fout.Close()
