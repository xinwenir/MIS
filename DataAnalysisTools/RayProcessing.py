# import array
# import os
# import sys
#
# from ROOT import TH2F, TTree, TFile, TVector3, TMath, gROOT
#
# from DataFactory import Utils
# from DataFactory.Utils import find_key_in_file
# from DataAnalysisTools import dic_hists_names
#
# import logging
#
# logger = logging.getLogger(__name__)
#
#
# class RayProcessing:
#     def __init__(self):
#         self._tree = None
#         self._tree_pro = None
#
#     def write_tree_into_file(self, fout_name, opt=''):
#         logger.info('write tree into file: %s' % fout_name)
#
#         fout = TFile(fout_name, 'recreate')
#
#         if opt.lower().__contains__('pro'):
#             self._tree_pro.Write()
#         else:
#             self._tree.Write()
#
#         fout.Save()
#         fout.Close()
#
#     def write_obs_file(self, fin_name, fout_name, leff_name='meas_Leff'):
#         logger.info('write %s into %s' % (fin_name, fout_name))
#
#         fin = TFile(fin_name, 'read')
#         t_tree = fin.Get('t_tree')
#
#         fout = open(fout_name, "w")
#
#         t_arr = array.array('d', [-10.])
#
#         i_print = 0
#         for entry in t_tree:
#             det_id = entry.det_id
#             pos_x = entry.pos_x
#             pos_y = entry.pos_y
#             pos_z = entry.pos_z
#             theta = entry.theta
#             phi = entry.phi
#             sim_Ltopo = entry.sim_Ltopo
#
#             str_cmd = "t_arr[0] = entry.%s" % leff_name
#             exec(str_cmd)
#             Leff = t_arr[0]
#             str_cmd = "t_arr[0] = entry.%s_err" % leff_name
#             exec(str_cmd)
#             Leff_err = t_arr[0]
#
#             if sim_Ltopo > 0. or sim_Ltopo < -9.:
#                 str_print = '%i %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f' % (det_id, pos_x, pos_y, pos_z, theta, phi, sim_Ltopo, Leff, Leff_err)
#                 print(str_print, file=fout)
#                 i_print = i_print + 1
#
#         fout.close()
#
#         str_cmd = 'sed -i "1i %d" %s' % (i_print, fout_name)
#         os.system(str_cmd)
#
#         return 0
#
#     def cvt_hist_to_tree(self, v_dic_ray_hist2tree: list):
#         # dic_ray_hist2tree = {'det_name': _h2tr_det_name,
#         #                      'meas_fnames': _h2tr_meas_fnames,
#         #                      'sim_fnames': _h2tr_sim_fnames,
#         #                      'meas_hist_prefixes': _h2tr_meas_hist_prefixes,
#         #                      'sim_hist_prefixes': _h2tr_sim_hist_prefixes}
#         # _tree_vars = ['DetID', 'PosX', 'PosY', 'PosZ', 'MeasRatio', 'MeasLeff', 'SimRatio', 'SimLeff', 'Ltopo']
#
#         self._tree = TTree('t_tree', 't_tree')
#
#         det_id = array.array('i', [0])
#         pos_x = array.array('d', [0.])
#         pos_y = array.array('d', [0.])
#         pos_z = array.array('d', [0.])
#         theta = array.array('d', [0.])
#         phi = array.array('d', [0.])
#         meas_ratio = array.array('d', [-10.])
#         sim_ratio = array.array('d', [-10.])
#         meas_Leff = array.array('d', [-10.])
#         sim_Leff = array.array('d', [-10.])
#         sim_Ltopo = array.array('d', [-10.])
#         meas_ratio_err = array.array('d', [-10.])
#         sim_ratio_err = array.array('d', [-10.])
#         meas_Leff_err = array.array('d', [-10.])
#         sim_Leff_err = array.array('d', [-10.])
#         meas_counts = array.array('d', [-10.])
#         sim_counts = array.array('d', [-10.])
#         meas_counts_err = array.array('d', [-10.])
#         sim_counts_err = array.array('d', [-10.])
#
#         self._tree.Branch('det_id', det_id, 'det_id/I')
#         self._tree.Branch('pos_x', pos_x, 'pos_x/D')
#         self._tree.Branch('pos_y', pos_y, 'pos_y/D')
#         self._tree.Branch('pos_z', pos_z, 'pos_z/D')
#         self._tree.Branch('theta', theta, 'theta/D')
#         self._tree.Branch('phi', phi, 'phi/D')
#         self._tree.Branch('meas_ratio', meas_ratio, 'meas_ratio/D')
#         self._tree.Branch('sim_ratio', sim_ratio, 'sim_ratio/D')
#         self._tree.Branch('meas_Leff', meas_Leff, 'meas_Leff/D')
#         self._tree.Branch('sim_Leff', sim_Leff, 'sim_Leff/D')
#         self._tree.Branch('sim_Ltopo', sim_Ltopo, 'sim_Ltopo/D')
#         self._tree.Branch('meas_ratio_err', meas_ratio_err, 'meas_ratio_err/D')
#         self._tree.Branch('sim_ratio_err', sim_ratio_err, 'sim_ratio_err/D')
#         self._tree.Branch('meas_Leff_err', meas_Leff_err, 'meas_Leff_err/D')
#         self._tree.Branch('sim_Leff_err', sim_Leff_err, 'sim_Leff_err/D')
#         self._tree.Branch('meas_counts', meas_counts, 'meas_counts/D')
#         self._tree.Branch('sim_counts', sim_counts, 'sim_counts/D')
#         self._tree.Branch('meas_counts_err', meas_counts_err, 'meas_counts_err/D')
#         self._tree.Branch('sim_counts_err', sim_counts_err, 'sim_counts_err/D')
#
#         for dic_ray_hist2tree in v_dic_ray_hist2tree:
#             logger.info('Convert hists to tree for det %s', dic_ray_hist2tree['det_name'])
#
#             # initial
#             meas_ratio[0] = sim_ratio[0] = meas_Leff[0] = sim_Leff[0] = sim_Ltopo[0] = -10.
#
#             dic_hists = self._fill_dic_hist(dic_ray_hist2tree)
#             if len(dic_hists) < 1:
#                 raise Exception('Did not get any hists for det %s' % dic_ray_hist2tree['det_name'])
#
#             det_id[0] = dic_ray_hist2tree['det'].get_det_id()
#             pos_x[0] = dic_ray_hist2tree['det'].get_positions()[0]
#             pos_y[0] = dic_ray_hist2tree['det'].get_positions()[1]
#             pos_z[0] = dic_ray_hist2tree['det'].get_positions()[2]
#             rot_phi = -1. * dic_ray_hist2tree['det'].get_rotAngle()[2]
#
#             t_hist = list(dic_hists.values())[0]
#
#             if t_hist.GetName().lower().__contains__('phi'):
#                 mode = 0
#             elif t_hist.GetName().lower().__contains__('thetax'):
#                 mode = 1
#             else:
#                 raise ValueError('Given hists should be in thetaphi or thetaxy')
#
#             lis_variables = list(dic_hists.keys())
#             lis_hists = list(dic_hists.values())
#
#             if len(lis_hists) > 1:
#                 for i in range(1, len(lis_hists)):
#                     if_same_size = Utils.check_two_hists_the_same_size(lis_hists[0], lis_hists[i])
#                     if not if_same_size:
#                         raise Exception('Give hists with same size here.')
#
#             xbins = t_hist.GetNbinsX()
#             ybins = t_hist.GetNbinsY()
#
#             for i_xbin in range(1, xbins + 1):
#                 xx = t_hist.GetXaxis().GetBinCenter(i_xbin)
#                 for i_ybin in range(1, ybins + 1):
#                     yy = t_hist.GetYaxis().GetBinCenter(i_ybin)
#
#                     ibin = t_hist.GetBin(i_xbin, i_ybin)
#
#                     if mode == 0:
#                         theta[0] = yy
#                         phi[0] = xx + rot_phi
#                     else:
#                         vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
#                         theta[0] = vec.Theta()
#                         phi[0] = vec.Phi() + rot_phi
#
#                     for i in range(len(lis_hists)):
#                         str_cmd = "%s[0] = lis_hists[i].GetBinContent(ibin)" % lis_variables[i]
#                         exec(str_cmd)
#                         if not lis_variables[i].__contains__('topo'):
#                             str_cmd = "%s_err[0] = lis_hists[i].GetBinError(ibin)" % lis_variables[i]
#                             exec(str_cmd)
#
#                     self._tree.Fill()
#
#         return 0
#
#     def _fill_dic_hist(self, dic_ray_hist2tree):
#         _det_name = dic_ray_hist2tree['det_name']
#         _meas_fnames = dic_ray_hist2tree['meas_fnames']
#         _sim_fnames = dic_ray_hist2tree['sim_fnames']
#         _meas_hist_names = dic_ray_hist2tree['meas_hist_names']
#         _sim_hist_names = dic_ray_hist2tree['sim_hist_names']
#
#         dic_hists = {}
#
#         if len(_meas_fnames) > 0 and os.path.isfile(_meas_fnames[0]):
#             dic_hists_meas = self._get_hists_from_file(_meas_fnames, _meas_hist_names, 'meas')
#             dic_hists.update(dic_hists_meas)
#
#         if len(_sim_fnames) > 0 and os.path.isfile(_sim_fnames[0]):
#             dic_hists_sim = self._get_hists_from_file(_sim_fnames, _sim_hist_names, 'sim')
#             dic_hists.update(dic_hists_sim)
#
#         return dic_hists
#
#     @classmethod
#     def _get_hist_type(cls, hist_name: str):
#         key_list = list(dic_hists_names.keys())
#         val_list = list(dic_hists_names.values())
#         for i in range(len(key_list)):
#             i_value = val_list[i]
#             for i_type_name in i_value:
#                 if hist_name.__contains__(i_type_name):
#                     return key_list[i]
#         logger.warning('Cannot find the type of the hist: %s' % hist_name)
#         return ''
#
#     def _get_hists_from_file(self, fin_names: list, hist_names: list, opt: str):
#         if opt.lower().__contains__('meas'):
#             opt = 'meas'
#         elif opt.lower().__contains__('sim'):
#             opt = 'sim'
#         else:
#             raise ValueError('opt must be meas or sim. opt = %s' % opt)
#
#         if len(fin_names) != len(hist_names):
#             raise Exception('[_get_hists_from_file] len(fin_names) == len(hist_names):')
#
#         dic_hists = {}
#         for i in range(len(fin_names)):
#             fin_name = fin_names[i]
#             hist_name = find_key_in_file(fin_name, hist_names[i])
#
#             if hist_name == '':
#                 continue
#
#             hist_type = self._get_hist_type(hist_name)
#             if hist_type.__contains__('count'):
#                 hist_name = find_key_in_file(fin_name, hist_names[i], 'whole')
#
#             if hist_name.count('_') > 2 \
#                     and not (hist_name.__contains__('Leff')
#                              or hist_name.__contains__('Ltopo')
#                              or hist_name.__contains__('ratio')):
#                 continue
#
#             hist_key = '%s_%s' % (opt, hist_type)
#             if dic_hists.get(hist_key) is None:
#                 fin = TFile(fin_name, 'read')
#                 t_hist = fin.Get(hist_name)
#                 t_hist.SetDirectory(gROOT)
#                 fin.Close()
#                 dic_hists[hist_key] = t_hist
#                 logger.debug('Get hist_key: %s from file %s' % (t_hist.GetName(), fin_name))
#
#         return dic_hists
#
# # def write_hist_into_tree(h_sim_ratio, h_sim_leff, det, h_meas_ratio, h_meas_leff, h_meas_count, h_sim_count):
# #     t_tree = TTree('t_tree', 't_tree')
# #     _tree_vars = ['DetID', 'PosX', 'PosY', 'PosZ', 'MeasRatio', 'MeasLeff', 'SimRatio', 'SimLeff', 'Ltopo']
# #     for i_var in _tree_vars:
# #         define_var = "%s = array('d', [0])" % i_var
# #         tree_branch = "t_tree.Branch('%s', %s, '%s/D')" % (i_var, i_var, i_var)
# #         if i_var.__contains__('ID'):
# #             define_var = define_var.replace("'d'", "'i'")
# #             tree_branch = tree_branch.replace('/D', '/I')
# #         exec(define_var)
# #         exec(tree_branch)
# #
# #     exec('PosX = det.get_positions()[0]')
# #     exec('PosY = det.get_positions()[1]')
# #     exec('PosZ = det.get_positions()[2]')
# #     exec('DetID = det.get_det_id()')
# #     rot_phi = det.get_rotAngle()[-1]
# #
# #     if h_meas_ratio is not None and not Utils.check_two_hists_the_same_size(h_meas_ratio, h_meas_leff):
# #         raise Exception('[write_hist_into_tree] h_meas_ratio, h_meas_leff should be with the same size.')
# #     if h_sim_ratio is not None and not Utils.check_two_hists_the_same_size(h_sim_ratio, h_sim_leff):
# #         raise Exception('[write_hist_into_tree] h_sim_ratio, h_sim_leff should be with the same size.')
# #
# #     if h_meas_ratio.GetName().lower().__contains__('phi'):
# #         mode = 0
# #     elif h_meas_ratio.GetName().lower().__contains__('thetax'):
# #         mode = 1
# #     else:
# #         raise ValueError('hist should be xxx_thetaphi or xxx_thetaxy')
# #
# #     xbins = h_meas_ratio.GetNbinsX()
# #     ybins = h_meas_ratio.GetNbinsY()
# #
# #     for i_xbin in range(1, xbins + 1):
# #         xx = h_meas_ratio.GetXaxis().GetBinCenter(i_xbin)
# #         for i_ybin in range(1, ybins + 1):
# #             yy = h_meas_ratio.GetYaxis().GetBinCenter(i_ybin)
# #
# #             ibin = h_meas_ratio.GetBin(i_xbin, i_ybin)
# #             meas_ratio = h_meas_ratio.GetBinContent(ibin)
# #             sim_ratio = h_sim_ratio.GetBinContent(ibin)
# #             meas_leff = h_meas_leff.GetBinContent(ibin)
# #             sim_leff = h_sim_leff.GetBinContent(ibin)
# #
# #             if mode == 0:
# #                 theta = yy
# #                 phi = xx + rot_phi
# #             else:
# #                 vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
# #                 theta = vec.Theta()
# #                 phi = vec.Phi() + rot_phi
# #
# #
# # def find_hist_from_file(fname, hist_prefixes):
# #     h_ratio = h_leff = None
# #     if fname.__contains__('.root') and not hist_prefixes == '':
# #         t_prefix_ratio = hist_prefixes + '_ratio'
# #         t_prefix_leff = hist_prefixes + '_Leff'
# #
# #         h_ratio_name = find_key_in_file(fname, t_prefix_ratio)
# #         h_leff_name = find_key_in_file(fname, t_prefix_leff)
# #
# #         if h_ratio_name == '' or h_leff_name == '':
# #             logger.warning('Cannot find hist %s_xxx in file %s' % (hist_prefixes, fname))
# #         else:
# #             t_fin = TFile(fname, 'read')
# #             h_ratio = t_fin.Get(h_ratio_name)
# #             h_leff = t_fin.Get(h_leff_name)
# #     return h_ratio, h_leff
