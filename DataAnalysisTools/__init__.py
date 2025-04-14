_hthetaxy_name_prefix = 'hthetax_thetay'
_hthetaphi_name_prefix = 'htheta_phi'

dic_tree_branch_name = {'det_id': 'DetID',
                        'det_x': 'DetX',
                        'det_y': 'DetY',
                        'det_z': 'DetZ',
                        'theta': 'Theta',
                        'phi': 'Phi',
                        'thetax': 'Thetax',
                        'thetay': 'Thetay',
                        'ltopo': 'Ltopo',
                        'meas_ratio': 'MeasRatio',
                        'meas_ratio_err': 'MeasRatioErr',
                        'sim_ratio': 'SimRatio',
                        'sim_ratio_err': 'SimRatioErr',
                        'meas_leff': 'MeasLeff',
                        'meas_leff_err': 'MeasLeffErr',
                        'sim_leff': 'SimLeff',
                        'sim_leff_err': 'SimLeffErr',
                        'meas_target_count': 'MeasTargetCount',
                        'meas_background_count': 'MeasBackgroundCount',
                        'sim_target_count': 'SimTargetCount',
                        'sim_background_count': 'SimBackgroundCount',
                        'rnd_meas_ratio': 'RndMeasRatio',
                        'rnd_meas_ratio_err': 'RndMeasRatioErr',
                        'rnd_sim_ratio': 'RndSimRatio',
                        'rnd_sim_ratio_err': 'RndSimRatioErr',
                        'rnd_meas_leff': 'RndMeasLeff',
                        'rnd_meas_leff_err': 'RndMeasLeffErr',
                        'rnd_sim_leff': 'RndSimLeff',
                        'rnd_sim_leff_err': 'RndSimLeffErr',
                        'sig': 'Sig',
                        'smooth_sig': 'SmoothSig'}

# """
# 包含一些通用的参数
# """
# # -----------数据索引-----------
# _Index = []
# # 0-35
# for i in range(36):
#     _Index.append('chn_{:2>d}'.format(i))
# _Index.append('temperature')  # 36
# _Index.append('triggerID')  # 37
# _Index.append('boardID')  # 38
# _Index.append("timeTag")  # 39
#
# _barIndex = []
# for i in range(36):
#     _barIndex.append('bar_{:2>d}'.format(i))
#
# # -----------指示数据大小的枚举-------
# from enum import Enum
#
#
# class sizeUnit_binary(Enum):
#     B = 1
#     KB = 1024
#     MB = 1024 * 1024
#
#
# # ---------共享数据地址---------
# _address = ("127.0.0.1", 50001)
# _authkey = b'jimbook'
#
#
# dic_hists_names = {'ratio': ['ratio'],
#                    'Ltopo': ['Ltopo'],
#                    'Leff': ['Leff'],
#                    'counts': ['counts', 'TOT', 'MM', 'SS', 'MS'],}
