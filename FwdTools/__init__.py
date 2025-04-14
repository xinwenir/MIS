htheta_phi_hist_name = 'htheta_phi'
hthetax_thetay_hist_name = 'hthetax_thetay'

dic_hist_name = {'hthetaphi_ltopo': '%s_Ltopo' % htheta_phi_hist_name,
                 'hthetaxy_ltopo': '%s_Ltopo' % hthetax_thetay_hist_name,
                 'hthetaphi_leff': '%s_Leff' % htheta_phi_hist_name,
                 'hthetaxy_leff': '%s_Leff' % hthetax_thetay_hist_name,
                 'hthetaphi_ratio': '%s_ratio' % htheta_phi_hist_name,
                 'hthetaxy_ratio': '%s_ratio' % hthetax_thetay_hist_name,
                 'hthetaphi_flux': '%s_flux' % htheta_phi_hist_name,
                 'hthetaxy_flux': '%s_flux' % hthetax_thetay_hist_name,
                 'hthetaphi_flux_x_SldAgl': '%s_flux_x_SldAgl' % htheta_phi_hist_name,
                 'hthetaxy_flux_x_SldAgl': '%s_flux_x_SldAgl' % hthetax_thetay_hist_name,
                 'hthetaphi_counting_rate': '%s_counting_rate' % htheta_phi_hist_name,
                 'hthetaxy_counting_rate': '%s_counting_rate' % hthetax_thetay_hist_name,
                 'hthetaphi_counts': '%s_counts' % htheta_phi_hist_name,
                 'hthetaxy_counts': '%s_counts' % hthetax_thetay_hist_name,
                 'thetaphi_axis_title': 'phi [rad]; theta [rad]',
                 'thetaxy_axis_title': 'thetay [rad]; thetax [rad]'
                 }

points_columns = ['pt1_x', 'pt1_y', 'pt1_z', 'pt2_x', 'pt2_y', 'pt2_z', 'pt3_x', 'pt3_y', 'pt3_z']
centers_columns = ['ctr_x', 'ctr_y', 'ctr_z']
normal_columns = ['vn_x', 'vn_y', 'vn_z']
plane_equation_columns = ['pl_a0', 'pl_a1', 'pl_a2']    # z = a0 + a1x + a2y
intersection_columns = ['int_x', 'int_y', 'int_z']
path_length_columns = ['path_length']
in_out_columns = ['int_out']
