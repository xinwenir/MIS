    # """城墙"""
    # # start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\optimize_res\new_methods_ref.den"
    # # mesh_file=r"data\seed_study\optimize_res\17_58MaMian.msh"
    # # matrix_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\%s"
    # # out_put_dir=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result"
    # # bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\optimize_res\new_methods_bounds.den"
    # """仿真"""
start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5_ref_0.den"
# # start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result\wall"

mesh_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5.msh"
matrix_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\%s"
out_put_dir=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result"
bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5_bnd_seed_extra1.den"
moudles_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5_ref_seed_extra1.den"

max_count=60
model = 1
# mesh_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG3.msh"
# matrix_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\%s"
# out_put_dir=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\output\4_obs\seed_res"
# bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG3_bnd_1.den"
# moudles_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG3_ref_ore.den"
# start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG3_ref_1.den"



#268408城墙小房子 185437 CRB5_中心值 179895 CRB5_非中心
start_seed=185437


# zzg实测
punishment_diff_value_multiple=1
punishment_Search_distancevalue_multiple=0.06*1.5
prompt_value_multiple=0.001 #促使当前的体素密度接近原始seed 不宜太大除非异常体密度比较均匀
misfit_all_multiple=1
smooth_all_multiple=0.01

threshold=0.3
seed_Minimum_density_variation=0.1#种子最小密度变化,小于该值将不能成为种子,该值必须小于threshold



# # CRB5仿真
# punishment_diff_value_multiple=1
# punishment_Search_distancevalue_multiple=0.06
# prompt_value_multiple=0.08 #促使当前的体素密度接近原始seed 不宜太大除非异常体密度比较均匀
# misfit_all_multiple=1
# smooth_all_multiple=0.01

# threshold=0.4
# seed_Minimum_density_variation=0.2#种子最小密度变化,小于该值将不能成为种子,该值必须小于threshold
