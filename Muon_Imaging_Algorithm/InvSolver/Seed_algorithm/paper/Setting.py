"""论文2仿真"""
start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6\refs"
refs_file=start_data_file
mesh_file=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\mesh.txt"
matrix_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\%s"
out_put_dir=r"G:\重要\论文和材料\第二篇\实验数据\BCD-NES\随机噪声10_05"
bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6\bounds"
moudles_file=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6\moudles"

"""城墙修正墙皮"""
# start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\optimize_res\new_methods_ref.den"
# refs_file=start_data_file
# mesh_file=r"data\seed_study\optimize_res\17_58MaMian.msh"
# matrix_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\%s"
# out_put_dir=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result"
# bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\optimize_res\new_methods_bounds.den"
"""城墙原始反演"""
# start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\wall\17_58MaMian_ref.den"
# refs_file=start_data_file 
# mesh_file=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\wall\17_58MaMian.msh"
# matrix_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\%s"
# out_put_dir=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\wall\res\seed_res"
# bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\wall\17_58MaMian_bnd.den"

# # zzg
# start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5_ref_0.den"
# # # start_data_file=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result\wall"

# mesh_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5.msh"
# matrix_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\%s"
# out_put_dir=r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result"
# bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5_bnd_seed_extra1.den"
# moudles_file=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\20_muto3d.den"

max_count=7

#268408城墙小房子 185437 CRB5_中心值 179895 CRB5_非中心 #第二篇仿真 88896
start_seed=88896

# paper2
punishment_diff_value_multiple=1 
punishment_Search_distancevalue_multiple=0.00*1.3 #限制搜索距离
prompt_value_multiple=0.05 #促使当前的体素密度接近原始seed 不宜太大除非异常体密度比较均匀
misfit_all_multiple=1
smooth_all_multiple=0.5 

punishment_factor=1.0 #外点罚函数的因子

threshold=1.0
seed_Minimum_density_variation=0.4 #种子最小密度变化,小于该值将不能成为种子,该值必须小于threshold

# wall--
# punishment_diff_value_multiple=1 
# punishment_Search_distancevalue_multiple=0.006*1.7 #限制搜索距离
# prompt_value_multiple=0.1 #促使当前的体素密度接近原始seed 不宜太大除非异常体密度比较均匀
# misfit_all_multiple=1
# smooth_all_multiple=0.5 

# punishment_factor=1.0 #外点罚函数的因子

# threshold=0.7
# seed_Minimum_density_variation=0.6 #种子最小密度变化,小于该值将不能成为种子,该值必须小于threshold

print(out_put_dir)
# print("高斯噪声1  0.5")