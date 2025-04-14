# author:高金磊
# datetime:2022/4/3 10:36

"""
不允许继承现有的任何本地类
"""
from ast import If
from hashlib import new
import json
from math import fabs
import os
import threading
import time
from datetime import datetime
from tkinter import N

import yaml
import configparser

from InvSysTools.MyTools import myPrint

class Setting(object):
    """
    严格单例模式
    """
    _instance = None
    @staticmethod
    def get_instance(need_new=False, new_setting_file_path=None):
        if Setting._instance is None or need_new:
            Setting._instance = Setting(new_setting_file_path)
        return Setting._instance

    """全局懒加载模式"""
    def __init__(self,setting_file_path) -> None:  # type: ignore
        self.configuration_path=setting_file_path
        if setting_file_path is None or not os.path.exists(setting_file_path):
            raise Exception("配置文件路径不正确或者不存在:%s"%(setting_file_path))
        f = open(setting_file_path, 'r', encoding='utf-8')
        config = yaml.unsafe_load(f)
        self.config=config
        Inversion_procedure_config = config["Inversion_procedure"]
        """env"""
        self._log_model_dic = Inversion_procedure_config["_log_model_dic"]
        self._solver_meth = Inversion_procedure_config["_solver_meth"]
        self.STP_dic = Inversion_procedure_config["STP_dic"]
        """参数"""
        _parameter = Inversion_procedure_config["_parameter"]
        self._beta = _parameter["_beta"]



        self._alpha_s = _parameter["_alpha_s"]
        self._alpha_x = _parameter["_alpha_x"]
        self._alpha_y = _parameter["_alpha_y"]
        self._alpha_z = _parameter["_alpha_z"]
        """输入文件"""
        _all_input = Inversion_procedure_config["_All_input"]
        self._obs_file = _all_input["_obs_file"]
        self._bnd_file = _all_input["_bnd_file"]
        self._refs_file = _all_input["_refs_file"]
        self._mesh_file = _all_input["_mesh_file"]
        self._topo_file = _all_input["_topo_file"]
        self._sparse_matrix_file = _all_input["_sparse_matrix_file"]

        """输出文件"""
        _Output_file = Inversion_procedure_config["_Output_file"]
        self._res_file = _Output_file["_res_file"]
        self._res_refs_file = _Output_file["_res_refs_file"]
        self._res_smooth_file = _Output_file["_res_smooth_file"]
        """日志文件"""
        _log_config = Inversion_procedure_config["_log_config"]
        self._golable_log_file = _log_config["_golable_log_file"]
        self._local_log_path = _log_config["_local_log_path"]
        self._log_model = _log_config["_log_model"]  # 1 or 2
        self._log_detailed=_log_config["_log_detailed"]
        """中间文件"""
        _Temp_file = Inversion_procedure_config["_Temp_file"]
        self._all_new = _Temp_file["_all_new"]  # 刷新所有文件,建议仅在当更换obs或者密度矩阵后为True
        self._jxyz_file = _Temp_file["_jxyz_file"]  # 此文件被淘汰,如果需要生成请单独生成或者指定需要生成
        self._need_jxyz_file = _Temp_file["_need_jxyz_file"]
        self._Assist_j_file = _Temp_file["_Assist_j_file"]
        self._G_file = _Temp_file["_G_file"]
        self._d_file = _Temp_file["_d_file"]
        self._air_j_file = _Temp_file["_air_j_file"]
        self._pred_obsd_derr_file = _Temp_file["_pred_obsd_derr_file"]
        self._all_misfit_ms_smooth_file = _Temp_file["_all_misfit_ms_smooth_file"]
        self._beta_misfit_norms_file = _Temp_file["_beta_misfit_norms_file"]
        self._Ray_way_j_file = _Temp_file["_Ray_way_j_file"]
        """求解器"""
        _Solve = Inversion_procedure_config["_Solve"]
        self._solve_model = _Solve["_solve_model"]
        self._max_iter = _Solve["_max_iter"]  # 这个仅仅指目标函数调用次数,具体内部迭代次数根据方法特性写死的
        self._max_time = _Solve["_max_time"]  # 允许求解器运行的最大时间
        self._rm_unable_ray = _Solve["_rm_unable_ray"]
        self._rm_air_cells = _Solve["_rm_air_cells"]  # 排除属于空气的格子,不参与运算,默认不改变方程组的维度
        self._rm_unray_cells = _Solve["_rm_unray_cells"]  # 排除没有射线穿过的格子,不参与的运算,默认不改变方程组的维度
        self._rm_unable_cells = _Solve["_rm_unable_cells"]  # _rm_air_cells和_rm_unray_cells排除掉的格子是否完全去掉(改变方程组的维度)---此项可以大幅度加快运算
        self._show_detailed = _Solve["_show_detailed"]
        """其他配置"""
        _Other_setting = Inversion_procedure_config["_Other_setting"]
        self._show_res_picture = _Other_setting["_show_res_picture"]
        self._store_res_picture = _Other_setting["_store_res_picture"]  # 没必要让用户指定文件名称
        self._res_recover_refs = _Other_setting["_res_recover_refs"]
        self._res_recover_air = _Other_setting["_res_recover_air"]
        self._air_value = _Other_setting["_air_value"]
        """检查输入文件是否存在"""
        for file in _all_input.values():
            if file is None:
                #部分用不到的文件跳过检查
                continue
            if not os.path.exists(file):
                myPrint.myPrint_Wran(file, "---Not Found!!!!")
                # 因为有很多小工具不需要全部数据,所以这里暂时不加异常?????
        if self._log_model not in self._log_model_dic.keys():
            myPrint.myPrint_Err("日志模式选择错误,请在以下内容中选择:", self._log_model_dic)
            raise Exception("日志模式不合法")
        if self._solve_model not in self._solver_meth.keys():
            myPrint.myPrint_Err("求解器选择错误,请在以下内容中选择:", self._solver_meth)
            raise Exception("求解器不合法")
        #????????????????????????????????
        # if not os.path.exists(self._res_file):
        #     os.makedirs(self._res_file, exist_ok=True)
        # if not os.path.exists(self._G_file):
        #     os.makedirs(self._G_file, exist_ok=True)
        self._log = None
        
    def set_beta(self, beta):  # 临时开放
        self._beta = beta
    def get_loger(self) -> myPrint.Loger_Print:
        """
        获取全局的Loger_Print对象
        对loger和print进行封装
        可以像使用loger一样使用它
        :return:
        :rtype: Tools.MyTools.myPrint.Loger_Print
        """
        if self._log is None:
            # 根据模式加载
            if self._log_model == 1:
                loger = open(self._golable_log_file, 'a')
                printer = myPrint.myPrint_Success
                self._log = myPrint.Loger_Print(loger, printer)
            if self._log_model == 2:
                import datetime
                # 获得当前时间
                now = datetime.datetime.now()
                # 转换为指定的格式:
                otherStyleTime = now.strftime("%Y_%m_%d_%H_%M_%S")
                log = open(self._local_log_path + otherStyleTime + ".txt", 'w')
                printer = myPrint.myPrint_Success
                self._log = myPrint.Loger_Print(log, printer)
        return self._log

    def get_beta(self):
        return self._beta

    def get_alpha(self):
        return self._alpha_s, self._alpha_x, self._alpha_y, self._alpha_z

    def get_rm_unable_ray(self):
        return self._rm_unable_ray

    def get_rm_unable_cells(self):
        return self._rm_unable_cells

    def get_rm_unray_cells(self):
        return self._rm_unray_cells

    def get_rm_air_cells(self):
        return self._rm_air_cells

    def get_max_time(self):
        return self._max_time

    def get_log_detailed(self,level):
        """
        用来判断是否需要记录日志,等级低于设置的等级将会被打印
        :param level:要记录的日志信息的等级
        :return:
        """
        return level<=self._log_detailed

    def record_env_detailed(self):
        """
        记录当前的详细配置
        :return:
        """
        log = self.get_loger()
        # print_s = MyTools.MyTools.myPrint.myPrint_Success()
        log.write("使用配置文件:%s" % self.configuration_path)
        log.write(json.dumps(self.config, sort_keys=False, indent=2, ensure_ascii=False))
        log.write("日志模式:%s" % (self._log_model_dic[self._log_model]), printer=myPrint.myPrint_Hint)
        log.write("系统环境记录.....%s" % (datetime.now()), printer=myPrint.myPrint_Hint)
        log.write("配置数据%s" % (""), printer=myPrint.myPrint_Hint)
        log.write("%s:%s" % ('_beta', self._beta))
        log.write("%s:%s" % (
            '_alpha_s, _alpha_x, _alpha_y, _alpha_z', (self._alpha_s, self._alpha_x, self._alpha_y, self._alpha_z)))
        log.write("%s:%s" % ('求解器类型:', self._solver_meth[self._solve_model]))
        log.write("%s:%s" % ('_obs_file', self._obs_file))
        log.write("%s:%s" % ('_bnd_file', self._bnd_file))
        log.write("%s:%s" % ('_refs_file', self._refs_file))
        log.write("%s:%s" % ('_mesh_file', self._mesh_file))
        log.write("%s:%s" % ('_topo_file', self._topo_file))
        log.write("%s:%s" % ('_sparse_matrix_file', self._sparse_matrix_file))
        log.write("输出文件记录.....%s" % (""), printer=myPrint.myPrint_Hint)
        log.write("%s:%s" % ('_res_file', self._res_file))
        log.write("%s:%s" % ('_res_refs_file', self._res_refs_file))
        log.write("中间文件记录.....%s" % (""), printer=myPrint.myPrint_Hint)
        log.write("%s:%s" % ('_jxyz', self._jxyz_file))
        log.write("%s:%s" % ('_G', self._G_file))
        log.write("%s:%s" % ('_d', self._d_file))
        log.write("%s:%s" % ('_air_j', self._air_j_file))
        
    def get_solve_method(self):
        return self._solver_meth[self._solve_model]

    def get_air_value(self):
        return self._air_value

    def get_obs_file(self):
        return self._obs_file

    def get_bnd_file(self):
        return self._bnd_file

    def get_refs_file(self):
        return self._refs_file

    def get_mesh_file(self):
        return self._mesh_file

    def get_topo_file(self):
        return self._topo_file

    def get_sparse_matrix_file(self):
        return self._sparse_matrix_file

    def get_res_file(self):
        return self._res_file

    def get_res_refs_file(self):
        return self._res_refs_file

    def get_jxyz_file(self):
        return self._jxyz_file

    def get_G_file(self):
        return self._G_file

    def get_d_file(self):
        return self._d_file

    def get_air_j_file(self):
        return self._air_j_file

    def get_max_iter(self):
        return self._max_iter

    def get_show_detailed(self):
        return self._show_detailed

    def get_res_recover_refs(self):
        return self._res_recover_refs

    def get_res_recover_air(self):
        return self._res_recover_air

    def get_all_new(self):
        return self._all_new

    def get_need_jxyz_file(self):
        return self._need_jxyz_file

    def get_res_smooth_file(self):
        return self._res_smooth_file

    def get_pred_obsd_derr_file(self):
        return self._pred_obsd_derr_file

    def get_all_misfit_ms_smooth_file(self):
        return self._all_misfit_ms_smooth_file

    def get_Ray_way_j_file(self):
        return self._Ray_way_j_file

    def get_beta_misfit_norms_file(self):
        return self._beta_misfit_norms_file
    @property
    def show_res_picture(self):
        return self._show_res_picture

    @property
    def Assist_j_file(self):
        return self._Assist_j_file

    # _instance_lock = threading.Lock()
    # _singleton=None
    # def __new__(cls, *args, **kwargs):
    #     if not getattr(cls, '_singleton', None):
    #         cls._singleton = super().__new__(cls, *args, **kwargs)
    #     return cls._singleton

    _sys_setting = None

    @property
    def sys_setting(self):
        if self._sys_setting is None:
            self._sys_setting = self._Sys_setting()
        return self._sys_setting

    _tools_setting = None

    @property
    def tools_setting(self):
        if self._tools_setting is None:
            self._tools_setting = self._Tools_setting()
        return self._tools_setting

    class _Sys_setting(object):

        def __init__(self):
            setting = Setting.get_instance()
            self.config = setting.config["Sys_setting"]
            setting.get_loger().write("初始化系统设置")

            self._process_monitoring_model=self.config["process_monitoring"]["model"]
            self._process_monitoring_interval=self.config["process_monitoring"]["interval"]
            process_monitoring_dict = self.config["process_monitoring"]["model_dict"]
            if self._process_monitoring_model not in process_monitoring_dict.keys():
                self._process_monitoring_model=2
            if self._process_monitoring_interval<1:
                self._process_monitoring_interval=5
            setting.get_loger().write("性能监视器模式: %s 记录间隔: %s" % (process_monitoring_dict[self._process_monitoring_model],self._process_monitoring_interval))
            try:
                if self._process_monitoring_model==2:#####无界面版本等待工具箱升级
                    from InvSysTools.MyTools.Monitor.Process_monitor import process_monitoring
                    process_monitoring(interval=self._process_monitoring_interval)

            except Exception as e:
                setting.get_loger().err("性能监视启动失败:\n%s"%(str(e)))

            self._bound_min_value = self.config["Bonds_tool"]["bound_min_value"]
            self._bound_max_value = self.config["Bonds_tool"]["bound_max_value"]


            setting.get_loger().write("Bounds限制最小最大值:(%s,%s)" % (self._bound_min_value, self._bound_max_value))
            ByterChooser_config = self.config["ByterChooser"]
            self._beta_choose_mode = ByterChooser_config["beta_choose_mode"]
            if type(self._beta_choose_mode) is not int:
                setting.get_loger().write("beta_choose_mode必须是数字而不是%s" % (self._beta_choose_mode),
                                          myPrint.myPrint_Err)
                self._beta_choose_mode = ByterChooser_config["beta_choose_mode_dict"]["default"]
            self._beta_choose_mode_dict = ByterChooser_config["beta_choose_mode_dict"]
            if self._beta_choose_mode not in self._beta_choose_mode_dict.keys():
                setting.get_loger().write(
                    "beta的选择策略%d不存在,现在使用默认的%s:" % (self._beta_choose_mode, self._beta_choose_mode_dict["default"]))
                setting.get_loger().write("建议使用的值有:")
                for key in self._beta_choose_mode_dict.keys():
                    setting.get_loger().write("%s--%s" % (key, self._beta_choose_mode_dict[key]))
            else:
                setting.get_loger().write(
                    "beta的选择策略:%s" % (self._beta_choose_mode_dict[self._beta_choose_mode]))
            beta_by_X2_config = ByterChooser_config["beta_by_X2"]
            self._beta_accuracy = beta_by_X2_config["beta_accuracy"]
            self._search_mode = beta_by_X2_config["search_mode"]
            search_mode_dic = beta_by_X2_config["search_mode_dic"]
            self._beta_coefficient = beta_by_X2_config["beta_coefficient"]
            if type(self._search_mode) is not int:
                setting.get_loger().write("beta_choose_mode必须是数字而不是%s" % (self._beta_choose_mode),
                                          myPrint.myPrint_Err)
            if self._beta_choose_mode == 2:
                # 用户选择使用卡方校验,打印设置详情
                setting.get_loger().write("搜索策略:%s" % (search_mode_dic[self._search_mode]))
                setting.get_loger().write("beta_coefficient:%f" % (self._beta_coefficient))
            Air_tools = self.config["Air_tools"]
            self._Air_tools_strategy = Air_tools["strategy"]
            Ref_tool = self.config["Ref_tool"]
            self._air_value = Ref_tool["air_value"]
            self._wall_value = Ref_tool["wall_value"]
            self._inner_value = Ref_tool["inner_value"]

        @property
        def bound_min_value(self):
            return self._bound_min_value

        @property
        def bound_max_value(self):
            return self._bound_max_value

        @property
        def beta_choose_mode(self):
            return self._beta_choose_mode

        @property
        def beta_accuracy(self):
            return self._beta_accuracy

        @property
        def search_mode(self):
            return self._search_mode

        @property
        def beta_coefficient(self):
            return self._beta_coefficient

        @property
        def air_tools_strategy(self):
            return self._Air_tools_strategy

        @property
        def air_value(self):
            return self._air_value

        @property
        def wall_value(self):
            return self._wall_value

        @property
        def inner_value(self):
            return self._inner_value

        _instance_lock = threading.Lock()

        def __new__(cls, *args, **kwargs):
            if not hasattr(cls, '_instance'):
                with Setting._Sys_setting._instance_lock:
                    if not hasattr(cls, '_instance'):
                        Setting._Sys_setting._instance = super().__new__(cls)
                return Setting._Sys_setting._instance

    class _Tools_setting(object):

        def __init__(self):
            pass

        _instance_lock = threading.Lock()

        def __new__(cls, *args, **kwargs):
            if not hasattr(cls, '_instance'):
                with Setting._Tools_setting._instance_lock:
                    if not hasattr(cls, '_instance'):
                        Setting._Tools_setting._instance_lock._instance = super().__new__(cls)
                return Setting._Tools_setting._instance_lock._instance
