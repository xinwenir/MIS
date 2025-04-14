# author:高金磊
# datetime:2022/4/5 10:10

"""
根据输入数据生成反演所需的所有数据文件
"""
import os
import time

from scipy.sparse import csc_matrix

from InvDataTools.Cell_Weight_Tools import Cell_Weight_Tools
from InvSysTools.MyTools import myPrint
from InvDataTools import MeshTools, Gij_tools
from InvDataTools.Air_j import Air_j
from InvDataTools.Assist_j_tool import Assist_j_tool
from InvDataTools.Bonds_tools import Bonds_tool
from InvDataTools.ByterChooser import ByterChooser
from InvDataTools.Jxyz_Tools import Make_jxyz
from InvDataTools.Visibility_res_tools import Show_pred_obsd_derr, Norm_tools, Beta_log, Show_data_distribution
from InvDataTools.calculate.Calcu_sensitivity import Calcsensitivity 
from InvDataTools.d_tools import d_tools
from InvDataTools.new_res_tools import new_res_tools
from InvDataTools.obs_tools import obs_tools
from InvDataTools.ref_tools import Ref_tools
from InvDataTools.res_tools import restore_res
from InvDataFactory import Setting


class DataManager:
    '''
    严格单例模式,应该避免其他类对数据的所有操作
    '''
    _instance = None
    identification=time.time()
    @staticmethod
    def get_instance(new=False):
        """
        生成全局单例的DataManager对象
        :param new: 弃用之前的结果创建一个新的对象为全局对象
        :return: 全局的DataManager对象
        """
        if DataManager._instance is None or new:
            DataManager.identification=time.time()
            DataManager._instance = DataManager(DataManager.identification)
        return DataManager._instance

    def __init__(self,identification=1.0):
        """生成中间文件"""
        global setting
        if identification !=DataManager.identification:
            setting.get_loger().write("请通过get_instance获取对象!!",myPrint.myPrint_Err)
        setting = Setting.Setting.get_instance()
        setting.record_env_detailed()

        self._meshObj = MeshTools.MeshTools(setting.get_mesh_file())
        self.obs_tool = obs_tools(setting.get_obs_file())
        self.air_j_tool = Air_j()
        self.bonds_tool = Bonds_tool(setting.get_bnd_file(), setting.sys_setting.bound_min_value,
                                     setting.sys_setting.bound_max_value)
        self.refs_tool = Ref_tools(setting.get_refs_file())
        if os.path.exists(setting.get_sparse_matrix_file()) and not setting.get_all_new():
            setting.get_loger().write(setting.get_sparse_matrix_file() + "存在,将继续使用",
                                      printer=myPrint.myPrint_Wran)
        else:
            setting.get_loger().write(setting.get_sparse_matrix_file() + "将自动生成",
                                      printer=myPrint.myPrint_Hint)
            calcsensitivity=Calcsensitivity(mesh_tool=self.mesh,loger=setting.get_loger(),obs_file=setting.get_obs_file())
            calcsensitivity.calc_all_rays_from_obs_file(setting.get_sparse_matrix_file(),setting.get_Ray_way_j_file())
            
            
            #ryj            
            # calcsensitivity = Calcsensitivity(obsf=setting.get_obs_file(), meshf=setting.get_mesh_file())
            # calcsensitivity.calcsensitivity(setting.get_sparse_matrix_file(), setting.get_Ray_way_j_file(),
            #                                 setting.get_loger(), isprint=setting.get_log_detailed(2))
        if os.path.exists(setting.get_G_file()) and not setting.get_all_new():
            setting.get_loger().write(setting.get_G_file() + "存在,将继续使用", printer=myPrint.myPrint_Wran)
        else:
            setting.get_loger().write(setting.get_G_file() + "将自动生成", printer=myPrint.myPrint_Hint)
            Gij_tools.get_g(setting.get_sparse_matrix_file(), setting.get_G_file())
        if os.path.exists(setting.get_d_file()) and not setting.get_all_new():
            setting.get_loger().write(setting.get_d_file() + "存在,将继续使用", printer=myPrint.myPrint_Wran)
        else:
            setting.get_loger().write(setting.get_d_file() + "将自动生成", printer=myPrint.myPrint_Hint)
            self.obs_tool.make_d_from_obs(setting.get_d_file())
        if setting.get_need_jxyz_file():
            setting.get_loger().write("注意您选择了生成jxyz文件,此文件将在未来版本不受支持!!!!", printer=myPrint.myPrint_Err)
            jxyz = Make_jxyz()
            jxyz.jxyz_form_shape(self.mesh.get_shape(), setting.get_jxyz_file())
            del jxyz
        else:
            setting.get_loger().write("不生成jxyz文件", printer=myPrint.myPrint_Success())
        if setting.get_topo_file() is None:
            setting.get_loger().write("没有提供topo文件,将不再标记空气,对应的%s也将失去作用"%(setting.get_air_j_file()), printer=myPrint.myPrint_Wran)
            #在具体使用处采取行动   get_unneed_j(self):
        else:
            if os.path.exists(setting.get_air_j_file()) and not setting.get_all_new():
                setting.get_loger().write(setting.get_air_j_file() + "存在,将继续使用", printer=myPrint.myPrint_Wran)
            else:
                setting.get_loger().write(setting.get_air_j_file() + "将自动生成", printer=myPrint.myPrint_Hint)
                if setting.get_need_jxyz_file():
                    self.air_j_tool.get_air_j_to_file_Jxyzfile_cache(self.mesh,
                                                                     setting.get_topo_file(),
                                                                     setting.get_air_j_file(),
                                                                     setting.get_jxyz_file(),
                                                                     strategy=setting.sys_setting.air_tools_strategy)

                else:
                    self.air_j_tool.get_air_j_to_file(self.mesh,
                                                      setting.get_topo_file(),
                                                      setting.get_air_j_file(), strategy=setting.sys_setting.air_tools_strategy)
        if os.path.exists(setting.Assist_j_file) and not setting.get_all_new():
            setting.get_loger().write(setting.Assist_j_file + "存在,将继续使用", printer=myPrint.myPrint_Wran)
        else:
            setting.get_loger().write(setting.Assist_j_file + "将自动生成", printer=myPrint.myPrint_Hint)
            Assist_j_tool.make_assist_j_file(self.mesh.get_shape(), setting.Assist_j_file)

            # 过时的求解空气编号的策略
            # xyzcells = self.mesh.get_xyz_start()
            # self.air_j_tool.get_air_j(xyzcells, self.mesh.get_xs()[0], self.mesh.get_ys()[0], self.mesh.get_zs()[0],
            #                           setting.get_topo_file(),
            #                           setting.get_jxyz_file(),
            #                           setting.get_air_j_file()
            #                           )
        self._g_ij = Gij_tools.G_Data(setting.get_G_file(), setting.get_sparse_matrix_file())
        self.d_tool = d_tools(setting.get_d_file())
        self.cell_weight_tool=Cell_Weight_Tools(path=setting.get_Ray_way_j_file())

    _A_others = None

    def Make_A(self,new=False):
        """
        根据配置中的信息，构造方程组
        单例模式优先使用缓存
        :param new:  是否丢弃之前结果，创建新的单例对象
        :return: 方程组信息
        """
        if self._A_others is not None or new:
            return self._A_others
        vi, vj = self._get_ij()
        value = self._get_g_values()
        d = self._get_d_values()
        d_err = self._get_d_err()
        refs = self._get_refs()
        unneed_j = self.get_unneed_j()
        bounds = self._get_bounds()
        if len(vi) != len(vj):
            raise Exception("i和j的数量不相等!!")
        if len(vi) != len(value):
            raise Exception("ij和真是长度g的数量不相等!!")
        if len(refs) != len(bounds):
            raise Exception("约束和refs的数量不相等!!")

        # # #对探测器附近的点收紧上下限约束
        # ray_near_cells=self.cell_weight_tool.get_cell_weight_by_rays(amount=15)
        # setting.get_loger().waring("将会有%s个格子的约束被收紧"%(len(ray_near_cells)))
        # for i in list(ray_near_cells):
        #
        #     """以下代码可以用来显示都有那些被标记了"""
        #     # refs[i-1]=12.34
        #     # bounds[i - 1][0]=12.33
        #     # bounds[i-1][1]=12.35
        #     """收紧约束"""
        #     bounds[i-1][0]=max(bounds[i-1][0],refs[i-1]-0.1)
        #     # bounds[i-1][1]=min(bounds[i-1][1],refs[i-1]+0.3)
        # #     """减小权重"""


        ## 进一步精简方程组--合并相同或者相似的射线信息
        # from InvDataTools.Sparse_Matrix_Tools import Find_similar_is
        # similar_is,nozero_js_is_map = Find_similar_is([vi, vj])
        # #稀疏矩阵降维
        # from InvDataTools.Sparse_Matrix_Tools import Reducing_dimensions
        # vi,vj,value,d,d_err=Reducing_dimensions(vi,vj,value,d,d_err,nozero_js_is_map,similar_is)


        i_new, j_new, oldj_newj, new_value = compress_jvalue(vi,
                                                             vj, unneed_j, value,
                                                             setting.get_rm_air_cells(),
                                                             setting.get_rm_unray_cells(),
                                                             self.mesh.cells_count())
        i_new, oldi_newi = compress_ivalue(i_new)
        if setting.get_rm_unable_ray():
            m = max(i_new) + 1
            d_err_new = [0] * m
            d_new = [0] * m
            # for i in range(len(d)):
            for i in oldi_newi.keys():
                d_new[oldi_newi[i]] = d[i]
                d_err_new[oldi_newi[i]] = d_err[i]
        else:
            m = len(d)
            d_err_new = d_err
            d_new = d
            oldi_newi = {}
            for i in range(len(d)):
                oldi_newi[i] = i
        if setting.get_rm_unable_cells():
            n = max(oldj_newj.values()) + 1
            refs_new = [0] * n
            bounds_new = [[0, 0.1] for i in range(n)]
            # for oldj in range(len(bounds)):
            for oldj in oldj_newj.keys():
                try:
                    refs_new[oldj_newj[oldj]] = refs[oldj]
                except:
                    print(oldj)
                bounds_new[oldj_newj[oldj]] = bounds[oldj]
        else:
            n = self.mesh.cells_count()
            refs_new = refs
            bounds_new = bounds
            oldj_newj = {}
            for j in range(len(bounds)):
                oldj_newj[j] = j
        # """降低权重"""
        # vj_count={}
        # #统计
        # for i in range(len(vj)):
        #     if vj[i] not in vj_count.keys():
        #         vj_count[vj[i]]=0
        #     vj_count[vj[i]]=vj_count[vj[i]]+value[i]
        # #加权
        # for i in range(len(vj)):
        #     value[i]/=vj_count[vj[i]]

        self.A = csc_matrix((new_value, (i_new, j_new)), shape=(m, n))


        self._A_others = (self.A, d_new, d_err_new, bounds_new, refs_new, oldi_newi, oldj_newj,new_value,(i_new, j_new))
        return self._A_others
    
    def Make_A_by_TF(self):
        #不可以根据稀疏矩阵创建偶尔会出现数据不相等的情况
        import tensorflow as tf
        ij_new=self._A_others[8]
        indices=[]
        for i in range(len(ij_new[0])):
            indices.append([ij_new[0][i],ij_new[1][i]])
        values=self._A_others[7]
        return tf.sparse.SparseTensor(indices=indices, values=values, dense_shape=self.A.shape)
    _unneed_j = None

    def get_unneed_j(self):
        """
        获取空气对应的格子的编号
        !未验证，2023年8月8日 增加通过上下限约束阀值去除空气（不需要优化的体素）。降低内存占用，避免罚函数等因为可行域太小出现的问题
        全局单例
        :return: 格子编号（n，1）
        """
        self._unneed_j = None
        if self._unneed_j is None:
            self._unneed_j=set()
            if setting.get_topo_file() is None:
                #不标记空气
                pass
            else:
                self._unneed_j = self.air_j_tool.get_air_j_from_file(setting.get_air_j_file())
            ####!上下限约束小于某个值时，标记为空气（没有topo文件的补充版本）
            bounds=self._get_bounds()
            for i in range(len(bounds)):
                j=i+1
                #todo 取值需要来自配置文件
                if(bounds[i][1]-bounds[i][0]<=0.4 and self._refs[i]==0):#!增加参考模型不为0的检查避免有的非空气被标记
                    self._unneed_j.add(j)
        return self._unneed_j

    def record_data_detail(self):
        """
        记录当前的数据详情，包含方程组的维度，稀疏矩阵的压缩信息等
        :return:
        """
        setting.get_loger().write("方程组详情:", printer=myPrint.myPrint_Hint)
        if self._A_others is None:
            self._A_others = self.Make_A()
        A, d_new, d_err_new, bounds_new, refs_new, oldi_newi, oldj_newj,values,ijnews = self._A_others
        setting.get_loger().write("mesh shape:%s" % (str(self.mesh.get_shape())),
                                  printer=myPrint.myPrint_Hint)
        setting.get_loger().write("格子总数:%d" % (self.mesh.cells_count()), printer=myPrint.myPrint_Hint)
        setting.get_loger().write("射线总数:%d" % (len(self._get_d_values())), printer=myPrint.myPrint_Hint)
        setting.get_loger().write("有效射线总数:%d" % (len(d_new)), printer=myPrint.myPrint_Hint)
        setting.get_loger().write("方程组shape:%d  %d" % A.shape, printer=myPrint.myPrint_Hint)
        setting.get_loger().write(
            "压缩比:%f %%" % (A.count_nonzero() * 3 / self.mesh.cells_count() / len(self._get_d_values()) * 100)
            , printer=myPrint.myPrint_Hint)

    _beta_chooser = None

    def choose_beta(self)->ByterChooser:
        """
        根据设置生成β选择器
        全局单例
        :return: 一个β选择服务提供者
        """
        if self._beta_chooser is not None:
            return self._beta_chooser
        if self._A_others is None:
            self._A_others = self.Make_A()
        A, d_new, d_err_new, bounds_new, refs_new, oldi_newi, oldj_newj,values,ijnews = self._A_others
        self._beta_chooser = ByterChooser(d_new, d_err_new, setting.get_beta(), setting.sys_setting.beta_coefficient)
        return self._beta_chooser

    def Postprocessor(self, x):
        """
        做结果的后置处理
        比如可视化，结果加工等
        :param x: 一维数组
        :return:
        """
        self._output_res_data_process(x)
        if setting.show_res_picture:
            self._Visible_computer_process()
            self._Visible_res(x)
            self._show_beta_history()
            self._show_res_distribution(x)

    _misfit_refs_smooth = []

    def collect_misfit_refs_smooth(self, all, misfit, refs_x, smooth):
        """
        用来收集反演中的目标函数值,misfit,结果偏离refs的值,平滑度误差
        都是float
        :param all:
        :param misfit:
        :param refs_x:
        :param smooth:
        :return:
        """
        self._misfit_refs_smooth.append((all, misfit, refs_x, smooth))
        return True

    _beta_misfit_norms = []

    def collect_beta_misfit_norms(self, beta, misfit, norm):
        """
        用来收集 beta, misfit, norm的关系
        :param beta:
        :param misfit:
        :param norm:
        :return:
        """
        self._beta_misfit_norms.append([beta, misfit, norm])

    def _output_res_data_process(self, res):
        """
        按照用户要求对结果进行加工
        :param res:
        :return:
        """
        """恢复结果"""
        A, d_new, d_err_new, bounds_new, refs_new, oldi_newi, oldj_newj,values,ijnews = self.Make_A()
        d_pred = A * res
        refs = self._get_refs()
        shape = self.mesh.get_shape()
        res = restore_res(res, oldj_newj, shape)

        if setting.get_res_recover_refs():
            res = self.refs_tool.recover_resj_by_refj(res)
        if setting.get_res_recover_air():
            if setting.get_topo_file() is not None:
                res = self.air_j_tool.recover_resj_by_airj(res, default=setting.get_air_value())
            else:
                setting.get_loger().write("由于没有topo文件,空气不能被去除,res_recover_air被视为无效",printer=myPrint.myPrint_Err)
        tool = new_res_tools(res_list=res)
        tool.mode0(target_file=setting.get_res_file())
        """恢复差异对比"""
        ref_res = self.refs_tool.make_refs_ps_err_list(refs, res)
        # ref_res = restore_res(ref_res, oldj_newj, shape[0] * shape[1] * shape[2])
        tool = new_res_tools(res_list=ref_res)
        tool.mode0(target_file=setting.get_res_refs_file())
        """对结果进行一次平滑过度处理"""
        smooth_cells_res = tool.smooth_cells(ref_res, refs, self._unneed_j, oldj_newj.keys(), shape, self._get_bounds())
        tool = new_res_tools(res_list=smooth_cells_res)
        tool.mode0(target_file=setting.get_res_smooth_file())
        #############获取差异比较大的射线#####################
        d_obs_pred_diff=[]
        for item in oldi_newi.items():
            d_obs_pred_diff.append([item[0],(d_pred[item[1]]-d_new[item[1]])/d_err_new[item[1]]])
        d_obs_pred_diff = sorted(d_obs_pred_diff, key=lambda x: x[1], reverse=True)
        # _d_obs_pred_diff_file=r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\d_obs_pred_diff"
        # file=open(_d_obs_pred_diff_file,'w')
        # for d_diff in d_obs_pred_diff:
        #     file.write(str(d_diff[0]))
        #     file.write(' ')
        #     file.write(str(d_diff[1]))
        #     file.write('\n')
        # file.close()

    def _Visible_res(self, x):
        """
        后置处理,组织Visibility_res_tools,显示结果
        :param x: 反演结果向量
        :return:
        """
        if self._A_others is None:
            self._A_others = self.Make_A()
        A, d_new, d_err_new, bounds_new, refs_new, oldi_newi, oldj_newj,values,ijnews = self._A_others
        pred = A * x
        pred_obsd_derr = []
        for i in range(len(d_new)):
            pred_obsd_derr.append([pred[i], d_new[i], d_err_new[i]])
        tool = Show_pred_obsd_derr(data=pred_obsd_derr)
        tool.store_res(setting.get_pred_obsd_derr_file())
        tool.show_pred_obsd()
        tool.show_pred_obsd_diff(interval=0.04, max_ignore=5)

    def _show_res_distribution(self, x):
        """
        可视化结果分布
        :param x:
        :return:
        """
        tool = Show_data_distribution()
        #显示res_ps的结果
        res_ps=[]
        for i in range(len(x)):
            res_ps.append(x[i]-self._refs[i])
        tool.show_obsx_num(res_ps, 200)

    def _Visible_computer_process(self):
        """
        可视化计算过程
        :return:
        """
        tool = Norm_tools(all_misfit_ms_smooth_list=self._misfit_refs_smooth)
        tool.store_res(file=setting.get_all_misfit_ms_smooth_file())
        tool.show_norm()

    def _show_beta_history(self):
        """
        可视化β的历史记录
        :return:
        """
        tool = Beta_log(setting.get_beta_misfit_norms_file())
        for beta, misfit, norm in self._beta_misfit_norms:
            tool.record(beta, misfit, norm)
        tool.show()
    @property
    def mesh(self):
        """
        提供mesh对象
        :return:
        """
        return self._meshObj
    """
    以下数据恐怕在内部发生改变而导致歧义故不共享
    """
    _g_values = None

    def _get_g_values(self):
        if self._g_values is None:
            self._g_values = self._g_ij.get_GV()
        return self._g_values

    _ij = None

    def _get_ij(self):
        if self._ij is None:
            ij = self._g_ij.get_ij()
            self._ij = [[i - 1 for i in ij[0]], [i - 1 for i in ij[1]]]
        return self._ij

    _d_value = None

    def _get_d_values(self):
        if self._d_value is None:
            self._d_value = self.d_tool.get_d()
        return self._d_value

    _d_err = None

    def _get_d_err(self):
        if self._d_err is None:
            self._d_err = self.d_tool.get_d_err()
        return self._d_err

    _bounds = None

    def _get_bounds(self):
        if self._bounds is None:
            self._bounds = self.bonds_tool.get_bonds_min_max()
            ###################临时代码,假设没有上下限约束的场景###############################
            # for i in range(len(self._bounds)):
            #     # self._bounds[i][0]=min(0,self._bounds[i][0])
            #     if self._bounds[i][1]==1.73:
            #         self._bounds[i][1]=1.8
            #     if self._bounds[i][0]==1.68:
            #         self._bounds[i][0]=1.3


        return self._bounds

    _refs = None

    def _get_refs(self):
        if self._refs is None:
            self._refs = self.refs_tool.get_data()
            ###################临时代码###############################
            # for i in range(len(self._refs)):
            #     if self._refs[i]>1:
            #         self._refs[i]=1.7
        return self._refs

    # _instance_lock = threading.Lock()
    # _instance=None
    # def __new__(cls, *args, **kwargs):
    #     if not hasattr(cls, '_instance'):
    #         with DataManager._instance_lock:
    #             if not hasattr(cls, '_instance'):
    #                 DataManager._instance = super().__new__(cls)
    #         return DataManager._instance


def compress_jvalue(vi, vj, unneed_j, value, _rm_air_cells, rm_unray_cells, all_cells_count):
    """
    压缩矩阵,去掉没有用的j
    :param vi: 稀疏矩阵行向量
    :param vj: 稀疏矩阵的列向量
    :param unneed_j: 属于空气的j
    :param value: (i,j)对应的value
    :param _rm_air_cells:  是否移除属于空气的格子
    :param rm_unray_cells:  是否移除没有射线穿过的格子
    :param all_cells_count:  对于没有射线穿过的格子为了保证其位置存在,对方程组填充使得自变量维度达到all_cells_count
    :return:
    """
    oldj_newj = {}
    if not _rm_air_cells and not rm_unray_cells:
        # for j in vj:
        #     oldj_newj[j]=j
        oldj_newj = None
        return vi, vj, oldj_newj, value
    count = 0
    new_value = []
    if type(unneed_j) != set:
        unneed_j = set(unneed_j)

    if _rm_air_cells:
        # 需要去掉空气
        if not rm_unray_cells:
            # 保留没有射线穿过的格子
            for j in range(0, all_cells_count - 1):
                if j + 1 not in unneed_j:
                    oldj_newj[j] = count
                    count += 1
            # if max(oldj_newj.keys()) != count - 1:
            #     # 最后一个或者很多个格子因为都是0,在构建矩阵的时候会出现被误删除
            #     vi.append(vi[0])
            #     vj.append(count - 1)
            #     value.append(0)
        else:
            # 需要去掉没有射线穿过的格子
            for j in vj:
                if j + 1 not in unneed_j and j not in oldj_newj.keys():
                    oldj_newj[j] = count
                    count += 1
    else:
        pass
    j_new = []
    i_new = []
    for i in range(len(vi)):
        if vj[i] in oldj_newj.keys():
            j_new.append(oldj_newj[vj[i]])
            i_new.append(vi[i])
            new_value.append(value[i])
    return i_new, j_new, oldj_newj, new_value


def compress_ivalue(i_old):
    """
    压缩方程组,去掉没有意义的等式
    :param i_old:
    :return:
    """
    effective = Setting.Setting.get_instance().get_rm_unable_ray()
    oldi_newi = {}
    if not effective:
        # 不去无用的射线

        # for i in i_old:
        #     i_new.append(i)
        #     oldi_newi[i]=i
        oldi_newi = None
        return i_old, oldi_newi
    i_new = []
    count = 0
    for i in i_old:
        if i not in oldi_newi.keys():
            oldi_newi[i] = count
            count += 1
        i_new.append(oldi_newi[i])
    return i_new, oldi_newi
