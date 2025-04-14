from DatabaseTool import *
from DataFactory import CommonUtils, ROOTUtils


class DataframeTools:
    def __init__(self):
        pass

    @staticmethod
    def get_all_caching(h5_name: str, directory: str, pageCapacity: int = 100000):
        mtx_data = dic_data = None
        pageIndex = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, pageCapacity, directory):
            t_mtx_data, t_dic_data = CommonUtils.cvt_dataframe_to_numpy(i)

            if pageIndex == 0:
                dic_data = t_dic_data
                mtx_data = t_mtx_data
            else:
                mtx_data = np.vstack((mtx_data, t_mtx_data))

            pageIndex += 1

        return mtx_data, dic_data
