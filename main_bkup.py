# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


from DataAnalysisTools import Statistics
from DataFactory import Setting, DataManager
import logging
import sys
import yaml
import logging.config
import os


def setup_logging(default_path="loggingCfg.yaml", default_level=logging.DEBUG, env_key="LOG_CFG"):
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, "r") as f:
            config = yaml.unsafe_load(f)
            logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)


flt_dictionary = {"forward": ['Topo2Ltopo',
                              'Ltopo2Leff',
                              'Leff2Res'],
                  "dataAnalysis": ['RawProcessing',
                                   'GetRatioLeff'],
                  "inversion": ['Inversion',
                                'SetRefBnd']}


if __name__ == '__main__':
    setup_logging(default_path="loggingCfg.yaml")

    # setting = Setting.Setting(r'/data/guorui/SynologyDrive/liugreen/muography/prj.ZaoziGou/step.ini')
    print("argv1 = ", sys.argv[1])
    setting = Setting.Setting(sys.argv[1])
    data_manager = DataManager.DataManager(setting)

    lis_steps = setting.get_step_lists()
    lis_steps_info = []

    for step in lis_steps:
        t_dic = setting.get_step_info(step)
        lis_steps_info.append(t_dic)

        logging.info(step)
        logging.info(t_dic.__str__())

    for i in range(len(lis_steps)):
        step_name = lis_steps[i]
        t_dic = lis_steps_info[i]

        dir_step = t_dic.get('OutDir', '')
        t_kwargs = t_dic.get('kwargs', {})
        t_kwargs['OutDir'] = dir_step

        t_flt = t_dic.get('Filter')
        if t_flt is None:
            raise Exception('Input Filter name for %s.\nFilter dictionary:\n%s'
                            % (step_name, flt_dictionary.__str__()))

        # logging.info('Run %s: %s\n%s' % (step_name, t_flt, t_kwargs.__str__()))
        getattr(data_manager, t_flt)(t_kwargs)
