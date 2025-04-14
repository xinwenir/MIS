import pandas as pd
import numpy as np

__version__ = "0.1.03"

UserName = "ReaderConnect"
Host = "210.26.49.10"
Port = 4561
PassWord = "IAMREADER"

ProjectName = "zaozigou"

from .onlineReader import Cormis_DataReader_online
from .offlineReader import Cormis_DataReader_offline

class Cormis_DataReader_hdf(Cormis_DataReader_offline):
    def __init__(self,*args,**kwargs):
        super(Cormis_DataReader_hdf, self).__init__(*args,**kwargs)


