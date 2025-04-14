DatabaseTool 是田恒开发的数据库接口，统一了新旧两种 Hdf5 文件格式，为 Hdf5 文件的处理和分析奠定基础，如需查看旧的 Hdf5d 格式可参考 baseCore，新的 Hdf5 格式是利用 pandas 自带功能 to_hfd5()以
key=”data”存储数据表， Hdf5 文件中不包含其他探测器本身性质的信息；
