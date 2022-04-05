import os
from datetime import datetime, timedelta
import netCDF4
import pandas as pd
import CONST
from format_trans import EcToSounding

def get_time(delta):
    start = datetime(1900, 1, 1)
    delta = timedelta(hours=int(delta))
    return start + delta

def get_loc_meta(dataset):
    """获取nc文件位置描述信息"""
    longitudes = dataset['longitude'][:]
    latitudes = dataset['latitude'][:]
    meta = [
        longitudes[1] - longitudes[0],
        latitudes[1] - latitudes[0],
        longitudes[0],
        longitudes[-1],
        latitudes[0],
        latitudes[-1]
    ]
    return meta

def get_heights(prs_list):
    heights = []
    for prs in prs_list:
        heights.append(EcToSounding.prs2gph(prs))
    return heights

def save_file(obs_time, station, df):
    filename = str(station) + '_' + datetime.strftime(obs_time, '%Y%m%d%H%M%S') + 'SURP.txt'
    print(filename)
    year = obs_time.year
    month = obs_time.month
    if not os.path.exists(f"./Sounding/{station}/{year}/{month}"):
        os.makedirs(f"./Sounding/{station}/{year}/{month}")

    df.to_csv(f"./Sounding/{station}/{year}/{month}/{filename}", sep=' ', index=False)
    with open(f"./Sounding/{station}/{year}/{month}/{filename}", "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write("ECthin Query Data\n" + content)


if __name__ == '__main__':
    for root, _, files in os.walk(r"D:\Data\microwave radiometer\Origin EC\NC"):
        for file in files:
            path = os.path.join(root, file)
            dataset = netCDF4.Dataset(path)
            loc_meta = get_loc_meta(dataset)
            print(file)
            # 遍历各个时间的各个气压层， 提取其中的温度、相对湿度数据  一个时间的一个站点对应一个ec文件
            for ti, hour in enumerate(dataset['time'][:]):  # ti: time_index
                # hour为1900年至观测时间的小时数， 需转化为标准时间
                time = get_time(hour)
                data = {
                    'PRS_HWC': dataset['level'][:].tolist()[::-1],
                    'GPH': get_heights(dataset['level'][:].tolist())[::-1]
                }
                temperatures = []
                humiditys = []
                # 对每个站点进行处理
                for station in CONST.locations.keys():
                    print(f'正在处理{time} {station}')
                    nearest_loc = EcToSounding.get_nearest(CONST.locations[station], loc_meta)
                    for li, level in enumerate(dataset['level'][:]):  # li: level_index
                        # 温度湿度平面(二维矩阵)
                        t_grid = dataset['t'][ti][li]
                        r_grid = dataset['r'][ti][li]
                        temperatures.append(round(t_grid[nearest_loc[1]][nearest_loc[0]] - 273.15, 2))
                        humiditys.append(round(r_grid[nearest_loc[1]][nearest_loc[0]], 2))
                    data['TEM'] = temperatures[::-1]
                    data['RHU'] = humiditys[::-1]
                    df = pd.DataFrame(data, columns=['TEM', 'PRS_HWC', 'RHU', 'GPH'])
                    df.insert(0, 'Second', 0)
                    print(df)
                    save_file(time, station, df)