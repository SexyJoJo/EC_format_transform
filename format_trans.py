import math
import shutil

import pandas as pd
import CONST
import os.path
from log import trans_log
from scipy import interpolate

INPUT_DIR = 'Sounding'
TEMP_DIR = 'Temp'
TEMPLATE_PATH = './template'

HIGH93 = [409, 433, 458, 483, 508, 533, 558, 583, 608, 633, 658, 683, 708, 733, 758, 783, 808, 833, 858, 883, 908, 958,
          1008, 1058, 1108, 1158, 1208, 1258, 1308, 1358, 1408, 1458, 1508, 1558, 1608, 1658, 1708, 1758, 1808, 1858,
          1908, 1958, 2008, 2058, 2108, 2158, 2208, 2258, 2308, 2358, 2408, 2658, 2908, 3158, 3408, 3658, 3908, 4158,
          4408, 4658, 4908, 5158, 5408, 5658, 5908, 6158, 6408, 6658, 6908, 7158, 7408, 7658, 7908, 8158, 8408, 8658,
          8908, 9158, 9408, 9658, 9908, 10158, 10408]


class FileUtils:
    @staticmethod
    def get_file_header(filepath):
        """读取文件头部信息"""
        fullpath = None
        for root, _, files in os.walk(filepath):
            flag = False
            for file in files:
                if file.endswith('000'):
                    flag = True
                    fullpath = os.path.join(root, file)
                    break
                else:
                    continue
            if flag:
                break

        if fullpath:
            with open(fullpath) as f:
                for i in range(3):
                    line = f.readline()
                    while line.isspace():
                        line = f.readline()
                    header = line.strip().split()
            return header

    @staticmethod
    def save_df(filename, station, df):
        """保存df为txt文件"""
        year = '20' + filename[:2]
        month = filename[2:4]
        if not os.path.exists(f"./Sounding/{station}/{year}"):
            os.mkdir(f"./Sounding/{station}/{year}")

        if not os.path.exists(f"./Sounding/{station}/{year}/{month}"):
            os.mkdir(f"./Sounding/{station}/{year}/{month}")

        trans_log.logger.info(f"保存结果文件{station}_20{filename[:-4]}0000SURP.txt")
        df.to_csv(f"./Sounding/{station}/{year}/{month}/{station}_20{filename[:-4]}0000SURP.txt", sep=' ', index=False)

        with open(f"./Sounding/{station}/{year}/{month}/{station}_20{filename[:-4]}0000SURP.txt", "r+") as f:
            content = f.read()
            f.seek(0, 0)
            f.write("ECthin Query Data\n" + content)


class EcToSounding:
    @staticmethod
    def get_nearest(location, loc_meta):
        """获取站点距离最近的格点->(经度索引坐标， 纬度索引坐标)"""
        start_longitude = float(loc_meta[2])
        end_longitude = float(loc_meta[3])
        start_latitude = float(loc_meta[4])
        end_latitude = float(loc_meta[5])
        delta_longitude = float(loc_meta[0])
        delta_latitude = float(loc_meta[1])

        x = start_longitude
        if delta_longitude >= 0:
            while x < end_longitude:
                if x > location[0]:
                    break
                else:
                    x += delta_longitude
        else:
            while x > end_latitude:
                if x < location[0]:
                    break
                else:
                    x += delta_longitude

        y = start_latitude
        if delta_latitude >= 0:
            while y < end_latitude:
                if y > location[1]:
                    break
                else:
                    y += delta_latitude
        else:
            while y > end_latitude:
                if y < location[1]:
                    break
                else:
                    y += delta_latitude

        # 站点周边的4个格点
        loc1 = (x - delta_longitude, y)
        loc2 = (x - delta_longitude, y - delta_latitude)
        loc3 = (x, y - delta_latitude)
        loc4 = (x, y)

        def cal_distance(location1, location2):
            """计算两点间距离"""
            dis = math.sqrt(
                (location1[0] - location2[0]) ** 2 + (location1[1] - location2[1]) ** 2
            )
            return dis

        dis1 = cal_distance(location, loc1)
        dis2 = cal_distance(location, loc2)
        dis3 = cal_distance(location, loc3)
        dis4 = cal_distance(location, loc4)

        min_dis = min(dis1, dis2, dis3, dis4)
        if min_dis == dis1:
            nearest = loc1
        elif min_dis == dis2:
            nearest = loc2
        elif min_dis == dis3:
            nearest = loc3
        elif min_dis == dis4:
            nearest = loc4
        else:
            nearest = None

        i_nearest = (int((nearest[0] - start_longitude) / delta_longitude),
                     int((nearest[1] - start_latitude) / delta_latitude))
        return i_nearest

    @staticmethod
    def init_df(un_df):
        """
        传入的dataframe的列名会用第一行的数据作为dataframe
        :param un_df: 待初始化的df
        :return: 初始化后的df
        """
        temp = pd.DataFrame(un_df.columns).T
        un_df.columns = temp.columns
        df = pd.concat([temp, un_df], axis=0, ignore_index=True)
        return df

    @staticmethod
    def prs2gph(prs):
        """气压转化海拔。  海拔每上升9m，大气压降低1百帕。"""
        gph = round(44331 * (1 - (prs / 1013.25) ** 0.1903), 2)  # 海拔
        return gph

    @staticmethod
    def parse_and_save(station, location, grid_point_cnt):
        """
        解析原始EC文件
        :param station: 处理的站号
        :param location: 站点距离最近的格点坐标（2元组）
        :param grid_point_cnt: 纬线条数
        """
        pressures = os.listdir(os.path.join(CONST.EC_path, 'T'))
        pressures = [int(x) for x in pressures]
        pressures.sort(reverse=True)

        root = os.path.join(CONST.EC_path, 'T', str(pressures[0]))  # T目录下1000气压目录
        # 遍历T目录下1000气压目录下的文件，每个文件对应一个结果dataframe
        for file in os.listdir(root):
            trans_log.logger.info(f"解析原始文件{file}")
            T_full_path = os.path.join(root, file)
            R_full_path = os.path.join(CONST.EC_path, 'R', str(pressures[0]), file)

            # 数据文件可能存在多种格式，分别处理
            T_line_count = len(open(T_full_path, 'r').readlines())
            R_line_count = len(open(R_full_path, 'r').readlines())

            if T_line_count == grid_point_cnt + 3:
                T_df = pd.read_table(T_full_path, skiprows=3, sep=r'\s+', header=0)
            elif T_line_count == 2 * grid_point_cnt + 5:
                T_df = pd.read_table(T_full_path, skiprows=3, sep=r'\s+', header=1)
            else:
                continue

            if R_line_count == grid_point_cnt + 3:
                R_df = pd.read_table(R_full_path, skiprows=3, sep=r'\s+', header=0)
            elif R_line_count == 2 * grid_point_cnt + 5:
                R_df = pd.read_table(T_full_path, skiprows=3, sep=r'\s+', header=1)
            else:
                continue

            # df初始化
            T_df = EcToSounding.init_df(T_df)
            R_df = EcToSounding.init_df(R_df)
            item = {
                "Second": [0],
                "TEM": [T_df.iloc[location[1], location[0]]],
                "PRS_HWC": [pressures[0]],
                "RHU": [R_df.iloc[location[1], location[0]]],
                "GPH": [EcToSounding.prs2gph(pressures[0])]
            }
            final_df = pd.DataFrame(item)

            # 追加不同气压下同名文件中的要素数据
            for pressure in pressures:
                if pressure == pressures[0]:  # 第一层气压已经过初始化
                    continue
                T_full_path = os.path.join(CONST.EC_path, 'T', str(pressure), file)
                R_full_path = os.path.join(CONST.EC_path, 'R', str(pressure), file)

                # 数据文件可能存在多种格式，分别处理
                try:
                    T_line_count = len(open(T_full_path, 'r', encoding='gb2312').readlines())
                    R_line_count = len(open(R_full_path, 'r', encoding='gb2312').readlines())
                except FileNotFoundError:
                    continue
                except UnicodeError:
                    continue

                if T_line_count == grid_point_cnt + 3:
                    T_df = pd.read_table(T_full_path, skiprows=3, sep=r'\s+', header=0)
                elif T_line_count == 2 * grid_point_cnt + 5:
                    T_df = pd.read_table(T_full_path, skiprows=3, sep=r'\s+', header=1)
                else:
                    continue

                if R_line_count == grid_point_cnt + 3:
                    R_df = pd.read_table(R_full_path, skiprows=3, sep=r'\s+', header=0)
                elif R_line_count == 2 * grid_point_cnt + 5:
                    R_df = pd.read_table(R_full_path, skiprows=3, sep=r'\s+', header=1)
                else:
                    continue

                # df初始化
                T_df = EcToSounding.init_df(T_df)
                R_df = EcToSounding.init_df(R_df)
                item = {
                    "Second": [0],
                    "TEM": [T_df.iloc[location[1], location[0]]],
                    "PRS_HWC": [pressure],
                    "RHU": [R_df.iloc[location[1], location[0]]],
                    "GPH": [EcToSounding.prs2gph(pressure)]
                }
                line = pd.DataFrame(item)

                # 合并dataframe
                final_df = pd.concat([final_df, line], axis=0)

            # 保存
            print(file)
            print(final_df)
            FileUtils.save_df(file, station, final_df)

    @staticmethod
    def ec2sounding():
        """ec原始文件转化为探空格式文件"""
        # 获取文件头的元数据信息
        trans_log.logger.info("读取文件头部信息")
        loc_info = FileUtils.get_file_header(CONST.EC_path)

        # 创建目录
        if not os.path.exists('./Sounding'):
            os.mkdir("./Sounding")

        # 处理每个站点
        for station_id in CONST.locations.keys():
            trans_log.logger.info(f"处理站点{station_id}")
            if not os.path.exists(f'./Sounding/{station_id}'):
                os.mkdir(f"./Sounding/{station_id}")
            # 找到当前站点距离最近的格点
            trans_log.logger.info("寻找最近格点")
            nearest_loc = EcToSounding.get_nearest(CONST.locations[station_id], loc_info)

            # 解析并保存
            EcToSounding.parse_and_save(station_id, nearest_loc, int(loc_info[7]))


class SoundingToMono:
    @staticmethod
    def sounding2mono():

        shutil.rmtree(TEMP_DIR, ignore_errors=True)
        os.makedirs(TEMP_DIR, exist_ok=True)

        # pool_size = os.cpu_count()
        # print('pool size: %s' % pool_size)
        # pool = Pool(pool_size)
        f = os.walk(INPUT_DIR)
        for root, dirs, files in f:
            for file in files:

                upar_file = os.path.join(root, file)
                if upar_file.upper().endswith("TXT"):
                    split_path = root.split("\\")
                    out_dir = os.path.join('Temp', split_path[1], split_path[2], split_path[3])
                    os.makedirs(out_dir, exist_ok=True)

                    # 计算83层高度数组
                    alt = CONST.DEV_ALT[int(split_path[1])]  # 海拔高度
                    base_height = [1, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425,
                                   450, 475,
                                   500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250,
                                   1300, 1350,
                                   1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2250,
                                   2500, 2750,
                                   3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250,
                                   6500, 6750,
                                   7000, 7250, 7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750, 10000]
                    height_83 = [(i + round(alt)) for i in base_height]
                    height_1000 = [i / 1000 for i in height_83]

                    print('Process input: %s' % upar_file)
                    # pool.apply_async(SoundingToMono.do_process, (upar_file, file, height_1000))
                    SoundingToMono.do_process(upar_file, file, height_1000, out_dir, alt)

        # print('Waiting for all subprocesses done...')
        # pool.close()
        # pool.join()
        # print('All subprocesses done.')

    # @staticmethod
    # def process(filedir, input_file, height):
    #     try:
    #         SoundingToMono.do_process(filedir, input_file, height)
    #     except:
    #         print('process error!')

    @staticmethod
    def do_process(filedir, input_file, height, out_dir, alt):
        height_list_m = [i * 1000 for i in height]
        df = pd.read_table(filedir, skiprows=1, header=None).iloc[:, 1:]
        df.columns = [0, 1, 2, 3]
        SoundingToMono.convert_numberic(df, [0, 1, 2, 3])
        base_height = df[3].min() - 1
        # df[3] -= base_height
        max_height = df[3].max()
        if max_height < 10000 + base_height:
            print('Skip %s, max height: %s' % (input_file, max_height))
            return

        # 对缺测列进行插值
        df.interpolate(inplace=True)

        # 生成线性插值函数
        func_temperature = interpolate.interp1d(df[3], df[0], fill_value='extrapolate')
        func_pressure = interpolate.interp1d(df[3], df[1], fill_value='extrapolate')
        func_humidity = interpolate.interp1d(df[3], df[2], fill_value='extrapolate')

        SoundingToMono.rewrite_monoRTM_template(alt)

        # 代入插值函数生成83层文件
        no_scale_file = input_file[6:20] + '.IN_NOSCALE_IATM1_dn'
        SoundingToMono.write_no_scale_file(func_humidity, func_pressure, func_temperature, max_height, no_scale_file,
                                           height=height_list_m, out_dir=out_dir)

    @staticmethod
    def convert_numberic(df, cols):
        for col in cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    @ staticmethod
    def rewrite_monoRTM_template(equip_record_alt):
        # format：83层高度层及template模板文件
        base_height = [1, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475,
                       500,
                       550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
                       1400,
                       1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2250, 2500, 2750, 3000,
                       3250,
                       3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250,
                       7500,
                       7750, 8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750, 10000]
        height_83 = [(i + round(equip_record_alt)) for i in base_height]
        height_1000 = [i / 1000 for i in height_83]

        for_h = map(lambda x: format(x, '.3f'), height_1000)
        re = ''
        count = 0
        for one in for_h:
            if count == 8:
                count = 0
                re = re + '\n'
            re = re + '     ' + one
            count = count + 1
        print('template高度模板')
        print(re)

        fre_43 = [0.74165, 0.74168, 0.74185, 0.75052, 0.76833, 0.76836, 0.76853, 0.79502, 0.79505,
                  0.79522, 0.83391, 0.84859, 0.87507, 0.87511, 0.87527, 0.92864, 0.93398, 1.00069,
                  1.04739, 1.70945, 1.70952, 1.70985, 1.72653, 1.74387, 1.76122, 1.76135, 1.77910,
                  1.77923, 1.79618, 1.79624, 1.79658, 1.81459, 1.83260, 1.85128, 1.86863, 1.88997,
                  1.91092, 1.91099, 1.91132, 1.93334, 1.93347, 1.93467, 1.96136
                  ]

        # frequence = [fre_43[int(i)] for i in equip_record.bands.split(',')]

        if os.path.exists('template'):
            os.remove('template')
        with open('template', 'w') as file:
            # 写入模板头
            with open('temphead', 'r') as temphead:
                for line in temphead:
                    file.write(line)

            file.write('\n')
            file.write(str(len(fre_43)))
            for i in fre_43:
                file.write('\n')
                file.write(str(i))

            # 设备海拔高度
            equip_alt = lambda x: format(x / 1000, '.3f')
            file.writelines(
                '\n0.000     0.000\n    0    2   83    1    0    7    1                                      \n')
            file.write('     ' + str(equip_alt(equip_record_alt + 1)) + '    ' + str(
                equip_alt(equip_record_alt + 10001)) + '     0.000')
            file.write('\n')
            file.writelines(re)
            file.writelines('\n   83 TEST\n')

        return height_1000

    @staticmethod
    def write_no_scale_file(func_humidity, func_pressure, func_temperature, max_height, no_scale_file, height, out_dir):
        with open(out_dir + "/" + no_scale_file, 'w') as file_out:
            # 写入模板头
            with open(TEMPLATE_PATH, 'r') as template:
                for line in template:
                    file_out.write(line)

            # 写入温湿压数据
            last_pressure = -1
            last_temperature = -1
            for high in height:
                if high > max_height:
                    pass

                pressure = func_pressure(high)
                if pressure == last_pressure:
                    pressure -= 0.001
                last_pressure = pressure

                temperature = 273.16 + func_temperature(high)
                if temperature == last_temperature:
                    temperature -= 0.001
                last_temperature = temperature

                humidity = func_humidity(high)

                file_out.write('%10.3f %9.4f %9.3f     %s   %s\n%10.4e %9.3f  %.6f  %.6f  %.6f  %.6f      %.2f\n' %
                               (high / 1000, pressure, temperature, 'AA', 'H222222', humidity, 0, 0, 0, 0, 0, 0))
            file_out.write('-1.\n%%%\n\n')

            print('Write %s finished' % (TEMP_DIR + no_scale_file))
        return file_out


if __name__ == '__main__':
    EcToSounding.ec2sounding()
    # SoundingToMono.sounding2mono()
