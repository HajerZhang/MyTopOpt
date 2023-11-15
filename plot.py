import os
import re
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from tqdm import tqdm  # 引入 tqdm

# 包含数据文件的子目录路径
directory_path = "data"  # 请用实际的子目录路径替换

# 获取目录中所有文件的列表
file_list = sorted([f for f in os.listdir(directory_path) if f.endswith('.txt')], key=lambda x: int(re.search(r'\d+', x).group()))

# 创建一个列表以存储图像
images = []

# 使用 tqdm 遍历每个文件，添加进度条
for file_name in tqdm(file_list, desc="Processing Files", unit="file"):
    # 构造完整的文件路径
    file_path = os.path.join(directory_path, file_name)

    # 从文件中读取数据
    data_array = np.loadtxt(file_path)

    # 获取数组的形状
    rows, cols = data_array.shape

    # 绘制灰度图
    plt.imshow(data_array, cmap='gray_r', vmin=0, vmax=1, aspect='equal', extent=[0, cols, 0, rows])
    plt.colorbar()

    # 从文件名中提取数字
    file_number = int(re.search(r'\d+', file_name).group())

    # 设置标题和标签
    plt.title(f'Grayscale Image - {file_number}')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    # 设置坐标轴刻度，间隔为2
    plt.xticks(np.arange(0, cols+1, 5))
    plt.yticks(np.arange(0, rows+1, 5))

    # 保存图像并添加到列表中
    output_file = f"{file_number}_plot.png"
    plt.savefig(output_file)
    images.append(Image.open(output_file))
    
    # 关闭当前图形
    plt.close()

    # 删除中间图像文件
    os.remove(output_file)

# 将图像保存为 GIF 动画，将 duration 设置为 100 毫秒
gif_output_file = "./plot/result.gif"
images[0].save(gif_output_file, save_all=True, append_images=images[1:], duration=100, loop=0)

# 单独保存最后一张图像
last_image_file = "./plot/final.png"
images[-1].save(last_image_file)

# 显示最后一张图像
images[-1].show()
