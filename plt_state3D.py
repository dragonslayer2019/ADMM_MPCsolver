import json
import matplotlib.pyplot as plt
import numpy as np

file_open = open('test.out')
file_read = file_open.read()
data = json.loads(file_read)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import pathpatch_2d_to_3d
import mpl_toolkits.mplot3d.art3d as art3d


# 创建一个Figure对象
fig = plt.figure()

# 在Figure对象中创建一个3D坐标轴
ax = fig.add_subplot(111, projection='3d')

t = [i * 0.2 for i in range(50)]
plt.plot([0.7, 0.7], [1, 1], [2.2, 3.8], c = 'black');
plt.plot([1.2, 1.2], [1, 1], [2.2, 3.8], c = 'black');
plt.plot([-0.5, -0.5],[1,1],[6.2, 7.8],  c = 'black');
plt.plot([-1.2, -1.2],[1,1],[6.2, 7.8], c = 'black');


# 计算对应的x,y,z值
X = data['x']
Y = data['y']
Z = t

# 绘制3D曲线
ax.plot(X, Y, Z)

# 定义椭圆的参数
#xc = [2]  # x坐标中心
#yc = [1]  # y坐标中心
#a = [2]  # 长半轴
#b = [1]  # 短半轴
#angles = [np.pi/4]  # 椭圆与x轴正向的夹角

# 生成t的取值范围
ti = np.linspace(2.2, 3.8, 50)
xc = np.full(ti.shape, 2)
yc = np.full(ti.shape, 1)
a = np.full(ti.shape, 2)
b = np.full(ti.shape, 1)
angles = np.full(ti.shape, np.pi/4)

def get_rotated_ellipse(xc, yc, a, b, angle):
    """ 生成旋转后的椭圆点坐标 """
    theta = np.linspace(0, 2*np.pi, 100)
    x = a * np.cos(theta)
    y = b * np.sin(theta)
    
    # 构造旋转矩阵
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, -s], [s, c]])
    
    # 旋转点坐标
    pts = np.column_stack((x, y)) @ R
    pts += np.array([xc, yc])
    
    return pts


# 绘制椭圆
for i in range(len(xc)):
    ellipse_pts = get_rotated_ellipse(xc[i], yc[i], a[i], b[i], angles[i])
    ellipse = plt.Polygon(ellipse_pts, fill=False, color='r')
    ax.add_patch(ellipse)
    art3d.pathpatch_2d_to_3d(ellipse, z=ti[i], zdir="z")

# 设置坐标轴标签
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('T')

# 显示图形
ax.auto_scale_xyz([-5, 5], [-5, 5], [0, 10])
plt.show()
