import json
import matplotlib.pyplot as plt
import numpy as np

file_open = open('test1.out')
file_read = file_open.read()
data = json.loads(file_read)

import numpy as np
import matplotlib.pyplot as plt

# 创建一个Figure对象
#fig = plt.figure()
fig, ax = plt.subplots(figsize=(8, 6))

t = [i * 0.2 for i in range(50)]
plt.plot(data['x'], data['y'])
plt.plot([-2, 4], [0.9, 0.9], c = 'b');
plt.plot([-2, 4], [1.3, 1.3], c = 'b');
plt.plot([1.7, 1.7], [1.1, 1.9], c = 'b');
plt.plot([2.2, 2.2], [1.1, 1.9], c = 'b');

plt.plot([-2, 4], [3.1, 3.1], c = 'b');
plt.plot([-2, 4], [3.9, 3.9], c = 'b');
plt.plot([-1.2, -1.2], [3.1, 3.9], c = 'b');
plt.plot([-0.5, -0.5], [3.1, 3.9], c = 'b');

ti = np.linspace(2.2, 3.8, 50)
xc = np.full(ti.shape, 2)
yc = np.full(ti.shape, 1)
a = np.full(ti.shape, 1)
b = np.full(ti.shape, 0.5)
angles = np.full(ti.shape, np.pi/3)

def get_rotated_ellipse(xc, yc, a, b, angle):
    """ 生成旋转后的椭圆点坐标 """
    theta = np.linspace(0, 2*np.pi, 100)
    x = a * np.cos(theta)
    y = b * np.sin(theta)
    
    # 构造旋转矩阵
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, s], [-s, c]])
    
    # 旋转点坐标
    pts = np.column_stack((x, y)) @ R
    pts += np.array([xc, yc])
    
    return pts


# 绘制椭圆
for i in range(len(ti)):
    pts = get_rotated_ellipse(xc[i], yc[i], a[i], b[i], angles[i])
    ax.plot(pts[:, 0], pts[:, 1], 'r')

plt.show()