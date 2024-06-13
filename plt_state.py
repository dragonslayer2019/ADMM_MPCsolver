import json
import matplotlib.pyplot as plt
import numpy as np

file_open = open('test.out')
file_read = file_open.read()
data = json.loads(file_read)
plt.xlabel('Time(s)')

# t = [i * 0.2 for i in range(25)] + [i * 0.5 + 5 for i in range(5)]
t = [i * 0.2 for i in range(50)]
plt.plot(t, data['x'], label = 'x')
plt.plot(t, data['y'], label = 'y')
plt.plot(t, data['Xref'], label = 'Xref')
plt.plot(t, data['v'], label = 'v')
plt.plot(t, data['a'], label = 'a')
plt.plot(t, data['u'], label = 'u')
plt.plot([2.2, 3.8], [1.7, 1.7], c = 'black');
plt.plot([2.2, 3.8], [2.2, 2.2], c = 'black');
plt.plot([6.2, 7.8], [-0.5, -0.5], c = 'black');
plt.plot([6.2, 7.8], [-1.2, -1.2], c = 'black');
# plt.axhline(y = 0.6, c = 'r', ls = '--');
# plt.axhline(y = 1.4, c = 'r', ls = '--');
# plt.axhline(y = -0.4, c = 'g', ls = '--');
# plt.axhline(y = -1.4, c = 'g', ls = '--');
plt.legend()
# plt.title('x0=0, K=5000, Complicated Case2, wei = 10, weig = 50, t=0.021s')
plt.savefig('tmp.pdf')
plt.show()
