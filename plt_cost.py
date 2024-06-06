import json
import matplotlib.pyplot as plt
import numpy as np
import math

file_open = open('testori.out')
file_read = file_open.read()
data = json.loads(file_read)
plt.xlabel('iteration')

#file_open = open('testacc.out')
#file_read = file_open.read()
#dataacc = json.loads(file_read)
#plt.xlabel('Time(s)')

#file_open = open('testori.out')
#file_read = file_open.read()
#dataori = json.loads(file_read)
#plt.xlabel('Time(s)')

data_log = [math.log2(x) for x in data['cost']]
#data_log_acc = [math.log2(x) for x in dataacc['cost']]
#data_log_ori = [math.log2(x) for x in dataori['cost']]
# t = [i * 0.2 for i in range(25)] + [i * 0.5 + 5 for i in range(5)]
t = [i for i in range(len(data_log))]
#plt.plot(t, data_log_ori, label = 'cost, no acceleration')
plt.plot(t, data_log, label = 'cost, sigma = 1')
#plt.plot(t, data_log_acc, label = 'cost, sigma = 0.1')
# plt.axhline(y = 0.6, c = 'r', ls = '--');
# plt.axhline(y = 1.4, c = 'r', ls = '--');
# plt.axhline(y = -0.4, c = 'g', ls = '--');
# plt.axhline(y = -1.4, c = 'g', ls = '--');
plt.legend()
# plt.title('x0=0, K=5000, Complicated Case2, wei = 10, weig = 50, t=0.021s')
plt.savefig('tmp.pdf')
plt.show()
