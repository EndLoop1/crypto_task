import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

histogram = pd.read_csv('output/histogram.csv')
gap_hist = pd.read_csv('output/gap_hist.csv')
normal_data = np.loadtxt('output/data_normal.csv')
uniform_data = np.loadtxt('output/data_uniform.csv')

# ===== 频率分布图 =====
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.bar(histogram['index'], histogram['value'], color='skyblue')
plt.title('nouce(rand() % N)')
plt.xlabel('(0 ~ N-1)')
plt.ylabel('times')
plt.grid(True, linestyle='--', alpha=0.5)

# ===== 间隔分布图 =====
plt.subplot(1, 2, 2)
plt.plot(gap_hist['gap'], gap_hist['count'], marker='o', linestyle='-', color='orange')
plt.title('Interval Distribution')
plt.xlabel('Length')
plt.ylabel('Times')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

# 绘制正态分布直方图
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.hist(normal_data, bins=50, color='skyblue', edgecolor='black', density=True)
plt.title("Normal Distribution")
plt.xlabel("Value")
plt.ylabel("Density")
plt.grid(True, linestyle='--', alpha=0.5)

# 均匀分布图（适配整数或浮点）
plt.subplot(1, 2, 2)
n = int(max(uniform_data)) + 1
bins = np.arange(-0.5, n + 0.5, 1)
plt.hist(uniform_data, bins=bins, color='lightgreen', edgecolor='black', density=True)
plt.title("Uniform Distribution")
plt.xlabel("Value")
plt.ylabel("Density")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
