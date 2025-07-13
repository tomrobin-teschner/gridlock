import matplotlib.pyplot as plt

with open("temp.dat", "r") as f:
  data = f.read().split()

x = []
for i in range(0, len(data)):
  x.append(float(data[i]))

plt.plot(x)
plt.show()
