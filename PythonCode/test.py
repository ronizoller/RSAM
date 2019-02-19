import matplotlib.pyplot as plt

times_effi = {200: 0.058, 300: 0.1227, 250: 0.0976, 100: 0.012, 150: 0.0338}
times_naive = {200: 0.0679, 300: 0.2384, 250: 0.1364, 100: 0.0053, 150: 0.026}


xs = []
effi = []
naive = []
names = []
for number_of_leafs,time in times_effi.items():
    xs.append(number_of_leafs)
    effi.append(time)
for number_of_leafs,time in times_naive.items():
    naive.append(time)
plt.xlabel('Number of Leafs', fontsize=10)
plt.ylabel('Running Time', fontsize=10)
fig, ax = plt.subplots()
plt.plot(xs, effi, 'ro')
plt.plot(xs, naive, 'g^')
plt.axis([0, max(xs)+50, 0, max(effi+naive)+0.5])
plt.show()