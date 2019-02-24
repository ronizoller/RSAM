import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt

path = '/Users/ronizoller/Documents/school/Master/מחקר/DATA/'

times_effi = {200: 0.058, 300: 0.1227, 250: 0.0976, 100: 0.012, 150: 0.0338,400: 0.2293, 450: 0.2757, 500: 0.328, 550: 0.4594, 350: 0.1812}
times_naive = {200: 0.0679, 300: 0.2384, 250: 0.1364, 100: 0.0053, 150: 0.026,400: 0.5899, 450: 0.8389, 500: 1.1318, 550: 1.5775, 350: 0.3938}


xs = []
effi = []
naive = []
names = []
for number_of_leafs,time in times_effi.items():
    xs.append(number_of_leafs)
    effi.append(time)
for number_of_leafs,time in times_naive.items():
    naive.append(time)
fig, ax = plt.subplots()
plt.xlabel('Number of Leafs', fontsize=10)
plt.ylabel('Running Time', fontsize=10)
plt.plot(xs, effi, 'ro')
plt.plot(xs, naive, 'g^')
plt.axis([0, max(xs)+50, 0, max(effi+naive)+0.5])
fig = ax.get_figure()
fig.savefig(path + '/plot_noise.png')