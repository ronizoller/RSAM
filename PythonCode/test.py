import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
import draw
import seaborn as sns

times_effi ={200: 0.058, 300: 0.1227, 250: 0.0976, 100: 0.012, 150: 0.0338,950: 1.4324, 1000: 1.6553,400: 0.2293, 450: 0.2757, 500: 0.328, 350: 0.1812,600: 0.5669, 550: 0.4321,800: 1.0328, 850: 1.1174, 900: 1.2949, 650: 0.6406, 700: 0.7542, 750: 0.9349}
times_naive ={200: 0.0679, 300: 0.2384, 250: 0.1364, 100: 0.0053, 150: 0.026,400: 0.5899, 450: 0.8389, 500: 1.1318, 350: 0.3938,600: 2.0923, 550: 1.571,800: 5.0184, 850: 6.1313, 900: 7.1809, 650: 2.5381, 700: 3.2661, 750: 4.2031}

xs = []
xs_naive = []
effi = []
naive = []
names = []
for number_of_leafs,time in times_effi.items():
    xs.append(number_of_leafs)
    effi.append(time)
for number_of_leafs,time in times_naive.items():
    xs_naive.append(number_of_leafs)
    naive.append(time)
plt.xlabel('Number of Leafs', fontsize=10)
plt.ylabel('Running Time', fontsize=10)
fig, ax = plt.subplots()
plt.plot(xs, effi, 'g^')
plt.plot(xs_naive, naive, 'ro')
plt.axis([50, max(xs)+50, 0, max(effi+naive)+0.5])
plt.savefig('/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/running_time.png')
