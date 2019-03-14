import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
import draw
import seaborn as sns

times_effi ={800: 1.0713, 850: 1.0365, 900: 1.2679, 950: 1.3324, 700: 0.7143, 750: 0.7823,400: 0.2293, 450: 0.2757, 500: 0.328, 350: 0.1812,600: 0.5669, 550: 0.4321,650: 0.6606}
times_naive ={ 800: 0,850:0,900:0,950:0,750: 4.2031,400: 0.5899, 450: 0.8389, 500: 1.1318, 350: 0.3938,600: 2.0923, 550: 1.571,650: 2.5381, 700: 3.2661}

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
plt.plot(xs, effi, 'g^')
plt.plot(xs, naive, 'ro')
plt.axis([300, max(xs)+50, 0, max(effi+naive)+0.5])
plt.savefig('/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/running_time.png')
