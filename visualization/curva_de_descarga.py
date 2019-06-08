import numpy as np
import matplotlib.pyplot as plt

PATHS = [
    '../data/2.0_140.0_1e-5_false_exitTimes_{}.csv',
    '../data/6.0_140.0_1e-5_false_exitTimes_{}.csv',
    '../data/10.0_140.0_1e-5_false_exitTimes_{}.csv'
]
LABELS = [
    'V = 2 m/s',
    'V = 6 m/s',
    'V = 10 m/s'
]

for i in range(len(PATHS)):
    times = []
    for j in range(10):
        run_times = np.genfromtxt(PATHS[i].format(j+1))
        times.append(run_times)
    times = np.vstack(times)
    times_mean = times.mean(axis=0)
    times_std = times.std(axis=0)
    plt.errorbar(range(1, times.shape[1] + 1), times_mean, yerr=times_std, fmt='-', label=LABELS[i])

plt.xlabel('Particulas Evacuadas')
plt.ylabel('Tiempo [s]')
plt.legend()
plt.show()

