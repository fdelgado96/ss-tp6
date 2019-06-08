import numpy as np
import matplotlib.pyplot as plt

PATHS = [
    '../data/2.0_140.0_1e-5_false_exitTimes_{}.csv',
    '../data/6.0_140.0_1e-5_false_exitTimes_{}.csv',
    '../data/10.0_140.0_1e-5_false_exitTimes_{}.csv'
]
VELOCITIES = [
    2,
    6,
    10
]

times = []
for i in range(len(PATHS)):
    vel_times = []
    for j in range(10):
        run_times = np.genfromtxt(PATHS[i].format(j+1))
        vel_times.append(run_times[-1])
    times.append(vel_times)

times = np.vstack(times)
times_mean = times.mean(axis=1)
times_std = times.std(axis=1)
plt.errorbar(VELOCITIES, times_mean, yerr=times_std, fmt='-')

plt.xlabel('Velocidad Deseada [m/s]')
plt.ylabel('Tiempo total de evacuación [s]')
plt.legend()
plt.show()