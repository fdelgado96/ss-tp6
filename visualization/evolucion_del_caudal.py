#Asume un archivo de entrada con el tiempo de salida de cada particula
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math

PATHS = [
    '../data/6.0_140.0_1e-5_false_exitTimes_{}.csv'
]
LABELS = [
    'V = 6 m/s'
]

BUCKET_SIZE = 5


def get_sliding_window_measures(path):
    times = []
    deltas = []

    for i in range(1, 11):
        file_times = np.genfromtxt(path.format(i))
        # convertir de tiempo acumulado a diferencias
        file_deltas = np.diff(file_times)

        times.append(file_times[1:])
        deltas.append(file_deltas)

    times = np.concatenate(times)
    deltas = np.concatenate(deltas)

    bins = np.array(range(0, math.ceil(times.max()) + BUCKET_SIZE, BUCKET_SIZE))
    digitized = np.digitize(times, bins)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bin_means = [1 / deltas[digitized == i].mean() for i in range(1, len(bins))]
    bin_stds = [1 / deltas[digitized == i].std() for i in range(1, len(bins))]

    return bin_centers, bin_means, bin_stds


total_measures = []
for i in range(len(PATHS)):
    times, means, stds = get_sliding_window_measures(PATHS[i])
    m, b, r, _, _ = stats.linregress(times, means)
    plt.errorbar(times, means, fmt='-o', label=LABELS[i])
    #plt.plot([times[0], times[-1]], [m*times[0] + b, m*times[-1] + b])
    print('m: {:.2E} - b: {:.2E} - r {:.2E}'.format(m, b, r))
    print('mean: {:.2E} - std: {:.2E}'.format(np.mean(means), np.std(means)))

plt.ylabel('Caudal en partículas por segundo')
plt.xlabel('Tiempo [s]')
plt.ylim(0, 10)
plt.legend()
plt.show()
