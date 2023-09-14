from radiodata import RD
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import random
from matplotlib import cm

global t_0, time_range, time_resolution, a
t_0 = 1  # s
time_range = 40
time_resolution = 6 # number of symbols after koma
a = 2


class ToyModel():
    def __init__(self, nucleus_name):
        try:
            self.nucleus_of_interest = nucleus_name
        except NameError:
            print(f'There is no nucleus named {nucleus_name} in this radiodata. :c')

    def timemoment_of_decay(self, half_life_time):
        y = random.random()
        tau = half_life_time / np.log(2)
        t = (-1) * tau * np.log(1 - y)
        return round(t, time_resolution)


    def energy_smearing(self, energy):
            # sigma = np.sqrt(1 + 1 * energy) / 2.355
            sigma = 25 / 2.355
            return round(random.gauss(energy, sigma), 3)

    def data_generating(self, primary_nucleus, gen):
        time, energy = [], []
        for _ in range(gen):
            nucleus = primary_nucleus
            t = 0
            while nucleus in RD.keys():
                nucleus_data = RD.get(nucleus)
                half_life_time = nucleus_data[0]
                t += self.timemoment_of_decay(half_life_time)
                if t == 0:
                    continue

                channel_probabilities = [el[0] for el in RD.get(nucleus)[1]]
                if sum(channel_probabilities) != 1:
                    channel_probabilities.append(1-sum(channel_probabilities))
                    decay_channel = random.choices(list(range(len(RD.get(nucleus)[1]))), weights=channel_probabilities, k=1)[0]
                    if decay_channel == len(channel_probabilities) - 1:
                        continue
                else:
                    decay_channel = random.choices(list(range(len(RD.get(nucleus)[1]))), weights=channel_probabilities, k=1)[0]
                decay_energy = nucleus_data[1][decay_channel][2]
                decay_energy = self.energy_smearing(decay_energy)

                time.append(t)
                energy.append(decay_energy)

                nucleus = nucleus_data[1][decay_channel][3]
        return time, energy


if __name__ == "__main__":
    primary_nucleus = '215Th'
    num_of_gens = 10000
    model = ToyModel(primary_nucleus)

    time, energy = model.data_generating(primary_nucleus, num_of_gens)
    with open(f'{primary_nucleus}_{num_of_gens}_ms.txt', 'w') as f:
        f.write(','.join([str(i) for i in time]) + '\n' + ','.join([str(i) for i in energy]))

    exit()

    with open(f'{primary_nucleus}_{num_of_gens}.txt', 'r') as f:
        file_data = f.read()
    time, energy = file_data.split('\n')
    time = [round(float(el), 6) * 1e3 for el in time.split(',')]
    energy = [float(el) for el in energy.split(',')]

    t = [0] + [10 * 4 ** i for i in range(-3, 10)]

    # for j in range(len(t) - 1):
    #     # temp_time = [np.log(0.00001 + el) for el in time].copy()
    #     temp_time = []
    #     temp_energy = []
    #     for i in range(len(time)):
    #         # if time[i] >= t[j] and time[i] <= t[j+1]:
    #         #     temp_time.append((time[i]))
    #         #     # temp_time.append((time[i] - t[j]))
    #         #     temp_energy.append(energy[i])
    #         if time[i] >= t[j]:
    #             # temp_time.append((time[i]))
    #             temp_time.append((time[i] - t[j]))
    #             temp_energy.append(energy[i])
    #     temp_time = [np.log(0.0000001 + el) for el in temp_time]
    #
    #     plt.hist2d(temp_energy, temp_time, bins=(100, 100), cmap=plt.cm.jet, cmin=1)
    #     plt.colorbar()
    #     # plt.ylim(-4, 20)
    #     # plt.xlim(6800, 7500)

        # plt.savefig(f'./pics/log(t{j+1}-t{j})_vs_E.png')
        # plt.clf()

    time = [np.log(0.00001 + el) for el in time]
    plt.hist2d(energy, time, bins=(300, 300), cmap=plt.cm.jet, cmin=1)
    plt.colorbar()
    plt.show()
    plt.savefig('./pics/log(t)_vs_E.png')
    plt.clf()

    exit()

    fig = plt.figure()  # create a canvas, tell matplotlib it's 3d
    ax = fig.add_subplot(111, projection='3d')

    # make histogram stuff - set bins - I choose 20x20 because I have a lot of data
    hist, xedges, yedges = np.histogram2d(energy, time, bins=(200, 200))
    xpos, ypos = np.meshgrid(xedges[:-1] + xedges[1:], yedges[:-1] + yedges[1:])

    xpos = xpos.flatten() / 2.
    ypos = ypos.flatten() / 2.
    zpos = np.zeros_like(xpos)

    dx = xedges[1] - xedges[0]
    dy = yedges[1] - yedges[0]
    dz = hist.flatten()

    cmap = cm.get_cmap('jet')  # Get desired colormap - you can change this!
    max_height = np.max(dz)  # get range of colorbars so we can normalize
    min_height = np.min(dz)
    # scale each z to [0,1], and get their rgb values
    rgba = [cmap((k - min_height) / max_height) for k in dz]

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
    plt.title("X vs. Y Amplitudes for ____ Data")
    plt.xlabel("My X data source")
    plt.ylabel("My Y data source")
    plt.savefig("Your_title_goes_here")
    plt.show()