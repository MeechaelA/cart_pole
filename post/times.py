#!/usr/bin/python3
import json
import pathlib
import numpy as np
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt


def main():
    print("Plotting Times")

    fig_0 = plt.figure()
    fig_1 = plt.figure()
    fig_2 = plt.figure()

    ax_0 = fig_0.gca(projection='3d')
    # ax_0.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    ax_1 = fig_1.gca(projection='3d')
    # ax_1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    ax_2 = fig_2.gca(projection='3d')
    # ax_2.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    sims = pathlib.Path("../build/study").glob("*.sim")
    times = pathlib.Path("../build/study").glob("*.time")

    threads = []
    num_simulations = []
    total_times = []
    outers_max = []
    inners_max = []
    
    for (i, time) in enumerate(times):
        with open(time) as time_data:
            data = time_data.read()
            json_data = json.loads(data)

            thread = json_data.get("num_threads")
            threads.append(thread)

            num_simulation = json_data.get("num_simulations")
            num_simulations.append(num_simulation)

            total_time = json_data.get("total_time")
            total_times.append(total_time)

            # for outer_time in outer_times:
            outer_times = json_data.get("outer_times")
            outers_max.append(np.median(outer_times))

            # for inner_time in inner_times:
            inner_times = json_data.get("inner_times")
            inners_max.append(np.median(inner_times))

            ax_0.scatter(thread, num_simulation, total_time)
            ax_1.scatter(thread, num_simulation, np.median(outer_times))
            ax_2.scatter(thread, num_simulation, np.median(inner_times))

    # ax_0.plot_trisurf(np.array(threads), np.array(num_simulations), np.array(total_times), cmap=cm.jet)
    ax_0.tricontourf(np.array(threads), np.array(num_simulations), np.array(total_times), zdir='z', offset=0, cmap=cm.coolwarm)

    # ax_1.plot_trisurf(np.array(threads), np.array(num_simulations), np.array(outers_max), cmap=cm.jet)
    ax_1.tricontourf(np.array(threads), np.array(num_simulations), np.array(outers_max), zdir='z', offset=0, cmap=cm.coolwarm)

    # ax_2.plot_trisurf(np.array(threads), np.array(num_simulations), np.array(inners_max), cmap=cm.jet)
    ax_2.tricontourf(np.array(threads), np.array(num_simulations), np.array(inners_max), zdir='z', offset=0, cmap=cm.coolwarm)

    ax_0.set_xlabel("Threads [#]")
    ax_0.set_ylabel("Simulations [#]")
    ax_0.set_zlabel("Time [s]")
    ax_0.dist = 15
    ax_0.view_init(30, -45)

    ax_1.set_xlabel("Threads [#]")
    ax_1.set_ylabel("Simulations [#]")
    ax_1.set_zlabel("Time [s]")
    ax_1.dist = 15
    ax_1.view_init(30, -45)

    ax_2.set_xlabel("Threads [#]")
    ax_2.set_ylabel("Simulations [#]")
    ax_2.set_zlabel("Time [s]")
    ax_2.dist = 15
    ax_2.view_init(30, -45)


    fig_0.savefig('total_times.png')
    fig_1.savefig('outer_times.png')
    fig_2.savefig('inner_times.png')

if __name__ == '__main__':
    main()