#!/usr/bin/python3
import json
import pathlib
import numpy as np
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt

def amdahl(thread_num, num_simulation, total_time, outer_time, inner_time):
    print()


def main():
    print("Plotting Times")

    fig_0 = plt.figure(figsize=(7.0, 7.0))
    fig_1 = plt.figure(figsize=(7.0, 7.0))
    fig_2 = plt.figure(figsize=(7.0, 7.0))

    ax_0 = fig_0.add_subplot(projection='3d')
    # ax_0.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    ax_1 = fig_1.add_subplot(projection='3d')
    # ax_1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    ax_2 = fig_2.add_subplot(projection='3d')
    # ax_2.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    sims = pathlib.Path("../studies").glob("*.sim")
    times = pathlib.Path("../studies").glob("*.time")

    sim_num = 60 # Change to compare which # of simulations to compare per processor count
    one_total_times = []
    one_outer_times = []
    one_inners_times = []

    five_total_times = []
    five_outer_times = []
    five_inners_times = []

    ten_total_times = []
    ten_outer_times = []
    ten_inners_times = []

    fifteen_total_times = []
    fifteen_outer_times = []
    fifteen_inners_times = []

    twenty_total_times = []
    twenty_outer_times = []
    twenty_inners_times = []

    twenty_five_total_times = []
    twenty_five_outer_times = []
    twenty_five_inners_times = []

    thirty_total_times = []
    thirty_outer_times = []
    thirty_inners_times = []

    sixty_total_times = []
    sixty_outer_times = []
    sixty_inners_times = []

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
            num_simulation = json_data.get("num_simulations")

            if (thread == 1 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                one_total_times.append(total_time)


            if (thread == 5 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                five_total_times.append(total_time)
                # # for outer_time in outer_times:
                outer_times = json_data.get("outer_times")
                one_outer_times.append(np.max(outer_times))

                # for inner_time in inner_times:
                inner_times = json_data.get("inner_times")
                one_inners_times.append(np.max(inner_times))
               
            if (thread == 10 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                ten_total_times.append(total_time)

                # for outer_time in outer_times:
                outer_times = json_data.get("outer_times")
                ten_outer_times.append(np.max(outer_times))

                # for inner_time in inner_times:
                inner_times = json_data.get("inner_times")
                ten_inners_times.append(np.max(inner_times))

            if (thread == 15 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                fifteen_total_times.append(total_time)

                # for outer_time in outer_times:
                outer_times = json_data.get("outer_times")
                fifteen_outer_times.append(np.max(outer_times))

                # for inner_time in inner_times:
                inner_times = json_data.get("inner_times")
                fifteen_inners_times.append(np.max(inner_times))

            if (thread == 20 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                twenty_total_times.append(total_time)

                # for outer_time in outer_times:
                outer_times = json_data.get("outer_times")
                twenty_outer_times.append(np.max(outer_times))

                # for inner_time in inner_times:
                inner_times = json_data.get("inner_times")
                twenty_inners_times.append(np.max(inner_times))

            if (thread == 25 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                twenty_five_total_times.append(total_time)

                # for outer_time in outer_times:
                outer_times = json_data.get("outer_times")
                twenty_five_outer_times.append(np.max(outer_times))

                # for inner_time in inner_times:
                inner_times = json_data.get("inner_times")
                twenty_five_inners_times.append(np.max(inner_times))


            if (thread == 30 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                thirty_total_times.append(total_time)

                # for outer_time in outer_times:
                outer_times = json_data.get("outer_times")
                thirty_outer_times.append(np.max(outer_times))

                # for inner_time in inner_times:
                inner_times = json_data.get("inner_times")
                thirty_inners_times.append(np.max(inner_times))

            if (thread == 60 and num_simulation == sim_num):
                total_time = json_data.get("total_time")
                sixty_total_times.append(total_time)

                # for outer_time in outer_times:
                outer_times = json_data.get("outer_times")
                sixty_outer_times.append(np.max(outer_times))

                # for inner_time in inner_times:
                inner_times = json_data.get("inner_times")
                sixty_inners_times.append(np.max(inner_times))


            threads.append(thread)

            num_simulations.append(num_simulation)

            total_time = json_data.get("total_time")
            total_times.append(total_time)

            # for outer_time in outer_times:
            outer_times = json_data.get("outer_times")
            outers_max.append(max(outer_times))

            # for inner_time in inner_times:
            inner_times = json_data.get("inner_times")
            inners_max.append(max(inner_times))

            ax_0.scatter(thread, num_simulation, total_time)
            ax_1.scatter(thread, num_simulation, max(outer_times))
            ax_2.scatter(thread, num_simulation, max(inner_times))

    # ax_0.plot_trisurf(np.array(threads), np.array(num_simulations), np.array(total_times), cmap=cm.coolwarm)
    # print(threads, num_simulations, total_times)
    tc_0 = ax_0.tricontourf(np.array(threads), np.array(num_simulations), np.array(total_times), zdir='z', offset=0, cmap=cm.coolwarm)
    fig_0.colorbar(tc_0, orientation="horizontal", pad=0.01, label="(s)")

    # ax_1.plot_trisurf(np.array(threads), np.array(num_simulations), np.array(outers_max), cmap=cm.coolwarm)
    tc_1 = ax_1.tricontourf(np.array(threads), np.array(num_simulations), np.array(outers_max), zdir='z', offset=0, cmap=cm.coolwarm)
    fig_1.colorbar(tc_1, orientation="horizontal", pad=0.01, label="(s)")


    # ax_2.plot_trisurf(np.array(threads), np.array(num_simulations), np.array(inners_max), cmap=cm.coolwarm)
    tc_2 = ax_2.tricontourf(np.array(threads), np.array(num_simulations), np.array(inners_max), zdir='z', offset=0, cmap=cm.coolwarm)
    fig_2.colorbar(tc_2, orientation="horizontal", pad=0.01, label="(s)")

    ax_0.set_xlabel("Threads [#]")
    ax_0.set_ylabel("Simulations [#]")
    ax_0.set_zlabel("Time [s]")
    # ax_0.dist = 15
    ax_0.view_init(30, -45)

    ax_1.set_xlabel("Threads [#]")
    ax_1.set_ylabel("Simulations [#]")
    ax_1.set_zlabel("Time [s]")
    # ax_1.dist = 15
    ax_1.view_init(30, -45)

    ax_2.set_xlabel("Threads [#]")
    ax_2.set_ylabel("Simulations [#]")
    ax_2.set_zlabel("Time [s]")
    # ax_2.dist = 15
    ax_2.view_init(30, -45)

    fig_0.savefig('total_times.png', dpi=500)
    fig_1.savefig('outer_times.png', dpi=500)
    fig_2.savefig('inner_times.png', dpi=500)

    cost = 1*10e-6
    print("5 Max Factor Speedup:",  compare_speedup_times(one_total_times, five_total_times),           ",", compare_efficiency_times(one_total_times, five_total_times, 5),            ",",cost*five_total_times[0], ",", five_total_times)
    print("10 Max Factor Speedup:", compare_speedup_times(one_total_times, ten_total_times),            ",", compare_efficiency_times(one_total_times, ten_total_times, 10),            ",",cost*ten_total_times[0], ",", ten_total_times)
    print("15 Max Factor Speedup:", compare_speedup_times(one_total_times, fifteen_total_times),        ",", compare_efficiency_times(one_total_times, fifteen_total_times, 15),        ",",cost*fifteen_total_times[0], ",", fifteen_total_times)
    print("20 Max Factor Speedup:", compare_speedup_times(one_total_times, twenty_total_times),         ",", compare_efficiency_times(one_total_times, twenty_total_times, 20),         ",",cost*twenty_total_times[0], ",", twenty_total_times)
    print("25 Max Factor Speedup:", compare_speedup_times(one_total_times, twenty_five_total_times),    ",", compare_efficiency_times(one_total_times, twenty_five_total_times, 25),    ",",cost*twenty_five_total_times[0], ",", twenty_five_total_times)
    print("30 Max Factor Speedup:", compare_speedup_times(one_total_times, thirty_total_times),         ",", compare_efficiency_times(one_total_times, thirty_total_times, 30),         ",",cost*thirty_total_times[0], ",", thirty_total_times)
    print("60 Max Factor Speedup:", compare_speedup_times(one_total_times, sixty_total_times),          ",", compare_efficiency_times(one_total_times, sixty_total_times, 60),          ",",cost*sixty_total_times[0], ",", sixty_total_times)


def compare_speedup_times(one_time, end_time):
    factor = one_time[0] / end_time[0]
    return factor

def compare_efficiency_times(one_time, end_time, num_threads):
    factor = (one_time[0] / (num_threads * end_time[0]))
    return factor

if __name__ == '__main__':
    main()


# 3.7595052959726165 , 0.7519010591945232 , 0.0002274653 , [22.74653]
# 8.805587422039808 , 0.8805587422039808 , 9.711527000000001e-05 , [9.711527]
# 10.718990639464064 , 0.7145993759642709 , 7.977962e-05 , [7.977962]
# 10.242681839501849 , 0.5121340919750924 , 8.348956e-05 , [8.348956]
# 10.125630678414662 , 0.4050252271365865 , 8.445469e-05 , [8.445469]
# 8.719541503617384 , 0.29065138345391284 , 9.807362e-05 , [9.807362]
# 7.085566326953351 , 0.11809277211588917 , 0.00012069000000000001 , [12.069]    