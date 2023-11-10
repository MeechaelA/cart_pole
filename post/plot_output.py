#!/usr/bin/python3
import sys
from plotter import Plotter

def main():
    args = sys.argv[1:]
    plotter = Plotter()
    plotter.plot_monte_carlo(args[0])

if __name__ == "__main__":
    main()