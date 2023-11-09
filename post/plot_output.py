#!/usr/bin/python3
import sys
from plotter import Plotter

def main():
    args = sys.argv[1:]
    plotter = Plotter()
    plotter.open_file(args[0])
    plotter.plot()

if __name__ == "__main__":
    main()