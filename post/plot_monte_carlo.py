#!/usr/bin/python3
import sys
from plotter import MultiPlotter

def main():
    args = sys.argv[1:]
    plotter = MultiPlotter()
    plotter.plot(args[0])

if __name__ == "__main__":
    main()