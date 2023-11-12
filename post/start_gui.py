#!/usr/bin/python3
import sys
from plotter import GUIPlotter

def main():
    args = sys.argv[1:]
    plotter = GUIPlotter()
    plotter.start()

if __name__ == "__main__":
    main()