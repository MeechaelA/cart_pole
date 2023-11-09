import csv
import matplotlib.pyplot as plt

class Plotter:
    def __init__(self):
        self.headers = []
        self.data = {}

    def open_file(self, file_path):
        with open(file_path, newline='\n') as csvfile:
            reader = csv.DictReader(csvfile)
            for (i,row) in enumerate(reader):
                if (i != 0):
                    for (i, item) in enumerate(row):
                        print(i, item)
                else:
                    for item in row:
                        self.headers.append(item)
                        self.data[item] = []

    def plot(self):
        fig, ax = plt.subplots(len(self.headers)-1)
        print(self.data[self.headers[0]][1])
        for (i, header) in enumerate(self.headers[1:]):
            ax[i].plot(self.data[self.headers[0]])
            ax[i].plot(self.data[self.headers[i]])
            ax[i].set_xlabel(self.headers[0])
            ax[i].set_ylabel(self.headers[i])
            ax[i].set_title(self.headers[i])

        fig.savefig('output.png')