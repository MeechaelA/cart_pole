import os
import csv
import pyqtgraph as pg
import pyqtgraph.exporters



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
                        self.data[item].append(float(row[item]))
                else:
                    for item in row:
                        self.headers.append(item)
                        self.data[item] = []

    def plot(self, name, file_path):
        self.open_file(file_path)
        for (i, header) in enumerate(self.headers):
            if i != 0:

                plot_widget = pg.plot(title = self.headers[i])
                vb = plot_widget.getViewBox()
                vb.setAutoVisible(y=True)

                # generate something to export
                plot_widget.plot(self.data[self.headers[0]], self.data[self.headers[i]])

                # create an exporter instance, as an argument give it
                # the item you wish to export
                exporter = pg.exporters.ImageExporter(plot_widget.plotItem)

                # set export parameters if needed
                exporter.parameters()['width'] = 1000   # (note this also affects height parameter)

                # save to file
                exporter.export(name + "_" + self.headers[i] + ".png")

    def plot_monte_carlo(self, directory):
        csv_files = []
        with os.scandir(directory) as it:
            for entry in it:
                if entry.name.endswith('.csv') and entry.is_file():
                    csv_files.append(entry.name)
        
        for file in csv_files:
            self.plot(file.replace('.csv',''), directory+ "/" + file)

class MultiPlotter:
    def __init__(self):
        self.headers = []
        self.data = {}