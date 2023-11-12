import os
import sys
import csv
import math
import random
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
        self.study_dir = ''
        self.headers = []
        self.data = {}

    def plot(self, study_dir):
        self.open_study(study_dir)
        key_list = list(self.data.keys())
        col_headers = list(self.data[key_list[0]].keys())
        for key in col_headers:
            plot_widget = pg.plot()
            for key_list_key in key_list:
                vb = plot_widget.getViewBox()
                vb.setAutoVisible(y=True)
                
                # generate something to export
                plot_widget.plot(self.data[key_list_key][col_headers[0]], self.data[key_list_key][key])

            # create an exporter instance, as an argument give it
            # the item you wish to export
            exporter = pg.exporters.ImageExporter(plot_widget.plotItem)

            # set export parameters if needed
            exporter.parameters()['width'] = 1000   # (note this also affects height parameter)

            # save to file
            exporter.export(key + ".png")


    def open_study(self, study_dir):
        csv_files = []
        with os.scandir(study_dir) as it:
            for entry in it:
                if entry.name.endswith('.csv') and entry.is_file():
                    csv_files.append(entry.name)
        
        for file in csv_files:
            key = file.replace('.csv','')
            self.data[key] = {}
            self.open_file(key, study_dir + "/" + file)
    

    def open_file(self, key, file_path):
        with open(file_path, newline='\n') as csvfile:
            reader = csv.DictReader(csvfile)
            for (i,row) in enumerate(reader):
                if (i != 0):
                    for (i, item) in enumerate(row):
                        self.data[key][item].append(float(row[item]))
                else:
                    for item in row:
                        self.headers.append(item)
                        self.data[key][item] = []



from PyQt6.QtWidgets import QGraphicsTransform, QPushButton, QSlider, QFileDialog, QMainWindow, QApplication, QWidget, QHBoxLayout, QGridLayout, QGraphicsItem, QGraphicsView, QGraphicsScene, QGraphicsRectItem, QGraphicsEllipseItem
from PyQt6.QtCore import Qt, QPointF, QTimer, QVariantAnimation
from PyQt6.QtGui import QColor, QBrush, QTransform

class GUIPlotter(MultiPlotter):
    def __init__(self):
        self.study_dir = ''
        self.headers = []
        self.data = {}
        self.plots = []
        self.carts = []
        self.poles = []
        self.balls = []
        self.times_ms = []
        self.norm_times = [0]
        
        self.right_layout = QGridLayout()
        self.plots_layout = QGridLayout()
        self.button_layout = QGridLayout()
        self.iteration = 0

        return

    def start_graphics_view(self):
        self.graphics_view = QGraphicsView()


    def load_data(self):
        directory = QFileDialog.getExistingDirectory(None, 'Simulation File Dialog', os.getcwd())
        if (directory):
            self.plots.clear()
            self.open_study(directory)
            return True
        return False

    def load_simulation(self):
        loaded = self.load_data()
        self.discrete_update()
        self.set_iteration(0)


        self.iteration_slider = QSlider()
        print(self.norm_times[-1])     
        self.iteration_slider.setMinimum(0)
        self.iteration_slider.setTickInterval(1)
        self.iteration_slider.setMaximum(len(self.times)-1)
        self.iteration_slider.sliderMoved.connect(self.set_iteration)       
        self.button_layout.addWidget(self.iteration_slider)


        # if (loaded):
            # key_list = list(self.data.keys())
            # col_headers = list(self.data[key_list[0]].keys())
            # for key in col_headers:
                # plot_widget = pg.plot()
                # for key_list_key in key_list:
                    # vb = plot_widget.getViewBox()
                    # vb.setAutoVisible(y=True)
                    # 
                    # generate something to export
                    # plot_widget.plot(self.data[key_list_key][col_headers[0]], self.data[key_list_key][key])
                    # self.plots.append(plot_widget)
            # 
        # else:
            # print("No Sim")

    def update_plots(self):
        for i in reversed(range(self.plots_layout.count())):
            self.plots_layout.itemAt(i).widget().setParent(None)

        for i in self.plots:
            self.plots_layout.addWidget(i)

    def generate_cartpole(self):
        cart_width = 20
        cart_height = 20
        mass_radius = 20
        pole_length = 3 * cart_height
        pole_width = cart_width / 2.0
        move_scale = 1.0

        self.scale_times(self.data[list(self.data.keys())[0]]["time"])

        for (i, simulation) in enumerate(self.data.keys()):

            x_offset = self.graphics_view.size().width() / 2
            y_offset = self.graphics_view.size().height() / 2

            cart = Cart(self.data[simulation]["cart_dis"][0], y_offset, cart_width, cart_height)
            cart.set_x_offset(x_offset)
            cart.set_y_offset(y_offset)
            cart.set_width(cart_width)
            cart.set_height(cart_height)
            cart.set_times(self.data[simulation]["time"])
            cart.set_move_scale(move_scale)
            cart.set_origin(x_offset, y_offset)
            cart.set_positions(self.data[simulation]["cart_dis_linear"])
            self.carts.append(cart)

            pole = Pole(self.data[simulation]["pole_dis"][0], y_offset,  pole_width, pole_length)
            pole.set_x_offset(x_offset)
            pole.set_y_offset(y_offset)
            pole.set_width(pole_width)
            pole.set_length(pole_length)
            pole.set_times(self.data[simulation]["time"])
            pole.set_move_scale(move_scale)
            pole.set_origin(x_offset, y_offset)
            pole.set_positions(self.data[simulation]["cart_dis_linear"])
            self.poles.append(pole)

            circle = Circle(self.data[simulation]["pole_dis"][0], y_offset, mass_radius)
            circle.set_x_offset(x_offset - mass_radius/ 2.0)
            circle.set_y_offset(y_offset)
            circle.set_radius(mass_radius)
            circle.set_times(self.data[simulation]["time"])
            circle.set_move_scale(move_scale)
            circle.set_origin(x_offset - mass_radius / 2.0, y_offset)
            circle.set_rotations(self.data[simulation]["pole_dis_linear"])
            circle.set_positions(self.data[simulation]["cart_dis_linear"])
            self.balls.append(circle)

            self.graphics_view.scene().addItem(cart)
            self.graphics_view.scene().addItem(pole)
            self.graphics_view.scene().addItem(circle)

    def set_iteration(self, iteration):
        self.iteration = iteration
        for (i, cart) in enumerate(self.carts):
            self.carts[i].set_iteration(self.iteration)
            self.poles[i].set_iteration(self.iteration)
            self.balls[i].set_iteration(self.iteration)


    def discrete_update(self):
        self.generate_cartpole()
        # self.update_plots()


    def start(self):       
        app = QApplication(sys.argv)
        window = QMainWindow()
        self.start_graphics_view()

        rhs = QWidget()
        plots = QWidget()
        buttons = QWidget()
        load_sim_file_button = QPushButton("Load Simulation")
        load_sim_file_button.clicked.connect(self.load_simulation)
        load_sim_file_button.clicked.connect(self.update_plots)
        self.button_layout.addWidget(load_sim_file_button)
        

    
        rhs.setLayout(self.right_layout)
        self.right_layout.addWidget(plots)
        self.right_layout.addWidget(buttons)

        ## Create a grid layout to manage the widgets size and position
        plots.setLayout(self.plots_layout)
        buttons.setLayout(self.button_layout)

        animation_layout = QGridLayout()

        scene = QGraphicsScene(self.graphics_view)
        self.graphics_view.setScene(scene)
        self.graphics_view.setSceneRect(0, 0, 600, 600)

        animation_layout.addWidget(self.graphics_view)

        window_layout = QHBoxLayout()
        widget = QWidget()
        window_layout.addWidget(self.graphics_view)
        window_layout.addWidget(rhs)
        widget.setLayout(window_layout)
        window.setCentralWidget(widget)


        # Create a Qt widget, which will be our window.
        window.show()  # IMPORTANT!!!!! Windows are hidden by default.

        # Start the event loop.
        app.exec()

    def scale_times(self, times):
        self.times = times
        time_min = min(self.times)
        time_max = max(self.times)
        self.times_ms = [time * 100 for time in self.times]
        for (i, value) in enumerate(self.times):
            self.norm_times.append((value - time_min) / (time_max - time_min))
        print(self.norm_times)

class Cart(QGraphicsRectItem):
    def __init__(self, x, y, width, height):
        super().__init__(0, 0, width, height)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges, True)

        self.setPos(x, y)
        brush = QBrush()
        brush.setColor(Qt.GlobalColor.black)
        self.setBrush(brush)

        self._animation = QVariantAnimation(duration=1000)
        self._animation.valueChanged.connect(self.setPos)
        self._animation.finished.connect(self.create_random_point)
        
        self.current_iteration = 0
        self.times=[]
        self.positions = []

        self.x_offset = 0
        self.y_offset = 0 
        self.width = 0
        self.height = 0
        self.move_scale = 1


    def create_random_point(self):
        self._animation.setStartValue(self.pos())
        self._animation.setEndValue(QPointF(self.positions[self.iteration], self.pos().y()))
        self._animation.start()
        self.iteration += 1

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            brush = QBrush()
            if self.collidingItems():
                brush.setColor(Qt.GlobalColor.blue)
            else:
                brush.setColor(Qt.GlobalColor.black)
            self.setBrush(brush)
        return super().itemChange(change, value)

    def set_positions(self, positions):
        self.positions = positions

    def set_times(self, times):
        self.times = times
        time_min = min(self.times)
        time_max = max(self.times)
        self.times_ms = [time * 100 for time in self.times]
        for (i, value) in enumerate(self.times):
            self._animation.setDuration(math.ceil(self.times_ms[i]))
            norm_time = (value - time_min) / (time_max - time_min)
            self._animation.setKeyValueAt(norm_time, self)

    def set_iteration(self, iteration):
        self.setPos(self.move_scale * self.positions[iteration] - self.width / 2.0, 0)

    def set_x_offset(self, x_offset):
        self.x_offset = x_offset
    
    def set_y_offset(self, y_offset):
        self.y_offset = y_offset

    def set_width(self, width):
        self.width = width

    def set_height(self, height):
        self.height = height

    def set_move_scale(self, move_scale):
        self.move_scale = move_scale

    def set_origin(self, x, y):
        transform = QTransform()
        transform.translate(x, y)
        self.setTransform(transform)        


class Pole(QGraphicsRectItem):
    def __init__(self, x, y, width, height):
        super().__init__(0, 0, width, height)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges, True)

        self.setPos(x, y)
        brush = QBrush()
        brush.setColor(Qt.GlobalColor.black)
        self.setBrush(brush)

        self._animation = QVariantAnimation(duration=1000)
        self._animation.valueChanged.connect(self.setPos)
        self._animation.finished.connect(self.create_random_point)
        
        self.current_iteration = 0
        self.times=[]
        self.positions = []

        self.x_offset = 0
        self.y_offset = 0 
        self.width = 0
        self.length = 0
        self.move_scale = 1


    def create_random_point(self):
        self._animation.setStartValue(self.pos())
        self._animation.setEndValue(QPointF(self.positions[self.iteration], self.pos().y()))
        self._animation.start()
        self.iteration += 1

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            brush = QBrush()
            if self.collidingItems():
                brush.setColor(Qt.GlobalColor.blue)
            else:
                brush.setColor(Qt.GlobalColor.black)
            self.setBrush(brush)
        return super().itemChange(change, value)

    def set_positions(self, positions):
        self.positions = positions

    def set_times(self, times):
        self.times = times
        time_min = min(self.times)
        time_max = max(self.times)
        self.times_ms = [time * 100 for time in self.times]
        for (i, value) in enumerate(self.times):
            self._animation.setDuration(math.ceil(self.times_ms[i]))
            norm_time = (value - time_min) / (time_max - time_min)
            self._animation.setKeyValueAt(norm_time, self)

    def set_iteration(self, iteration):
        self.setPos(self.move_scale * self.positions[iteration] - self.width / 2.0, self.width)

    def set_x_offset(self, x_offset):
        self.x_offset = x_offset
    
    def set_y_offset(self, y_offset):
        self.y_offset = y_offset

    def set_width(self, width):
        self.width = width

    def set_length(self, length):
        self.length = length

    def set_move_scale(self, move_scale):
        self.move_scale = move_scale

    def set_origin(self, x, y):
        transform = QTransform()
        transform.translate(x, y)
        self.setTransform(transform)        


class Circle(QGraphicsEllipseItem):
    def __init__(self, x, y, r):
        super().__init__(0, 0, r, r)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges, True)

        self.setPos(x, y)
        brush = QBrush()
        brush.setColor(Qt.GlobalColor.black)
        self.setBrush(brush)

        self._animation = QVariantAnimation()

        self.times=[]
        self.times_ms = []
        self.positions = []
        self.rotations = []
        self.rotation_value = 0
        self.iteration = 0

        self.x_offset = 0
        self.y_offset = 0 

        self.x = 0
        self.y = 0

        self.radius = 0
        self.move_scale = 1
        self.transform = QTransform()

    def create_random_point(self):
        self._animation.setStartValue(self.pos())
        self._animation.setEndValue(QPointF(self.positions[self.iteration], self.pos().y()))
        self._animation.start()
        self.iteration += 1

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            brush = QBrush()
            if self.collidingItems():
                brush.setColor(Qt.GlobalColor.blue)
            else:
                brush.setColor(Qt.GlobalColor.black)
            self.setBrush(brush)
        return super().itemChange(change, value)

    def set_positions(self, positions):
        self.positions = positions

    def set_rotations(self, rotations):
        self.rotations = rotations

    def set_times(self, times):
        self.times = times
        time_min = min(self.times)
        time_max = max(self.times)
        self.times_ms = [time * 100 for time in self.times]
        for (i, value) in enumerate(self.times):
            self._animation.setDuration(math.ceil(self.times_ms[i]))
            norm_time = (value - time_min) / (time_max - time_min)
            self._animation.setKeyValueAt(norm_time, self)

    def set_radius(self, radius):
        self.radius = radius

    def set_x_offset(self, x_offset):
        self.x_offset = x_offset
    
    def set_y_offset(self, y_offset):
        self.y_offset = y_offset

    def set_iteration(self, iteration):
        # + x_offset - mass_radius / 2.0, y_offset + pole_length, mass_radius)
        self.set_rotation(self.x_origin+3*self.radius, self.y_origin+3*self.radius, self.rotations[iteration] + math.pi/2 + math.pi - math.pi/4)
        self.x = self.x_origin + self.move_scale * self.positions[iteration] + self.x_rotation 
        self.y = self.y_origin + self.y_rotation
        self.setPos(self.x, self.y)
       

    def set_move_scale(self, move_scale):
        self.move_scale = move_scale

    def set_origin(self, x, y):
        self.x_origin = x/2.0
        self.y_origin = y/2.0
        # self.transform.translate(self.x_origin, self.y_origin)
        # self.setTransform(self.transform)        

    def set_rotation(self, x_location, y_location, rotation):
        self.x_rotation = self.x_origin + (x_location - self.x_origin) * math.cos(rotation) - (y_location - self.y_origin) * math.sin(rotation)
        self.y_rotation = self.y_origin + (y_location - self.x_origin) * math.sin(rotation) + (y_location - self.y_origin) * math.cos(rotation)
        