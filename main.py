#Import necessary packages
import sys

#Import PyQt6 for graphics

from PyQt6.QtCore import QSize
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *

#System path imports
import os.path
from os import path

# Math imports
import numpy as np

#Import other files
import GUI
import GoodwinModel as GM
from Plot import Plot


class App(QMainWindow):
    # Define variables
    settings = {

    }

    def __init__(self):
        super().__init__()
        self.title = "Goodwin Model Simulator"
        # Dimensions for App
        self.left = 10
        self.top = 10
        self.width = 1920
        self.height = 1020
        self.initUI()
        self.show()

    def initUI(self):
        '''Initialize necessary UI elements'''
        # print(self.xRayPath)
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Center Self
        qtRectangle = self.frameGeometry()
        centerPoint = self.screen().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())

        # Make plot of Phillips Curve
        phillipsPlot = Plot(self, 5, 4, 100)
        unemployment = np.linspace(0, 1, 100)
        alpha = 1
        beta = 0.1
        inflation = -alpha + beta * (1 - unemployment)
        print(len(unemployment))
        print(len(inflation))
        phillipsPlot.setX(unemployment)
        phillipsPlot.setY(inflation)
        phillipsPlot.plot()
        phillipsPlot.show()
        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QVBoxLayout(self._main)
        layout.addWidget(phillipsPlot.makeToolbar())
        layout.addWidget(phillipsPlot)

    # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec())