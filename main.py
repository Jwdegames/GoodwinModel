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
from GoodwinModel import GoodwinModel
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

        # Begin making layout
        self._main = QWidget()
        self.setCentralWidget(self._main)
        mainLayout = QGridLayout(self._main)
        graphLayout = QVBoxLayout()
        mainLayout.addLayout(graphLayout, 1, 1)

        # Make plot of Phillips Curve
        gM = GoodwinModel()
        phillipsPlot = gM.makePhillipsPlot(self, False)
        phillipsPlot.plot()


        graphLayout.addWidget(phillipsPlot.makeToolbar())
        graphLayout.addWidget(phillipsPlot)
        phillipsPlot.show()

        # Make parametric plots
        parametricPlot = gM.makeParametricPlot()

        graphLayout.addWidget(parametricPlot.makeToolbar())
        graphLayout.addWidget(parametricPlot)
        parametricPlot.show()
        parametricPlot.showLegend()

        # Make Goodwin Cycle Plot
        goodwinPlot = gM.makeGoodwinPlot()

        graphLayout.addWidget(goodwinPlot.makeToolbar())
        graphLayout.addWidget(goodwinPlot)
        goodwinPlot.show()

    # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec())