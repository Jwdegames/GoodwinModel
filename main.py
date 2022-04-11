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
        self.width = 1800
        self.height = 1020
        self.initUI()
        self.show()

    def updatePlots(self):
        '''Updates the plots'''
        # print("Called update")
        # print(self.vField)
        v = self.vField.text()
        # print("V success")
        theta = self.thetaField.text()
        n = self.nField.text()
        alpha = self.alphaField.text()
        beta = self.betaField.text()
        k = self.kField.text()
        # print("Got Values")
        self.gM.updateVals(float(v), float(theta), float(n), float(alpha), float(beta), float(k))
        # print("Updated Values")
        self.gM.updatePlots()
        # print("Updated Plots")
        self.gM.phillipsPlot.draw()

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
        self.showMaximized()

        # Begin making layout
        self._main = QWidget()
        self.setCentralWidget(self._main)
        mainLayout = QGridLayout(self._main)
        leftLayout = QGridLayout()
        updateLayout = QFormLayout()
        mainLayout.addLayout(leftLayout, 0, 0)
        leftLayout.addLayout(updateLayout, 0, 0)
        graphLayout = QVBoxLayout()
        mainLayout.addLayout(graphLayout, 0, 1)


        # Add Input fields to the left side
        self.vField = QLineEdit()
        self.thetaField = QLineEdit()
        self.nField = QLineEdit()
        self.alphaField = QLineEdit()
        self.betaField = QLineEdit()
        self.kField = QLineEdit()
        updateLayout.addRow("v (Capital-Output Ratio):", self.vField)
        updateLayout.addRow("\u03B8 (Natural Growth Rate of Labor Productivity):", self.thetaField)
        updateLayout.addRow("n (Natural Growth Rate of Amount of Workers):", self.nField)
        updateLayout.addRow("\u03B1 (Negative Inflation At Full Unemployment):", self.alphaField)
        updateLayout.addRow("\u03B2 (Negative Inflation Slope):", self.betaField)
        updateLayout.addRow("k (Inflation Exponent):", self.kField)

        # Make plot of Phillips Curve
        gM = GoodwinModel()
        self.gM = gM
        phillipsPlot = gM.makePhillipsPlot(self, False)
        phillipsPlot.plot()

        # Update the fields with the initial values of the Goodwin Model
        self.vField.setText(str(gM.v))
        self.thetaField.setText(str(gM.theta))
        self.nField.setText(str(gM.n))
        self.alphaField.setText(str(gM.alpha))
        self.betaField.setText(str(gM.beta))
        self.kField.setText(str(gM.k))


        graphLayout.addWidget(phillipsPlot.makeToolbar())
        graphLayout.addWidget(phillipsPlot)
        phillipsPlot.show()

        # Make parametric plots
        parametricPlot = gM.makeParametricPlot()

        graphLayout.addWidget(parametricPlot.makeToolbar())
        graphLayout.addWidget(parametricPlot)
        parametricPlot.show()
        parametricPlot.showLegend()
        # Shift plot over
        paraPos = parametricPlot.axes.get_position()
        parametricPlot.setPos(paraPos.x0 - 0.07, paraPos.y0, paraPos.width * 0.85, paraPos.height)

        # Make Goodwin Cycle Plot
        goodwinPlot = gM.makeGoodwinPlot()

        graphLayout.addWidget(goodwinPlot.makeToolbar())
        graphLayout.addWidget(goodwinPlot)
        goodwinPlot.show()

        # Make the rest of the left side
        # Make the update button
        updateButton = GUI.makeButton(self, "Update", 0, 0, 100, 100)
        updateButton.setFont(QFont('Times', 16))
        updateButton.clicked.connect(lambda: self.updatePlots())
        leftLayout.addWidget(updateButton, 0, 1)

    # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec())