#Import necessary packages
import sys

#Import PyQt6 for graphics

from PyQt6.QtCore import *
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
        #print("Called update")
        #print(self.vField)
        v = self.vField.text()
        #print("V success")
        theta = self.thetaField.text()
        n = self.nField.text()
        alpha = self.alphaField.text()
        beta = self.betaField.text()
        k = self.kField.text()
        #print("Got Values")
        self.gM.updateVals(float(v), float(theta), float(n), float(alpha), float(beta), float(k))
        #print("Updated Values")
        self.gM.updatePlots()
        #print("Updated Plots")
        self.updateStats()

    def togglePhillips(self):
        '''Switches the mode of the Phillips Curve'''
        if self.gM.exponential == True:
            self.phillipsButton.setText("Switch To Exponential Phillips Curve")
            self.alphaField.setText(str(0.6))
            self.betaField.setText(str(1))
            # Disable exponent field
            self.kField.setReadOnly(True)
            GUI.setBGColor(self.kField, "128, 128, 128")
        else:
            self.phillipsButton.setText("Switch To Linear Phillips Curve")
            self.alphaField.setText(str(1))
            self.betaField.setText(str(0.1))
            self.kField.setText(str(4.5))
            # Enable exponent field
            self.kField.setReadOnly(False)
            GUI.setBGColor(self.kField, "255, 255, 255")
        self.vField.setText(str(5))
        self.thetaField.setText(str(0.009))
        self.nField.setText(str(0.075))
        self.gM.exponential = not self.gM.exponential

        self.updatePlots()

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

        graphLayout = QVBoxLayout()
        mainLayout.addLayout(graphLayout, 0, 1, 1, 1)
        leftLayout.addLayout(updateLayout, 0, 0, 1, 2)
        mainLayout.addLayout(leftLayout, 0, 0, 1, 1)



        # Add Input fields to the left side
        self.vField = QLineEdit()
        self.thetaField = QLineEdit()
        self.nField = QLineEdit()
        self.alphaField = QLineEdit()
        self.betaField = QLineEdit()
        self.kField = QLineEdit()
        self.kField.setReadOnly(True)
        GUI.setBGColor(self.kField, "128, 128, 128")
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
        goodwinPlot.showLegend(bbox_to_anchor = (1, 0.5))
        # Shift plot over
        gPos = goodwinPlot.axes.get_position()
        goodwinPlot.setPos(gPos.x0 - 0.07, gPos.y0, gPos.width * 0.85, gPos.height)

        # Make the rest of the left side
        # Make the update and phillips switch button
        updateButton = GUI.makeButton(self, "ㅤㅤㅤㅤㅤㅤㅤUpdate Plotsㅤㅤㅤㅤㅤㅤㅤ", 0, 0, 100, 100)
        updateButton.setFont(QFont('Times', 16))
        updateButton.clicked.connect(lambda: self.updatePlots())
        phillipsButton = GUI.makeButton(self, "Switch To Exponential Phillips Curve", 0, 0, 100, 100)
        self.phillipsButton = phillipsButton
        phillipsButton.setFont(QFont('Times', 16))
        phillipsButton.clicked.connect(lambda: self.togglePhillips())
        leftLayout.addWidget(updateButton, 1, 0)
        leftLayout.addWidget(phillipsButton, 1, 1)

        # Add stats so the user can understand data better
        self.statLabel = QLabel()
        self.statLabel.setFont(QFont('Times', 16))
        self.updateStats()
        leftLayout.addWidget(self.statLabel, 2, 0, 1, 2)
        #leftLayout.addWidget(self.wsEquiLabel, 3, 0)
        #leftLayout.addWidget(self.erUpperLabel, 4, 0)
        #leftLayout.addWidget(self.erLowerLabel, 5, 0)
        #leftLayout.addWidget(self.wsUpperLabel, 6, 0)
        #leftLayout.addWidget(self.wsLowerLabel, 7, 0)

    def updateStats(self):
        '''Updates the stats label with data from the Goodwin Model'''
        #stats = "Employment Rate Equilibrium: " + str(self.gM.erEqui) + "\n"
        #stats += "Wage Share Equilibrium: " + str(self.gM.wsEqui) + "\n"
        #stats += "Employment Rate Upper Bound: " + str(self.gM.erUpperB) + "\n"
        #stats += "Employment Rate Lower Bound: " + str(self.gM.erLowerB) + "\n"
        #stats += "Wage Share Upper Bound: " + str(self.gM.wsUpperB) + "\n"
        #stats += "Wage Share Lower Bound: " + str(self.gM.wsLowerB)
        stats = "Employment Rate Equilibrium:\t {erEqui: 7.5f}\n"
        stats += "Wage Share Equilibrium:\t\t {wsEqui: 7.5f}\n"
        stats += "Employment Rate Upper Bound:\t {erUpperB: 7.5f}\n"
        stats += "Employment Rate Lower Bound:\t {erLowerB: 7.5f}\n"
        stats += "Wage Share Upper Bound:\t {wsUpperB: 7.5f}\n"
        stats += "Wage Share Lower Bound:\t {wsLowerB: 7.5f}"
        stats = stats.format(
            erEqui = self.gM.erEqui,
            wsEqui = self.gM.wsEqui,
            erUpperB = self.gM.erUpperB,
            erLowerB = self.gM.erLowerB,
            wsUpperB = self.gM.wsUpperB,
            wsLowerB = self.gM.wsLowerB
        )
        self.statLabel.setText(stats)



    # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec())