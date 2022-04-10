# Computation imports
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, FigureCanvasQT
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import scipy as sp
from scipy.integrate import odeint
from scipy.special import lambertw

#Graphical Imports
from PyQt6.QtCore import QSize
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *

class Plot(FigureCanvasQTAgg):
    '''Makes a plot for graphical use that can be updated'''

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig = fig
        canvas = FigureCanvas(fig)
        self.canvas = canvas
        self.x = []
        self.y = []
        self.x2 = []
        self.y2 = []
        self.label = []
        self.label2 = []
        self.axes = fig.add_subplot(111)
        super(Plot, self).__init__(fig)
        self.parent = parent

    def setX(self, array):
        '''Set x coordinates of plot'''
        self.x = array

    def setX2(self, array):
        '''Set 2nd x coordinates of plot'''
        self.x2 = array

    def setY(self, array):
        '''Set y coordinates of plot'''
        self.y = array

    def setY2(self, array):
        '''Set 2nd y coordinates of plot'''
        self.y2 = array

    def plot(self):
        '''Plots the coordinates'''
        self.axes.plot(self.x, self.y, label = self.label)
        self.axes.plot(self.x2, self.y2, label = self.label2)

    def makeToolbar(self):
        '''Makes a toolbar to interact with the graph'''
        self.toolbar = NavigationToolbar(self, self.parent)
        return self.toolbar

    def clear(self):
        '''Clears the plot'''
        self.axes.cla()

    def setTitle(self, title):
        self.title = title
        self.axes.set_title(title)

    def setXAxisTitle(self, title):
        self.xAxisTitle = title
        self.axes.set_xlabel(title)

    def setYAxisTitle(self, title):
        self.yAxisTitle = title
        self.axes.set_ylabel(title)

    def setAxesTitles(self, titleX, titleY):
        self.setXAxisTitle(titleX)
        self.setYAxisTitle(titleY)

    def setLabel(self, label):
        '''Set the label to be used for the first plot in a legend'''
        self.label = label

    def setLabel2(self, label):
        '''Set the label to be used for the second plot used in a legend'''
        self.label2 = label

    def makeHLine(self, y, c, linestyle = "solid", label = ""):
        '''Makes a horizontal line on the plot'''
        self.axes.axhline(y = y, c = c, linestyle = linestyle, label=label)

    def showLegend(self):
        '''Shows the legend in the plot'''
        self.axes.legend()

