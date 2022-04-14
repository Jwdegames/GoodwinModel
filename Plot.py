# Computation imports
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQT, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

class Plot(FigureCanvasQTAgg):
    '''Makes a plot for graphical use that can be updated'''

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig = fig
        canvas = FigureCanvasQT(fig)
        self.canvas = canvas
        self.x = []
        self.y = []
        self.x2 = []
        self.y2 = []
        self.label = "_nolegend_"
        self.label2 = "_nolegend_"
        self.color = 'b'
        self.color2 = 'r'
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
        self.axes.plot(self.x, self.y, label = self.label, c = self.color)
        self.axes.plot(self.x2, self.y2, label = self.label2, c = self.color2)

    def makeToolbar(self):
        '''Makes a toolbar to interact with the graph'''
        self.toolbar = NavigationToolbar(self, self.parent)
        return self.toolbar

    def clear(self):
        '''Clears the plot'''
        self.axes.clear()
        self.axes.cla()

    def setTitle(self, title):
        '''Set's the plot's title'''
        self.title = title
        self.axes.set_title(title)

    def setXAxisTitle(self, title):
        '''Sets the x-axis of the plot'''
        self.xAxisTitle = title
        self.axes.set_xlabel(title)

    def setYAxisTitle(self, title):
        '''Sets the y-axis of the plot'''
        self.yAxisTitle = title
        self.axes.set_ylabel(title)

    def setAxesTitles(self, titleX, titleY):
        '''Sets the axis titles of the plot'''
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

    def makeVLine(self, x, c, linestyle="solid", label=""):
        '''Makes a vertical line on the plot'''
        self.axes.axvline(x=x, c=c, linestyle=linestyle, label=label)

    def setXTicks(self, ind = False, labels = False):
        '''Makes the x ticks on a set of labels'''
        if ind != False:
            self.axes.set_xTicks(ind)
        if labels != False:
            self.axes.set_xticklabels(labels)

    def showLegend(self, loc = 'center left', bbox_to_anchor=(0.93, 0.5)):
        '''Shows the legend in the plot'''
        self.axes.legend(loc = loc, bbox_to_anchor=bbox_to_anchor)
        # print("Showing legend")

    def setPos(self, x0, y0, width, height):
        '''Sets the plot's position'''
        self.axes.set_position([x0, y0, width, height])

