#Import necessary packages
import sys

#Import PyQt6 for graphics

from PyQt6.QtCore import QSize
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *

#System path imports
import os.path
from os import path

#Import other files
import GUI;


class App(QWidget):
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
        self.height = 1080
        self.initUI()
        self.show()

    def initUI(self):
        # print(self.xRayPath)
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

    # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec())