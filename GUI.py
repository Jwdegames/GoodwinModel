#Import PyQt6 for graphics

from PyQt6.QtCore import QSize
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *

#Make an image
def makeImage(label, src, left, top, width, height):
    im = QPixmap(src)
    im = im.scaled(width, height)
    label.setPixmap(im)
    label.move(left, top)
    return im

#Make a button
def makeButton(widget, text, left, top, width, height):
    button = QPushButton(widget)
    button.setText(text)
    button.setGeometry(left, top, width, height)
    return button
