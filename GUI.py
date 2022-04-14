#Import PyQt6 for graphics

from PyQt6.QtWidgets import QPushButton
from PyQt6.QtGui import QPixmap

#Make an image
def makeImage(label, src, left, top, width, height):
    '''Makes an image'''
    im = QPixmap(src)
    im = im.scaled(width, height)
    label.setPixmap(im)
    label.move(left, top)
    return im

#Make a button
def makeButton(widget, text, left, top, width, height):
    '''Makes a button'''
    button = QPushButton(widget)
    button.setText(text)
    button.setGeometry(left, top, width, height)
    return button

def setBGColor(widget, color):
    '''Sets the background color of a widget'''
    widget.setStyleSheet("""
        background-color: rgb(""" + color + """);
    """)
