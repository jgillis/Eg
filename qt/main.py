#!/bin/python

import sys
from PyQt4.QtCore import QLocale, QTranslator
from PyQt4.QtGui import QApplication
from mainwindow import MainWindow

__version__ = "1.0"

if __name__ == "__main__":

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
