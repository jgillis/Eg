from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtNetwork import *
from ui_mainwindow import Ui_MainWindow

class MainWindow(QMainWindow, Ui_MainWindow):
  def __init__(self, parent = None):
    QMainWindow.__init__(self, parent)
    self.setupUi(self)

    # For a list of SIGNALs, see http://doc.qt.nokia.com/stable/qabstractbutton.html
    self.connect(self.pushButton, SIGNAL("clicked(bool)"),self.buttonclick)
    
  def buttonclick(self,pressed):
    print pressed
