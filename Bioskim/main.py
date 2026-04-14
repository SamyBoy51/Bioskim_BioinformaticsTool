import sys
from PyQt5 import QtWidgets
from BiosKimGUI import Ui_MainWindow
from analysis1 import *
def startUp():
    app = QtWidgets.QApplication(sys.argv)
    ventana = Ui_MainWindow()
    ventana.show()

    app.exec_()
##la ventana.ruta es mi path para poder aceder la varaible en los demás archivos.
    if ventana.ruta:
        print(f"El script main recibió la ruta: {ventana.ruta}")
        file = ventana.ruta

    else:
        print("No se seleccionó ningún archivo.")
    return file


if __name__=="__main__":
    startUp()
