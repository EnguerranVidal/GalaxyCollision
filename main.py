import os
import sys

import qdarktheme
from PyQt5.QtWidgets import QApplication

from src.gui.mainWindow import MainWindow


def main():
    import engine
    print(engine.__file__)
    print(hasattr(engine, "SimulatorParameters"))
    qdarktheme.enable_hi_dpi()
    app = QApplication(sys.argv)
    qdarktheme.setup_theme('dark', additional_qss='QToolTip {color: black;}')
    currentDirectory = os.path.dirname(os.path.realpath(__file__))
    mainWindow = MainWindow(currentDirectory)
    mainWindow.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
