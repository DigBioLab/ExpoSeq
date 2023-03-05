import sys
import os
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QGridLayout, QToolBar, QAction, QLabel
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MyApp(QWidget):

    def __init__(self):
        super().__init__()

        self.title = 'NGS Pipeline Application'
        self.left = 50
        self.top = 50
        self.width = 800
        self.height = 600

        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Create toolbar
        toolbar = QToolBar(self)
        self.addToolBar(toolbar)

        # Add file upload action
        fileUploadAction = QAction('Upload Data', self)
        fileUploadAction.triggered.connect(self.uploadData)
        toolbar.addAction(fileUploadAction)

        # Add settings action
        settingsAction = QAction('Settings', self)
        settingsAction.triggered.connect(self.openSettings)
        toolbar.addAction(settingsAction)

        # Create plot buttons
        self.plotButtons = {}
        plotButtonLayout = QGridLayout()
        plotButtonLayout.setSpacing(10)
        plotButtonLayout.setContentsMargins(20, 20, 20, 20)

        # Add a button for each plot
        for i, plot in enumerate(['Plot 1', 'Plot 2', 'Plot 3']):
            button = QPushButton(plot, self)
            if plot == 'Plot 1':
                button.clicked.connect(self.plot1)
            elif plot == 'Plot 2':
                button.clicked.connect(self.plot2)
            elif plot == 'Plot 3':
                button.clicked.connect(self.plot3)
            self.plotButtons[plot] = button
            plotButtonLayout.addWidget(button, 0, i)

        # Create figure canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)

        # Create layout for figure canvas and plot buttons
        layout = QGridLayout()
        layout.addWidget(self.canvas, 0, 0)
        layout.addLayout(plotButtonLayout, 1, 0)

        self.setLayout(layout)

        self.show()

    def uploadData(self):
        fileName, _ = QFileDialog.getOpenFileName(self, 'Upload Data')
        if fileName:
            print(f"Uploading file: {fileName}")
            # Call your data upload function here

    def openSettings(self):
        fileName, _ = QFileDialog.getOpenFileName(self, 'Open Settings File', '', 'Settings files (*.txt)')
        if fileName:
            print(f"Opening settings file: {fileName}")
            # Call your settings function here

    def plotData(self, plot):
        print(f"Plotting data for: {plot}")
        # Call your plotting function here
        # Update the figure canvas with the plot

    def updatePlot(self):
        # Update the figure canvas with the latest plot
        pass


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    sys.exit(app.exec_())
