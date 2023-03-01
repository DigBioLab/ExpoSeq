from PyQt5.QtWidgets import QApplication, QFileDialog
def get_file_path():
    app = QApplication([])

    # Create a file chooser dialog
    filename, _ = QFileDialog.getOpenFileName(
        None,
        "Select a file",
        "",
        "Text Files (*.txt);;All Files (*)"
    )

    # Close the QApplication object
    app.quit()

    # Return the selected filename
    return filename

def get_directory_path():
    app = QApplication([])

    # Create a directory chooser dialog
    directory_path = QFileDialog.getExistingDirectory(
        None,
        "Select a directory",
        ""
    )

    # Close the QApplication object
    app.quit()

    # Return the selected directory path
    return directory_path