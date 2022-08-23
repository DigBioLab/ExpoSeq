from tkinter import *
t = Tk()
t.geometry("600x600")
t.title("ExpoSeq")
fileOptions = ["New", "open", "Save", "Save as"]
fileOptionsAfterseparator = ["Import", "Export", "Exit"]
viewOptions = ["Transform", "Edit", "Create"]
menuBar = Menu(t)
file = Menu(menuBar, tearoff=0)
for i in fileOptions:
    file.add_command(label=i, command=None)
file.add_separator()
for i in fileOptionsAfterseparator:
    file.add_command(label=i, command=None)
menuBar.add_cascade(label="File", menu=file)
View = Menu(menuBar, tearoff=0)
for i in viewOptions:
    View.add_command(label=i, command=None)
menuBar.add_cascade(label="View", menu=View)
t.config(menu=menuBar)
t.mainloop()