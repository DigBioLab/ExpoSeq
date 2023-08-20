import tkinter as tk
from tkinter import simpledialog

class SequenceManipulator:
    def __init__(self, sequence):
        self.root = tk.Tk()
        self.root.title('Sequence Manipulator')

        # Set the root window's height and width
        self.root.geometry("800x600")

        # Create a frame for canvas and scrollbar
        frame = tk.Frame(self.root, height=150)  # Define the frame's height
        frame.pack(fill="x", padx=20, pady=20)

        self.canvas = tk.Canvas(frame, height=150)  # Define canvas height
        self.canvas.pack(fill="x", expand=True, side="top")

        scrollbar = tk.Scrollbar(frame, command=self.canvas.xview, orient="horizontal")
        scrollbar.pack(side="bottom", fill="x")

        # Create the canvas frame and attach to the canvas
        self.canvas_frame = tk.Frame(self.canvas)
        self.canvas.create_window((0, 0), window=self.canvas_frame, anchor="nw")

        self.sequence = tk.StringVar()
        self.sequence.set(sequence)

        spacing = 20  # consistent spacing
        for idx, char in enumerate(sequence):
            self.canvas.create_text(10 + idx * spacing, 10, anchor="nw", text=char, font=('Courier', 10))

        # Add positions below the sequence
        positions = [str(i) if i % 5 == 0 else '' for i in range(1, len(sequence) + 1)]
        for idx, position in enumerate(positions):
            pos_width = (idx + 1) * spacing
            self.canvas.create_text(10 + idx * spacing, 30, anchor="nw", text=position, font=('Courier', 10))

        # Update scroll region after items are added to canvas
        self.canvas.config(xscrollcommand=scrollbar.set)
        self.canvas.update_idletasks()  # ensure the canvas items are created
        self.canvas.config(scrollregion=self.canvas.bbox("all"))

        # Frame for buttons
        btn_frame = tk.Frame(self.root)
        btn_frame.pack(fill="both", pady=20, expand=True)

        btn_insert = tk.Button(btn_frame, text="Insert AA", command=self.insert_amino_acid)
        btn_insert.pack(side="left", fill="both", expand=True, padx=5, pady=10)

        btn_delete = tk.Button(btn_frame, text="Delete AA", command=self.delete_amino_acid)
        btn_delete.pack(side="left", fill="both", expand=True, padx=5, pady=10)

        btn_replace = tk.Button(btn_frame, text="Replace AA", command=self.replace_amino_acid)
        btn_replace.pack(side="left", fill="both", expand=True, padx=5, pady=10)

    def insert_amino_acid(self):
        position = simpledialog.askinteger("Input", "At which position?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get())+1)
        if position is None: return
        aa = simpledialog.askstring("Input", "Which amino acid to insert?", parent=self.root)
        if not aa: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:position-1] + aa + current_seq[position-1:]
        self.sequence.set(new_seq)

    def delete_amino_acid(self):
        position = simpledialog.askinteger("Input", "At which position?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get()))
        if position is None: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:position-1] + current_seq[position:]
        self.sequence.set(new_seq)

    def replace_amino_acid(self):
        position = simpledialog.askinteger("Input", "At which position?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get()))
        if position is None: return
        aa = simpledialog.askstring("Input", "Replace with which amino acid?", parent=self.root)
        if not aa: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:position-1] + aa + current_seq[position:]
        self.sequence.set(new_seq)

    def run(self):
        self.root.mainloop()

