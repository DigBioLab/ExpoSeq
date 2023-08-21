from tkinter import TclError
import tkinter.font as tkfont
import tkinter as tk
from tkinter import simpledialog, messagebox
from tkinter.colorchooser import askcolor
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logomaker
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
AMINO_ACID_COLORS = {
    "S": "green", "T": "green","Q": "green","N": "green",
    "A": "grey", "V": "grey", "L": "grey", "I": "grey", "F": "grey",  "W": "grey", "M": "grey", "Y": "grey",
    "C": "yellow","P": "yellow","G": "yellow",
    "D": "red", "E": "red",
    "H": "cyan","K": "blue","R": "blue",
    '-': 'black', '*': 'black','end': "black"
}
AMINO_ACID_COLORS = {
    "S": "#388E3C",  # Green 700
    "T": "#388E3C",  # Green 700
    "Q": "#7CB342",  # Light Green 600
    "N": "#7CB342",  # Light Green 600
    "A": "#616161",  # Grey 700
    "V": "#616161",  # Grey 700
    "L": "#424242",  # Grey 800
    "I": "#424242",  # Grey 800
    "F": "#212121",  # Grey 900
    "W": "#212121",  # Grey 900
    "M": "#455A64",  # Blue Grey 700
    "Y": "#455A64",  # Blue Grey 700
    "C": "#F57C00",  # Orange 700 (Instead of Yellow)
    "P": "#F57C00",  # Orange 700 (Instead of Yellow)
    "G": "#FFA000",  # Amber 700 (Brighter than Orange)
    "D": "#D32F2F",  # Red 700
    "E": "#C62828",  # Red 800
    "H": "#00838F",  # Cyan 700
    "K": "#1565C0",  # Blue 700
    "R": "#0D47A1",  # Blue 900
    '-': "#1B1B1B",  # Almost Black
    '*': "#1B1B1B",  # Almost Black
    'end': "#FFFFFF", # White,
    '.': "#1B1B1B",
    ':': "#1B1B1B"
}

def parse_clustal_format(filename):
    sequences = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Remove any initial non-sequence lines (e.g., metadata or blank lines)
        while lines and (lines[0].strip() == '' or lines[0].startswith('CLUSTAL') or lines[0].startswith(' ')):
            lines.pop(0)
        while lines:
            seq_block = {}
            while lines and lines[0].strip() != '':
                line = lines.pop(0).strip()
                if not line.startswith(" "):  # Avoiding lines starting with space which can be alignment indicators
                    parts = line.split()
                    seq_id = parts[0]
                    seq_segment = parts[1]
                    if seq_id not in seq_block:
                        seq_block[seq_id] = seq_segment
                    else:
                        seq_block[seq_id] += seq_segment
            # Update the master sequence dictionary
            for seq_id, seq_segment in seq_block.items():
                if seq_id not in sequences:
                    sequences[seq_id] = seq_segment
                else:
                    sequences[seq_id] += seq_segment
            # Skip any spacer lines between sequence blocks
            while lines and lines[0].strip() == '':
                lines.pop(0)
    return sequences



class SequenceManipulator:
    def __init__(self, sequence, spacing = 20):
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

        self.spacing = spacing  # consistent spacing
        for idx, char in enumerate(sequence):
            self.canvas.create_text(10 + idx * self.spacing, 10, anchor="nw", text=char, font=('Courier', 10),
                                    tags="sequence_char")

        # Add positions below the sequence
        positions = [str(i) if i % 5 == 0 else '' for i in range(1, len(sequence) + 1)]
        for idx, position in enumerate(positions):
            pos_width = (idx + 1) * self.spacing
            self.canvas.create_text(10 + idx * self.spacing, 30, anchor="nw", text=position, font=('Courier', 10))
        self.indicator_position = 0
        self.clipboard = ""
        # Update scroll region after items are added to canvas
        self.canvas.config(xscrollcommand=scrollbar.set)
        self.canvas.update_idletasks()  # ensure the canvas items are created
        self.canvas.config(scrollregion=self.canvas.bbox("all"))
        self.canvas.bind("<Button-1>", self.on_canvas_click)
        self.clipboard = ""

        self.canvas.bind("<Button-1>", self.focus_canvas)

        self.canvas.bind("<Left>", self.move_indicator_left)
        self.canvas.bind("<Right>", self.move_indicator_right)
        self.canvas.bind("<Key>", self.on_key_press)

        self.canvas.bind("<BackSpace>", self.on_backspace)
        self.canvas.bind('<Control-v>', self.on_paste)
        self.canvas.bind('<Control-c>', self.on_copy)

        # If you're on macOS, also add:
        self.canvas.bind("<Command-c>", self.on_copy)
        self.canvas.bind("<Command-v>", self.on_paste)

        self.canvas.focus_set()

        # Frame for buttons
        btn_frame = tk.Frame(self.root)
        btn_frame.pack(fill="both", pady=20, expand=True)

        btn_insert = tk.Button(btn_frame, text="Insert AA", command=self.insert_amino_acid)
        btn_insert.pack(side="left", fill="both", expand=True, padx=5, pady=10)

        btn_delete = tk.Button(btn_frame, text="Delete AA", command=self.delete_amino_acid)
        btn_delete.pack(side="left", fill="both", expand=True, padx=5, pady=10)

        btn_replace = tk.Button(btn_frame, text="Replace AA", command=self.replace_amino_acid)
        btn_replace.pack(side="left", fill="both", expand=True, padx=5, pady=10)

        btn_insert_region = tk.Button(btn_frame, text="Insert Region", command=self.insert_region)
        btn_insert_region.pack(side="left", fill="both", expand=True, padx=5, pady=10)

        btn_replace_region = tk.Button(btn_frame, text="Replace Region", command=self.replace_region)
        btn_replace_region.pack(side="left", fill="both", expand=True, padx=5, pady=10)

    def on_copy(self, event):
        idx = (self.indicator_position - 10) // self.spacing
        if 0 <= idx < len(self.sequence.get()):
            char_to_copy = self.sequence.get()[idx]
            self.clipboard_clear()
            self.clipboard_append(char_to_copy)
            print(f"Copied: {char_to_copy}")
        else:
            print("Nothing to copy")

    def on_paste(self, event):
        try:
            clipboard_text = self.text_widget.clipboard_get()
            # Now, you can process or use the clipboard_text as you wish
            # For example, if you want to insert it into the Text widget:
            self.text_widget.insert(tk.END, clipboard_text)
        except:
            print("error")

    #     except TclError:  # No data in clipboard
          #  print("Clipboard is empty")

    def focus_canvas(self, event):
        self.canvas.focus_set()
    def on_key_press(self, event):
        char = event.char.upper()
        if char in 'ACDEFGHIKLMNPQRSTVWY':  # This allows any character, modify as needed.
            idx = (self.indicator_position - 10) // self.spacing
            current_seq = self.sequence.get()
            updated_seq = current_seq[:idx] + char + current_seq[idx:]
            self.sequence.set(updated_seq)
            self.move_indicator_right(None)  # Move the caret to the right after insertion
            self.draw_sequence()

    def on_backspace(self, event):
        idx = (self.indicator_position - 10) // self.spacing
        if idx >= 0:  # This will allow the first character to be deleted
            current_seq = self.sequence.get()
            updated_seq = current_seq[:idx] + current_seq[idx + 1:]
            self.sequence.set(updated_seq)
            self.move_indicator_left(None)
            self.draw_sequence()

    def insert_sequence_at_position(self, idx):
        new_sequence = self.entry.get().upper()  # Get the content of the entry widget
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids
        if all(char in amino_acids for char in new_sequence):
            # Check if the entered characters are valid nucleotides
            current_seq = self.sequence.get()
            updated_seq = current_seq[:idx] + new_sequence + current_seq[idx:]
            self.sequence.set(updated_seq)
            self.draw_sequence()
            self.entry.destroy()  # Destroy the entry widget after inserting
        else:
            messagebox.showwarning("Invalid Input", "Please enter a valid nucleotide sequence (A, C, G, T).")

    def move_indicator_left(self, event):
        # Decrease the indicator position by the spacing only if it won't go beyond the starting position
        if self.indicator_position > 10:
            self.indicator_position -= self.spacing
            self.draw_position_indicator()

    def move_indicator_right(self, event):
        # Increase the indicator position by the spacing
        if self.indicator_position < (len(self.sequence.get()) * self.spacing):
            self.indicator_position += self.spacing
            self.draw_position_indicator()

    def draw_position_indicator(self):
        self.canvas.delete("position_indicator")
        self.canvas.create_line(self.indicator_position, 0, self.indicator_position, 50, tags="position_indicator",
                                fill="black")

    def on_key_press(self, event):
        char = event.char.upper()
        if char in 'ACDEFGHIKLMNPQRSTVWY':  # If it's an amino acid
            idx = (self.indicator_position - 10) // self.spacing
            current_seq = self.sequence.get()
            updated_seq = current_seq[:idx] + char + current_seq[idx:]
            self.sequence.set(updated_seq)
            self.move_indicator_right(None)
            self.draw_sequence()
        elif event.keysym == "BackSpace":
            self.on_backspace(event)

    def on_canvas_click(self, event):
        # Remove any existing vertical lines (with the tag "position_indicator")
        self.canvas.delete("position_indicator")

        # Calculate the closest position index based on the click
        idx = round((event.x - 10) / self.spacing)
        if idx >= len(self.sequence.get()):  # if clicked beyond the last character
            idx = len(self.sequence.get()) - 1
        self.indicator_position = 10 + idx * self.spacing

        # Draw the position indicator at the calculated position
        self.draw_position_indicator()

        # Calculate approximate sequence index from x position
        idx = (event.x - 10) // self.spacing  # Assuming 10 is the starting x-coordinate, adjust if different.

        # Position the Entry widget at the clicked position
        self.entry = tk.Entry(self.canvas_frame, width=5)
        self.entry.place(x=event.x, y=10)
        self.entry.focus_set()  # Set focus to Entry widget

        # Bind events to the entry
        self.entry.bind("<Return>", lambda e: self.insert_sequence_at_position(idx))
        self.entry.bind("<Escape>", lambda e: self.entry.destroy())  # Destroy entry widget if escape key pressed

    def draw_positions(self):
        self.canvas.delete("position_number")  # Clear previous positions
        positions = [str(i) if i % 5 == 0 else '' for i in range(1, len(self.sequence.get()) + 1)]
        for idx, position in enumerate(positions):
            self.canvas.create_text(10 + idx * self.spacing, 30, anchor="nw", text=position, font=('Courier', 10),
                                    tags="position_number")
    def draw_sequence(self):
        self.canvas.delete("all")
        self.draw_positions()
  # Clear previous sequence

        sequence = self.sequence.get()
        for idx, char in enumerate(sequence):
            self.canvas.create_text(10 + idx * self.spacing, 10, anchor="nw", text=char, font=('Courier', 10),
                                    tags="sequence_char")
        self.draw_position_indicator()

    def insert_amino_acid(self):
        position = simpledialog.askinteger("Input", "At which position?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get())+1)
        if position is None: return
        aa = simpledialog.askstring("Input", "Which amino acid to insert?", parent=self.root)
        if not aa: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:position-1] + aa + current_seq[position-1:]
        self.sequence.set(new_seq)
        self.draw_sequence()
    def delete_amino_acid(self):
        position = simpledialog.askinteger("Input", "At which position?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get()))
        if position is None: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:position-1] + current_seq[position:]
        self.sequence.set(new_seq)
        self.draw_sequence()

    def replace_amino_acid(self):
        position = simpledialog.askinteger("Input", "At which position?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get()))
        if position is None: return
        aa = simpledialog.askstring("Input", "Replace with which amino acid?", parent=self.root)
        if not aa: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:position-1] + aa + current_seq[position:]
        self.sequence.set(new_seq)
        self.draw_sequence()

    def insert_region(self):
        start_position = simpledialog.askinteger("Input", "Starting position?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get())+1)
        end_position = simpledialog.askinteger("Input", "Ending position?", parent=self.root, minvalue=start_position, maxvalue=len(self.sequence.get())+1)
        if start_position is None or end_position is None: return
        region = simpledialog.askstring("Input", "Which region to insert?", parent=self.root)
        if not region: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:start_position-1] + region + current_seq[end_position-1:]
        self.sequence.set(new_seq)
        self.draw_sequence()

    def replace_region(self):
        start_position = simpledialog.askinteger("Input", "Starting position of region to replace?", parent=self.root, minvalue=1, maxvalue=len(self.sequence.get())+1)
        end_position = simpledialog.askinteger("Input", "Ending position of region to replace?", parent=self.root, minvalue=start_position, maxvalue=len(self.sequence.get())+1)
        if start_position is None or end_position is None: return
        region = simpledialog.askstring("Input", "Replace with which region?", parent=self.root)
        if not region: return
        current_seq = self.sequence.get()
        new_seq = current_seq[:start_position-1] + region + current_seq[end_position:]
        self.sequence.set(new_seq)
        self.draw_sequence()

    def run(self):

        self.root.mainloop()

import re

def extract_number(seq_id):
    match = re.search(r'(\d+)', seq_id)
    return int(match.group(1)) if match else 999999  # Large default number


class RegionVisualizer:
    def __init__(self, parent, sequence_length, regions, char_width, x_offset, canvas_height=40):
        self.canvas_width = char_width * sequence_length
        self.canvas = tk.Canvas(parent, width=self.canvas_width, height=canvas_height)
        self.canvas.pack(fill="x", expand=True)
        self.canvas.config(scrollregion=(0, 0, self.canvas_width, canvas_height))

        for region in regions:
            start_pos = region['start']
            end_pos = region['end']
            color = region['color']
            label = region['label']

            x_start = x_offset + start_pos * char_width
            x_end = x_offset + end_pos * char_width

            self.canvas.create_line(x_start, canvas_height / 2, x_end, canvas_height / 2, fill=color, width=5)
            self.canvas.create_text((x_start + x_end) / 2, canvas_height / 4, text=label, anchor="center", font=('Courier', 8))

    def xview_moveto(self, fraction):
        self.canvas.xview_moveto(fraction)


class LimitedScrollText(tk.Text):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, wrap=tk.NONE, **kwargs)
        self.bind("<Configure>", self.calculate_scroll_width)
        self.max_x_scroll = kwargs.pop('max_x_scroll', 1.0)

    def calculate_scroll_width(self, event=None):
        # Determine the maximum width of content in the widget
        content = self.get(1.0, tk.END).splitlines()
        longest_line = max(content, key=len, default="")

        # Measure the pixel width of the longest line
        font_name = self.cget('font')
        font = tkfont.Font(font=font_name)
        max_pixel_width = font.measure(longest_line)

        # Limit the xscrollcommand to the measured width
        self.xscrollcommand = lambda f, l: self.limit_xview_scroll(f, l, max_pixel_width)

    def limited_xview(self, *args):
        if args[0] == 'moveto':
            fraction = float(args[1])
            if 0 <= fraction <= self.max_x_scroll:
                self.xview(*args)
        else:
            # The scroll "units" scenario
            shift = int(args[1])
            current_x = self.xview()[0]  # Current x fraction
            target_fraction = current_x + (shift * self.unit_fraction)
            if 0 <= target_fraction <= self.max_x_scroll:
                self.xview(*args)
    def limit_xview_scroll(self, first, last, max_pixel_width):
        if float(last) <= 1.0:
            # Allow normal scrolling until the end of content
            self.xview_moveto(float(first))
        else:
            # Calculate the over-scroll distance in terms of fraction
            excess_fraction = max_pixel_width * (float(last) - 1.0) / (self.winfo_width() * float(last))
            # Limit the maximum scrolling position
            self.xview_moveto(float(first) - excess_fraction)

    # Override the default xview method to restrict its behavior
    def xview(self, *args):
        if args[0] == "moveto":
            fraction = float(args[1])
            self.xview_moveto(fraction)
        elif args[0] == "scroll":
            self.xview_scroll(int(args[1]), args[2])


class ExtendedSequenceManipulator(SequenceManipulator):
    def plot_amino_acid_composition(self):
        # Calculate amino acid composition for all sequences
        all_seqs = list(self.sequences.values()) # Change here
        print(all_seqs)
        unique_aas = list(set(all_seqs))
        chosen_sequence_length = len(all_seqs[0])
        aminoacids = "ACDEFGHIKLMNPQRSTVWY"
        compDict = {aa: chosen_sequence_length * [0] for aa in aminoacids}
        sequences = all_seqs
        for seq in sequences:
            for aa_position in range(len(seq)):
                aminoacid = seq[aa_position]
                if aminoacid in ["*", "-", "", ".", ":"]:
                    pass
                else:
                    compDict[aminoacid][aa_position] += 1
        aa_distribution = pd.DataFrame.from_dict(compDict)
        # aa_distribution = aa_distribution.divide(aa_distribution.shape[1])
        aa_distribution = aa_distribution.divide(len(all_seqs), axis=0)



        # Create the bar plot
        # Create the bar plot
        fig_width = chosen_sequence_length * 0.05  # adjust 0.05 as needed
  # dynamically set width based on sequence length
        fig, ax = plt.subplots(figsize=(fig_width, 4))
        logo_plot = logomaker.Logo(aa_distribution,
                                   shade_below=.5,
                                   fade_below=.5,
                                   font_name='Arial Rounded MT Bold',
                                   color_scheme="skylign_protein",
                                   show_spines=False,
                                   ax=ax,
                                   figsize= (fig_width, 4)
                                   )
        logo_plot.style_xticks(anchor=1,
                               spacing=1,
                               rotation=0)
        labels_true = list(range(0, chosen_sequence_length))
        numbers_true = list(range(1, chosen_sequence_length + 1))
        step = 5  # display every 5th label
        plt.xticks(labels_true[::step], numbers_true[::step])

        plot_frame = tk.Frame(self.root)
        plot_frame.pack(side=tk.TOP, fill=tk.BOTH) #expand 1  # pack it so that it fills available space

        # Embed the plot in the plot frame
        canvas = FigureCanvasTkAgg(fig, master=plot_frame)

        canvas_widget = canvas.get_tk_widget()


        canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)


        canvas.draw()  # Force drawing
    def __init__(self, sequence, fasta_file, master=None):
        super().__init__(sequence)

        # Regions and visual parameters
        regions = [
            {'start': 10, 'end': 50, 'color': 'red', 'label': 'Region A'},
            {'start': 55, 'end': 80, 'color': 'blue', 'label': 'Region B'}
        ]
        x_offset = 10
        self.char_width = 20

        # RegionVisualizer
       # self.region_visualizer = RegionVisualizer(self.root, len(sequence), regions, char_width, x_offset)

        # Parsing and sorting sequences
        self.sequences = parse_clustal_format(fasta_file)
        sorted_seq_ids = sorted(self.sequences.keys(), key=extract_number)
        self.sequences = {seq_id: self.sequences[seq_id] for seq_id in sorted_seq_ids}

        self.max_length = max(map(len, self.sequences.values()))

        # Multi-Sequence Alignment Display
        msa_frame = tk.Frame(self.root, height=400)
        msa_frame.pack(fill="both", expand=True)

        self.msa_display = LimitedScrollText(msa_frame, font=('Courier', 10))


        vertical_scroll = tk.Scrollbar(msa_frame, command=self.msa_display.yview, orient="vertical")
        horizontal_scroll = tk.Scrollbar(msa_frame, command=self.msa_display.limited_xview, orient="horizontal")

        # Configurations
        self.msa_display.config(yscrollcommand=vertical_scroll.set)
        horizontal_scroll.config(command=self.msa_display.xview)
        self.msa_display.config(xscrollcommand=horizontal_scroll.set)


        # Pack components
        horizontal_scroll.pack(side="bottom", fill="x")
        vertical_scroll.pack(side="right", fill="y")
        self.msa_display.pack(side="left", fill="both", expand=True)

        # Populate the MSA display with sequences
        # Populate the MSA display with sequences
        for seq_id, seq in self.sequences.items():
            formatted_seq = ' '.join([char for char in seq]).ljust(self.char_width * self.max_length)
            self.msa_display.insert("end", f">{seq_id}\n", "header")

            for char in seq:
                tag_name = f"{AMINO_ACID_COLORS[char]}_tag"
                self.msa_display.insert("end", char, tag_name)
                self.msa_display.insert("end", " ")
            self.msa_display.insert("end", "\n")
            spacing = 2
            # Insert positions below the sequence
            positions = [str(i) if i % 5 == 0 else '' for i in range(1, len(seq) + 1)]
            for idx, position in enumerate(positions):
                self.msa_display.insert("end", position.ljust(spacing))  # ljust will ensure each position takes up 'spacing' spaces.

            self.msa_display.insert("end", "\n\n")
        # Apply tag configurations
        for color in AMINO_ACID_COLORS.values():
            tag_name = f"{color}_tag"
            self.msa_display.tag_configure(tag_name, foreground=color)
        self.msa_display.tag_configure("header", font=("Courier", 10, "bold"))

        self.change_color_button = tk.Button(self.root, text="Change Color", command=self.on_change_color_button_pressed)
        self.change_color_button.pack(pady=5)
        self.zoom_scale = tk.Scale(self.root, from_=5, to=20, orient=tk.HORIZONTAL, label="Zoom", command=self.on_zoom)
        self.zoom_scale.set(10)  # set initial value
        self.zoom_scale.pack(pady=5)
       # self.plot_amino_acid_composition()
    # You can use this method to sync other visual elements with the scroll of the MSA.
    # For now, it's just a placeholder that only interacts with the MSA's own scroll.
    def on_zoom(self, value):
        font_size = int(value)
        self.msa_display.config(font=('Courier', font_size))

    def sync_with_msa_scroll(self, *args):
        fraction = float(args[0])
        self.region_visualizer.xview_moveto(fraction)
        self.msa_display.xview_moveto(fraction)

    def refresh_display(self):
        """Refresh the sequences displayed in the Text widget."""

        # Clear the existing content
        self.msa_display.delete(1.0, tk.END)

        # Insert the updated sequences
        for seq_id, seq in self.sequences.items():
            formatted_seq = ' '.join([char for char in seq]).ljust(self.char_width * self.max_length)
            self.msa_display.insert("end", f">{seq_id}\n", "header")

            for char in seq:
                tag_name = f"{AMINO_ACID_COLORS[char]}_tag"
                self.msa_display.insert("end", char, tag_name)
                self.msa_display.insert("end", " ")
            spacing = 2
            # Insert positions below the sequence
            positions = [str(i) if i % 5 == 0 else '' for i in range(1, len(seq) + 1)]
            for idx, position in enumerate(positions):
                self.msa_display.insert("end", position.ljust(spacing))  # ljust will ensure each position takes up 'spacing' spaces.

            self.msa_display.insert("end", "\n\n")

        for color in AMINO_ACID_COLORS.values():
            tag_name = f"{color}_tag"
            self.msa_display.tag_configure(tag_name, foreground=color)
        self.msa_display.tag_configure("header", font=("Courier", 10, "bold"))
        #self.draw_sequence()

    def change_color(self, amino_acid):
        """
        Allow the user to select a color for the provided amino acid and update the display.
        """
        color = askcolor()[1]  # Ask the user to pick a color
        if color:  # Check if a color was selected
            AMINO_ACID_COLORS[amino_acid] = color
            tag_name = f"{color}_tag"
            self.msa_display.tag_configure(tag_name, foreground=color)
            # You might also want to redraw or refresh your sequence display after this
        self.refresh_display()

    def on_change_color_button_pressed(self):
        amino_acid = simpledialog.askstring("Input", "Enter the amino acid you want to change the color for:")
        if amino_acid and amino_acid in AMINO_ACID_COLORS:
            self.change_color(amino_acid)


