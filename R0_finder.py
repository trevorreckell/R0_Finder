import sys
import os
import re
import tkinter as tk
from tkinter import ttk
from tkinter import simpledialog, messagebox, Toplevel, Checkbutton, IntVar, Label, Text, END, Frame, Scrollbar, \
    LEFT, RIGHT, Y, BOTH, GROOVE, filedialog, Entry, Menu
import math
import traceback
import json
import sympy as sp
from PIL import Image, EpsImagePlugin
# Matplotlib imports
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from R0_logic import VariableWeightPetriNet, ODESystem, StochasticPetriNet, decode_latex

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

LATEX_RENDER_MAP = {
    'alpha': 'α', 'beta': 'β', 'gamma': 'γ', 'delta': 'δ', 'epsilon': 'ε',
    'zeta': 'ζ', 'eta': 'η', 'theta': 'θ', 'iota': 'ι', 'kappa': 'κ',
    'lambda': 'λ', 'mu': 'μ', 'nu': 'ν', 'xi': 'ξ', 'omicron': 'ο',
    'pi': 'π', 'rho': 'ρ', 'sigma': 'σ', 'tau': 'τ', 'upsilon': 'υ',
    'phi': 'φ', 'chi': 'χ', 'psi': 'ψ', 'omega': 'ω',
    'Alpha': 'Α', 'Beta': 'Β', 'Gamma': 'Γ', 'Delta': 'Δ', 'Epsilon': 'Ε',
    'Zeta': 'Ζ', 'Eta': 'Η', 'Theta': 'Θ', 'Iota': 'Ι', 'Kappa': 'Κ',
    'Lambda': 'Λ', 'Mu': 'Μ', 'Nu': 'Ν', 'Xi': 'Ξ', 'Omicron': 'Ο',
    'Pi': 'Π', 'Rho': 'Ρ', 'Sigma': 'Σ', 'Tau': 'Τ', 'Upsilon': 'Υ',
    'Phi': 'Φ', 'Chi': 'Χ', 'Psi': 'Ψ', 'Omega': 'Ω'
}

def render_latex(text):
    """Converts \beta to β strictly for visual display on Canvas."""
    if not text: return ""
    # Regex to find \word
    def repl(match):
        key = match.group(1)
        return LATEX_RENDER_MAP.get(key, match.group(0)) # Return symbol or original \word
    return re.sub(r'\\([a-zA-Z]+)', repl, text)

#font settings
FONT_UI = ("Arial", 12)
FONT_BOLD = ("Arial", 12, "bold")
FONT_HEADER = ("Arial", 16, "bold")
FONT_TITLE = ("Arial", 36, "bold")
FONT_CANVAS = ("Arial", 14, "bold")
FONT_PLACE_NAME = ("Arial", 18, "bold")

WATERMARK_TEXT = "Reckell, Trevor, Beckett Sterner, and Petar Jevtić. \"The basic reproduction number for petri net models: A next-generation matrix approach.\" Applied Sciences 15.23 (2025): 12827."
FONT_WATERMARK = ("Arial", 9, "italic")

#   custom warning dialog for non int
class SPNWeightWarning(Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Warning: Non-Integer/Variable Weight")
        self.geometry("600x250")
        self.result = None

        Label(self,
              text="Current SPN modeling software cannot handle non-integer rational numbers\n or variables as arc weights. A Petri net with an arc with a non-integer rational\n number or a variable as the arc weight makes it a Continuous Petri net or\n Hybrid Petri net. Consider switching to a VAPN.",
              font=("Arial", 12), fg="red", pady=20).pack()

        btn_frame = Frame(self)
        btn_frame.pack(fill=tk.X, pady=10)

        tk.Button(btn_frame, text="Switch to VAPN builder", command=lambda: self.set_res("switch")).pack(side=LEFT,
                                                                                                         padx=10,
                                                                                                         expand=True)
        tk.Button(btn_frame, text="Edit Arc Weight", command=lambda: self.set_res("edit")).pack(side=LEFT, padx=10,
                                                                                                expand=True)
        tk.Button(btn_frame, text="Confirm Arc Weight", command=lambda: self.set_res("confirm")).pack(side=LEFT,
                                                                                                      padx=10,
                                                                                                      expand=True)

        self.transient(parent)
        self.grab_set()
        self.wait_window()

    def set_res(self, res):
        self.result = res
        self.destroy()

#warning dialog:ODe mapping
class ODEMappingWarning(Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Warning: Invalid SPN Weights Detected")
        self.geometry("650x250")
        self.result = None

        msg = ("Current SPN modeling software cannot handle non-integer rational numbers\n"
               "or variables as arc weights. A Petri net with an arc with a non-integer rational\n"
               "number or a variable as the arc weight makes it a Continuous Petri net or\n"
               "Hybrid Petri net. Consider switching to a VAPN.")

        Label(self, text=msg, font=("Arial", 11), fg="red", pady=20).pack()

        btn_frame = Frame(self)
        btn_frame.pack(fill=tk.X, pady=10)

        tk.Button(btn_frame, text="Edit ODE Model", command=lambda: self.set_res("edit")).pack(side=LEFT, padx=10,
                                                                                               expand=True)
        tk.Button(btn_frame, text="Switch to VAPN Model", command=lambda: self.set_res("vapn")).pack(side=LEFT, padx=10,
                                                                                                     expand=True)
        tk.Button(btn_frame, text="Continue Mapping to SPN Model", command=lambda: self.set_res("continue")).pack(
            side=LEFT, padx=10, expand=True)

        self.transient(parent)
        self.grab_set()
        self.wait_window()

    def set_res(self, res):
        self.result = res
        self.destroy()


#Tab 1: Start/home screen
class StartScreen(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.configure(bg="#f0f0f0")

        lbl_title = tk.Label(self, text="Select Model Type", font=FONT_TITLE, bg="#f0f0f0")
        lbl_title.pack(side=tk.TOP, pady=(40, 20), fill=tk.X)

        btn_frame = tk.Frame(self, bg="#f0f0f0")
        btn_frame.pack(side=tk.TOP, expand=True, fill=tk.BOTH, padx=50, pady=20)

        btn_frame.columnconfigure(0, weight=1)
        btn_frame.rowconfigure(0, weight=1)
        btn_frame.rowconfigure(1, weight=1)
        btn_frame.rowconfigure(2, weight=1)

        btn_font = ("Arial", 20, "bold")

        btn_vapn = tk.Button(btn_frame, text="VAPN\n(Variable Arc Weight Petri Net)", font=btn_font,
                             bg="#d0e0ff", activebackground="#b0c0ff",
                             command=lambda: controller.add_tab("VAPN"))
        btn_vapn.grid(row=0, column=0, sticky="nsew", pady=10)

        btn_spn = tk.Button(btn_frame, text="SPN\n(Stochastic Petri Net)", font=btn_font,
                            bg="#ffd0d0", activebackground="#ffb0b0",
                            command=lambda: controller.add_tab("SPN"))
        btn_spn.grid(row=1, column=0, sticky="nsew", pady=10)

        btn_ode = tk.Button(btn_frame, text="ODE\n(Ordinary Differential Equations)", font=btn_font,
                            bg="#d0ffd0", activebackground="#b0ffb0",
                            command=lambda: controller.add_tab("ODE"))
        btn_ode.grid(row=2, column=0, sticky="nsew", pady=10)

        # Watermark
        tk.Label(self, text=WATERMARK_TEXT, font=FONT_WATERMARK, fg="gray50", bg="#f0f0f0", wraplength=1000).pack(
            side=tk.BOTTOM, pady=10)


# PN editor tab
class PetriNetEditor(tk.Frame):
    def __init__(self, parent, controller, mode_type="VAPN"):
        super().__init__(parent)
        self.controller = controller
        self.mode_type = mode_type

        if mode_type == "SPN":
            self.pn_logic = StochasticPetriNet()
        else:
            self.pn_logic = VariableWeightPetriNet()

        self.mode = "SELECT"
        self.nodes_ui = []
        self.arcs_ui = []
        self.temp_arc_start = None
        self.spawn_x = 50
        self.spawn_y = 100
        self.drag_data = {"x": 0, "y": 0, "item": None}

        toolbar = tk.Frame(self, bg="#ddd", height=60)
        toolbar.pack(side=tk.TOP, fill=tk.X)

        def tb_btn(txt, cmd, bg=None, side=tk.LEFT):
            b = tk.Button(toolbar, text=txt, command=cmd, font=FONT_UI, bg=bg, height=2)
            b.pack(side=side, padx=5, pady=5)
            return b

        tb_btn("Close Tab", lambda: controller.close_tab(self), bg="#ffcccc")

        self.lbl_mode = tk.Label(toolbar, text=f" | Mode: {mode_type} | ", font=FONT_HEADER, bg="#ddd")
        self.lbl_mode.pack(side=tk.LEFT)

        tb_btn("Save", self.save_model)
        tb_btn("Load", self.load_model)
        tb_btn("Export Diagram", self.export_diagram, bg="#e6e6fa")
        tk.Label(toolbar, text="|", bg="#ddd", font=FONT_UI).pack(side=tk.LEFT)
        tb_btn("Place (O)", self.trigger_add_place)
        tb_btn("Transition (|)", self.trigger_add_trans)
        tb_btn("Arc (->)", lambda: self.set_mode("ARC"))
        tk.Label(toolbar, text="|", bg="#ddd", font=FONT_UI).pack(side=tk.LEFT)
        tb_btn("Select / Move", self.trigger_select_mode, bg="#e0e0e0")
        tb_btn("Clear", self.trigger_clear_all, bg="#ffdddd")
        tk.Label(toolbar, text="|", bg="#ddd", font=FONT_UI).pack(side=tk.LEFT)

        if mode_type == "SPN":
            tk.Button(toolbar, text="CALCULATE R₀", bg="#d0f0c0", font=("Arial", 14, "bold"), height=2,
                      command=self.open_calculation_window).pack(side=tk.RIGHT, padx=20, pady=5)
        else:
            tk.Button(toolbar, text="CALCULATE R₀", bg="#d0f0c0", font=("Arial", 14, "bold"), height=2,
                      command=self.open_calculation_window).pack(side=tk.RIGHT, padx=20, pady=5)
        # Watermark at bottom
        tk.Label(self, text=WATERMARK_TEXT, font=FONT_WATERMARK, fg="gray50", wraplength=1000).pack(side=tk.BOTTOM,
                                                                                                    fill=tk.X, pady=2)

        self.canvas = tk.Canvas(self, bg="white")
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.canvas.bind("<ButtonPress-1>", self.on_left_click_down)
        self.canvas.bind("<B1-Motion>", self.on_left_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_left_click_up)

        self.canvas.bind("<Double-Button-1>", self.on_double_click)
        self.canvas.bind("<Button-3>", self.on_right_click)
        self.canvas.bind("<Button-2>", self.on_right_click)
        self.canvas.bind("<Control-Button-1>", self.on_right_click)

    def export_diagram(self):
        #get the coordinates of all drawn items
        bbox = self.canvas.bbox("all")
        if not bbox:
            messagebox.showinfo("Info", "Canvas is empty, nothing to export.")
            return

        #add Padding
        pad = 20  # pixels
        x1, y1, x2, y2 = bbox

        x1 -= pad
        y1 -= pad
        x2 += pad
        y2 += pad

        #calculate width/height for the export
        w = x2 - x1
        h = y2 - y1

        #ask for filename
        file_types = [("PNG Image", "*.png"), ("EPS Vector", "*.eps")]
        filepath = filedialog.asksaveasfilename(defaultextension=".png", filetypes=file_types)
        if not filepath:
            return

        #generate PostScript ("master" vector source), temporarily save the PS data to memory or a file
        #'pageheight' and 'pagewidth' ensure the PS bounding box matches our crop
        try:
            #native Tkinter method to generate EPS data for a specific region
            #use 'colormode' to ensure color preservation
            ps_data = self.canvas.postscript(colormode='color', x=x1, y=y1, width=w, height=h, pagewidth=w,
                                             pageheight=h)

            if filepath.lower().endswith(".eps"):
                # Save EPS directly
                with open(filepath, "w") as f:
                    f.write(ps_data)
                messagebox.showinfo("Success", f"Diagram exported to {filepath}")

            else:
                #convert to PNG using Pillow
                #note: This requires Ghostscript to be installed on the system for high-quality rendering.
                # if Ghostscript is missing, this might fail or fallback to low-res.
                try:
                    # have to save PS to a temp file for Pillow to read it reliably on some OSs
                    temp_ps = "temp_export.eps"
                    with open(temp_ps, "w") as f:
                        f.write(ps_data)

                    with Image.open(temp_ps) as img:
                        img.load()  # Force load
                        # Save as PNG
                        img.save(filepath, "png")

                    if os.path.exists(temp_ps):
                        os.remove(temp_ps)

                    messagebox.showinfo("Success", f"Diagram exported to {filepath}")

                except Exception as e:
                    messagebox.showerror("PNG Error",
                                         f"Could not convert to PNG.\n(Do you have Ghostscript installed?)\n\nError: {e}\n\nTry saving as .eps instead.")

        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export diagram: {e}")

    def set_mode(self, m):
        self.mode = m

    def trigger_select_mode(self):
        self.set_mode("SELECT")

        # Check if there is an active highlighting (red node from Arc mode)
        if self.temp_arc_start:
            node = self.temp_arc_start
            # Restore original look based on type
            if node['type'] == 'place':
                self.canvas.itemconfig(node['id'], outline="black", width=2)
            else:
                self.canvas.itemconfig(node['id'], outline="black", width=1)

            # Clear the reference
            self.temp_arc_start = None

        # Also clear any drag data to be safe
        self.drag_data = {"x": 0, "y": 0, "item": None}

    def get_spawn_pos(self):
        self.canvas.update_idletasks()
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        if w < 100: w = 800
        if h < 100: h = 600

        x, y = self.spawn_x, self.spawn_y
        self.spawn_x += 80

        if self.spawn_x > w - 50:
            self.spawn_x = 50
            self.spawn_y += 80

        if self.spawn_y > h - 50:
            self.spawn_y = 50

        return x, y

    def trigger_add_place(self):
        name = simpledialog.askstring("Input", "Place Name:")
        if name:
            for n in name.split(','):
                n = n.strip()
                if n: self.add_place_ui(*self.get_spawn_pos(), n); self.pn_logic.add_place(n)

    def trigger_add_trans(self):
        if self.mode_type == "SPN":
            name = simpledialog.askstring("Input", "Transition Name:")
            if name:
                rate = simpledialog.askstring("Rate", f"Propensity Function (Firing Rate) for {name}:",
                                              initialvalue="1.0")
                if rate is None: rate = "1.0"
                self.add_trans_ui(*self.get_spawn_pos(), name, rate=rate)
                self.pn_logic.add_transition(name, rate)
        else:
            name = simpledialog.askstring("Input", "Transition Name:")
            if name:
                for n in name.split(','):
                    n = n.strip()
                    if n: self.add_trans_ui(*self.get_spawn_pos(), n); self.pn_logic.add_transition(n)

    def add_place_ui(self, x, y, name):
        r = 25
        uid = self.canvas.create_oval(x - r, y - r, x + r, y + r, fill="white", outline="black", width=2)
        font_use = ("Arial", 22, "bold") if self.mode_type == "SPN" else ("Arial", 18, "bold")
        lid = self.canvas.create_text(x, y, text=render_latex(name), font=font_use)
        self.nodes_ui.append(
            {'id': uid, 'label_id': lid, 'type': 'place', 'name': name, 'x': x, 'y': y, 'w': r, 'h': r})

    def add_trans_ui(self, x, y, name, rate="1.0"):
        if self.mode_type == "SPN":
            # Dynamic width based on name and rate length
            w_name = len(name) * 7  # Increased multiplier for larger font
            w_rate = len(rate) * 6
            w = max(40, w_name, w_rate) + 10
            h = 25

            uid = self.canvas.create_rectangle(x - w, y - h, x + w, y + h, fill="black", outline="black")

            lid = self.canvas.create_text(x, y - 10, text=render_latex(name), font=("Arial", 16, "bold"), fill="white")
            rid = self.canvas.create_text(x, y + 10, text=render_latex(rate), font=("Arial", 12, "italic"), fill="white")

            self.nodes_ui.append({
                'id': uid, 'label_id': lid, 'rate_id': rid,
                'type': 'trans', 'name': name, 'rate': rate,
                'x': x, 'y': y, 'w': w, 'h': h
            })
        else:
            w, h = 30, 20
            uid = self.canvas.create_rectangle(x - w, y - h, x + w, y + h, fill="black")

            lid = self.canvas.create_text(x, y, text=render_latex(name), font=("Arial", 16, "bold"), fill="white")

            self.nodes_ui.append({
                'id': uid, 'label_id': lid, 'type': 'trans', 'name': name, 'x': x, 'y': y, 'w': w, 'h': h})

    def on_left_click_down(self, event):
        x, y = event.x, event.y
        clicked = self.find_node_at(x, y)
        if self.mode == "SELECT" and clicked:
            self.drag_data = {"x": x, "y": y, "item": clicked}
        elif self.mode == "ARC" and clicked:
            if not self.temp_arc_start:
                self.temp_arc_start = clicked
                self.canvas.itemconfig(clicked['id'], outline="red", width=3)
            else:
                self.finish_arc(self.temp_arc_start, clicked)
                self.canvas.itemconfig(self.temp_arc_start['id'], outline="black", width=2)
                self.temp_arc_start = None

    def on_right_click(self, event):
        clicked_node = self.find_node_at(event.x, event.y)
        if clicked_node:
            self.show_node_context_menu(event, clicked_node)
            return

        clicked_arc = self.find_arc_at(event.x, event.y)
        if clicked_arc:
            self.show_arc_context_menu(event, clicked_arc)

    def find_arc_at(self, x, y):
        items = self.canvas.find_overlapping(x - 5, y - 5, x + 5, y + 5)
        for item_id in items:
            for arc in self.arcs_ui:
                if arc['line_id'] == item_id:
                    return arc
        return None

    def show_node_context_menu(self, event, node):
        menu = Menu(self, tearoff=0)
        menu.add_command(label=f"Rename '{node['name']}'", command=lambda: self.edit_node(node, "name"))
        if self.mode_type == "SPN" and node['type'] == 'trans':
            menu.add_command(label="Edit Propensity Function", command=lambda: self.edit_node(node, "rate"))
        menu.add_separator()
        menu.add_command(label="Delete", command=lambda: self.delete_node(node), foreground="red")
        menu.post(event.x_root, event.y_root)

    def show_arc_context_menu(self, event, arc):
        menu = Menu(self, tearoff=0)
        menu.add_command(label="Edit Weight", command=lambda: self.edit_arc(arc))
        menu.add_separator()
        menu.add_command(label="Delete Arc", command=lambda: self.delete_arc(arc), foreground="red")
        menu.post(event.x_root, event.y_root)

    def on_double_click(self, event):
        clicked = self.find_node_at(event.x, event.y)
        if clicked:
            self.edit_node(clicked)
        else:
            arc = self.find_arc_at(event.x, event.y)
            if arc: self.edit_arc(arc)

    def edit_node(self, node, target="all"):
        if node['type'] == 'place':
            new_name = simpledialog.askstring("Edit", "Place Name:", initialvalue=node['name'])
            if new_name and new_name != node['name']:
                self.pn_logic.remove_place(node['name'])
                self.pn_logic.add_place(new_name)
                node['name'] = new_name
                self.canvas.itemconfig(node['label_id'], text=render_latex(new_name))

        elif node['type'] == 'trans':
            if target == "name" or target == "all":
                new_name = simpledialog.askstring("Edit", "Transition Name:", initialvalue=node['name'])
                if new_name and new_name != node['name']:
                    node['name'] = new_name
                    self.canvas.itemconfig(node['label_id'], text=render_latex(new_name))
                    if self.mode_type == "SPN": self.resize_trans_node(node)

            if self.mode_type == "SPN" and (target == "rate" or target == "all"):
                new_rate = simpledialog.askstring("Edit", f"Propensity Function (Firing Rate) for {node['name']}:",
                                                  initialvalue=node.get('rate', '1.0'))
                if new_rate:
                    node['rate'] = new_rate
                    if 'rate_id' in node:
                        self.canvas.itemconfig(node['rate_id'], text=render_latex(new_rate))
                        self.resize_trans_node(node)
                    self.pn_logic.update_transition_rate(node['name'], new_rate)

    def resize_trans_node(self, node):
        if node['type'] != 'trans': return

        name_len = len(node['name'])
        rate_len = len(node.get('rate', ''))

        w_name = name_len * 5
        w_rate = rate_len * 4
        new_w = max(40, w_name, w_rate) + 10

        node['w'] = new_w
        x, y = node['x'], node['y']
        h = node['h']

        self.canvas.coords(node['id'], x - new_w, y - h, x + new_w, y + h)
        self.redraw_arcs(node)

    def edit_arc(self, arc):
        w = simpledialog.askstring("Weight", "Weight:", initialvalue=arc['weight_str'])
        if w:
            if self.mode_type == "SPN" and not w.isdigit():
                warn = SPNWeightWarning(self)
                if warn.result == "switch":
                    self.controller.add_tab("VAPN")
                    return
                elif warn.result == "edit":
                    self.edit_arc(arc)
                    return
                # if confirm, fall through to update

            self.pn_logic.remove_arc(arc['start_node']['name'], arc['end_node']['name'])
            self.pn_logic.add_arc(arc['start_node']['name'], arc['end_node']['name'], w)
            arc['weight_str'] = w
            txt = w
            if self.mode_type == "SPN" and txt == "1": txt = ""
            self.canvas.itemconfig(arc['text_id'], text=render_latex(txt))

    def delete_node(self, node):
        if messagebox.askyesno("Delete", f"Delete {node['name']}?"):
            self.nodes_ui.remove(node)
            self.canvas.delete(node['id'])
            self.canvas.delete(node['label_id'])
            if 'rate_id' in node: self.canvas.delete(node['rate_id'])

            if node['type'] == 'place':
                self.pn_logic.remove_place(node['name'])
            else:
                self.pn_logic.remove_transition(node['name'])

            to_remove = [a for a in self.arcs_ui if a['start_node'] == node or a['end_node'] == node]
            for a in to_remove:
                self.delete_arc(a, logic_already_cleared=True)

    def delete_arc(self, arc, logic_already_cleared=False):
        if arc in self.arcs_ui:
            self.arcs_ui.remove(arc)
        self.canvas.delete(arc['line_id'])
        self.canvas.delete(arc['text_id'])
        if not logic_already_cleared:
            self.pn_logic.remove_arc(arc['start_node']['name'], arc['end_node']['name'])

    def finish_arc(self, start, end):
        if start == end:
            return

        w = simpledialog.askstring("Weight", f"{start['name']}->{end['name']} weight:", initialvalue="1")
        if w:
            if self.mode_type == "SPN" and not w.isdigit():
                warn = SPNWeightWarning(self)
                if warn.result == "switch":
                    self.controller.add_tab("VAPN")
                    return
                elif warn.result == "edit":
                    self.finish_arc(start, end)
                    return
                #if confirm, fall through

            try:
                self.pn_logic.add_arc(start['name'], end['name'], w)
                self.add_arc_ui(start, end, w)
            except Exception as e:
                messagebox.showerror("Error", str(e))

    def add_arc_ui(self, start, end, w_str):
        line_id = self.canvas.create_line(0, 0, 0, 0, arrow=tk.LAST, width=2, smooth=True)
        text_val = w_str
        if self.mode_type == "SPN" and w_str == "1": text_val = ""

        text_id = self.canvas.create_text(0, 0, text=render_latex(text_val), fill="blue", font=("Arial", 16, "bold"))

        arc = {'line_id': line_id, 'text_id': text_id, 'start_node': start, 'end_node': end, 'weight_str': w_str}
        self.arcs_ui.append(arc)
        self.update_arc_coords(arc)

    def update_arc_coords(self, arc):
        start = arc['start_node']
        end = arc['end_node']

        sx, sy = self.get_edge_coords(start, end['x'], end['y'])
        ex, ey = self.get_edge_coords(end, start['x'], start['y'])

        mx, my = (sx + ex) / 2, (sy + ey) / 2
        dx, dy = ex - sx, ey - sy
        dist = math.hypot(dx, dy)

        if dist == 0: return

        nx, ny = -dy / dist, dx / dist
        offset = 30

        cx, cy = mx + nx * offset, my + ny * offset

        self.canvas.coords(arc['line_id'], sx, sy, cx, cy, ex, ey)
        self.canvas.coords(arc['text_id'], cx, cy - 10)

    def get_edge_coords(self, node_ui, target_x, target_y):
        cx, cy = node_ui['x'], node_ui['y']
        dx = target_x - cx;
        dy = target_y - cy
        if dx == 0 and dy == 0: return cx, cy
        if node_ui['type'] == 'place':
            angle = math.atan2(dy, dx);
            r = node_ui['w']
            return cx + r * math.cos(angle), cy + r * math.sin(angle)
        elif node_ui['type'] == 'trans':
            w = node_ui['w'];
            h = node_ui['h']
            slope_safe_dx = dx if dx != 0 else 0.0001
            slope_safe_dy = dy if dy != 0 else 0.0001
            tx = (w if dx > 0 else -w) / slope_safe_dx;
            ty = (h if dy > 0 else -h) / slope_safe_dy
            if abs(tx) <= abs(ty):
                return cx + tx * slope_safe_dx, cy + tx * slope_safe_dy
            else:
                return cx + ty * slope_safe_dx, cy + ty * slope_safe_dy

    def find_node_at(self, x, y):
        for n in self.nodes_ui:
            if math.hypot(x - n['x'], y - n['y']) < 40: return n
        return None

    def on_left_drag(self, event):
        if self.mode == "SELECT" and self.drag_data["item"]:
            n = self.drag_data["item"]

            # boundary constraint logic
            w = self.canvas.winfo_width()
            h = self.canvas.winfo_height()

            padding = 40
            new_x = max(padding, min(w - padding, event.x))
            new_y = max(padding, min(h - padding, event.y))

            dx = new_x - n['x']
            dy = new_y - n['y']

            self.canvas.move(n['id'], dx, dy)
            self.canvas.move(n['label_id'], dx, dy)
            if 'rate_id' in n:
                self.canvas.move(n['rate_id'], dx, dy)

            n['x'] = new_x
            n['y'] = new_y
            self.drag_data['x'] = new_x
            self.drag_data['y'] = new_y
            self.redraw_arcs(n)

    def on_left_click_up(self, e):
        self.drag_data['item'] = None

    def redraw_arcs(self, node):
        for a in self.arcs_ui:
            if a['start_node'] == node or a['end_node'] == node:
                self.update_arc_coords(a)

    def trigger_clear_all(self):
        if not self.nodes_ui: return
        if messagebox.askyesno("Clear", "Delete everything?"):
            self.canvas.delete("all")
            self.nodes_ui = [];
            self.arcs_ui = [];
            if self.mode_type == "SPN":
                self.pn_logic = StochasticPetriNet()
            else:
                self.pn_logic = VariableWeightPetriNet()

    def save_model(self):
        if not self.nodes_ui: messagebox.showwarning("Warning", "Nothing to save!"); return
        filepath = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON Files", "*.json")])
        if not filepath: return
        data = {"places": [], "transitions": [], "arcs": []}
        for n in self.nodes_ui:
            d = {"name": n["name"], "type": n["type"], "x": n["x"], "y": n["y"]}
            if self.mode_type == "SPN" and n["type"] == "trans":
                d["rate"] = n.get("rate", "1.0")

            if n["type"] == "place":
                data["places"].append(d)
            else:
                data["transitions"].append(d)
        for a in self.arcs_ui:
            data["arcs"].append(
                {"source": a["start_node"]["name"], "target": a["end_node"]["name"], "weight": a["weight_str"]})
        try:
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=4)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def load_model(self):
        filepath = filedialog.askopenfilename(filetypes=[("JSON Files", "*.json")])
        if not filepath: return
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            self.load_from_data(data)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def load_from_data(self, data):
        self.canvas.delete("all")
        self.nodes_ui = []
        self.arcs_ui = []

        # check for spn mapping
        # if in SPN mode, check if incoming data has non-integer weights like "e+1"
        if self.mode_type == "SPN":
            has_non_int = False
            for a in data.get("arcs", []):
                w = str(a.get("weight", "1"))
                if not w.isdigit():
                    has_non_int = True
                    break

            if has_non_int:
                warn = ODEMappingWarning(self)
                if warn.result == "edit":
                    # User wants to fix ODE so close this invalid SPN tab immediately
                    self.controller.close_tab(self)
                    return
                elif warn.result == "vapn":
                    # User wants to switch to VAPN so update mode and continue loading
                    self.mode_type = "VAPN"
                    self.lbl_mode.config(text=" | Mode: VAPN | ")
                # if "continue", do nothing and let it load as broken SPN (per users choice)

        # Initialize the correct logic engine based on (potentially updated) mode
        if self.mode_type == "SPN":
            self.pn_logic = StochasticPetriNet()
        else:
            self.pn_logic = VariableWeightPetriNet()

        name_map = {}

        for p in data.get("places", []):
            x = p.get("x", 0)
            y = p.get("y", 0)
            self.add_place_ui(x, y, p["name"])
            self.pn_logic.add_place(p["name"])
            name_map[p["name"]] = self.nodes_ui[-1]

        for t in data.get("transitions", []):
            x = t.get("x", 0)
            y = t.get("y", 0)
            if self.mode_type == "SPN":
                rate = t.get("rate", "1.0")
                self.add_trans_ui(x, y, t["name"], rate)
                self.pn_logic.add_transition(t["name"], rate)
            else:
                self.add_trans_ui(x, y, t["name"])
                self.pn_logic.add_transition(t["name"])
            name_map[t["name"]] = self.nodes_ui[-1]

        for a in data.get("arcs", []):
            if a["source"] in name_map and a["target"] in name_map:
                w_str = str(a.get("weight", "1"))
                self.pn_logic.add_arc(a["source"], a["target"], w_str)
                self.add_arc_ui(name_map[a["source"]], name_map[a["target"]], w_str)

        if any(n['x'] == 0 for n in self.nodes_ui):
            self.apply_circle_layout()

    def apply_circle_layout(self):
        n = len(self.nodes_ui)
        if n == 0: return
        self.canvas.update_idletasks()
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        if w < 100: w = 1000
        if h < 100: h = 800
        cx, cy = w / 2, h / 2
        padding = 80
        max_r_x = (w / 2) - padding
        max_r_y = (h / 2) - padding
        radius = min(max_r_x, max_r_y)
        if radius < 50: radius = 50
        angle_step = 2 * math.pi / n

        for i, node in enumerate(self.nodes_ui):
            angle = i * angle_step
            node['x'] = cx + radius * math.cos(angle)
            node['y'] = cy + radius * math.sin(angle)
            if node['type'] == 'place':
                r = node['w']
                self.canvas.coords(node['id'], node['x'] - r, node['y'] - r, node['x'] + r, node['y'] + r)
            else:
                w_node, h_node = node['w'], node['h']
                self.canvas.coords(node['id'], node['x'] - w_node, node['y'] - h_node, node['x'] + w_node,
                                   node['y'] + h_node)
            self.canvas.coords(node['label_id'], node['x'], node['y'])
            if 'rate_id' in node:
                self.canvas.coords(node['rate_id'], node['x'], node['y'] + 10)
                self.canvas.coords(node['label_id'], node['x'], node['y'] - 10)
            self.redraw_arcs(node)

    def open_calculation_window(self):
        CalcWindow(self, self.pn_logic, "PETRI")


# ODE editor tab
class ODEEditor(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.ode_logic = ODESystem()

        # UI Toolbar
        toolbar = tk.Frame(self, bg="#ddd", height=60)
        toolbar.pack(side=tk.TOP, fill=tk.X)

        tk.Button(toolbar, text="Close Tab", command=lambda: controller.close_tab(self), font=FONT_UI, height=2,
                  bg="#ffcccc").pack(side=tk.LEFT, padx=5, pady=5)

        tk.Label(toolbar, text=" | Mode: ODE Editor | ", font=FONT_HEADER, bg="#ddd").pack(side=tk.LEFT)

        tk.Button(toolbar, text="Save", command=self.save_model, font=FONT_UI, height=2).pack(side=tk.LEFT, padx=5)
        tk.Button(toolbar, text="Load", command=self.load_model, font=FONT_UI, height=2).pack(side=tk.LEFT, padx=5)

        tk.Button(toolbar, text="CALCULATE R₀", bg="#d0f0c0", font=("Arial", 14, "bold"), height=2,
                  command=self.open_calc).pack(side=tk.RIGHT, padx=20, pady=5)

        tk.Button(toolbar, text="Map to VAPN", command=lambda: self.map_to_petri("VAPN"), font=FONT_UI, height=2,
                  bg="#d0e0ff").pack(side=tk.RIGHT, padx=5, pady=5)

        tk.Button(toolbar, text="Map to SPN", command=lambda: self.map_to_petri("SPN"), font=FONT_UI, height=2,
                  bg="#ffd0d0").pack(side=tk.RIGHT, padx=5, pady=5)

        # Watermark
        tk.Label(self, text=WATERMARK_TEXT, font=FONT_WATERMARK, fg="gray50", bg="#f0f0f0", wraplength=1000).pack(
            side=tk.BOTTOM, pady=5)


        sys_frame = tk.LabelFrame(self, text="Full System (Auto-Updates)", font=FONT_BOLD, bg="#f9f9f9")
        sys_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=False, padx=10, pady=5)

        self.txt_system = Text(sys_frame, height=8, font=("Courier", 12), bg="#fffff0")  # Monospace for code
        self.txt_system.pack(side=LEFT, fill=BOTH, expand=True, padx=5, pady=5)

        sb_sys = Scrollbar(sys_frame, command=self.txt_system.yview)
        sb_sys.pack(side=RIGHT, fill=Y)
        self.txt_system.config(yscrollcommand=sb_sys.set)

        instruct_frame = tk.Frame(self, bg="#fff8e1", bd=1, relief="solid")
        instruct_frame.pack(side="top", fill="x", padx=10, pady=(10, 0))

        lbl_tip = tk.Label(instruct_frame, justify="left", bg="#fff8e1", text=(
            "⚠️  MAPPING TIP: Separate terms = Separate Arrows\n"
            "• To create distinct transitions (e.g. Death vs Migration), write them separately:  alpha*I + beta*I\n"
            "• To create a single transition with a complex rate, group them:  (1-p)*gamma*I"
        ), font=("Arial", 10), fg="#e65100", padx=10, pady=5)
        lbl_tip.pack(side="left")

        # Scrollable Area for Entries (Top)
        self.canvas_area = tk.Canvas(self, borderwidth=0, background="#ffffff")
        self.scrollbar = Scrollbar(self, orient="vertical", command=self.canvas_area.yview)
        self.scroll_frame = tk.Frame(self.canvas_area, background="#ffffff")

        self.scroll_frame.bind("<Configure>",
                               lambda e: self.canvas_area.configure(scrollregion=self.canvas_area.bbox("all")))
        self.canvas_area.create_window((0, 0), window=self.scroll_frame, anchor="nw")
        self.canvas_area.configure(yscrollcommand=self.scrollbar.set)

        self.canvas_area.pack(side="top", fill="both", expand=True, padx=10, pady=10)
        self.scrollbar.pack(side="right", fill="y")

        # Headers
        tk.Label(self.scroll_frame, text="Please Select\nHuman Compartments", font=FONT_BOLD, bg="white").grid(row=0,
                                                                                                               column=0,
                                                                                                               padx=10,
                                                                                                               pady=5)
        tk.Label(self.scroll_frame, text="Variable", font=FONT_BOLD, bg="white").grid(row=0, column=2, padx=10, pady=5)
        tk.Label(self.scroll_frame, text="Equation (dx/dt = ...)", font=FONT_BOLD, bg="white").grid(row=0, column=4,
                                                                                                    padx=10, pady=5,
                                                                                                    sticky="w")

        self.row_frames = []
        self.btn_add_row = tk.Button(self.scroll_frame, text="+ Add Variable", command=self.add_variable_row,
                                     font=FONT_UI, bg="#e0e0e0")

        self.add_variable_row()

    def add_variable_row(self):
        r = len(self.row_frames) + 1
        self.btn_add_row.grid_forget()

        var_human = IntVar()
        chk_hum = Checkbutton(self.scroll_frame, variable=var_human, bg="white")
        chk_hum.grid(row=r, column=0, padx=10, pady=5)

        lbl_d = tk.Label(self.scroll_frame, text="d", font=FONT_UI, bg="white")
        lbl_d.grid(row=r, column=1, sticky="e")

        ent_var = Entry(self.scroll_frame, width=8, font=FONT_UI, justify="center")
        ent_var.grid(row=r, column=2, padx=2)

        # binding for live update
        ent_var.bind("<KeyRelease>", lambda e: self.update_system_text())

        lbl_dt = tk.Label(self.scroll_frame, text="/ dt  =", font=FONT_UI, bg="white")
        lbl_dt.grid(row=r, column=3, sticky="w")

        ent_eq = Entry(self.scroll_frame, width=50, font=FONT_UI)
        ent_eq.grid(row=r, column=4, padx=5, sticky="ew")

        ent_eq.bind("<KeyRelease>", lambda e: self.update_system_text())

        widgets = {
            "chk": chk_hum, "var": var_human,
            "d": lbl_d, "ent_var": ent_var, "dt": lbl_dt, "ent_eq": ent_eq
        }

        def delete_this_row():
            for w in [widgets["chk"], widgets["d"], widgets["ent_var"], widgets["dt"], widgets["ent_eq"],
                      widgets["btn"]]:
                w.destroy()
            if widgets in self.row_frames:
                self.row_frames.remove(widgets)
            self.update_system_text()  # Update on delete

        btn_del = tk.Button(self.scroll_frame, text="X", font=("Arial", 10, "bold"), fg="red", command=delete_this_row)
        btn_del.grid(row=r, column=5, padx=10)
        widgets["btn"] = btn_del

        self.row_frames.append(widgets)
        self.btn_add_row.grid(row=r + 1, column=0, columnspan=5, pady=15, sticky="w", padx=20)

    def update_system_text(self):
        """Generates the full system text for display."""
        lines = []
        for row in self.row_frames:
            v = row["ent_var"].get().strip()
            eq = row["ent_eq"].get().strip()
            if v:
                # Format: dX/dt = equation
                lines.append(f"d{v}/dt = {eq}")

        full_text = "\n".join(lines)

        self.txt_system.config(state="normal")
        self.txt_system.delete("1.0", END)
        self.txt_system.insert("1.0", full_text)
        # self.txt_system.config(state="disabled") # Keep normal so they can copy/paste easily

    def save_model(self):
        filepath = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON Files", "*.json")])
        if not filepath: return

        data = []
        for row in self.row_frames:
            # We save the state of each row
            item = {
                "var": row["ent_var"].get(),
                "eq": row["ent_eq"].get(),
                "human": row["var"].get()
            }
            data.append(item)

        try:
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=4)
            messagebox.showinfo("Saved", f"Model saved to {filepath}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))

    def load_model(self):
        filepath = filedialog.askopenfilename(filetypes=[("JSON Files", "*.json")])
        if not filepath: return

        try:
            with open(filepath, 'r') as f:
                data = json.load(f)

            # Clear existing UI rows
            for row in self.row_frames:
                for w in [row["chk"], row["d"], row["ent_var"], row["dt"], row["ent_eq"], row["btn"]]:
                    w.destroy()
            self.row_frames = []

            # Rebuild UI
            for item in data:
                self.add_variable_row()
                row = self.row_frames[-1]
                row["ent_var"].insert(0, item.get("var", ""))
                row["ent_eq"].insert(0, item.get("eq", ""))
                row["var"].set(item.get("human", 0))

            self.update_system_text()
            messagebox.showinfo("Loaded", f"Model loaded from {filepath}")

        except Exception as e:
            messagebox.showerror("Load Error", str(e))

    def sync_logic(self):
        self.ode_logic = ODESystem()
        for row in self.row_frames:
            var_name = row["ent_var"].get().strip()
            eq_str = row["ent_eq"].get().strip()
            is_hum = (row["var"].get() == 1)

            if var_name:
                self.ode_logic.add_variable(var_name)
                if eq_str: self.ode_logic.set_equation(var_name, eq_str)
                self.ode_logic.set_properties(var_name, False, is_hum)

    def open_calc(self):
        self.sync_logic()
        if not self.ode_logic.variables:
            messagebox.showerror("Error", "Please add at least one variable.")
            return
        CalcWindow(self, self.ode_logic, "ODE")

    def map_to_petri(self, target_type):
        self.sync_logic()
        if not self.ode_logic.variables:
            messagebox.showerror("Error", "No ODE system defined.")
            return

        try:
            # Step 1) Smart scaling detection, filter terms to decide if need to ask user anything

            has_risky_scaling = False

            # 1. Extract all flows temporarily to analyze them
            state_syms = {self.ode_logic.symbols[v] for v in self.ode_logic.variables}
            all_terms = []

            # Use the same "Smart Expansion" logic here to ensure consistent detection
            for var, eq_str in self.ode_logic.equations.items():
                if not eq_str: continue
                try:
                    # Smart Parse
                    raw_expr = sp.sympify(eq_str, locals=self.ode_logic.symbols)
                    expanded_terms = []
                    if isinstance(raw_expr, sp.Add):
                        base_terms = raw_expr.args
                    else:
                        base_terms = [raw_expr]

                    for term in base_terms:
                        has_state_add = False
                        for atom in term.atoms(sp.Add):
                            if atom.free_symbols & state_syms:
                                has_state_add = True;
                                break
                        if has_state_add:
                            exp = term.expand()
                            if isinstance(exp, sp.Add):
                                expanded_terms.extend(exp.args)
                            else:
                                expanded_terms.append(exp)
                        else:
                            expanded_terms.append(term)

                    for term in expanded_terms:
                        coeff = term.as_coeff_Mul()[0]
                        direction = -1 if coeff < 0 else 1
                        rate_struct = sp.simplify(term) if coeff > 0 else sp.simplify(-term)
                        all_terms.append({'place': var, 'expr': rate_struct, 'direction': direction})
                except:
                    continue

            outflows = [t for t in all_terms if t['direction'] < 0]
            inflows = [t for t in all_terms if t['direction'] > 0]

            # 2. Check for "Risky" Proportionality
            for in_item in inflows:
                for out_item in outflows:
                    # Is there a mathematical relationship?
                    try:
                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                    except:
                        continue

                    if not (ratio.free_symbols & state_syms):
                        #it is proportional (ratio is constant parameters/numbers)

                        #case A: Identical (Ratio == 1) is safer, Auto-Merge, don't ask
                        if ratio == 1:
                            continue

                        #case B: Different (Ratio != 1) riskier
                        # Could be Efficiency (Merge) OR Migration/Death (Separate).
                        #Need to  ask the user.
                        has_risky_scaling = True
                        break
                if has_risky_scaling: break

            use_scaling_choice = True  #default behavior (Merge Identicals)

            if has_risky_scaling:
                msg = ("Ambiguous Proportional Flows Detected!\n\n"
                       "The system detected terms where Inflow is proportional to Outflow (e.g. In = 0.5 * Out).\n"
                       "Please select a DEFAULT STRATEGY for the entire model:\n\n"
                       "YES (Merge - Recommended for most models):\n"
                       "   • Treat as single processes with loss/gain (e.g. Efficiency, Stoichiometry).\n"
                       "   • Preserves connections between places.\n\n"
                       "NO (Separate):\n"
                       "   • Treat as unrelated events (e.g. Independent Source & Sink).\n"
                       "   • May break flow connections if terms are actually linked.")
                use_scaling_choice = messagebox.askyesno("Mapping Strategy", msg)

            # Step 2) prep loop variables
            current_target_type = target_type
            merge_decisions = {}
            final_data = None

            # Outer loop handles swithc to VAPN restart
            while True:

                # The inner loop does ambiguity detection
                while True:
                    if current_target_type == "SPN":
                        result = self.ode_logic.convert_to_spn_structure(use_scaling=use_scaling_choice,
                                                                         merge_strategies=merge_decisions)
                    else:
                        result = self.ode_logic.convert_to_vapn_structure(use_scaling=use_scaling_choice,
                                                                          merge_strategies=merge_decisions)

                    if result.get("status") == "success":
                        final_data = result
                        break

                    elif result.get("status") == "ambiguous":
                        issue = result["details"][0]
                        decision = self.ask_ambiguity_strategy(issue, current_target_type)
                        if not decision: return  # User cancelled
                        merge_decisions[issue["key"]] = decision

                # Step 3) check for non valid spn weights
                restart_needed = False
                if current_target_type == "SPN":
                    has_non_int = False
                    for a in final_data.get("arcs", []):
                        w = str(a.get("weight", "1"))
                        if not w.isdigit():
                            has_non_int = True
                            break

                    if has_non_int:
                        warn = ODEMappingWarning(self)
                        if warn.result == "edit":
                            return
                        elif warn.result == "vapn":
                            current_target_type = "VAPN"
                            restart_needed = True

                if not restart_needed:
                    break

            # Step 4) create new tab
            if not final_data.get('transitions'):
                messagebox.showinfo("Info", "No valid transitions found.")
                return

            self.controller.map_from_ode(current_target_type, final_data)

        except Exception as e:
            messagebox.showerror("Mapping Error", str(e))

    def draw_ambiguity_viz(self, canvas, mode, place_name, catalyst_name, target_type):
        """Draws the topology based on Target Type (SPN vs VAPN)."""
        canvas.delete("all")

        r = 20
        tw, th = 30, 15
        y_mid = 100
        p_str = render_latex(place_name)
        c_str = render_latex(catalyst_name)

        #VAPN Visualization
        if target_type == "VAPN":
            x_N = 60
            x_P = 290

            canvas.create_oval(x_N - r, y_mid - r, x_N + r, y_mid + r, outline="black", width=2, fill="white")
            canvas.create_text(x_N, y_mid, text=p_str, font=("Arial", 12, "bold"))
            canvas.create_oval(x_P - r, y_mid - r, x_P + r, y_mid + r, outline="black", width=2, fill="white")
            canvas.create_text(x_P, y_mid, text=c_str, font=("Arial", 12, "bold"))

            if mode == "merge":
                x_t = 175
                y_t = 100
                canvas.create_rectangle(x_t - tw, y_t - th, x_t + tw, y_t + th, fill="black")
                canvas.create_text(x_t, y_t - 25, text="t1 (Net)", font=("Arial", 10, "bold"))

                canvas.create_line(x_N + r, y_mid, x_t - tw, y_t, arrow=tk.LAST, width=2)
                canvas.create_line(x_t + tw, y_t, x_P - r, y_mid, arrow=tk.LAST, width=2)
                canvas.create_text(175, 160, text=f"Structure: {p_str} converts to {c_str} (Net Flow)", fill="blue",
                                   font=("Arial", 10, "italic"))

            else:
                x_t = 175
                y_t1, y_t2 = 150, 50
                # t1 (Removal)
                canvas.create_rectangle(x_t - tw, y_t1 - th, x_t + tw, y_t1 + th, fill="black")
                canvas.create_text(x_t, y_t1 + 25, text="t1", font=("Arial", 9, "bold"))
                canvas.create_line(x_N + r, y_mid + 10, x_t - tw, y_t1, arrow=tk.LAST, width=2, smooth=True)
                canvas.create_line(x_t + tw, y_t1, x_P - r, y_mid + 10, arrow=tk.LAST, width=2, smooth=True)
                # t2 (Creation)
                canvas.create_rectangle(x_t - tw, y_t2 - th, x_t + tw, y_t2 + th, fill="black")
                canvas.create_text(x_t, y_t2 - 25, text="t2", font=("Arial", 9, "bold"))
                canvas.create_line(x_P - r, y_mid - 10, x_t + tw, y_t2, arrow=tk.LAST, width=2, smooth=True)
                canvas.create_line(x_t - tw, y_t2, x_N + r, y_mid - 10, arrow=tk.LAST, width=2, smooth=True)

                canvas.create_text(175, 100,
                                   text=f"Structure: Cycle ({c_str} creates {p_str}, {p_str} creates {c_str})",
                                   fill="#1b5e20", font=("Arial", 9, "italic"))

        # SPN Visualization
        else:
            x_N = 60
            x_P = 290
            canvas.create_oval(x_N - r, y_mid - r, x_N + r, y_mid + r, outline="black", width=2, fill="white")
            canvas.create_text(x_N, y_mid, text=p_str, font=("Arial", 14, "bold"))
            canvas.create_oval(x_P - r, y_mid - r, x_P + r, y_mid + r, outline="black", width=2, fill="white")
            canvas.create_text(x_P, y_mid, text=c_str, font=("Arial", 14, "bold"))

            if mode == "merge":
                x_t = 175
                y_t = 100
                canvas.create_rectangle(x_t - tw, y_t - th, x_t + tw, y_t + th, fill="black")
                canvas.create_text(x_t, y_t - 25, text="t1", font=("Arial", 10, "bold"))
                canvas.create_line(x_N + r, y_mid, x_t - tw, y_t, arrow=tk.LAST, width=2)
                # P Loop
                canvas.create_line(x_P - r, y_mid - 5, x_t + tw, y_t - 5, arrow=tk.LAST, width=2)
                canvas.create_line(x_t + tw, y_t + 5, x_P - r, y_mid + 5, arrow=tk.LAST, width=2)

                canvas.create_text(175, 160, text="Structure: Net Interaction (Catalyzed)", fill="blue",
                                   font=("Arial", 10, "italic"))

            else:
                x_t = 175
                y_t1, y_t2 = 150, 50
                # t1 Removal
                canvas.create_rectangle(x_t - tw, y_t1 - th, x_t + tw, y_t1 + th, fill="black")
                canvas.create_text(x_t, y_t1 + 25, text="t1", font=("Arial", 9, "bold"))
                canvas.create_line(x_N + r, y_mid + 10, x_t - tw, y_t1, arrow=tk.LAST, width=2, smooth=True)
                canvas.create_line(x_P - r, y_mid + 10, x_t + tw, y_t1 - 5, arrow=tk.LAST, width=2, smooth=True)
                canvas.create_line(x_t + tw, y_t1 + 5, x_P - r, y_mid + 20, arrow=tk.LAST, width=2, smooth=True)
                # t2 Creation
                canvas.create_rectangle(x_t - tw, y_t2 - th, x_t + tw, y_t2 + th, fill="black")
                canvas.create_text(x_t, y_t2 - 25, text="t2", font=("Arial", 9, "bold"))
                canvas.create_line(x_P - r, y_mid - 20, x_t + tw, y_t2 + 5, arrow=tk.LAST, width=2, smooth=True)
                canvas.create_line(x_t + tw, y_t2 - 5, x_P - r, y_mid - 10, arrow=tk.LAST, width=2, smooth=True)
                canvas.create_line(x_t - tw, y_t2, x_N + r, y_mid - 10, arrow=tk.LAST, width=2, smooth=True)

                # CHANGED COLOR FROM RED TO GREEN (#1b5e20)
                canvas.create_text(175, 100, text=f"Structure: Separate {p_str} Removal & Creation", fill="#1b5e20",
                                   font=("Arial", 9, "italic"))

    def ask_ambiguity_strategy(self, issue, target_type):
        win = Toplevel(self)
        win.title("Biological Ambiguity Detected")
        win.geometry("900x600")
        win.transient(self)
        win.grab_set()

        self.user_choice = None

        # Header
        tk.Label(win, text="Ambiguous Biological Interaction Detected", font=("Arial", 16, "bold"), fg="#d32f2f").pack(
            pady=(15, 5))

        msg = (
            f"Variable '{render_latex(issue['place'])}' has a SOURCE term that depends on '{render_latex(issue['catalyst'])}'.\n"
            f"However, '{render_latex(issue['catalyst'])}' also REMOVES '{render_latex(issue['place'])}' in another term.")
        tk.Label(win, text=msg, font=("Arial", 12), wraplength=800, justify="center").pack(pady=5)

        # Main Content Frame
        frame_vis = tk.Frame(win, bd=2, relief="groove", bg="white")
        frame_vis.pack(fill="both", expand=True, padx=20, pady=10)

        # left side: merge
        f1 = tk.Frame(frame_vis, bg="#e3f2fd")  # Light Blue
        f1.pack(side="left", fill="both", expand=True, padx=5, pady=5)

        tk.Label(f1, text="Option 1: CORRECTION\n(Merge Terms)", font=("Arial", 14, "bold"), bg="#e3f2fd",
                 fg="#0d47a1").pack(pady=10)

        c1 = tk.Canvas(f1, width=350, height=200, bg="#e3f2fd", highlightthickness=0)
        c1.pack(pady=5)
        self.draw_ambiguity_viz(c1, "merge", issue['place'], issue['catalyst'], target_type)

        lbl_rate1 = f"Net Rate = {render_latex(issue['outflow_expr'])} - {render_latex(issue['inflow_expr'])}"
        tk.Label(f1, text=lbl_rate1, font=("Courier", 10, "bold"), bg="white", relief="solid", padx=5, pady=5,
                 wraplength=300).pack(pady=10)

        # right side: separate
        f2 = tk.Frame(frame_vis, bg="#e8f5e9")  # Light Green (Green 50)
        f2.pack(side="right", fill="both", expand=True, padx=5, pady=5)

        tk.Label(f2, text="Option 2: NEW PROCESS\n(Keep Separate)", font=("Arial", 14, "bold"), bg="#e8f5e9",
                 fg="#1b5e20").pack(pady=10)

        c2 = tk.Canvas(f2, width=350, height=200, bg="#e8f5e9", highlightthickness=0)
        c2.pack(pady=5)
        self.draw_ambiguity_viz(c2, "separate", issue['place'], issue['catalyst'], target_type)

        lbl_rate2 = f"Rate 1: {render_latex(issue['outflow_expr'])}\nRate 2: {render_latex(issue['inflow_expr'])}"

        tk.Label(f2, text=lbl_rate2, font=("Courier", 10, "bold"), bg="white", fg="#1b5e20", relief="solid", padx=5,
                 pady=5, wraplength=300).pack(pady=10)

        # Buttons
        btn_frame = tk.Frame(win)
        btn_frame.pack(fill="x", pady=20)

        def set_choice(c):
            self.user_choice = c
            win.destroy()

        tk.Button(btn_frame, text="Use Option 1 (MERGE)", bg="#bbdefb", font=("Arial", 12, "bold"), height=2,
                  command=lambda: set_choice("merge")).pack(side="left", expand=True, padx=40, ipadx=10)

        tk.Button(btn_frame, text="Use Option 2 (SEPARATE)", bg="#a5d6a7", font=("Arial", 12, "bold"), height=2,
                  command=lambda: set_choice("separate")).pack(side="right", expand=True, padx=40, ipadx=10)

        self.wait_window(win)
        return self.user_choice


#shared calculation popup
class CalcWindow(Toplevel):
    def __init__(self, parent, logic_obj, mode_str):
        super().__init__(parent)
        self.logic = logic_obj
        self.title(f"R₀ Finder - Calculation ({mode_str})")
        self.geometry("1300x850")  # Slightly taller for the new button

        # Store latest results for sensitivity
        self.latest_sensitivity_data = None
        self.latest_r0_list = None

        #1. Data prep
        if mode_str == "ODE":
            self.places = logic_obj.variables
            self.infected_preselect = logic_obj.infected_vars
        else:
            self.places = logic_obj.places
            self.infected_preselect = []

            # 2. Layout splitting left/right
        main_split = Frame(self)
        main_split.pack(fill=BOTH, expand=True, padx=10, pady=10)

        # Left panel (Existing Controls)
        left_panel = Frame(main_split)
        left_panel.pack(side=LEFT, fill=BOTH, expand=True)

        # Right panel (Symbol Table)
        right_panel = Frame(main_split, width=300, relief=GROOVE, bd=1)
        right_panel.pack(side=RIGHT, fill=Y, padx=(15, 0))

        # 3)build right panel (Table)
        self.create_symbol_table(right_panel)

        # 4. build left panel

        # Step 1
        Label(left_panel, text="Step 1: Confirm INFECTED Variables", font=FONT_HEADER).pack(pady=10)
        self.check_vars_inf = {}
        frame_inf = Frame(left_panel)
        frame_inf.pack(fill=tk.X, padx=10)

        # Step 2
        Label(left_panel, text="Optional Step 2: Free Parameters (Assumed Non-Zero in DFE Calculation)", font=FONT_HEADER).pack(
            pady=(15, 5))
        self.check_vars_free = {}
        frame_free = Frame(left_panel)
        frame_free.pack(fill=tk.X, padx=10)

        for p in self.places:
            var_inf = IntVar()
            self.check_vars_inf[p] = var_inf
            var_free = IntVar()
            self.check_vars_free[p] = var_free

            if p in self.infected_preselect:
                var_inf.set(1)

            def on_check(place=p):
                if self.check_vars_inf[place].get() == 1:
                    self.check_vars_free[place].set(0)

            def on_check_free(place=p):
                if self.check_vars_free[place].get() == 1:
                    self.check_vars_inf[place].set(0)

            Checkbutton(frame_inf, text=p, variable=var_inf, command=on_check, font=FONT_UI).pack(side=LEFT, padx=5)
            Checkbutton(frame_free, text=p, variable=var_free, command=on_check_free, font=FONT_UI).pack(side=LEFT,
                                                                                                         padx=5)

        # Step 3
        Label(left_panel, text="Optional Step 3: Parameter Constraints for DFE Configuration", font=FONT_HEADER).pack(pady=15)

        frame_cons = Frame(left_panel)
        frame_cons.pack(fill=tk.X, padx=20, pady=5)

        Label(frame_cons, text="Optional Constraints (e.g. q+p=1):", font=FONT_UI).pack(side=LEFT)
        self.ent_cons = Entry(frame_cons, width=20, font=FONT_UI)
        self.ent_cons.pack(side=LEFT, padx=5)

        self.cons_list = []
        self.lbl_cons_display = Label(left_panel, text="", font=("Arial", 10, "italic"), fg="blue")
        self.lbl_cons_display.pack()

        def add_constraint():
            c = self.ent_cons.get().strip()
            if c:
                self.cons_list.append(c)
                self.update_cons_display()
                self.ent_cons.delete(0, END)

        tk.Button(frame_cons, text="Add", command=add_constraint).pack(side=LEFT)
        tk.Button(frame_cons, text="Clear All", command=lambda: [self.cons_list.clear(), self.update_cons_display()],
                  fg="red").pack(side=LEFT, padx=5)

        tk.Button(left_panel, text="Auto-Detect DFE", command=self.auto_detect_dfe, font=FONT_UI).pack(pady=5)
        tk.Label(left_panel, text="DFE can be edited directly. Inaccurate DFE can lead to an inaccurate R₀",
                 font=("Arial", 10, "italic"), fg="red").pack(pady=(0, 5))

        # DFE Grid Container
        container = Frame(left_panel, height=200, bd=1, relief=GROOVE)
        container.pack(fill=tk.BOTH, padx=10, pady=10, expand=False)

        canvas = tk.Canvas(container)
        scr_y = Scrollbar(container, orient="vertical", command=canvas.yview)
        scr_x = Scrollbar(container, orient="horizontal", command=canvas.xview)
        self.dfe_frame = Frame(canvas)
        self.dfe_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.dfe_frame, anchor="nw")
        canvas.configure(yscrollcommand=scr_y.set, xscrollcommand=scr_x.set)
        scr_y.pack(side="right", fill="y")
        scr_x.pack(side="bottom", fill="x")
        canvas.pack(side="left", fill="both", expand=True)

        self.dfe_widgets = {}
        self.num_dfe_cols = 1
        self.setup_dfe_grid()

        # Watermark
        tk.Label(left_panel, text=WATERMARK_TEXT, font=FONT_WATERMARK, fg="gray50", wraplength=700).pack(pady=(15, 0))

        # Calculate Button
        self.btn_calc = tk.Button(left_panel, text="CALCULATE R₀", bg="#d0f0c0", font=("Arial", 16, "bold"),
                                  command=self.run_calc)
        self.btn_calc.pack(pady=5, ipadx=20, ipady=5)

        # Sensitivity button
        self.btn_sens = tk.Button(left_panel, text="Sensitivity for R₀", bg="#c0e0ff", font=("Arial", 12, "bold"),
                                  state="disabled", command=self.run_sensitivity)
        self.btn_sens.pack(pady=5, ipadx=10)

        # Results Text
        self.txt_res = Text(left_panel, height=10, font=FONT_UI, spacing3=10)
        self.txt_res.pack(fill=BOTH, padx=10, pady=5, expand=True)

    def create_symbol_table(self, parent_frame):
        """Creates the right-hand table listing all detected symbols."""
        Label(parent_frame, text="Detected Symbols", font=FONT_BOLD, bg="#eee").pack(fill=tk.X, ipady=5)

        vars_set = set(self.places)
        all_syms_set = set(self.logic.symbols.keys())
        params_set = all_syms_set - vars_set
        if 'N' in params_set: params_set.remove('N')
        if 't' in params_set: params_set.remove('t')

        cols = ("Type", "Name")
        tree = ttk.Treeview(parent_frame, columns=cols, show='headings', height=25)
        tree.heading("Type", text="Type")
        tree.heading("Name", text="Symbol Name")
        tree.column("Type", width=100)
        tree.column("Name", width=150)

        sb = Scrollbar(parent_frame, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=sb.set)

        sorted_vars = sorted(list(vars_set))
        sorted_params = sorted(list(params_set))

        for v in sorted_vars:
            # Decode for display: greek_beta, β
            display_v = decode_latex(v)
            tree.insert("", "end", values=("Variable", display_v), tags=('var',))

        for p in sorted_params:
            display_p = decode_latex(p)
            tree.insert("", "end", values=("Parameter", display_p), tags=('param',))

        tree.tag_configure('var', background="#e0f7fa")
        tree.tag_configure('param', background="#fff3e0")

        lbl_hint = Label(parent_frame, text="Check this list for typos!\n(e.g. 'bta' vs 'beta')",
                         fg="red", font=("Arial", 10, "italic"), justify="center")
        lbl_hint.pack(side=tk.BOTTOM, pady=10)

        sb.pack(side=RIGHT, fill=Y)
        tree.pack(side=LEFT, fill=BOTH, expand=True)

    def update_cons_display(self):
        txt = "Active Constraints: " + ", ".join(self.cons_list)
        self.lbl_cons_display.config(text=txt)

    def setup_dfe_grid(self, solutions_list=None):
        for w in self.dfe_frame.winfo_children(): w.destroy()
        self.dfe_widgets = {}
        if solutions_list:
            self.num_dfe_cols = len(solutions_list)
        else:
            self.num_dfe_cols = 1;
            solutions_list = [{}]

        Label(self.dfe_frame, text="Var", width=10, relief="solid", font=FONT_BOLD).grid(row=0, column=0, sticky="nsew")
        for c in range(self.num_dfe_cols):
            Label(self.dfe_frame, text=f"DFE {c + 1}", width=35, relief="solid", font=FONT_BOLD).grid(row=0,
                                                                                                      column=c + 1,
                                                                                                      sticky="nsew")

        for r, p in enumerate(self.places):
            Label(self.dfe_frame, text=p, width=10, font=FONT_UI).grid(row=r + 1, column=0, sticky="nsew")
            for c in range(self.num_dfe_cols):
                txt = Text(self.dfe_frame, width=40, height=5, font=("Arial", 11), wrap="char", borderwidth=1,
                           relief="sunken")
                txt.grid(row=r + 1, column=c + 1, padx=2, pady=2, sticky="nsew")
                val = str(solutions_list[c].get(p, "0"))
                txt.insert("1.0", val)
                self.dfe_widgets[(p, c)] = txt

    def auto_detect_dfe(self):
        inf = [p for p, v in self.check_vars_inf.items() if v.get() == 1]
        free = [p for p, v in self.check_vars_free.items() if v.get() == 1]

        if not inf: messagebox.showerror("Error", "Select infected vars"); return

        sols, msg = self.logic.calculate_dfe(inf, free, self.cons_list)

        if not sols: messagebox.showerror("Error", msg); return
        self.setup_dfe_grid(sols)
        for c in range(self.num_dfe_cols):
            for p in inf:
                if (p, c) in self.dfe_widgets:
                    self.dfe_widgets[(p, c)].delete("1.0", END);
                    self.dfe_widgets[(p, c)].insert("1.0", "0")

    def run_calc(self):
        inf = [p for p, v in self.check_vars_inf.items() if v.get() == 1]
        if not inf: messagebox.showerror("Error", "Select infected vars"); return
        self.txt_res.delete("1.0", END)

        self.latest_sensitivity_data = []
        transition_decisions = {}

        for c in range(self.num_dfe_cols):
            dfe = {}
            for p in self.places:
                if (p, c) in self.dfe_widgets:
                    dfe[p] = self.dfe_widgets[(p, c)].get("1.0", "end-1c").strip()

            dfe_values_str = ", ".join([f"{k}={v}" for k, v in dfe.items()])
            self.txt_res.insert(END, f"  Results for DFE {c + 1}: ({dfe_values_str})  \n")

            while True:
                result = self.logic.calculate_r0(inf, dfe, source_classification=transition_decisions)

                if result[1] == "ambiguous":
                    ambiguous_list = result[2]
                    for item in ambiguous_list:
                        t_name = item['name']
                        rate = item['rate']
                        if t_name in transition_decisions: continue

                        msg = (f"Ambiguous Transition: '{t_name}'\nRate: {rate}\n\n"
                               f"Classify Source:\n"
                               f"YES = New Infection (F)\nNO = Migration (V)")
                        decision = messagebox.askyesno("Classify Source", msg)
                        transition_decisions[t_name] = 'F' if decision else 'V'
                    continue
                else:
                    # SAFE UNPACKING (Handles standard 6-item return)
                    if len(result) >= 6:
                        res, F_mat, V_mat, sens_data, r0_raw_list, n_def = result
                    else:
                        # Fallback for Error Messages
                        res = result[0]
                        F_mat = V_mat = n_def = None
                    break

            self.latest_sensitivity_data.append(sens_data)

            # Display N* Definition
            if n_def:
                self.txt_res.insert(END, f"Note: {n_def}\n\n")

            if F_mat: self.txt_res.insert(END, f"F Matrix (New Infections):\n{F_mat}\n\n")
            if V_mat: self.txt_res.insert(END, f"V Matrix (Transitions):\n{V_mat}\n\n")

            for r in res: self.txt_res.insert(END, f"R₀ = {r}\n")

            if len(res) > 1:
                self.txt_res.insert(END, "\n*** NOTE: MULTIPLE R₀ VALUES DETECTED ***\n")
                self.txt_res.insert(END, "System R₀ = max(R₀_1, R₀_2, ...)\n")

            self.txt_res.insert(END, "\n" + "=" * 40 + "\n")

        self.btn_sens.config(state="normal")

    def run_sensitivity(self):
        if not self.latest_sensitivity_data:
            messagebox.showinfo("Info", "No sensitivity data available. Run R0 calculation first.")
            return

        # Create Popup
        win = Toplevel(self)
        win.title("Sensitivity Analysis (Elasticity Indices)")
        win.geometry("1100x700")

        # Split Layout
        frame_text = Frame(win, width=400)
        frame_text.pack(side=LEFT, fill=BOTH, expand=False, padx=10, pady=10)

        frame_graph = Frame(win)
        frame_graph.pack(side=RIGHT, fill=BOTH, expand=True, padx=10, pady=10)

        # 1. Text Results (Left Side)
        lbl = Label(frame_text, text="Elasticity Indices (Symbolic)", font=FONT_BOLD)
        lbl.pack(anchor="w")
        txt = Text(frame_text, width=40, font=("Arial", 11))
        txt.pack(fill=BOTH, expand=True)

        # Populate Text
        for idx, dfe_data in enumerate(self.latest_sensitivity_data):
            txt.insert(END, f"  DFE {idx + 1} Indices  \n")
            for r_idx, sens_map in enumerate(dfe_data):
                txt.insert(END, f"[R0 Solution {r_idx + 1}]\n")
                if not sens_map:
                    txt.insert(END, "  (None)\n")
                for param, val_expr in sens_map.items():
                    txt.insert(END, f"  Upsilon_{param} = {val_expr}\n")
                txt.insert(END, "\n")

        # 2. Graph Input Area (Right Side)

        # Dropdown Selection
        select_frame = Frame(frame_graph)
        select_frame.pack(fill=tk.X, pady=(0, 10))
        Label(select_frame, text="Select R₀ Solution to Analyze:", font=FONT_BOLD).pack(side=LEFT)

        options = []
        map_registry = {}
        for dfe_i, dfe_list in enumerate(self.latest_sensitivity_data):
            for r0_i, _ in enumerate(dfe_list):
                label = f"DFE {dfe_i + 1} - R₀ Solution {r0_i + 1}"
                options.append(label)
                map_registry[label] = (dfe_i, r0_i)

        combo_select = ttk.Combobox(select_frame, values=options, state="readonly", width=30)
        if options: combo_select.current(0)
        combo_select.pack(side=LEFT, padx=5)

        Label(frame_graph, text="Enter Values for R₀ Parameters to View Normalized Sensitivity", font=FONT_BOLD).pack(pady=5)

        # Input Frame (Where entries will be rebuilt)
        input_frame = Frame(frame_graph)
        input_frame.pack(fill=tk.X, pady=5)

        # Storage for current entry widgets
        self.sens_entries = {}

        # this is helper to rebuild input grid based on selection
        def rebuild_inputs(target_map):
            # 1. Clear old inputs
            for widget in input_frame.winfo_children():
                widget.destroy()
            self.sens_entries.clear()

            if not target_map: return

            # 2. Identify all parameters needed
            # need the parameter itself (key) and any symbols used in its formula (value)
            active_symbols = set()
            for p_key, expr in target_map.items():
                active_symbols.add(str(p_key))
                try:
                    # Convert string expr back to sympy to find dependencies
                    sym_expr = sp.sympify(expr)
                    for s in sym_expr.free_symbols:
                        active_symbols.add(str(s))
                except:
                    pass

            # Filter out junk
            valid_params = sorted([p for p in active_symbols if p not in ["Error", "None"]])

            # 3. Build Grid
            r, c = 0, 0
            for p in valid_params:
                f = Frame(input_frame)
                f.grid(row=r, column=c, padx=5, pady=2, sticky="w")
                Label(f, text=f"{p} = ").pack(side=LEFT)
                e = Entry(f, width=8)
                e.insert(0, "1.0")  # Default value
                e.pack(side=LEFT)
                self.sens_entries[p] = e

                c += 1
                if c > 3:  # 4 columns max
                    c = 0;
                    r += 1

        # Canvas for Graph
        self.fig, self.ax = plt.subplots(figsize=(5, 4))
        self.canvas_widget = FigureCanvasTkAgg(self.fig, master=frame_graph)
        self.canvas_widget.get_tk_widget().pack(fill=BOTH, expand=True, pady=10)

        #  Scrollable Numerical Results
        res_frame = Frame(frame_graph, height=100)  # Fixed height for results
        res_frame.pack(fill=tk.X, pady=5, side=tk.BOTTOM)

        sb_res = Scrollbar(res_frame)
        sb_res.pack(side=RIGHT, fill=Y)

        txt_numeric_results = Text(res_frame, height=6, font=("Arial", 10), fg="blue",
                                   yscrollcommand=sb_res.set, relief=GROOVE, bd=1)
        txt_numeric_results.pack(side=LEFT, fill=BOTH, expand=True)
        sb_res.config(command=txt_numeric_results.yview)

        # update logic
        def update_graph(event=None):
            selection = combo_select.get()
            if selection not in map_registry: return

            dfe_idx, r0_idx = map_registry[selection]
            target_map = self.latest_sensitivity_data[dfe_idx][r0_idx]

            if not target_map: return

            # Get values
            vals = {}
            for p, ent in self.sens_entries.items():
                try:
                    vals[p] = float(ent.get())
                except:
                    vals[p] = 1.0

                    # Calculate & Plot
            labels = []
            values = []
            result_str = "Evaluated Elasticity Indices:\n"

            for p, expr in target_map.items():
                try:
                    sub_dict = {sp.Symbol(k, real=True): v for k, v in vals.items()}
                    sym_expr = sp.sympify(expr)
                    num_val = float(sym_expr.subs(sub_dict).evalf())

                    labels.append(p)
                    values.append(num_val)
                    result_str += f"  {p}: {num_val:.4f}\n"
                except:
                    pass

            # Update Scrollable Text
            txt_numeric_results.config(state="normal")
            txt_numeric_results.delete("1.0", END)
            txt_numeric_results.insert(END, result_str)
            txt_numeric_results.config(state="disabled")  # Read-only

            self.ax.clear()
            if values:
                bars = self.ax.bar(labels, values, color=['skyblue' if v >= 0 else 'salmon' for v in values])
                self.ax.axhline(0, color='black', linewidth=0.8)
                self.ax.set_title(f"R₀ Parameters Normalized Sensitivity ({selection})")
                self.ax.set_ylabel("Index Value")

            self.canvas_widget.draw()

        # event handler
        def on_selection_change(event):
            selection = combo_select.get()
            if selection in map_registry:
                dfe_idx, r0_idx = map_registry[selection]
                target_map = self.latest_sensitivity_data[dfe_idx][r0_idx]
                rebuild_inputs(target_map)
                update_graph()

        btn_plot = tk.Button(frame_graph, text="Update Graph", bg="#ddd", command=update_graph)
        btn_plot.pack(before=self.canvas_widget.get_tk_widget(), pady=5)

        # Bind Event
        combo_select.bind("<<ComboboxSelected>>", on_selection_change)

        # Initialize
        if options:
            on_selection_change(None)



#     Main controller (tabbed)

class MainApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("R₀ Finder")
        self.geometry("1300x800")

        try:
            icon_img = tk.PhotoImage(file="R_0applogo.png")
            self.iconphoto(True, icon_img)
        except Exception as e:
            print(f"Warning: Could not load app icon: {e}")

        # Setup Tabs (Notebook)
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True)

        # Create Home Tab (Dashboard)
        self.create_start_tab()

    def create_start_tab(self):
        start_frame = StartScreen(self.notebook, self)
        self.notebook.add(start_frame, text="Home")

    def add_tab(self, mode_type, data=None):
        # Limit to 6 model tabs (plus the 1 Home tab = 7 total)
        if len(self.notebook.tabs()) >= 7:
            messagebox.showwarning("Limit Reached",
                                   "You can only have 6 active model tabs open at once.\nPlease close a tab to create a new one.")
            return

        if mode_type == "ODE":
            editor = ODEEditor(self.notebook, self)
            title = "ODE Model"
        elif mode_type == "SPN":
            editor = PetriNetEditor(self.notebook, self, "SPN")
            title = "SPN Model"
        else:
            editor = PetriNetEditor(self.notebook, self, "VAPN")
            title = "VAPN Model"

        # Add and Select Tab
        self.notebook.add(editor, text=title)
        self.notebook.select(editor)

        # If mapping data exists, load it
        if data:
            editor.load_from_data(data)

    def close_tab(self, tab_frame):
        self.notebook.forget(tab_frame)
        tab_frame.destroy()

    def map_from_ode(self, target_type, data):
        # Maps ODE to SPN/VAPN by opening a new tab with the data
        self.add_tab(target_type, data)

    def report_callback_exception(self, exc_type, exc_value, exc_traceback):
        error_msg = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
        messagebox.showerror("Critical Error", f"An unexpected error occurred:\n\n{error_msg}")


if __name__ == "__main__":
    app = MainApp()
    app.mainloop()