import tkinter as tk
from tkinter import simpledialog, messagebox, Toplevel, Checkbutton, IntVar, Label, Text, END, Entry, Frame, Scrollbar, \
    LEFT, RIGHT, Y, BOTH
import math
import traceback  # Import traceback to get detailed error logs
import sympy as sp
from petrinet_logic import VariableWeightPetriNet


class PetriNetApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Petri Net R0 Calculator (VAPN)")

        self.pn_logic = VariableWeightPetriNet()

        self.mode = "SELECT"
        self.nodes_ui = []
        self.arcs_ui = []
        self.temp_arc_start = None
        self.spawn_x = 50
        self.spawn_y = 100
        self.drag_data = {"x": 0, "y": 0, "item": None}

        # --- Toolbar ---
        self.toolbar = tk.Frame(root, bg="#f0f0f0", height=40)
        self.toolbar.pack(side=tk.TOP, fill=tk.X, pady=2)

        # Consistent font for buttons to minimize macOS warnings
        btn_font = ("Arial", 11)

        # Tools
        tk.Button(self.toolbar, text="Select/Move (Esc)", font=btn_font, command=lambda: self.set_mode("SELECT")).pack(
            side=tk.LEFT, padx=2)
        tk.Button(self.toolbar, text="Add Place (O)", font=btn_font, command=self.trigger_add_place).pack(side=tk.LEFT,
                                                                                                          padx=2)
        tk.Button(self.toolbar, text="Add Transition (|)", font=btn_font, command=self.trigger_add_trans).pack(
            side=tk.LEFT, padx=2)
        tk.Button(self.toolbar, text="Add Arc (->)", font=btn_font, command=lambda: self.set_mode("ARC")).pack(
            side=tk.LEFT, padx=2)

        # Clear Button (Red)
        tk.Button(self.toolbar, text="Clear All", font=btn_font, bg="#ffdddd", command=self.trigger_clear_all).pack(
            side=tk.LEFT, padx=10)

        tk.Label(self.toolbar, text=" | ").pack(side=tk.LEFT)

        # Calculate Button (Green)
        btn_calc = tk.Button(self.toolbar, text="CALCULATE R0", font=("Arial", 11, "bold"), bg="#d0f0c0",
                             command=self.open_calculation_window)
        btn_calc.pack(side=tk.LEFT, padx=5)

        # --- Canvas ---
        self.canvas = tk.Canvas(root, bg="white", width=900, height=600)
        self.canvas.pack(fill=tk.BOTH, expand=True)

        # --- Status Bar ---
        self.status_var = tk.StringVar()
        self.status_var.set("Mode: SELECT - Click and drag to move nodes.")
        self.status_bar = tk.Label(root, textvariable=self.status_var, bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

        # --- Bindings ---
        self.canvas.bind("<ButtonPress-1>", self.on_left_click_down)
        self.canvas.bind("<B1-Motion>", self.on_left_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_left_click_up)

        # Handle Right Click for Windows/Linux (Button-3) and Mac (Button-2)
        self.canvas.bind("<Button-3>", self.on_right_click)
        if root.tk.call('tk', 'windowingsystem') == 'aqua':
            self.canvas.bind("<Button-2>", self.on_right_click)

        root.bind("<Escape>", lambda e: self.set_mode("SELECT"))

    def trigger_clear_all(self):
        """Clears the canvas and resets logic after confirmation."""
        if not self.nodes_ui and not self.arcs_ui:
            messagebox.showinfo("Info", "Canvas is already empty.")
            return

        confirm = messagebox.askyesno("Confirm Clear",
                                      "Are you sure you want to delete ALL places, transitions, and arcs?\n\nThis cannot be undone.")
        if confirm:
            self.nodes_ui = []
            self.arcs_ui = []
            self.temp_arc_start = None
            self.canvas.delete("all")
            self.pn_logic = VariableWeightPetriNet()
            self.spawn_x = 50
            self.spawn_y = 100
            self.status_var.set("Canvas cleared. Ready for new model.")

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
            tx = (w if dx > 0 else -w) / slope_safe_dx
            ty = (h if dy > 0 else -h) / slope_safe_dy
            if abs(tx) <= abs(ty):
                return cx + tx * slope_safe_dx, cy + tx * slope_safe_dy
            else:
                return cx + ty * slope_safe_dx, cy + ty * slope_safe_dy

    def get_spawn_pos(self):
        x, y = self.spawn_x, self.spawn_y
        self.spawn_x += 70
        if self.spawn_x > 800: self.spawn_x = 50; self.spawn_y += 70
        return x, y

    def trigger_add_place(self):
        prompt = "Enter Place Name (e.g., S)\nMultiple: S,I,R"
        name_input = simpledialog.askstring("Input", prompt, parent=self.root)
        if name_input:
            names = name_input.split(',')
            added = 0
            for name in names:
                name = name.strip()
                if not name or name in self.pn_logic.places: continue
                x, y = self.get_spawn_pos()
                self.add_place_ui(x, y, name)
                self.pn_logic.add_place(name)
                added += 1
            if added > 0: self.set_mode("SELECT")

    def trigger_add_trans(self):
        prompt = "Enter Transition Name (e.g., t1)\nMultiple: t1,t2"
        name_input = simpledialog.askstring("Input", prompt, parent=self.root)
        if name_input:
            names = name_input.split(',')
            added = 0
            for name in names:
                name = name.strip()
                if not name or name in self.pn_logic.transitions: continue
                x, y = self.get_spawn_pos()
                self.add_trans_ui(x, y, name)
                self.pn_logic.add_transition(name)
                added += 1
            if added > 0: self.set_mode("SELECT")

    def set_mode(self, mode):
        self.mode = mode
        if self.temp_arc_start:
            self.canvas.itemconfig(self.temp_arc_start['id'], outline="black", width=2)
            self.temp_arc_start = None
        self.status_var.set(f"Mode: {mode}")

    def on_left_click_down(self, event):
        x, y = event.x, event.y
        clicked = self.find_node_at(x, y)
        if self.mode == "SELECT" and clicked:
            self.drag_data["item"] = clicked
            self.drag_data["x"] = x
            self.drag_data["y"] = y
        elif self.mode == "ARC" and clicked:
            if not self.temp_arc_start:
                self.temp_arc_start = clicked
                self.canvas.itemconfig(clicked['id'], outline="red", width=3)
            else:
                start = self.temp_arc_start
                end = clicked
                if start != end:
                    w = simpledialog.askstring("Weight", f"Weight {start['name']}->{end['name']}:", initialvalue="1",
                                               parent=self.root)
                    if w:
                        try:
                            self.pn_logic.add_arc(start['name'], end['name'], w)
                            self.add_arc_ui(start, end, w)
                        except Exception as e:
                            messagebox.showerror("Error", str(e))
                self.canvas.itemconfig(start['id'], outline="black", width=2)
                self.temp_arc_start = None
                self.set_mode("SELECT")

    def on_left_drag(self, event):
        if self.mode == "SELECT" and self.drag_data["item"]:
            node = self.drag_data["item"]
            dx = event.x - self.drag_data["x"]
            dy = event.y - self.drag_data["y"]
            self.canvas.move(node['id'], dx, dy)
            self.canvas.move(node['label_id'], dx, dy)
            node['x'] += dx;
            node['y'] += dy
            self.drag_data["x"] = event.x;
            self.drag_data["y"] = event.y
            self.redraw_connected_arcs(node)

    def on_left_click_up(self, event):
        self.drag_data["item"] = None

    def on_right_click(self, event):
        item_id = self.canvas.find_closest(event.x, event.y)[0]
        tags = self.canvas.gettags(item_id)
        clicked_node = None
        clicked_arc = None
        for n in self.nodes_ui:
            if n['id'] == item_id or n['label_id'] == item_id: clicked_node = n; break
        if not clicked_node and "arc_text" in tags:
            for a in self.arcs_ui:
                if a['text_id'] == item_id: clicked_arc = a; break
        menu = tk.Menu(self.root, tearoff=0)
        if clicked_node:
            menu.add_command(label="Rename", command=lambda: self.rename_node(clicked_node))
            menu.add_command(label="Delete", command=lambda: self.delete_node(clicked_node))
            menu.tk_popup(event.x_root, event.y_root)
        elif clicked_arc:
            menu.add_command(label="Edit Weight", command=lambda: self.edit_arc_weight(clicked_arc))
            menu.add_command(label="Delete", command=lambda: self.delete_arc(clicked_arc))
            menu.tk_popup(event.x_root, event.y_root)

    def find_node_at(self, x, y):
        for n in self.nodes_ui:
            if math.sqrt((x - n['x']) ** 2 + (y - n['y']) ** 2) < max(n['w'], n['h']) + 5: return n
        return None

    def redraw_connected_arcs(self, node):
        for a in self.arcs_ui:
            if a['start_node'] == node or a['end_node'] == node:
                sx, sy = self.get_edge_coords(a['start_node'], a['end_node']['x'], a['end_node']['y'])
                ex, ey = self.get_edge_coords(a['end_node'], a['start_node']['x'], a['start_node']['y'])
                self.canvas.coords(a['line_id'], sx, sy, ex, ey)
                self.canvas.coords(a['text_id'], (sx + ex) / 2, (sy + ey) / 2 - 15)

    def add_place_ui(self, x, y, name):
        r = 20
        uid = self.canvas.create_oval(x - r, y - r, x + r, y + r, fill="white", outline="black", width=2, tags="node")
        lid = self.canvas.create_text(x, y, text=name, font=("Arial", 12, "bold"), tags="node")
        self.nodes_ui.append(
            {'id': uid, 'label_id': lid, 'type': 'place', 'name': name, 'x': x, 'y': y, 'w': r, 'h': r})

    def add_trans_ui(self, x, y, name):
        w, h = 25, 15
        uid = self.canvas.create_rectangle(x - w, y - h, x + w, y + h, fill="black", tags="node")
        lid = self.canvas.create_text(x, y, text=name, font=("Arial", 10, "bold"), fill="white", tags="node")
        self.nodes_ui.append(
            {'id': uid, 'label_id': lid, 'type': 'trans', 'name': name, 'x': x, 'y': y, 'w': w, 'h': h})

    def add_arc_ui(self, start, end, weight):
        sx, sy = self.get_edge_coords(start, end['x'], end['y'])
        ex, ey = self.get_edge_coords(end, start['x'], start['y'])
        lid = self.canvas.create_line(sx, sy, ex, ey, arrow=tk.LAST, width=2, smooth=True, tags="arc_line")
        tid = self.canvas.create_text((sx + ex) / 2, (sy + ey) / 2 - 15, text=weight, fill="blue", font=("Arial", 9),
                                      tags="arc_text")
        self.arcs_ui.append(
            {'line_id': lid, 'text_id': tid, 'start_node': start, 'end_node': end, 'weight_str': weight})

    def rename_node(self, node):
        new = simpledialog.askstring("Rename", f"Rename {node['name']}:", initialvalue=node['name'], parent=self.root)
        if new and new != node['name']:
            if node['type'] == 'place':
                self.pn_logic.remove_place(node['name']);
                self.pn_logic.add_place(new)
            else:
                self.pn_logic.remove_transition(node['name']);
                self.pn_logic.add_transition(new)
            node['name'] = new
            self.canvas.itemconfig(node['label_id'], text=new)

    def delete_node(self, node):
        self.canvas.delete(node['id']);
        self.canvas.delete(node['label_id'])
        to_del = [a for a in self.arcs_ui if a['start_node'] == node or a['end_node'] == node]
        for a in to_del: self.delete_arc(a)
        if node['type'] == 'place':
            self.pn_logic.remove_place(node['name'])
        else:
            self.pn_logic.remove_transition(node['name'])
        self.nodes_ui.remove(node)

    def delete_arc(self, arc):
        self.canvas.delete(arc['line_id']);
        self.canvas.delete(arc['text_id'])
        self.pn_logic.remove_arc(arc['start_node']['name'], arc['end_node']['name'])
        self.arcs_ui.remove(arc)

    def edit_arc_weight(self, arc):
        w = simpledialog.askstring("Edit", "Weight:", initialvalue=arc['weight_str'], parent=self.root)
        if w:
            try:
                s, t = arc['start_node']['name'], arc['end_node']['name']
                self.pn_logic.remove_arc(s, t);
                self.pn_logic.add_arc(s, t, w)
                arc['weight_str'] = w
                self.canvas.itemconfig(arc['text_id'], text=w)
            except Exception as e:
                messagebox.showerror("Error", str(e))

    # --- CALCULATION WINDOW ---
    def open_calculation_window(self):
        if not self.pn_logic.places:
            messagebox.showinfo("Error", "No places defined.")
            return

        win = Toplevel(self.root)
        win.title("R0 Calculation Configuration")
        win.geometry("500x650")

        # Step 1
        Label(win, text="Step 1: Select INFECTED Places", font=("Arial", 11, "bold")).pack(pady=5)
        self.check_vars = {}
        frame_checks = Frame(win)
        frame_checks.pack(fill=tk.X, padx=10)
        for p in self.pn_logic.places:
            var = IntVar()
            cb = Checkbutton(frame_checks, text=p, variable=var)
            cb.pack(side=LEFT, padx=5)
            self.check_vars[p] = var

        # Step 2
        Label(win, text="Step 2: DFE Configuration", font=("Arial", 11, "bold")).pack(pady=(15, 5))
        btn_auto = tk.Button(win, text="Auto-Detect DFE", command=self.auto_detect_dfe)
        btn_auto.pack()
        Label(win, text="Edit values below (e.g. N, 0, Pi/mu):").pack(pady=2)

        container = Frame(win, height=150)
        container.pack(fill=tk.X, padx=10, pady=5)
        canvas = tk.Canvas(container, height=150)
        scrollbar = Scrollbar(container, orient="vertical", command=canvas.yview)
        self.dfe_frame = Frame(canvas)
        self.dfe_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.dfe_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        self.dfe_entries = {}
        for p in self.pn_logic.places:
            row = Frame(self.dfe_frame)
            row.pack(fill=tk.X, pady=2)
            Label(row, text=f"{p}* =", width=10, anchor="e").pack(side=LEFT)
            ent = Entry(row)
            ent.insert(0, "0")
            ent.pack(side=LEFT, fill=tk.X, expand=True, padx=5)
            self.dfe_entries[p] = ent

        # Step 3
        Label(win, text="Step 3: Results", font=("Arial", 11, "bold")).pack(pady=(15, 5))
        btn_run = tk.Button(win, text="CALCULATE R0", bg="#d0f0c0", font=("Arial", 12, "bold"),
                            command=self.run_calculation_logic)
        btn_run.pack(pady=5)
        self.txt_result = Text(win, height=8, width=50)
        self.txt_result.pack(padx=10, pady=5)

    def auto_detect_dfe(self):
        infected = [p for p, var in self.check_vars.items() if var.get() == 1]
        if not infected:
            messagebox.showerror("Error", "Select infected places first.")
            return
        dfe_dict, status = self.pn_logic.calculate_dfe(infected)
        if dfe_dict is None:
            messagebox.showerror("DFE Error", status)
        else:
            for p, val_str in dfe_dict.items():
                if p in self.dfe_entries:
                    self.dfe_entries[p].delete(0, END)
                    self.dfe_entries[p].insert(0, str(val_str))
            for p in infected:
                if p in self.dfe_entries:
                    self.dfe_entries[p].delete(0, END)
                    self.dfe_entries[p].insert(0, "0")

    def run_calculation_logic(self):
        infected = [p for p, var in self.check_vars.items() if var.get() == 1]
        if not infected:
            messagebox.showerror("Error", "Please select at least one infected place.")
            return
        dfe_manual = {}
        for p, ent in self.dfe_entries.items():
            val = ent.get()
            dfe_manual[p] = val

        self.txt_result.delete("1.0", END)
        self.txt_result.insert(END, f"Infected: {', '.join(infected)}\n")
        self.txt_result.insert(END, "Using DFE:\n")
        for p, v in dfe_manual.items():
            self.txt_result.insert(END, f"  {p}* = {v}\n")

        self.txt_result.insert(END, "\nCalculating R0...\n")
        r0_results = self.pn_logic.calculate_r0(infected, dfe_manual)

        self.txt_result.insert(END, "R0 Result(s):\n")
        for res in r0_results:
            self.txt_result.insert(END, f"  {res}\n")


# --- GLOBAL ERROR HANDLER ---
def report_exception(exc_type, exc_value, exc_traceback):
    """
    Catches any crash in the GUI and shows a popup instead of silently closing.
    """
    error_msg = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
    messagebox.showerror("Critical Error", f"An unexpected error occurred:\n\n{error_msg}")


if __name__ == "__main__":
    root = tk.Tk()
    # Hook the error handler into Tkinter
    root.report_callback_exception = report_exception

    app = PetriNetApp(root)
    root.mainloop()