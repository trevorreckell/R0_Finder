# $R_0$ Finder

An interactive Python application designed to bridge the gap between continuous compartmental modeling (Ordinary Differential Equations) and stochastic network analysis (Petri Nets). The software automatically maps ODEs to Petri net topologies, resolves biological ambiguities via Abstract Syntax Tree (AST) parsing, and symbolically computes the basic reproduction number ($R_0$) and parameter elasticity indices directly from the network architecture.

## Cross-Platform Compatibility
$R_0$ Finder is designed to run natively on **macOS**, **Windows**, and **Linux**. The graphical user interface is built using `tkinter`, which requires no external installation on Mac or Windows.

##  Installation & Setup

**1. Clone the repository to your local machine:**
```bash
git clone [https://github.com/trevorreckell/R0_Finder.git](https://github.com/trevorreckell/R0_Finder.git)
cd R0_Finder
```

**2. Install the required dependencies:**
Ensure you have Python 3.8+ installed, then run:
```bash
pip install -r requirements.txt
```

**3. Linux Users Only:** While `tkinter` is bundled with Python on macOS and Windows, some Linux distributions require it to be installed separately via the system package manager:
```bash
sudo apt-get install python3-tk
```

##  How to Run the Program

To launch the $R_0$ Finder interface, open your terminal or command prompt, navigate to the folder containing the files, and execute the main script:

```bash
python R0_finder.py
```
*(Depending on your system aliases, you may need to type `python3` instead of `python`).*


## User Guide & Workflow Walkthrough

When you launch the program, the home screen allows you to select your starting model type: **VAPN** (Variable Arc Weight Petri Net), **SPN** (Stochastic Petri Net), or **ODE** (Ordinary Differential Equations).

### 1. Loading and Saving Models
$R_0$ Finder allows you to save and load complete model architectures using `.json` files. 
* To load a pre-built example, click the **Load** button and select a file from the `examples` folder.
* **Important:** Saved models are specific to their workspace. A model saved from the ODE screen can only be loaded as an ODE model.

### 2. Building Models Manually
* **ODE Mode:** Add human compartments and define the continuous rates in the `dx/dt = ...` equation boxes.
* **Petri Net Mode (SPN/VAPN):** Construct graphical topologies manually using the **Place**, **Transition**, and **Arc** buttons. Labels and weights can be edited directly on the canvas.

### 3. Mapping ODEs to Petri Nets
Once an ODE system is loaded or built, you can automatically translate it by clicking **Map to SPN** or **Map to VAPN**.
* **Resolving Ambiguities:** During mapping, the AST parser may trigger popups (e.g., "Ambiguous Proportional Flows" or "Ambiguous Biological Interactions"). Read these carefully. Choosing to "Merge" or "Separate" terms dictates the resulting biological effect and can drastically alter the final transitions, arc structures, and weights.
* **Untangling the Graph:** Newly mapped Petri nets automatically generate in a circular layout. **Tip:** To cleanly untangle the network, evenly space the circular "Places" across a horizontal line first, then move the rectangular "Transitions" into position.

### 4. Calculating $R_0$
$R_0$ Finder calculates the exact symbolic basic reproduction number natively from any of the three model types.
1. Click **Calculate $R_0$** to open the calculation window.
2. **Step 1:** Select the "Infected Variables" (compartments representing infected individuals).
3. **Step 2/3 (Optional):** Define free parameters or constraints, as these can fundamentally alter the equilibrium configuration.
4. Click **Auto-Detect DFE** (Disease-Free Equilibrium). The system will symbolically resolve the null space (e.g., `Pi/mu, 0, 0, 0`). You can manually edit these boxes if desired.
5. Click **Calculate $R_0$** again. 
6. **Classifying Sources:** If the engine detects an ambiguous infection term (e.g., $E \to I$), a popup will force you to classify it as a "New Infection" ($F$ matrix) or "Migration" ($V$ matrix).
7. The software will output the extracted $F$ matrix, $V$ matrix, and the final exact $R_0$ equation.

### 5. Parameter Sensitivity Analysis
Once $R_0$ is calculated, you can instantly extract elasticity indices to evaluate epidemiological interventions.
1. Click the **Sensitivity for $R_0$** button to open the analysis window.
2. Select the specific $R_0$ solution you wish to analyze.
3. Enter numerical values for the parameters. *(Note: These values must be fit to real-world data for the sensitivity analysis to be properly utilized)*.
4. Click **Update Graph** to generate a normalized elasticity bar chart visualizing the impact of each parameter on the epidemic threshold.

---

##  License
This project is open-source and released under the **BSD 3-Clause License**.

##  Citation
If you use $R_0$ Finder in your research, please cite our corresponding methodology paper:
> Reckell, T., Sterner, B., & Jevtić, P. (2026). $R_0$ Finder: An integrated tool for $R_0$ computation and ODE to Petri net model mapping. *PeerJ Computer Science* (Under Review).
