from petrinet_logic import VariableWeightPetriNet
import sympy as sp

# Initialize the Net
pn = VariableWeightPetriNet()

# 1. Define Places (S, I, R)
pn.add_place("S")
pn.add_place("I")
pn.add_place("R")

# 2. Define Transitions
pn.add_transition("t1") # Infection
pn.add_transition("t2") # Recovery
pn.add_transition("t3") # Loss of immunity

# 3. Define Arcs (From Figure 2a in your paper)
# Arc S -> t1 (weight beta * S * I)
pn.add_arc("S", "t1", "beta * S * I")
# Arc t1 -> I (weight beta * S * I)
pn.add_arc("t1", "I", "beta * S * I")

# Arc I -> t2 (weight gamma * I)
pn.add_arc("I", "t2", "gamma * I")
# Arc t2 -> R (weight gamma * I)
pn.add_arc("t2", "R", "gamma * I")

# Arc R -> t3 (weight delta * R)
pn.add_arc("R", "t3", "delta * R")
# Arc t3 -> S (weight delta * R)
pn.add_arc("t3", "S", "delta * R")

# 4. Define DFE and Infected Set
# At DFE: S = N, I = 0, R = 0
N, beta, gamma, delta = sp.symbols('N beta gamma delta')
dfe = {'S': N, 'I': 0, 'R': 0}
infected_places = ['I']

# 5. Calculate R0
r0 = pn.calculate_r0(infected_places, dfe)

print("Calculated R0:", r0)