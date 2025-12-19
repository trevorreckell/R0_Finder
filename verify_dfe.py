from petrinet_logic import VariableWeightPetriNet
import sympy as sp


def verify_model(name, setup_func, infected_list, expected_dfe_str):
    print(f"\n--- Verifying DFE for {name} ---")
    pn = VariableWeightPetriNet()
    setup_func(pn)

    dfe, status = pn.calculate_dfe(infected_list)
    print(f"Status: {status}")
    print(f"Calculated DFE: {dfe}")
    print(f"Expected DFE:   {expected_dfe_str}")


# --- Model Setups ---
def setup_sirs(pn):
    # SIRS (Section 4.1): S -> I -> R -> S
    # Note: Closed system (S+I+R=N). Solver might return S as free variable or N if we define it.
    # Without explicit birth/death, the balance eq for S is: -beta*S*I + delta*R = 0.
    # At I=0, R must be 0. S is free.
    pn.add_place("S");
    pn.add_place("I");
    pn.add_place("R")
    pn.add_transition("t1");
    pn.add_arc("S", "t1", "beta*S*I");
    pn.add_arc("t1", "I", "beta*S*I")
    pn.add_transition("t2");
    pn.add_arc("I", "t2", "gamma*I");
    pn.add_arc("t2", "R", "gamma*I")
    pn.add_transition("t3");
    pn.add_arc("R", "t3", "delta*R");
    pn.add_arc("t3", "S", "delta*R")


def setup_seir(pn):
    # SEIR (Section 4.2): Birth Pi -> S, Death mu*S
    pn.add_place("S");
    pn.add_place("E");
    pn.add_place("I");
    pn.add_place("R")
    pn.add_transition("birth");
    pn.add_arc("birth", "S", "Pi")  # Source transition
    pn.add_transition("death_s");
    pn.add_arc("S", "death_s", "mu*S")
    # Infection loop
    pn.add_transition("inf");
    pn.add_arc("S", "inf", "beta*S*I");
    pn.add_arc("inf", "E", "beta*S*I")
    # ... rest doesn't affect S DFE if I=0 ...


def setup_sveir(pn):
    # SVEIR (Section 4.4):
    # S -> V (phi*V ??? No, phi*V is S->V in paper eq 17? Wait.
    # Eq 16: dS/dt = ... + phi*V - xi*S
    # Eq 17: dV/dt = ... + xi*S - phi*V
    # We need to ensure the arcs match Eq 16/17 exactly for the solver to work.
    pn.add_place("S");
    pn.add_place("V");
    pn.add_place("E");
    pn.add_place("I");
    pn.add_place("R")

    # S <-> V exchange
    pn.add_transition("S_to_V");
    pn.add_arc("S", "S_to_V", "xi_t*S");
    pn.add_arc("S_to_V", "V", "xi_t*S")
    pn.add_transition("V_to_S");
    pn.add_arc("V", "V_to_S", "phi*V");
    pn.add_arc("V_to_S", "S", "phi*V")

    # For DFE solving, we need the "Conservation" N = S + V.
    # The solver might just give S in terms of V if we don't add birth/death or N constraint.
    # Let's add the N definition implicitly or interpret the result.
    # The paper result S* = N * phi / (phi + xi) relies on S+V=N.
    # Our generic solver might return S = V * phi/xi. Let's see.


def setup_nonlinear(pn):
    # Section 4.6: dS/dt = mu - ... - mu*S
    pn.add_place("S");
    pn.add_place("E");
    pn.add_place("I")
    pn.add_transition("birth");
    pn.add_arc("birth", "S", "mu")
    pn.add_transition("death");
    pn.add_arc("S", "death", "mu*S")
    # Infection arc (ignored at I=0)
    pn.add_transition("inf");
    pn.add_arc("S", "inf", "beta*S*I/(1+alpha*I**2)")


# --- Run Verifications ---
if __name__ == "__main__":
    # 1. SIRS
    verify_model("SIRS", setup_sirs, ['I'], "{S: S, I: 0, R: 0} (S is free, user sets S=N)")

    # 2. SEIR
    verify_model("SEIR", setup_seir, ['E', 'I'], "{S: Pi/mu, ...}")

    # 3. Nonlinear
    verify_model("Nonlinear", setup_nonlinear, ['E', 'I'], "{S: 1, ...}")

    # 4. SVEIR
    # Note: Because SVEIR relies on S+V=N (conservation), the solver will likely
    # return the ratio S/V.
    verify_model("SVEIR (Check Ratio)", setup_sveir, ['E', 'I'], "Expect S = V * phi/xi_t")