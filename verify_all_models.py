from petrinet_logic import VariableWeightPetriNet
import sympy as sp


def test_seir():
    print("\n--- Testing SEIR Model (Section 4.2) ---")
    pn = VariableWeightPetriNet()
    pn.add_place("S");
    pn.add_place("E");
    pn.add_place("I");
    pn.add_place("R")

    pn.add_transition("t_inf");
    pn.add_arc("S", "t_inf", "beta * S * I");
    pn.add_arc("t_inf", "E", "beta * S * I")
    pn.add_transition("t_inc");
    pn.add_arc("E", "t_inc", "eta * E");
    pn.add_arc("t_inc", "I", "eta * E")
    pn.add_transition("t_death_e");
    pn.add_arc("E", "t_death_e", "mu * E")
    pn.add_transition("t_rec");
    pn.add_arc("I", "t_rec", "alpha * I");
    pn.add_arc("t_rec", "R", "alpha * I")
    pn.add_transition("t_death_i");
    pn.add_arc("I", "t_death_i", "mu * I")

    Pi, mu, beta, eta, alpha = sp.symbols('Pi mu beta eta alpha')
    dfe = {'S': Pi / mu, 'E': 0, 'I': 0, 'R': 0}
    r0 = pn.calculate_r0(['E', 'I'], dfe)
    print(f"Calculated R0: {r0}")
    print("Look for:      [beta * Pi * eta / (mu * (alpha + mu) * (eta + mu))]")


def test_seeir():
    print("\n--- Testing SEEIR Model (Section 4.3) ---")
    pn = VariableWeightPetriNet()
    for p in ["S", "E1", "E2", "I", "R"]: pn.add_place(p)

    pn.add_transition("t_inf1");
    pn.add_arc("S", "t_inf1", "p * beta * S * I / N");
    pn.add_arc("t_inf1", "E1", "p * beta * S * I / N")
    pn.add_transition("t_inf2");
    pn.add_arc("S", "t_inf2", "(1-p) * beta * S * I / N");
    pn.add_arc("t_inf2", "E2", "(1-p) * beta * S * I / N")
    pn.add_transition("t_e1_i");
    pn.add_arc("E1", "t_e1_i", "nu1 * E1");
    pn.add_arc("t_e1_i", "I", "nu1 * E1")
    pn.add_transition("t_d1");
    pn.add_arc("E1", "t_d1", "mu * E1")
    pn.add_transition("t_e2_i");
    pn.add_arc("E2", "t_e2_i", "nu2 * E2");
    pn.add_arc("t_e2_i", "I", "nu2 * E2")
    pn.add_transition("t_d2");
    pn.add_arc("E2", "t_d2", "mu * E2")
    pn.add_transition("t_rec");
    pn.add_arc("I", "t_rec", "gamma * I");
    pn.add_arc("t_rec", "R", "gamma * I")
    pn.add_transition("t_di");
    pn.add_arc("I", "t_di", "mu * I")

    N, beta, mu, gamma, p, nu1, nu2 = sp.symbols('N beta mu gamma p nu1 nu2')
    dfe = {'S': N, 'E1': 0, 'E2': 0, 'I': 0, 'R': 0}
    r0 = pn.calculate_r0(['E1', 'E2', 'I'], dfe)
    print(f"Calculated R0: {r0}")
    print("Look for:      beta * (p*nu1/(mu+nu1) + (1-p)*nu2/(mu+nu2)) / (gamma+mu)")


def test_sveir():
    print("\n--- Testing SVEIR Model (Section 4.4) ---")
    pn = VariableWeightPetriNet()
    for p in ["S", "V", "E", "I", "R"]: pn.add_place(p)

    # Note: I -> Out rate is (alpha_rec + delta) based on the paper's V matrix
    pn.add_transition("t_inf_s");
    pn.add_arc("S", "t_inf_s", "beta_t * I * S / N");
    pn.add_arc("t_inf_s", "E", "beta_t * I * S / N")
    pn.add_transition("t_inf_v");
    pn.add_arc("V", "t_inf_v", "sigma * beta_t * I * V / N");
    pn.add_arc("t_inf_v", "E", "sigma * beta_t * I * V / N")
    pn.add_transition("t_inc");
    pn.add_arc("E", "t_inc", "alpha * E");
    pn.add_arc("t_inc", "I", "alpha * E")
    pn.add_transition("t_out_i");
    pn.add_arc("I", "t_out_i", "(alpha_rec + delta) * I")

    N, S_star, V_star, beta_t, sigma, alpha, delta, alpha_rec = sp.symbols(
        'N S_star V_star beta_t sigma alpha delta alpha_rec')
    dfe = {'S': S_star, 'V': V_star, 'E': 0, 'I': 0}

    r0 = pn.calculate_r0(['E', 'I'], dfe)
    print(f"Calculated R0: {sp.factor(r0[1])}")  # Index 1 usually holds the non-zero eigenvalue here
    print("Look for:      beta_t * (S_star + V_star*sigma) / (N * (alpha_rec + delta))")


def test_nonlinear():
    print("\n--- Testing Nonlinear Model (Section 4.6) ---")
    pn = VariableWeightPetriNet()
    pn.add_place("S");
    pn.add_place("E");
    pn.add_place("I")

    pn.add_transition("t_inf");
    pn.add_arc("S", "t_inf", "beta * S * I / (1 + alpha_nl * I**2)");
    pn.add_arc("t_inf", "E", "beta * S * I / (1 + alpha_nl * I**2)")
    pn.add_transition("t_inc");
    pn.add_arc("E", "t_inc", "sigma * E");
    pn.add_arc("t_inc", "I", "sigma * E")
    pn.add_transition("t_de");
    pn.add_arc("E", "t_de", "mu * E")
    pn.add_transition("t_out");
    pn.add_arc("I", "t_out", "gamma * I + mu * I")

    dfe = {'S': 1, 'E': 0, 'I': 0}
    r0 = pn.calculate_r0(['E', 'I'], dfe)
    print(f"Calculated R0: {sp.simplify(r0[1])}")
    print("Look for:      beta * sigma / ((gamma + mu) * (mu + sigma))")


def test_vector():
    print("\n--- Testing SIR Vector-Borne (Section 4.8) ---")
    pn = VariableWeightPetriNet()
    pn.add_place("Sh");
    pn.add_place("Ih");
    pn.add_place("Sv");
    pn.add_place("Iv")

    pn.add_transition("t_inf_h");
    pn.add_arc("Sh", "t_inf_h", "beta_hv * Sh * Iv");
    pn.add_arc("t_inf_h", "Ih", "beta_hv * Sh * Iv")
    pn.add_transition("t_inf_v");
    pn.add_arc("Sv", "t_inf_v", "beta_vh * Sv * Ih");
    pn.add_arc("t_inf_v", "Iv", "beta_vh * Sv * Ih")
    pn.add_transition("t_out_h");
    pn.add_arc("Ih", "t_out_h", "(alpha + mu_h + sigma - delta) * Ih")
    pn.add_transition("t_out_v");
    pn.add_arc("Iv", "t_out_v", "mu_v * Iv")

    Pi, mu_h, Lambda, mu_v = sp.symbols('Pi mu_h Lambda mu_v')
    dfe = {'Sh': Pi / mu_h, 'Sv': Lambda / mu_v, 'Ih': 0, 'Iv': 0}

    r0 = pn.calculate_r0(['Ih', 'Iv'], dfe)
    print(f"Calculated R0: {r0}")
    print("Look for:      sqrt((beta_hv*Pi)/(mu_h*mu_v) * (beta_vh*Lambda)/(mu_v*(alpha + mu_h + sigma - delta)))")


if __name__ == "__main__":
    test_seir()
    test_seeir()
    test_sveir()
    test_nonlinear()
    test_vector()

