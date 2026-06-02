import sympy as sp
import re
import itertools
from collections import defaultdict
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Greek Letters
LATEX_TO_UNICODE = {
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


def encode_latex(text):
    """Converts \beta to greek_beta for internal symbolic safety."""
    if not isinstance(text, str): return text
    return re.sub(r'\\([a-zA-Z]+)', r'greek_\1', text)


def decode_latex(text):
    """Converts greek_beta to β for display."""
    if not isinstance(text, str): return text

    def replace_match(match):
        word = match.group(1)
        return LATEX_TO_UNICODE.get(word, '\\' + word)

    return re.sub(r'greek_([a-zA-Z]+)', replace_match, text)


def format_expr_for_display(expr, variables):
    if expr is None: return ""
    expr_str = str(expr)
    targets = set(variables)
    targets.add('N')
    sorted_targets = sorted(list(targets), key=len, reverse=True)
    for var in sorted_targets:
        pattern = r'(?<![a-zA-Z0-9_])' + re.escape(var) + r'(?![a-zA-Z0-9_])'
        expr_str = re.sub(pattern, f"{var}^*", expr_str)
    return decode_latex(expr_str)


def sanitize_input_str(input_str):
    if not input_str: return "0"
    return input_str.replace("^*", "")


def calculate_sensitivity_indices(r0_expr):
    """
    Calculates the Normalized Forward Sensitivity Index for all free symbols in r0_expr.
    Formula: (p / R0) * (dR0 / dp)
    """
    if r0_expr is None or r0_expr == 0: return {}
    indices = {}
    params = r0_expr.free_symbols
    for p in params:
        try:
            derivative = sp.diff(r0_expr, p)
            elasticity = (p / r0_expr) * derivative
            indices[str(p)] = sp.simplify(elasticity)
        except:
            indices[str(p)] = sp.sympify(0)
    return indices

# Phase 1 and 2: SPN

class StochasticPetriNet:
    def __init__(self):
        self.places = []
        self.transitions = {}  # Maps name , rate_function (string)
        self.arcs = []
        self.symbols = {}

    def add_place(self, name):
        name = encode_latex(name)
        if name not in self.places:
            self.places.append(name)
            self.symbols[name] = sp.Symbol(name, real=True)

    def remove_place(self, name):
        name = encode_latex(name)
        if name in self.places:
            self.places.remove(name)
            if name in self.symbols: del self.symbols[name]
            self.arcs = [a for a in self.arcs if a['source'] != name and a['target'] != name]

    def add_transition(self, name, rate_function="1.0"):
        name = encode_latex(name)
        rate_function = encode_latex(rate_function)
        self.transitions[name] = rate_function
        self._register_new_symbols(rate_function)

    def update_transition_rate(self, name, new_rate):
        name = encode_latex(name)
        new_rate = encode_latex(new_rate)
        if name in self.transitions:
            self.transitions[name] = new_rate
            self._register_new_symbols(new_rate)

    def remove_transition(self, name):
        name = encode_latex(name)
        if name in self.transitions:
            del self.transitions[name]
            self.arcs = [a for a in self.arcs if a['source'] != name and a['target'] != name]

    def add_arc(self, source, target, weight=1):
        source = encode_latex(source)
        target = encode_latex(target)
        weight = encode_latex(weight)
        self._register_new_symbols(str(weight))
        self.arcs.append({'source': source, 'target': target, 'weight': weight})

    def remove_arc(self, source, target):
        source = encode_latex(source)
        target = encode_latex(target)
        self.arcs = [a for a in self.arcs if not (a['source'] == source and a['target'] == target)]

    def _register_new_symbols(self, expression_str):
        potential_symbols = re.findall(r'[a-zA-Z_][a-zA-Z0-9_]*', str(expression_str))
        for name in potential_symbols:
            if name not in self.symbols:
                self.symbols[name] = sp.Symbol(name, real=True)

    def get_input_arcs(self, node):
        return [a for a in self.arcs if a['target'] == node]

    def get_output_arcs(self, node):
        return [a for a in self.arcs if a['source'] == node]

    def get_net_flow(self, trans_name, place_name):
        flow = 0
        for arc in self.arcs:
            if arc['source'] == trans_name and arc['target'] == place_name:
                try:
                    w = sp.sympify(arc['weight'], locals=self.symbols)
                    flow += w
                except:
                    pass
        for arc in self.arcs:
            if arc['source'] == place_name and arc['target'] == trans_name:
                try:
                    w = sp.sympify(arc['weight'], locals=self.symbols)
                    flow -= w
                except:
                    pass
        return flow

    def calculate_dfe(self, infected_places, free_places=[], constraints=[]):
        P_N = [p for p in self.places if p not in infected_places]
        sym_N = self.symbols.get('N', sp.Symbol('N', real=True))

        for c in constraints: self._register_new_symbols(c)

        base_equations = []
        for place in P_N:
            net_rate_change = 0
            for t_name, rate_str in self.transitions.items():
                net_flow_stoich = self.get_net_flow(t_name, place)
                if net_flow_stoich != 0:
                    try:
                        rate_expr = sp.sympify(rate_str, locals=self.symbols)
                        net_rate_change += rate_expr * net_flow_stoich
                    except:
                        pass
            base_equations.append(net_rate_change)

        constraint_subs = {}
        remaining_constraints = []
        for c in constraints:
            if "=" in c:
                lhs_str, rhs_str = c.split("=")
                try:
                    lhs = sp.sympify(lhs_str, locals=self.symbols)
                    rhs = sp.sympify(rhs_str, locals=self.symbols)
                    eqn = lhs - rhs
                    syms = list(eqn.free_symbols)
                    if syms:
                        target = syms[0]
                        res = sp.solve(eqn, target)
                        if res:
                            constraint_subs[target] = res[0]
                        else:
                            remaining_constraints.append(eqn)
                    else:
                        remaining_constraints.append(eqn)
                except:
                    pass
            else:
                try:
                    remaining_constraints.append(sp.sympify(c, locals=self.symbols))
                except:
                    pass

        subs_I_zero = {self.symbols[p]: 0 for p in infected_places}
        eq_system = [eq.subs(constraint_subs).subs(subs_I_zero) for eq in base_equations]
        for eq in remaining_constraints:
            eq_system.append(eq.subs(constraint_subs).subs(subs_I_zero))

        vars_to_solve = [self.symbols[p] for p in P_N if p not in free_places]

        solutions = []
        try:
            solutions = sp.solve(eq_system, vars_to_solve, dict=True)
        except:
            solutions = []

        # validation logic
        valid = False
        if solutions:
            for sol in solutions:
                # 1. check for triviality (all 0s)
                is_zero = True
                for k, v in sol.items():
                    if v != 0: is_zero = False

                # 2. check for underdetermined System (aka dependency)
                # If a solution value depends on a variable we are trying to solve for (like S = k*V), the system is underdetermined and needs conservation law.
                is_dependent = False
                for val_expr in sol.values():
                    if hasattr(val_expr, 'free_symbols'):
                        # If the expression contains any symbol from vars_to_solve, it's dependent
                        if not set(vars_to_solve).isdisjoint(val_expr.free_symbols):
                            is_dependent = True
                            break

                # If solution implies free variables (missing keys), it's dependent
                if len(sol) < len(vars_to_solve):
                    is_dependent = True

                # Only accept if non-zero and fully resolved
                if not is_zero and not is_dependent:
                    valid = True
                    break

        if not valid:
            # Fallback: enforce conservation law (Sum(Places) = N)
            total_tokens = sum(self.symbols[p] for p in self.places)
            cons_eq = total_tokens.subs(subs_I_zero) - sym_N
            cons_eq = cons_eq.subs(constraint_subs)
            try:
                # Append conservation equation to the system
                solutions = sp.solve(eq_system + [cons_eq], vars_to_solve, dict=True)
            except:
                pass

        final_solutions = []
        if not solutions: solutions = [{}]

        for sol in solutions:
            full = sol.copy()
            full.update(subs_I_zero)
            for fp in free_places:
                if self.symbols[fp] not in full: full[self.symbols[fp]] = self.symbols[fp]

            str_sol = {}
            for k, v in full.items():
                try:
                    simplified_val = sp.simplify(v)
                    str_sol[k.name] = format_expr_for_display(simplified_val, self.places)
                except:
                    str_sol[k.name] = str(v)
            if str_sol and str_sol not in final_solutions: final_solutions.append(str_sol)

        return final_solutions, "Success"

    def calculate_r0(self, infected_places_names, dfe_substitutions_str, source_classification=None):
        if source_classification is None: source_classification = {}

        parse_locals = self.symbols.copy()
        parse_locals['N'] = sp.Symbol('N', real=True)
        dfe_substitutions = {}
        for k, v in dfe_substitutions_str.items():
            sym_key = self.symbols[k] if k in self.symbols else sp.Symbol(k, real=True)
            v_clean = sanitize_input_str(str(v))
            try:
                val_expr = sp.sympify(v_clean, locals=parse_locals)
            except:
                val_expr = sp.sympify(v_clean)
            dfe_substitutions[sym_key] = val_expr

        # Calculate N*
        N_val = sum(dfe_substitutions.values())
        sym_N = sp.Symbol('N', real=True)
        dfe_substitutions[sym_N] = N_val

        N_str_parts = [f"{p}^*" for p in self.places]
        N_def_lhs = " + ".join(N_str_parts)
        N_def_rhs = format_expr_for_display(N_val, self.places)
        N_def_str = f"N^* = {N_def_lhs} = {N_def_rhs}"

        P_I = infected_places_names
        P_N = [p for p in self.places if p not in P_I]

        T_inf = []
        ambiguous_transitions = []

        for t in self.transitions:
            inputs = [arc['source'] for arc in self.get_input_arcs(t)]
            outputs = [arc['target'] for arc in self.get_output_arcs(t)]
            has_susceptible_input = any(p in P_N for p in inputs)
            has_infected_output = any(p in P_I for p in outputs)
            is_source_like = (len(inputs) == 0) or all(p in P_I for p in inputs)

            if has_susceptible_input and has_infected_output:
                T_inf.append(t)
            elif is_source_like and has_infected_output:
                user_decision = source_classification.get(t)
                if user_decision == 'F':
                    T_inf.append(t)
                elif user_decision == 'V':
                    pass
                else:
                    rate_expr = self.transitions.get(t, "Unknown")
                    ambiguous_transitions.append({'name': t, 'rate': rate_expr})

        if ambiguous_transitions:
            return None, "ambiguous", ambiguous_transitions, [], [], None

        if not T_inf and not P_I:
            return ["Error: No Infection Transitions found."], None, None, [], [], None

        F_matrix_symbolic = []
        for i_place in P_I:
            rate_new_infection = 0
            for t in T_inf:
                t_rate_str = self.transitions.get(t, "0")
                try:
                    t_rate_expr = sp.sympify(t_rate_str, locals=self.symbols)
                except:
                    t_rate_expr = sp.Integer(0)

                w_in = 0
                for arc in self.get_input_arcs(t):
                    if arc['source'] == i_place:
                        try:
                            w = sp.sympify(arc['weight'], locals=self.symbols)
                        except:
                            w = 1
                        w_in += w

                w_out = 0
                for arc in self.get_output_arcs(t):
                    if arc['target'] == i_place:
                        try:
                            w = sp.sympify(arc['weight'], locals=self.symbols)
                        except:
                            w = 1
                        w_out += w

                net_gain = w_out - w_in
                if sp.simplify(net_gain) != 0:
                    try:
                        coeff = sp.simplify(net_gain).as_coeff_Mul()[0]
                        if coeff > 0 or not coeff.is_Number:
                            rate_new_infection += t_rate_expr * net_gain
                    except:
                        rate_new_infection += t_rate_expr * net_gain

            F_matrix_symbolic.append(rate_new_infection)

        V_matrix_symbolic = []
        for i_place in P_I:
            v_minus = 0
            v_plus = 0

            for t in self.transitions:
                if t in T_inf: continue

                t_rate = self.transitions.get(t, "0")
                try:
                    t_expr = sp.sympify(t_rate, locals=self.symbols)
                except:
                    t_expr = 0

                w_in = 0
                for arc in self.get_input_arcs(t):
                    if arc['source'] == i_place:
                        try:
                            w_in += sp.sympify(arc['weight'], locals=self.symbols)
                        except:
                            w_in += 1

                w_out = 0
                for arc in self.get_output_arcs(t):
                    if arc['target'] == i_place:
                        try:
                            w_out += sp.sympify(arc['weight'], locals=self.symbols)
                        except:
                            w_out += 1

                if sp.simplify(w_in) != 0: v_minus += t_expr * w_in
                if sp.simplify(w_out) != 0: v_plus += t_expr * w_out

            for t in T_inf:
                t_rate = self.transitions.get(t, "0")
                try:
                    t_expr = sp.sympify(t_rate, locals=self.symbols)
                except:
                    t_expr = 0
                w_in = 0
                for arc in self.get_input_arcs(t):
                    if arc['source'] == i_place:
                        try:
                            w_in += sp.sympify(arc['weight'], locals=self.symbols)
                        except:
                            w_in += 1
                w_out = 0
                for arc in self.get_output_arcs(t):
                    if arc['target'] == i_place:
                        try:
                            w_out += sp.sympify(arc['weight'], locals=self.symbols)
                        except:
                            w_out += 1
                net_change = w_out - w_in
                if sp.simplify(net_change) != 0:
                    try:
                        coeff = sp.simplify(net_change).as_coeff_Mul()[0]
                        if coeff < 0: v_minus += t_expr * abs(net_change)
                    except:
                        pass

            V_matrix_symbolic.append(v_minus - v_plus)

        try:
            infected_symbols = [self.symbols[name] for name in P_I]
            F_mat = sp.Matrix(F_matrix_symbolic).jacobian(infected_symbols).subs(dfe_substitutions)
            V_mat = sp.Matrix(V_matrix_symbolic).jacobian(infected_symbols).subs(dfe_substitutions)
            eigs = (F_mat * V_mat.inv()).eigenvals()
            r0_list = list(eigs.keys())
            r0_fmt = [format_expr_for_display(r, self.places) for r in r0_list]
            F_fmt = format_expr_for_display(F_mat, self.places)
            V_fmt = format_expr_for_display(V_mat, self.places)
            sensitivity_data = []
            for r_expr in r0_list:
                sensitivity_data.append(calculate_sensitivity_indices(r_expr))

            #  Return 6 items, for error
            return r0_fmt, F_fmt, V_fmt, sensitivity_data, r0_list, N_def_str
        except Exception as e:
            return [f"Math Error: {e}"], None, None, [], [], None


#    Petri net logic (VAPN)

class VariableWeightPetriNet:
    def __init__(self):
        self.places = []
        self.transitions = []
        self.arcs = []
        self.symbols = {}

    def add_place(self, name):
        name = encode_latex(name)
        if name not in self.places:
            self.places.append(name)
            self.symbols[name] = sp.Symbol(name, real=True)

    def remove_place(self, name):
        name = encode_latex(name)
        if name in self.places:
            self.places.remove(name)
            if name in self.symbols: del self.symbols[name]
            self.arcs = [arc for arc in self.arcs if arc['source'] != name and arc['target'] != name]

    def add_transition(self, name):
        name = encode_latex(name)
        if name not in self.transitions: self.transitions.append(name)

    def remove_transition(self, name):
        name = encode_latex(name)
        if name in self.transitions:
            self.transitions.remove(name)
            self.arcs = [arc for arc in self.arcs if arc['source'] != name and arc['target'] != name]

    def _register_new_symbols(self, expression_str):
        potential_symbols = re.findall(r'[a-zA-Z_][a-zA-Z0-9_]*', expression_str)
        for name in potential_symbols:
            if name not in self.symbols: self.symbols[name] = sp.Symbol(name, real=True)

    def add_arc(self, source, target, weight_expr):
        weight_expr = encode_latex(weight_expr)
        self._register_new_symbols(str(weight_expr))
        try:
            expr = sp.sympify(weight_expr, locals=self.symbols)
        except:
            expr = sp.sympify(weight_expr)
        self.arcs.append({'source': source, 'target': target, 'weight_str': str(weight_expr), 'weight_expr': expr})

    def remove_arc(self, source, target):
        source = encode_latex(source)
        target = encode_latex(target)
        self.arcs = [arc for arc in self.arcs if not (arc['source'] == source and arc['target'] == target)]

    def get_input_arcs(self, node):
        return [arc for arc in self.arcs if arc['target'] == node]

    def get_output_arcs(self, node):
        return [arc for arc in self.arcs if arc['source'] == node]

    def calculate_dfe(self, infected_places_names, free_places_names=[], constraints=[]):
        P_I = infected_places_names
        P_N = [p for p in self.places if p not in P_I]
        sym_N = sp.Symbol('N', real=True)
        subs_I_zero = {self.symbols[p]: 0 for p in P_I}
        equations = []
        vars_to_solve = [self.symbols[p] for p in P_N if p not in free_places_names]

        for c in constraints: self._register_new_symbols(c)

        for p in P_N:
            if p in free_places_names: continue
            inflow = sum(arc['weight_expr'].subs(subs_I_zero) for arc in self.get_input_arcs(p))
            outflow = sum(arc['weight_expr'].subs(subs_I_zero) for arc in self.get_output_arcs(p))
            equations.append(inflow - outflow)

        constraint_subs = {}
        remaining_constraints = []
        for c in constraints:
            if "=" in c:
                lhs_str, rhs_str = c.split("=")
                try:
                    lhs = sp.sympify(lhs_str, locals=self.symbols)
                    rhs = sp.sympify(rhs_str, locals=self.symbols)
                    eqn = lhs - rhs
                    syms = list(eqn.free_symbols)
                    if syms:
                        target = syms[0]
                        res = sp.solve(eqn, target)
                        if res:
                            constraint_subs[target] = res[0]
                        else:
                            remaining_constraints.append(eqn)
                    else:
                        remaining_constraints.append(eqn)
                except:
                    pass
            else:
                try:
                    remaining_constraints.append(sp.sympify(c, locals=self.symbols))
                except:
                    pass

        final_eqs = [eq.subs(constraint_subs) for eq in equations] + [eq.subs(constraint_subs) for eq in
                                                                      remaining_constraints]

        try:
            solutions = sp.solve(final_eqs, vars_to_solve, dict=True)
        except Exception as e:
            return None, f"Solver Error: {e}"

        final_solutions = []
        if not solutions: solutions = [{}]
        for sol in solutions:
            free_vars = [v for v in vars_to_solve if v not in sol]
            current_sol_list = [sol]
            if free_vars:
                sum_PN = sum(vars_to_solve) + sum(self.symbols[fp] for fp in free_places_names)
                constraint = sum_PN - sym_N
                try:
                    constrained_sols = sp.solve(final_eqs + [constraint], vars_to_solve, dict=True)
                    if constrained_sols: current_sol_list = constrained_sols
                except:
                    pass

            for s in current_sol_list:
                full_dfe_sym = s.copy()
                full_dfe_sym.update(subs_I_zero)
                for free_p in free_places_names:
                    sym = self.symbols[free_p]
                    if sym not in full_dfe_sym: full_dfe_sym[sym] = sym
                full_dfe_str = {}
                for k, v in full_dfe_sym.items():
                    simplified_val = sp.simplify(v)
                    full_dfe_str[k.name] = format_expr_for_display(simplified_val, self.places)
                if full_dfe_str not in final_solutions: final_solutions.append(full_dfe_str)

        if not final_solutions: return [], "No valid DFE found."
        return final_solutions, "Success"

    def calculate_r0(self, infected_places_names, dfe_substitutions_str, source_classification=None):
        if source_classification is None: source_classification = {}

        parse_locals = self.symbols.copy()
        parse_locals['N'] = sp.Symbol('N', real=True)
        dfe_substitutions = {}
        for k, v in dfe_substitutions_str.items():
            sym_key = self.symbols[k] if k in self.symbols else sp.Symbol(k, real=True)
            v_clean = sanitize_input_str(str(v))
            try:
                val_expr = sp.sympify(v_clean, locals=parse_locals)
            except:
                val_expr = sp.sympify(v_clean)
            dfe_substitutions[sym_key] = val_expr

        # Calculate N*
        N_val = sum(dfe_substitutions.values())
        sym_N = sp.Symbol('N', real=True)
        dfe_substitutions[sym_N] = N_val

        N_str_parts = [f"{p}^*" for p in self.places]
        N_def_lhs = " + ".join(N_str_parts)
        N_def_rhs = format_expr_for_display(N_val, self.places)
        N_def_str = f"N^* = {N_def_lhs} = {N_def_rhs}"

        P_I = infected_places_names
        P_N = [p for p in self.places if p not in P_I]

        T_inf = []
        ambiguous_transitions = []

        for t in self.transitions:
            inputs = [arc['source'] for arc in self.get_input_arcs(t)]
            outputs = [arc['target'] for arc in self.get_output_arcs(t)]
            has_susceptible_input = any(p in P_N for p in inputs)
            has_infected_output = any(p in P_I for p in outputs)
            is_source_like = (len(inputs) == 0) or all(p in P_I for p in inputs)

            if has_susceptible_input and has_infected_output:
                T_inf.append(t)
            elif is_source_like and has_infected_output:
                user_decision = source_classification.get(t)
                if user_decision == 'F':
                    T_inf.append(t)
                elif user_decision == 'V':
                    pass
                else:
                    rate_expr = "Unknown"
                    for arc in self.get_output_arcs(t):
                        if arc['target'] in P_I: rate_expr = arc['weight_str']; break
                    ambiguous_transitions.append({'name': t, 'rate': rate_expr})

        if ambiguous_transitions:
            return None, "ambiguous", ambiguous_transitions, [], [], None

        if not T_inf and not P_I:
            return ["Error: No Infection Transitions found."], None, None, [], [], None

        F_matrix_symbolic = []
        for i_place in P_I:
            rate_new_infection = 0
            for t in T_inf:
                w_in = 0
                for arc in self.get_input_arcs(t):
                    if arc['source'] == i_place: w_in += arc['weight_expr']
                w_out = 0
                for arc in self.get_output_arcs(t):
                    if arc['target'] == i_place: w_out += arc['weight_expr']
                net = w_out - w_in
                if sp.simplify(net) != 0: rate_new_infection += net
            F_matrix_symbolic.append(rate_new_infection)

        V_matrix_symbolic = []
        for i_place in P_I:
            outflow = 0
            for arc in self.get_output_arcs(i_place):
                outflow += arc['weight_expr']
            inflow = 0
            for t in self.transitions:
                if t in T_inf: continue
                for arc in self.get_output_arcs(t):
                    if arc['target'] == i_place: inflow += arc['weight_expr']
            V_matrix_symbolic.append(outflow - inflow)

        try:
            infected_symbols = [self.symbols[name] for name in P_I]
            F_mat = sp.Matrix(F_matrix_symbolic).jacobian(infected_symbols).subs(dfe_substitutions)
            V_mat = sp.Matrix(V_matrix_symbolic).jacobian(infected_symbols).subs(dfe_substitutions)
            eigs = (F_mat * V_mat.inv()).eigenvals()
            r0_list = list(eigs.keys())
            r0_fmt = [format_expr_for_display(r, self.places) for r in r0_list]
            F_fmt = format_expr_for_display(F_mat, self.places)
            V_fmt = format_expr_for_display(V_mat, self.places)
            sensitivity_data = []
            for r_expr in r0_list:
                sensitivity_data.append(calculate_sensitivity_indices(r_expr))

            #  Return 6 items
            return r0_fmt, F_fmt, V_fmt, sensitivity_data, r0_list, N_def_str
        except Exception as e:
            return [f"Math Error: {e}"], None, None, [], [], None

    # hybrid conversion VAPN
    def convert_to_vapn_structure(self, use_scaling=True, merge_strategies=None):
            if merge_strategies is None: merge_strategies = {}
            outflows = []
            inflows = []
            state_syms = {self.symbols[v] for v in self.variables}

            # 1. Parse Equations
            for var, eq_str in self.equations.items():
                if not eq_str.strip(): continue
                try:
                    expr = sp.sympify(eq_str, locals=self.symbols).expand()
                except:
                    continue
                for term in expr.as_ordered_terms():
                    coeff = term.as_coeff_Mul()[0]
                    if coeff < 0:
                        outflows.append({'place': var, 'expr': -term, 'matched': False})
                    elif term.is_Mul and term.args[0].is_Number and term.args[0] < 0:
                        outflows.append({'place': var, 'expr': -term, 'matched': False})
                    else:
                        inflows.append({'place': var, 'expr': term, 'matched': False})

            # Phase 0: apply strategies and fix coupled terms
            initial_strategies = list(merge_strategies.items())

            for strategy_key, decision in initial_strategies:
                parts = strategy_key.split("|")
                if len(parts) != 3: continue
                var_name, in_str, out_str = parts

                # A. Find Target Terms
                target_in_expr = None
                target_out_expr = None
                try:
                    target_in_expr = sp.sympify(in_str, locals=self.symbols)
                    target_out_expr = sp.sympify(out_str, locals=self.symbols)
                except:
                    pass

                t_in = None
                for i in inflows:
                    if i['place'] == var_name:
                        if str(i['expr']) == in_str: t_in = i; break
                        if target_in_expr is not None and sp.simplify(i['expr'] - target_in_expr) == 0: t_in = i; break

                t_out = None
                for o in outflows:
                    if o['place'] == var_name:
                        if str(o['expr']) == out_str: t_out = o; break
                        if target_out_expr is not None and sp.simplify(
                            o['expr'] - target_out_expr) == 0: t_out = o; break

                if not (t_in and t_out): continue

                # B. Apply Decision
                actual_t_in = t_in['expr']
                actual_t_out = t_out['expr']
                net_rate = sp.simplify(actual_t_out - actual_t_in)

                if decision == "merge":
                    t_out['expr'] = net_rate
                    t_in['matched'] = True

                # C. Auto-Fix Partner Variable
                for other_in in inflows:
                    if other_in['matched'] or other_in['place'] == var_name: continue

                    is_direct_match = False
                    if str(other_in['expr']) == str(actual_t_out):
                        is_direct_match = True
                    elif sp.simplify(other_in['expr'] - actual_t_out) == 0:
                        is_direct_match = True

                    ratio = sp.Integer(1)
                    if not is_direct_match:
                        ratio = sp.simplify(other_in['expr'] / actual_t_out)
                        if ratio.free_symbols & state_syms: continue

                    expected_out = sp.simplify(actual_t_in * ratio)
                    partner_out = None
                    for other_out in outflows:
                        if other_out['matched'] or other_out['place'] != other_in['place']: continue
                        if str(other_out['expr']) == str(expected_out): partner_out = other_out; break
                        if sp.simplify(other_out['expr'] - expected_out) == 0: partner_out = other_out; break

                    if partner_out:
                        p_key = f"{other_in['place']}|{str(other_in['expr'])}|{str(partner_out['expr'])}"

                        if decision == "merge":
                            other_in['expr'] = sp.simplify(ratio * net_rate)
                            partner_out['matched'] = True

                        elif decision == "separate":
                            if p_key not in merge_strategies:
                                merge_strategies[p_key] = "separate"
                        break


            #critical initialization block  (probably best to leave untouched)
            pn_places = [{"name": v, "type": "place", "x": 0, "y": 0} for v in self.variables]
            pn_transitions = []
            pn_arcs = []
            t_counter = 1

            #  Phase 0.5: ambiguity detection
            ambiguities = []
            for in_item in inflows:
                if in_item['matched']: continue

                # Check 1: Catalyst Shared Symbols
                in_syms = in_item['expr'].free_symbols & state_syms
                cross_syms = {s for s in in_syms if str(s) != in_item['place']}
                if cross_syms:
                    for out_item in outflows:
                        if out_item['place'] != in_item['place']: continue
                        if out_item['matched']: continue

                        out_syms = out_item['expr'].free_symbols & state_syms
                        shared = cross_syms.intersection(out_syms)
                        if shared:
                            in_expr = str(in_item['expr'])
                            out_expr = str(out_item['expr'])
                            key = f"{in_item['place']}|{in_expr}|{out_expr}"
                            if key not in merge_strategies:
                                ambiguities.append({
                                    "key": key, "place": in_item['place'], "type": "catalyst",
                                    "inflow_expr": in_expr, "outflow_expr": out_expr, "catalyst": str(list(shared)[0])
                                })

                # Check 2: Proportional Self-Loops
                if use_scaling:
                    for out_item in outflows:
                        if out_item['place'] != in_item['place']: continue  # Must be Self-Loop
                        if out_item['matched']: continue

                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                        # flag if not identical and proportional
                        if ratio != 1 and not (ratio.free_symbols & state_syms):
                            in_expr = str(in_item['expr'])
                            out_expr = str(out_item['expr'])
                            key = f"{in_item['place']}|{in_expr}|{out_expr}"

                            if key not in merge_strategies:
                                ambiguities.append({
                                    "key": key, "place": in_item['place'], "type": "loop",
                                    "inflow_expr": in_expr, "outflow_expr": out_expr, "catalyst": "Self"
                                })

            if ambiguities:
                return {"status": "ambiguous", "details": ambiguities}

            #  Step 1: generate transitions
            if use_scaling:
                # 1a. Identical Matches (Safe Auto-Merge)
                for out_item in outflows:
                    if out_item['matched']: continue
                    for in_item in inflows:
                        if in_item['matched']: continue

                        is_identical = False
                        if str(in_item['expr']) == str(out_item['expr']):
                            is_identical = True
                        elif sp.simplify(in_item['expr'] - out_item['expr']) == 0:
                            is_identical = True

                        if is_identical:
                            t_name = f"t{t_counter}"
                            t_counter += 1
                            pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})

                            w_in = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(out_item['expr']))
                            pn_arcs.append({"source": out_item['place'], "target": t_name, "weight": w_in})
                            out_item['matched'] = True

                            w_out = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(in_item['expr']))
                            pn_arcs.append({"source": t_name, "target": in_item['place'], "weight": w_out})
                            in_item['matched'] = True
                            break

                # 1b. Proportional Matches (Filtered)
                for out_item in outflows:
                    if out_item['matched']: continue
                    matched_inflows = []
                    for in_item in inflows:
                        if in_item['matched']: continue

                        # Strict strategy check (Separate)
                        key_check = f"{in_item['place']}|{str(in_item['expr'])}|{str(out_item['expr'])}"
                        if merge_strategies.get(key_check) == "separate": continue

                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                        if not (ratio.free_symbols & state_syms):
                            # Prevent proportional self-merging unless explicit
                            if in_item['place'] == out_item['place']:
                                if merge_strategies.get(key_check) != "merge":
                                    continue
                            matched_inflows.append({'item': in_item, 'ratio': ratio})

                    if matched_inflows:
                        t_name = f"t{t_counter}"
                        t_counter += 1
                        pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})

                        w_in = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(out_item['expr']))
                        pn_arcs.append({"source": out_item['place'], "target": t_name, "weight": w_in})
                        out_item['matched'] = True

                        for match in matched_inflows:
                            w_out_expr = match['item']['expr']
                            w_out = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(w_out_expr))
                            pn_arcs.append({"source": t_name, "target": match['item']['place'], "weight": w_out})
                            match['item']['matched'] = True
            else:
                # fallback logic
                for out_item in outflows:
                    if out_item['matched']: continue
                    candidates = []
                    for in_item in inflows:
                        if in_item['matched']: continue

                        key_check = f"{in_item['place']}|{str(in_item['expr'])}|{str(out_item['expr'])}"
                        if merge_strategies.get(key_check) == "separate": continue

                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                        if not (ratio.free_symbols & state_syms):
                            if in_item['place'] == out_item['place']: continue
                            candidates.append(in_item)

                    if not candidates: continue
                    match_found = None
                    total_in = sum(c['expr'] for c in candidates)
                    if sp.simplify(out_item['expr'] - total_in) == 0:
                        match_found = candidates

                    if match_found:
                        out_item['matched'] = True
                        rate = out_item['expr']
                        t_name = f"t{t_counter}"
                        t_counter += 1
                        pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})
                        w_src = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(rate))
                        pn_arcs.append({"source": out_item['place'], "target": t_name, "weight": w_src})
                        for m in match_found:
                            m['matched'] = True
                            w_dest = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(m['expr']))
                            pn_arcs.append({"source": t_name, "target": m['place'], "weight": w_dest})


            for item in inflows:
                if not item['matched']:
                    t_name = f"t{t_counter}"
                    t_counter += 1
                    pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})
                    w_item = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(item['expr']))
                    pn_arcs.append({"source": t_name, "target": item['place'], "weight": w_item})

            for item in outflows:
                if item['matched']: continue
                t_name = f"t{t_counter}"
                t_counter += 1
                pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})
                w_item = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(item['expr']))
                pn_arcs.append({"source": item['place'], "target": t_name, "weight": w_item})

            return {"places": pn_places, "transitions": pn_transitions, "arcs": pn_arcs, "status": "success"}

# ODE System Logic
class ODESystem:
    def __init__(self):
        self.variables = []
        self.symbols = {}
        self.equations = {}
        self.infected_vars = []
        self.human_vars = []

    def add_variable(self, name):
        name = encode_latex(name)
        if name not in self.variables:
            self.variables.append(name)
            self.symbols[name] = sp.Symbol(name, real=True)

    def remove_variable(self, name):
        if name in self.variables:
            self.variables.remove(name)
            if name in self.symbols: del self.symbols[name]
            if name in self.equations: del self.equations[name]
            if name in self.infected_vars: self.infected_vars.remove(name)
            if name in self.human_vars: self.human_vars.remove(name)

    def _register_symbols(self, expr_str):
        potential = re.findall(r'[a-zA-Z_][a-zA-Z0-9_]*', expr_str)
        for p in potential:
            if p not in self.symbols: self.symbols[p] = sp.Symbol(p, real=True)

    def set_equation(self, var, full_expr):
        var = encode_latex(var)
        full_expr = encode_latex(full_expr)
        self._register_symbols(full_expr)
        self.equations[var] = full_expr

    def set_properties(self, var, is_infected, is_human):
        if is_infected and var not in self.infected_vars:
            self.infected_vars.append(var)
        elif not is_infected and var in self.infected_vars:
            self.infected_vars.remove(var)
        if is_human and var not in self.human_vars:
            self.human_vars.append(var)
        elif not is_human and var in self.human_vars:
            self.human_vars.remove(var)

    def calculate_dfe(self, infected_vars, free_vars=[], constraints=[]):
        P_N = [v for v in self.variables if v not in infected_vars]
        subs_I_zero = {self.symbols[v]: 0 for v in infected_vars}
        for c in constraints: self._register_symbols(c)

        system_eqs = []
        vars_to_solve = [self.symbols[v] for v in P_N if v not in free_vars]

        for v in P_N:
            if v in free_vars: continue
            raw_expr = self.equations.get(v, "0")
            try:
                sym_expr = sp.sympify(raw_expr, locals=self.symbols)
                system_eqs.append(sym_expr.subs(subs_I_zero))
            except:
                return None, f"Parse Error in {v}"

        constraint_subs = {}
        remaining_constraints = []
        for c in constraints:
            if "=" in c:
                lhs_str, rhs_str = c.split("=")
                try:
                    lhs = sp.sympify(lhs_str, locals=self.symbols)
                    rhs = sp.sympify(rhs_str, locals=self.symbols)
                    eqn = lhs - rhs
                    syms = list(eqn.free_symbols)
                    if syms:
                        target = syms[0]
                        res = sp.solve(eqn, target)
                        if res:
                            constraint_subs[target] = res[0]
                        else:
                            remaining_constraints.append(eqn)
                    else:
                        remaining_constraints.append(eqn)
                except:
                    pass
            else:
                try:
                    remaining_constraints.append(sp.sympify(c, locals=self.symbols))
                except:
                    pass

        final_eqs = [eq.subs(constraint_subs) for eq in system_eqs] + [eq.subs(constraint_subs) for eq in
                                                                       remaining_constraints]

        try:
            solutions = sp.solve(final_eqs, vars_to_solve, dict=True)
        except Exception as e:
            return None, f"Solver Error: {e}"

        final_solutions = []
        if not solutions: solutions = [{}]
        sym_N = sp.Symbol('N', real=True)

        for sol in solutions:
            current_sol_list = [sol]
            missing_vars = [v for v in vars_to_solve if v not in sol]
            if missing_vars and self.human_vars:
                human_sum_expr = 0
                for h_var in self.human_vars:
                    h_sym = self.symbols[h_var]
                    if h_var in infected_vars:
                        human_sum_expr += 0
                    elif h_sym in sol:
                        human_sum_expr += sol[h_sym]
                    else:
                        human_sum_expr += h_sym
                constraint = human_sum_expr - sym_N
                try:
                    constrained_sols = sp.solve(final_eqs + [constraint], vars_to_solve, dict=True)
                    if constrained_sols: current_sol_list = constrained_sols
                except:
                    pass

            for s in current_sol_list:
                full_sol = s.copy()
                full_sol.update(subs_I_zero)
                for free_v in free_vars:
                    sym = self.symbols[free_v]
                    if sym not in full_sol: full_sol[sym] = sym
                full_str = {}
                for k, v in full_sol.items():
                    simplified_val = sp.simplify(v)
                    full_str[k.name] = format_expr_for_display(simplified_val, self.variables)
                if full_str not in final_solutions: final_solutions.append(full_str)

        if not final_solutions: return [], "No DFE found"
        return final_solutions, "Success"

    def calculate_r0(self, infected_vars, dfe_dict, source_classification=None):
        if source_classification is None: source_classification = {}

        subs_dfe = {}
        parse_locals = self.symbols.copy()
        parse_locals['N'] = sp.Symbol('N', real=True)
        for k, v in dfe_dict.items():
            sym = self.symbols[k] if k in self.symbols else sp.Symbol(k, real=True)
            v_clean = sanitize_input_str(str(v))
            try:
                expr = sp.sympify(v_clean, locals=parse_locals)
            except:
                expr = sp.sympify(v_clean)
            subs_dfe[sym] = expr

        #  Calculate N*
        N_val = sum(subs_dfe.values())
        sym_N = sp.Symbol('N', real=True)
        subs_dfe[sym_N] = N_val

        N_str_parts = [f"{v}^*" for v in self.variables]
        N_def_lhs = " + ".join(N_str_parts)
        N_def_rhs = format_expr_for_display(N_val, self.variables)
        N_def_str = f"N^* = {N_def_lhs} = {N_def_rhs}"

        F_list = []
        V_list = []
        infected_syms = [self.symbols[v] for v in infected_vars]
        non_inf_vars = [v for v in self.variables if v not in infected_vars]
        non_inf_syms = set()
        for v in non_inf_vars:
            if v in self.symbols: non_inf_syms.add(self.symbols[v])

        state_syms = {self.symbols[v] for v in self.variables}
        ambiguous_terms = []

        for v in infected_vars:
            raw_eq = self.equations.get(v, "0")
            terms = []
            try:
                raw_expr = sp.sympify(raw_eq, locals=self.symbols)
                if isinstance(raw_expr, sp.Add):
                    base_terms = raw_expr.args
                else:
                    base_terms = [raw_expr]

                for term in base_terms:
                    has_state_add = False
                    for atom in term.atoms(sp.Add):
                        if atom.free_symbols & state_syms: has_state_add = True; break
                    if has_state_add:
                        expanded = term.expand()
                        if isinstance(expanded, sp.Add):
                            terms.extend(expanded.args)
                        else:
                            terms.append(expanded)
                    else:
                        terms.append(term)
            except:
                return [f"Math Error in {v}"], None, None, [], [], None

            f_expr = 0
            v_expr = 0

            for term in terms:
                coeff = term.as_coeff_Mul()[0]
                is_positive = (coeff > 0)
                has_infected = any(term.has(i) for i in infected_syms)
                has_susceptible = any(term.has(s) for s in non_inf_syms)

                if has_infected and has_susceptible:
                    if is_positive:
                        f_expr += term
                    else:
                        v_expr += -term
                elif has_infected and not has_susceptible:
                    if is_positive:
                        term_key = f"{str(term)}_in_{v}"
                        decision = source_classification.get(term_key)
                        if decision == 'F':
                            f_expr += term
                        elif decision == 'V':
                            v_expr += -term
                        else:
                            ambiguous_terms.append(
                                {'name': term_key, 'rate': f"Term '{str(term)}' in equation for {v}"})
                    else:
                        v_expr += -term
                else:
                    v_expr += -term

            F_list.append(f_expr)
            V_list.append(v_expr)

        if ambiguous_terms:
            return None, "ambiguous", ambiguous_terms, [], [], None

        try:
            F_mat = sp.Matrix(F_list).jacobian(infected_syms).subs(subs_dfe)
            V_mat = sp.Matrix(V_list).jacobian(infected_syms).subs(subs_dfe)
            eigs = (F_mat * V_mat.inv()).eigenvals()
            r0_list = list(eigs.keys())
            r0_fmt = [format_expr_for_display(r, self.variables) for r in r0_list]
            F_fmt = format_expr_for_display(F_mat, self.variables)
            V_fmt = format_expr_for_display(V_mat, self.variables)
            sensitivity_data = []
            for r_expr in r0_list:
                sensitivity_data.append(calculate_sensitivity_indices(r_expr))

            #  Return 6 items
            return r0_fmt, F_fmt, V_fmt, sensitivity_data, r0_list, N_def_str
        except Exception as e:
            return [f"Math Error: {e}"], None, None, [], [], None

    #  Hybrid conversion:SPN
    def convert_to_spn_structure(self, use_scaling=True, merge_strategies=None):
            if merge_strategies is None: merge_strategies = {}
            state_syms = {self.symbols[v] for v in self.variables}
            all_terms = []

            #1.parse equations into rate structures (smart expansion v3 (key of conversion))
            for var, eq in self.equations.items():
                if not eq: continue
                try:
                    raw_expr = sp.sympify(eq, locals=self.symbols)

                    expanded_terms = []
                    if isinstance(raw_expr, sp.Add):
                        base_terms = raw_expr.args
                    else:
                        base_terms = [raw_expr]

                    for term in base_terms:
                        should_expand = False

                        if term.is_Mul:
                            add_part = None
                            multiplier_syms = set()

                            for arg in term.args:
                                if isinstance(arg, sp.Add):
                                    add_part = arg
                                else:
                                    multiplier_syms.update(arg.free_symbols)

                            if add_part:
                                inner_state_syms = add_part.free_symbols & state_syms

                                # Criteria 1: Multiple different state variables inside
                                if len(inner_state_syms) > 1:
                                    # EXCEPTION: Pure parameter multiplier to Keep Grouped
                                    if multiplier_syms.isdisjoint(state_syms):
                                        should_expand = False
                                    else:
                                        should_expand = True

                                # Criteria 2: Overlap to Expand
                                elif not inner_state_syms.isdisjoint(multiplier_syms & state_syms):
                                    should_expand = True

                        if should_expand:
                            expanded_part = term.expand()
                            if isinstance(expanded_part, sp.Add):
                                expanded_terms.extend(expanded_part.args)
                            else:
                                expanded_terms.append(expanded_part)
                        else:
                            expanded_terms.append(term)

                    for term in expanded_terms:
                        coeff = term.as_coeff_Mul()[0]
                        if coeff < 0:
                            rate_struct = sp.simplify(-term)
                            direction = -1
                        else:
                            rate_struct = sp.simplify(term)
                            direction = 1
                        all_terms.append(
                            {'place': var, 'rate_struct': rate_struct, 'direction': direction, 'matched': False})

                except:
                    continue

            outflows = [t for t in all_terms if t['direction'] < 0]
            inflows = [t for t in all_terms if t['direction'] > 0]

            #  Step 0.5: Ambiguity Detection
            ambiguities = []

            # Check 1: Proportional Self-Loops (The Vector Model Fix)
            if use_scaling:
                for in_item in inflows:
                    for out_item in outflows:
                        if out_item['place'] != in_item['place']: continue
                        if out_item['matched'] or in_item['matched']: continue

                        try:
                            ratio = sp.simplify(in_item['rate_struct'] / out_item['rate_struct'])
                        except:
                            continue

                        # If Ratio is constant but not 1, it's a proportional self-loop (Migration vs Death)
                        if ratio != 1 and not (ratio.free_symbols & state_syms):
                            in_expr = str(in_item['rate_struct'])
                            out_expr = str(out_item['rate_struct'])
                            key = f"{in_item['place']}|{in_expr}|{out_expr}"
                            if key not in merge_strategies:
                                ambiguities.append({
                                    "key": key, "place": in_item['place'], "type": "loop",
                                    "inflow_expr": in_expr, "outflow_expr": out_expr, "catalyst": "Self"
                                })

            # Check 2: Shared Catalyst Symbols
            for in_item in inflows:
                if in_item['matched']: continue
                in_syms = in_item['rate_struct'].free_symbols & state_syms
                cross_syms = {s for s in in_syms if str(s) != in_item['place']}
                if cross_syms:
                    for out_item in outflows:
                        if out_item['place'] != in_item['place']: continue

                        # it is a standard flow if terms are identical, auto-merge in step 1.
                        if sp.simplify(in_item['rate_struct'] - out_item['rate_struct']) == 0:
                            continue

                        out_syms = out_item['rate_struct'].free_symbols & state_syms
                        shared = cross_syms.intersection(out_syms)
                        if shared:
                            in_expr = str(in_item['rate_struct'])
                            out_expr = str(out_item['rate_struct'])
                            key = f"{in_item['place']}|{in_expr}|{out_expr}"
                            if key not in merge_strategies:
                                ambiguities.append({
                                    "key": key, "place": in_item['place'], "type": "catalyst",
                                    "inflow_expr": in_expr, "outflow_expr": out_expr, "catalyst": str(list(shared)[0])
                                })

            if ambiguities:
                return {"status": "ambiguous", "details": ambiguities}

            #  Step 1: generate transitions in unified loop
            temp_transitions = []

            if use_scaling:

                    # Pass 1: Exact matches (Priority)
                    # Handle Identity Flows (Conservation) first.
                    for out_item in outflows:
                        if out_item['matched']: continue

                        # Search for one exact match (Ratio == 1)
                        exact_match = None
                        for in_item in inflows:
                            if in_item['matched']: continue

                            try:
                                ratio = sp.simplify(in_item['rate_struct'] / out_item['rate_struct'])
                            except:
                                continue

                            if ratio == 1:
                                # Self-Loop Logic Check
                                if in_item['place'] == out_item['place']:
                                    key = f"{in_item['place']}|{str(in_item['rate_struct'])}|{str(out_item['rate_struct'])}"
                                    if merge_strategies.get(key) != "merge": continue

                                exact_match = in_item
                                break

                        if exact_match:
                            # Create Transition for Identity Flow
                            out_item['matched'] = True
                            exact_match['matched'] = True

                            rate = out_item['rate_struct']

                            # Calculate Inputs (Source + Catalysts)
                            inputs = defaultdict(int)
                            inputs[out_item['place']] += 1
                            for sym in rate.free_symbols & state_syms:
                                if str(sym) != out_item['place']: inputs[str(sym)] += 1

                            # Calculate Outputs (Target + Catalysts)
                            outputs = defaultdict(int)
                            outputs[exact_match['place']] += 1  # Ratio is 1

                            # Add Catalysts to Outputs (Preserve them)
                            for p, w in inputs.items():
                                if p != out_item['place']: outputs[p] += w

                            temp_transitions.append({'rate': rate, 'inputs': inputs, 'outputs': outputs})

                    # Pass 2: proportional matches (cleanup)
                    # Handle Splitting/Branching events (ex: E -> r*Ia + (1-r)*Is)
                    for out_item in outflows:
                        if out_item['matched']: continue

                        matched_inflows = []
                        for in_item in inflows:
                            if in_item['matched']: continue

                            try:
                                ratio = sp.simplify(in_item['rate_struct'] / out_item['rate_struct'])
                            except:
                                continue

                            # Check Independence & Scaling
                            if (ratio.free_symbols & state_syms): continue

                            # Self-Loop Check
                            if in_item['place'] == out_item['place']:
                                key = f"{in_item['place']}|{str(in_item['rate_struct'])}|{str(out_item['rate_struct'])}"
                                if merge_strategies.get(key) != "merge": continue

                            matched_inflows.append({'item': in_item, 'weight': ratio})

                        if matched_inflows:
                            out_item['matched'] = True
                            rate = out_item['rate_struct']

                            inputs = defaultdict(int)
                            inputs[out_item['place']] += 1
                            for sym in rate.free_symbols & state_syms:
                                if str(sym) != out_item['place']: inputs[str(sym)] += 1

                            outputs = defaultdict(int)
                            for m in matched_inflows:
                                m['item']['matched'] = True
                                outputs[m['item']['place']] += m['weight']

                            # Preserve Catalysts
                            for p, w in inputs.items():
                                if p != out_item['place']: outputs[p] += w

                            temp_transitions.append({'rate': rate, 'inputs': inputs, 'outputs': outputs})
            else:
                    # Fallback Logic (No Scaling)
                    for out_item in outflows:
                        if out_item['matched']: continue
                        candidates = []
                        for in_item in inflows:
                            if in_item['matched']: continue
                            try:
                                ratio = sp.simplify(in_item['rate_struct'] / out_item['rate_struct'])
                            except:
                                continue
                            if not (ratio.free_symbols & state_syms):
                                if in_item['place'] == out_item['place']: continue
                                candidates.append(in_item)

                        if not candidates: continue
                        match_found = None
                        total_in = sum(c['rate_struct'] for c in candidates)
                        if sp.simplify(out_item['rate_struct'] - total_in) == 0:
                            match_found = candidates

                        if match_found:
                            out_item['matched'] = True
                            for m in match_found: m['matched'] = True
                            rate = out_item['rate_struct']
                            inputs = defaultdict(int)
                            inputs[out_item['place']] += 1
                            for sym in rate.free_symbols & state_syms:
                                if str(sym) != out_item['place']: inputs[str(sym)] += 1
                            outputs = defaultdict(int)
                            for m in match_found:
                                w = sp.simplify(m['rate_struct'] / rate)
                                outputs[m['place']] += w
                            for p, w in inputs.items():
                                if p != out_item['place']: outputs[p] += w
                            temp_transitions.append({'rate': rate, 'inputs': inputs, 'outputs': outputs})

            #cleanup
            for out_item in outflows:
                if out_item['matched']: continue
                rate = out_item['rate_struct']
                inputs = defaultdict(int)
                inputs[out_item['place']] += 1
                for sym in rate.free_symbols & state_syms:
                    if str(sym) != out_item['place']: inputs[str(sym)] += 1
                outputs = defaultdict(int)
                for p, w in inputs.items():
                    if p != out_item['place']: outputs[p] += w
                temp_transitions.append({'rate': rate, 'inputs': inputs, 'outputs': outputs})

            remaining_inflows = [t for t in inflows if not t['matched']]
            while remaining_inflows:
                root = remaining_inflows[0]
                cluster = [root]
                remaining_inflows.remove(root)
                changed = True
                while changed:
                    changed = False
                    to_remove = []
                    for t in remaining_inflows:
                        for mem in cluster:
                            common = sp.gcd(mem['rate_struct'], t['rate_struct'])
                            param_syms = set(self.symbols.values()) - state_syms
                            if common.free_symbols & param_syms:
                                cluster.append(t);
                                to_remove.append(t);
                                changed = True;
                                break
                    for t in to_remove: remaining_inflows.remove(t)

                place_sums = defaultdict(int)
                for c in cluster: place_sums[c['place']] += c['rate_struct']
                sum_exprs = list(place_sums.values())
                common_rate = sum_exprs[0]
                for expr in sum_exprs[1:]: common_rate = sp.gcd(common_rate, expr)
                if common_rate == 1 and len(sum_exprs) > 0: common_rate = sum(sum_exprs)

                inputs = defaultdict(int)
                final_outputs = defaultdict(int)
                for p, p_sum in place_sums.items():
                    w = sp.simplify(p_sum / common_rate)
                    final_outputs[p] = w
                temp_transitions.append({'rate': common_rate, 'inputs': inputs, 'outputs': final_outputs})

            # Preserves Separate Arrows
            pn_trans = []
            pn_arcs = []
            t_id = 1

                #process transitions directly without merging
            for t in temp_transitions:
                    #format Rate
                raw_rate = str(sp.simplify(t['rate'])).replace("**", "^")
                rate_str = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', raw_rate)

                t_name = f"t{t_id}"
                pn_trans.append({"name": t_name, "type": "trans", "rate": rate_str, "x": 0, "y": 0})

                # input arcs
                for p, w in t['inputs'].items():
                     w_str = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(w))
                     pn_arcs.append({"source": p, "target": t_name, "weight": w_str})

                # create output arcs
                for p, w in t['outputs'].items():
                    w_str = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(w))
                    pn_arcs.append({"source": t_name, "target": p, "weight": w_str})

                t_id += 1

            pn_places = [{"name": v, "type": "place", "x": 0, "y": 0} for v in self.variables]
            return {"places": pn_places, "transitions": pn_trans, "arcs": pn_arcs, "status": "success"}
            # merged_map = {}
            # for t in temp_transitions:
            #     in_sig = tuple(sorted((k, str(v)) for k, v in t['inputs'].items()))
            #     out_sig = tuple(sorted((k, str(v)) for k, v in t['outputs'].items()))
            #     sig = (in_sig, out_sig)
            #     if sig not in merged_map: merged_map[sig] = []
            #     merged_map[sig].append(t['rate'])
            #
            # pn_trans = []
            # pn_arcs = []
            # t_id = 1
            # for sig, rates in merged_map.items():
            #     in_sig, out_sig = sig
            #     total_rate = sum(rates)
            #     raw_rate = str(sp.simplify(total_rate)).replace("**", "^")
            #     rate_str = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', raw_rate)
            #     t_name = f"t{t_id}"
            #     pn_trans.append({"name": t_name, "type": "trans", "rate": rate_str, "x": 0, "y": 0})
            #     for p, w in in_sig:
            #         w_str = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(w))
            #         pn_arcs.append({"source": p, "target": t_name, "weight": w_str})
            #     for p, w in out_sig:
            #         w_str = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(w))
            #         pn_arcs.append({"source": t_name, "target": p, "weight": w_str})
            #     t_id += 1
            #
            # pn_places = [{"name": v, "type": "place", "x": 0, "y": 0} for v in self.variables]
            # return {"places": pn_places, "transitions": pn_trans, "arcs": pn_arcs, "status": "success"}

            # Hybrid conversion of VAPN

    def convert_to_vapn_structure(self, use_scaling=True, merge_strategies=None):
        if merge_strategies is None: merge_strategies = {}
        outflows = []
        inflows = []
        state_syms = {self.symbols[v] for v in self.variables}

        #1.parse equations into rate structures (smart expansion v3 (key of conversion))
        for var, eq_str in self.equations.items():
            if not eq_str.strip(): continue
            try:
                raw_expr = sp.sympify(eq_str, locals=self.symbols)

                expanded_terms = []
                if isinstance(raw_expr, sp.Add):
                    base_terms = raw_expr.args
                else:
                    base_terms = [raw_expr]

                for term in base_terms:
                    should_expand = False

                    if term.is_Mul:
                        add_part = None
                        multiplier_syms = set()

                        for arg in term.args:
                            if isinstance(arg, sp.Add):
                                add_part = arg
                            else:
                                multiplier_syms.update(arg.free_symbols)

                        if add_part:
                            inner_state_syms = add_part.free_symbols & state_syms

                            # crit 1: Multiple different state variables inside ( S+I+R)
                            if len(inner_state_syms) > 1:
                                # exception if multiplier is pure parameters (no state vars),
                                # treat as a scalar group ( Total Pop N). Keep grouped.
                                if multiplier_syms.isdisjoint(state_syms):
                                    should_expand = False
                                else:
                                    should_expand = True

                            # crit 2: Multiplier overlaps with Inner (S*(1-S)) -> Expand
                            elif not inner_state_syms.isdisjoint(multiplier_syms & state_syms):
                                should_expand = True

                    if should_expand:
                        expanded_part = term.expand()
                        if isinstance(expanded_part, sp.Add):
                            expanded_terms.extend(expanded_part.args)
                        else:
                            expanded_terms.append(expanded_part)
                    else:
                        expanded_terms.append(term)

                # Process the refined terms
                for term in expanded_terms:
                    coeff = term.as_coeff_Mul()[0]
                    if coeff < 0:
                        outflows.append({'place': var, 'expr': -term, 'matched': False})
                    elif term.is_Mul and term.args[0].is_Number and term.args[0] < 0:
                        outflows.append({'place': var, 'expr': -term, 'matched': False})
                    else:
                        inflows.append({'place': var, 'expr': term, 'matched': False})
            except:
                continue

        # Step 0.5: AMBIGUITY DETECTION
        ambiguities = []

        # Check 1: Proportional Self-Loops
        if use_scaling:
            for in_item in inflows:
                for out_item in outflows:
                    if out_item['place'] != in_item['place']: continue
                    if out_item['matched'] or in_item['matched']: continue

                    try:
                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                    except:
                        continue

                    if ratio != 1 and not (ratio.free_symbols & state_syms):
                        in_expr = str(in_item['expr'])
                        out_expr = str(out_item['expr'])
                        key = f"{in_item['place']}|{in_expr}|{out_expr}"
                        if key not in merge_strategies:
                            ambiguities.append({
                                "key": key, "place": in_item['place'], "type": "loop",
                                "inflow_expr": in_expr, "outflow_expr": out_expr, "catalyst": "Self"
                            })

        # Check 2: Catalyst Shared Symbols
        for in_item in inflows:
            if in_item['matched']: continue
            in_syms = in_item['expr'].free_symbols & state_syms
            cross_syms = {s for s in in_syms if str(s) != in_item['place']}
            if cross_syms:
                for out_item in outflows:
                    if out_item['place'] != in_item['place']: continue

                    # If terms are identical, it is a standard flow. Auto-merge in Phase 1.
                    if sp.simplify(in_item['expr'] - out_item['expr']) == 0:
                        continue

                    out_syms = out_item['expr'].free_symbols & state_syms
                    shared = cross_syms.intersection(out_syms)
                    if shared:
                        in_expr = str(in_item['expr'])
                        out_expr = str(out_item['expr'])
                        key = f"{in_item['place']}|{in_expr}|{out_expr}"
                        if key not in merge_strategies:
                            ambiguities.append({
                                "key": key, "place": in_item['place'], "type": "catalyst",
                                "inflow_expr": in_expr, "outflow_expr": out_expr, "catalyst": str(list(shared)[0])
                            })

        if ambiguities:
            return {"status": "ambiguous", "details": ambiguities}

        # Step 1: GENERATE TRANSITIONS (UNIFIED)
        pn_places = [{"name": v, "type": "place", "x": 0, "y": 0} for v in self.variables]
        pn_transitions = []
        pn_arcs = []
        t_counter = 1

        if use_scaling:
                # Pass 1: EXACT MATCHES (Priority)
                # Handle Identity Flows (Conservation) first so they aren't stolen by proportional terms
                #
            for out_item in outflows:
                if out_item['matched']: continue

                # Search for one Exact Match (Ratio == 1)
                exact_match = None
                for in_item in inflows:
                    if in_item['matched']: continue

                    try:
                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                    except:
                        continue

                    if ratio == 1:
                        # Self-Loop Logic Check
                        if in_item['place'] == out_item['place']:
                            key = f"{in_item['place']}|{str(in_item['expr'])}|{str(out_item['expr'])}"
                            if merge_strategies.get(key) != "merge": continue

                        exact_match = in_item
                        break

                if exact_match:
                    # Create Transition for Identity Flow
                    t_name = f"t{t_counter}"
                    t_counter += 1
                    pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})

                    # Arc In
                    w_in = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(out_item['expr']))
                    pn_arcs.append({"source": out_item['place'], "target": t_name, "weight": w_in})
                    out_item['matched'] = True

                    # Arc Out
                    w_out = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(exact_match['expr']))
                    pn_arcs.append({"source": t_name, "target": exact_match['place'], "weight": w_out})
                    exact_match['matched'] = True


            # Pass 2: proportional matches
            # Handle Splitting/Branching events ex: E -> r*Ia + (1-r)*Is

            for out_item in outflows:
                if out_item['matched']: continue

                matched_inflows = []
                for in_item in inflows:
                    if in_item['matched']: continue

                    try:
                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                    except:
                        continue

                    # Check Independence & Scaling
                    if (ratio.free_symbols & state_syms): continue

                    # Self-Loop Check
                    if in_item['place'] == out_item['place']:
                        key = f"{in_item['place']}|{str(in_item['expr'])}|{str(out_item['expr'])}"
                        if merge_strategies.get(key) != "merge": continue

                    matched_inflows.append({'item': in_item, 'ratio': ratio})

                if matched_inflows:
                    #Create transition for branching flow
                    t_name = f"t{t_counter}"
                    t_counter += 1
                    pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})

                    w_in = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(out_item['expr']))
                    pn_arcs.append({"source": out_item['place'], "target": t_name, "weight": w_in})
                    out_item['matched'] = True

                    for match in matched_inflows:
                        match['item']['matched'] = True
                        w_out_expr = match['item']['expr']
                        w_out = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(w_out_expr))
                        pn_arcs.append({"source": t_name, "target": match['item']['place'], "weight": w_out})

        else:
            #No Scaling Mode
            for out_item in outflows:
                if out_item['matched']: continue
                candidates = []
                for in_item in inflows:
                    if in_item['matched']: continue
                    try:
                        ratio = sp.simplify(in_item['expr'] / out_item['expr'])
                    except:
                        continue
                    if not (ratio.free_symbols & state_syms):
                        if in_item['place'] == out_item['place']: continue
                        candidates.append(in_item)

                if not candidates: continue
                match_found = None
                total_in = sum(c['expr'] for c in candidates)
                if sp.simplify(out_item['expr'] - total_in) == 0:
                    match_found = candidates

                if match_found:
                    out_item['matched'] = True
                    rate = out_item['expr']
                    t_name = f"t{t_counter}"
                    t_counter += 1
                    pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})
                    w_src = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(rate))
                    pn_arcs.append({"source": out_item['place'], "target": t_name, "weight": w_src})
                    for m in match_found:
                        m['matched'] = True
                        w_dest = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(m['expr']))
                        pn_arcs.append({"source": t_name, "target": m['place'], "weight": w_dest})

        for item in inflows:
            if not item['matched']:
                t_name = f"t{t_counter}"
                t_counter += 1
                pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})
                w_item = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(item['expr']))
                pn_arcs.append({"source": t_name, "target": item['place'], "weight": w_item})

        for item in outflows:
            if item['matched']: continue
            t_name = f"t{t_counter}"
            t_counter += 1
            pn_transitions.append({"name": t_name, "type": "trans", "rate": "1.0", "x": 0, "y": 0})
            w_item = re.sub(r'greek_([a-zA-Z]+)', r'\\\1', str(item['expr']))
            pn_arcs.append({"source": item['place'], "target": t_name, "weight": w_item})

        return {"places": pn_places, "transitions": pn_transitions, "arcs": pn_arcs, "status": "success"}