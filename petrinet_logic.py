import sympy as sp
import re


class VariableWeightPetriNet:
    def __init__(self):
        self.places = []  # List of place names (strings)
        self.transitions = []  # List of transition names
        self.arcs = []  # List of dictionaries
        self.symbols = {}  # Maps names to SymPy symbols

    def add_place(self, name):
        if name not in self.places:
            self.places.append(name)
            self.symbols[name] = sp.Symbol(name, real=True)

    def remove_place(self, name):
        if name in self.places:
            self.places.remove(name)
            if name in self.symbols:
                del self.symbols[name]
            self.arcs = [arc for arc in self.arcs if arc['source'] != name and arc['target'] != name]

    def add_transition(self, name):
        if name not in self.transitions:
            self.transitions.append(name)

    def remove_transition(self, name):
        if name in self.transitions:
            self.transitions.remove(name)
            self.arcs = [arc for arc in self.arcs if arc['source'] != name and arc['target'] != name]

    def _register_new_symbols(self, expression_str):
        potential_symbols = re.findall(r'[a-zA-Z_][a-zA-Z0-9_]*', expression_str)
        for name in potential_symbols:
            if name not in self.symbols:
                self.symbols[name] = sp.Symbol(name, real=True)

    def add_arc(self, source, target, weight_expr):
        self._register_new_symbols(weight_expr)
        try:
            expr = sp.sympify(weight_expr, locals=self.symbols)
        except Exception:
            expr = sp.sympify(weight_expr)

        self.arcs.append({
            'source': source,
            'target': target,
            'weight_str': weight_expr,
            'weight_expr': expr
        })

    def remove_arc(self, source, target):
        self.arcs = [arc for arc in self.arcs if not (arc['source'] == source and arc['target'] == target)]

    def get_input_arcs(self, node):
        return [arc for arc in self.arcs if arc['target'] == node]

    def get_output_arcs(self, node):
        return [arc for arc in self.arcs if arc['source'] == node]

    # --- PHASE 3 LOGIC: DFE & R0 ---

    def calculate_dfe(self, infected_places_names):
        """
        Calculates DFE. Returns a dictionary with STRING keys { 'S': 'N' }
        so the GUI can display/edit them easily.
        """
        P_I = infected_places_names
        P_N = [p for p in self.places if p not in P_I]

        # We explicitly use N as a symbol here for the Closed System check
        sym_N = sp.Symbol('N', real=True)

        subs_I_zero = {self.symbols[p]: 0 for p in P_I}

        equations = []
        vars_to_solve = [self.symbols[p] for p in P_N]

        for p in P_N:
            inflow = sum(arc['weight_expr'].subs(subs_I_zero) for arc in self.get_input_arcs(p))
            outflow = sum(arc['weight_expr'].subs(subs_I_zero) for arc in self.get_output_arcs(p))
            equations.append(inflow - outflow)

        try:
            solutions = sp.solve(equations, vars_to_solve, dict=True)
        except Exception as e:
            return None, f"Solver Error: {e}"

        if not solutions:
            best_sol = {}
        else:
            best_sol = solutions[0]
            for sol in solutions:
                if len(sol) > 0 and not all(val == 0 for val in sol.values()):
                    best_sol = sol
                    break

        # Check for free variables (Closed System Logic)
        free_vars = [v for v in vars_to_solve if v not in best_sol]
        if free_vars:
            sum_PN = sum(vars_to_solve)
            constraint = sum_PN - sym_N
            try:
                constrained_sols = sp.solve(equations + [constraint], vars_to_solve, dict=True)
                if constrained_sols:
                    best_sol = constrained_sols[0]
            except Exception as e:
                return None, f"Constrained Solver Error: {e}"

        # Merge with I=0
        full_dfe_sym = best_sol.copy()
        full_dfe_sym.update(subs_I_zero)

        # CONVERT TO STRINGS for GUI
        full_dfe_str = {k.name: v for k, v in full_dfe_sym.items()}

        return full_dfe_str, "Success"

    def calculate_r0(self, infected_places_names, dfe_substitutions_str):
        """
        Calculates R0.
        dfe_substitutions_str: Dict { 'S': expression/string }
        """
        # 1. Prepare a parsing dictionary that includes all known symbols
        # CRITICAL FIX: Explicitly define 'N' as a Symbol to prevent "function N" error
        parse_locals = self.symbols.copy()
        parse_locals['N'] = sp.Symbol('N', real=True)

        # 2. Convert String keys back to Symbols using the safe dictionary
        dfe_substitutions = {}
        for k, v in dfe_substitutions_str.items():
            # Resolve Key
            if k in self.symbols:
                sym_key = self.symbols[k]
            else:
                sym_key = sp.Symbol(k, real=True)

            # Resolve Value (using parse_locals to safely handle "N")
            try:
                val_expr = sp.sympify(v, locals=parse_locals)
            except Exception:
                # Fallback
                val_expr = sp.sympify(v)

            dfe_substitutions[sym_key] = val_expr

        # 3. Standard R0 Logic
        P_I = infected_places_names
        P_N = [p for p in self.places if p not in P_I]

        T_inf = []
        for t in self.transitions:
            inputs = [arc['source'] for arc in self.get_input_arcs(t)]
            outputs = [arc['target'] for arc in self.get_output_arcs(t)]

            has_input_from_PN = any(p in P_N for p in inputs)
            has_output_to_PI = any(p in P_I for p in outputs)

            if has_input_from_PN and has_output_to_PI:
                T_inf.append(t)

        if not T_inf or not P_I:
            return ["Error: No Infection Transitions or Infected Places defined."]

        F_matrix_symbolic = []
        for i_place in P_I:
            rate_new_infection = 0
            input_arcs = self.get_input_arcs(i_place)
            for arc in input_arcs:
                if arc['source'] in T_inf:
                    rate_new_infection += arc['weight_expr']
            F_matrix_symbolic.append(rate_new_infection)

        V_matrix_symbolic = []
        for i_place in P_I:
            outflow = sum(arc['weight_expr'] for arc in self.get_output_arcs(i_place))
            inflow = sum(arc['weight_expr'] for arc in self.get_input_arcs(i_place) if arc['source'] not in T_inf)
            V_matrix_symbolic.append(outflow - inflow)

        infected_symbols = [self.symbols[name] for name in P_I]
        F_jacobian = sp.Matrix(F_matrix_symbolic).jacobian(infected_symbols)
        V_jacobian = sp.Matrix(V_matrix_symbolic).jacobian(infected_symbols)

        # Substitute DFE
        F_at_dfe = F_jacobian.subs(dfe_substitutions)
        V_at_dfe = V_jacobian.subs(dfe_substitutions)

        try:
            V_inv = V_at_dfe.inv()
            NextGenMatrix = F_at_dfe * V_inv
            eigenvalues = NextGenMatrix.eigenvals()
            return list(eigenvalues.keys())
        except Exception as e:
            return [f"Math Error: {e}"]