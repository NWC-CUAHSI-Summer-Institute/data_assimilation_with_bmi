from bmipy import Bmi
import numpy as np
import json
from new_EnKF import EnKF


class BMIEnKF(Bmi):
    def __init__(self, json_file):
        self.json_file = json_file  # store json file name for later
        with open(json_file, "r") as file:
            data = json.load(file)

        self.n = data['n']
        self.m = data['m']
        self.R = np.eye(self.m) * np.array(data['R'])
        self.Q = np.eye(self.n) * np.array(data['Q'])
        self.N = data['N']
        self.P = np.eye(self.n) * np.array(data['P'])
        self.smcmax = data['smcmax']
        self.D = data['D']
        self.enkf = EnKF(self.n, self.m, self.R, self.Q, self.N, self.P, self.smcmax, self.D)
        self.state_estimates = np.zeros((self.n))
        self.covariance_matrices = np.zeros((self.n, self.n))
        self.current_step = 0
        self.storage_max_m = self.smcmax * self.D
        self.storage_init = self.storage_max_m * 0.667
        self.F_results_init = np.full((self.n, self.N), self.storage_init)

    def initialize(self):
        with open(self.json_file, "r") as file:
            data = json.load(file)
        
        yi_init = np.array(0)
        H_results_init = np.array(0)
        F_results_init = self.F_results_init

        self.set_value('Observations', yi_init)
        self.set_value('Simulations', H_results_init)
        self.set_value('State Variables', F_results_init)

    def update(self):
        self.enkf.predict(self.F_results)
        self.enkf.update(self.yi, self.H_results)
        self.state_estimates = self.enkf.get_state_estimate()
        self.covariance_matrices = self.enkf.get_covariance()
        self.current_step += 1

    def finalize(self):
        pass

    def get_component_name(self):
        return "BMIEKFilter"

    def get_input_item_count(self):
        return 3

    def get_input_var_names(self):
        return ("Observations", "Simulations", "State Variables")

    def get_output_item_count(self):
        return 2

    def get_output_var_names(self):
        return ("state_estimates", "covariance_matrices")

    def get_var_type(self, name):
        if name in self.get_output_var_names():
            return "double"
        raise ValueError(f"Unknown variable name {name}")

    def get_value(self, name):
        if name == "state_estimates":
            return self.state_estimates
        elif name == "covariance_matrices":
            return self.covariance_matrices
        raise ValueError(f"Unknown variable name {name}")

    def set_value(self, name, value_array):
        if name == 'Observations':
            self.yi = value_array
        elif name == 'Simulations':
            self.H_results = value_array
        elif name == 'State Variables':
            self.F_results = value_array
        else:
            raise ValueError(f"Unknown variable name {name}")

    def get_current_time(self):
        return self.current_step

    def get_start_time(self):
        return 0

    def get_end_time(self):
        return self.n_steps

    def get_time_step(self):
        return 1

    def get_time_units(self):
        return "s"
    
    def get_grid_edge_count(self):
        raise NotImplementedError()

    def get_grid_edge_nodes(self, grid_edge_id):
        raise NotImplementedError()

    def get_grid_face_count(self):
        raise NotImplementedError()

    def get_grid_face_edges(self, grid_face_id):
        raise NotImplementedError()

    def get_grid_face_nodes(self, grid_face_id):
        raise NotImplementedError()

    def get_grid_node_count(self):
        raise NotImplementedError()

    def get_grid_nodes_per_face(self, grid_face_id):
        raise NotImplementedError()

    def get_grid_origin(self, grid_id):
        raise NotImplementedError()

    def get_grid_shape(self, grid_id):
        raise NotImplementedError()

    def get_grid_spacing(self, grid_id):
        raise NotImplementedError()

    def get_grid_x(self, grid_id):
        raise NotImplementedError()

    def get_grid_y(self, grid_id):
        raise NotImplementedError()

    def get_grid_z(self, grid_id):
        raise NotImplementedError()

    def get_var_location(self, name):
        raise NotImplementedError()
    
    def get_grid_rank(self, grid_id):
        raise NotImplementedError()

    def get_grid_size(self, grid_id):
        raise NotImplementedError()

    def get_grid_type(self, grid_id):
        raise NotImplementedError()

    def get_value_at_indices(self, name, inds):
        raise NotImplementedError()

    def get_value_ptr(self, name):
        raise NotImplementedError()

    def get_var_grid(self, name):
        raise NotImplementedError()

    def get_var_itemsize(self, name):
        raise NotImplementedError()

    def get_var_nbytes(self, name):
        raise NotImplementedError()

    def get_var_units(self, name):
        raise NotImplementedError()

    def set_value_at_indices(self, name, inds, src):
        raise NotImplementedError()

    def update_until(self, time):
        raise NotImplementedError()