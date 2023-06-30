from bmipy import Bmi
import numpy as np
import json
from EnKF import EnKF


class BMIEnKF(Bmi):

    def __init__(self, json_file):
        with open(json_file, "r") as file:
            data = json.load(file)

        self.n = data['n']
        self.m = data['m']
        self.R =  np.eye(self.m) * np.array(data['R'])
        self.Q = np.eye(self.n) * np.array(data['Q'])
        self.N = data['N']
        self.P = np.eye(self.n) * np.array(data['P'])
        self.x0 = np.zeros(self.n)

        self.enkf = EnKF(self.n, self.m, self.R, self.Q, self.N, self.x0, self.P)
        self.state_estimates = np.zeros((self.n))
        self.covariance_matrices = np.zeros((self.n, self.n))
        self.current_step = 0

    def initialize(self):
        pass

    def update(self, yi, F_results, H_results):
        self.enkf.predict(F_results)
        self.enkf.update(yi, H_results)
        self.state_estimates = self.enkf.get_state_estimate()
        self.covariance_matrices = self.enkf.get_covariance()
        self.current_step += 1

    def update_until(self, y, F_results, H_results):
        for yi in y.T:
            self.update(yi, F_results, H_results)

    def finalize(self):
        pass

    def get_component_name(self):
        return "BMIEKFilter"

    def get_input_item_count(self):
        return 1

    def get_input_var_names(self):
        return ("observations",)

    def get_output_item_count(self):
        return 2

    def get_output_var_names(self):
        return ("state_estimates", "covariance_matrices")

    def get_grid_rank(self):
        raise NotImplementedError()

    def get_grid_size(self):
        raise NotImplementedError()

    def get_grid_type(self):
        raise NotImplementedError()

    def get_var_grid(self, name):
        raise NotImplementedError()

    def get_var_type(self, name):
        if name in self.get_output_var_names():
            return "double"
        raise ValueError(f"Unknown variable name {name}")

    def get_var_itemsize(self, name):
        return np.dtype(self.get_var_type(name)).itemsize

    def get_var_nbytes(self, name):
        return self.get_var_itemsize(name) * np.prod(np.shape(self.get_value(name)))

    def get_var_units(self, name):
        if name in self.get_output_var_names():
            return "-"
        raise ValueError(f"Unknown variable name {name}")

    def get_value(self, name):
        if name == "state_estimates":
            return self.state_estimates
        elif name == "covariance_matrices":
            return self.covariance_matrices
        raise ValueError(f"Unknown variable name {name}")

    def get_value_at_indices(self, name, inds):
        raise NotImplementedError()

    def get_value_ptr(self, name):
        raise NotImplementedError()

    def set_value(self, name, value_array):
        raise NotImplementedError()

    def set_value_at_indices(self, name, inds, src):
        raise NotImplementedError()

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