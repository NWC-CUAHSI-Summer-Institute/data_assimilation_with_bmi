from bmipy import Bmi
import numpy as np
import json
from pf_method import particle_filter


class particle_filter_bmi:
    def __init__(self, json_file):
        with open(json_file, "r") as file:
            data = json.load(file)

        self.n = data['n']
        self.N = data['N']
        self.m = data['m']
        self.smcmax = data['smcmax']
        self.D = data['D']
        self.likelihoods = np.zeros(self.N)
        self.outputs = np.zeros((self.n))
   #     self.simulations = np.zeros((self.m, self.N))
        self.pf = particle_filter(self.N, self.n, self.m, self.smcmax, self.D, self.simulations)
        self.state_estimates = np.zeros((self.n))
        self.observation = np.zeros((self.n))
        self.current_step = 0

    def initialize(self):
        pass

    def update(self):
        particle_filter.add_noise(self)
        particle_filter.calculate_likelihood(self)
        particle_filter.execute_pf(self)
        self.current_step += 1

    def finalize(self):
        pass

    def get_component_name(self):
        return "PF_DA"

    def get_input_item_count(self):
        return 4

    def get_input_var_names(self):
        return ("Observations", "Simulations", "State Variables")

    def get_output_item_count(self):
        return 2

    def get_output_var_names(self):
        return ("state_estimates")

    def get_var_type(self, name):
        if name in self.get_output_var_names():
            return "double"
        raise ValueError(f"Unknown variable name {name}")

    def get_value(self, name):
        if name == "state_estimates":
            return updated_state

    def set_value(self, name, value_array):
        if name == 'Observations':
            self.observation = value_array
        elif name == 'Simulations':
            self.simulations = value_array
        elif name == 'State Variables':
            self.state_variables = value_array
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