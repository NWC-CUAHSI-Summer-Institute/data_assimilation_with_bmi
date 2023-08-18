import numpy as np
from scipy.stats import norm



class particle_filter:
    def __init__(self, N, n, m, smcmax, D):
        self.N = N  # number of ensembles
        self.n = n  # number of states
        self.current_step = 0
        self.smcmax = smcmax
        self.D = D
        self.m = m
        self.storage_max_m = self.smcmax * self.D
        self.storage_init = self.storage_max_m * 0.667
        self.outputs = np.zeros((self.m, self.N))
   #     self.state_estimates = np.zeros((self.n))
        self.likelihoods = np.zeros(self.N)
    #    self.observation = np.zeros((self.n))
        
    def add_noise (self, simulations):
        perturbation_factor_sim = np.random.randn()
        error_factor_sim = 0.2
        for i in range(self.N):
            self.outputs[:, i] = simulations[i]
            self.outputs[:, i] += (perturbation_factor_sim * simulations[i] * error_factor_sim)
        return self.outputs

    def calculate_likelihood(self):
        # Estimate the parameters of the Gaussian distribution
        mu, sigma = norm.fit(observation)
        # Calculate the PDF for each number
        self.likelihoods = [norm.pdf(number, loc=mu, scale=sigma) for number in self.outputs]
        return self.likelihoods


    def execute_pf (self):
        weights = self.likelihoods / self.likelihoods.sum()
        combined = list(zip(self.state_variables, self.outputs, weights))
        sorted_combined = sorted(combined, key=lambda x: x[2])
        sorted_states = [x[0] for x in sorted_combined]
        self.state_estimates = np.random.choice(sorted_states, size=1, replace=True, p=weights)
        return self.state_estimates