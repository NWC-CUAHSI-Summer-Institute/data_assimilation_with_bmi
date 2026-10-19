# EnKF.py
import numpy as np

class EnKF:
    def __init__(self, n, m, R, Q, N, P, smcmax, D):
        self.n = n  # number of states
        self.m = m  # number of measurements
        self.R = R  # measurement noise
        self.Q = Q  # process noise
        self.N = N  # number of ensembles
        self.P = P
        self.state_estimates = np.zeros((self.n))
        self.covariance_matrices = np.zeros((self.n, self.n))
        self.current_step = 0
        self.smcmax = smcmax
        self.D = D
        self.storage_max_m = self.smcmax * self.D
        self.storage_init = self.storage_max_m * 0.667
        self.ensembles = np.full((self.n, self.N), self.storage_init)

    def predict(self, F_results):
        error_factor_sim = 0.005
        for i in range(self.N):
            perturbation_factor_sim = np.random.standard_normal()
            self.ensembles[:, i] = F_results[i]
            self.ensembles[:, i] += (perturbation_factor_sim * F_results[i] * error_factor_sim)
        return self.ensembles

    def update(self, yi, H_results):
        if np.any(np.isnan(yi)):
            return np.mean(self.ensembles, axis=1)  # return current state estimate without updating

        y_ensembles = np.zeros((self.m, self.N))
        error_factor_sm = 0.03
        for i in range(self.N):
            perturbation_factor_sm = np.random.standard_normal()
            y_ensembles[:, i] = H_results[i]
            y_ensembles[:, i] += (perturbation_factor_sm * H_results[i] * error_factor_sm)

        ensemble_mean = np.mean(self.ensembles, axis=1)[:, np.newaxis]
        y_ensembles_mean = np.mean(y_ensembles, axis=1)[:, np.newaxis]
        Pxy = (self.ensembles - ensemble_mean).dot((y_ensembles - y_ensembles_mean).T) / (self.N - 1)
        Pyy = np.cov(y_ensembles, bias=True)
        K = Pxy.dot(np.linalg.pinv(Pyy + self.R))
        self.ensembles += K.dot(yi - y_ensembles_mean)
        return np.mean(self.ensembles, axis=1)

    def get_state_estimate(self):
        return np.mean(self.ensembles, axis=1)
    
    def get_covariance(self):
        ensemble_mean = np.mean(self.ensembles, axis=1)[:, np.newaxis]
        return (self.ensembles - ensemble_mean).dot((self.ensembles - ensemble_mean).T) / (self.N - 1)
