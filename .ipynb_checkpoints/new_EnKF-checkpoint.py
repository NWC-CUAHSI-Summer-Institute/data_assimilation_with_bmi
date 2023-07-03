import numpy as np

class EnKF:
    def __init__(self, n, m, R, Q, N, P, x0=None):
        self.n = n  # number of states
        self.m = m  # number of measurements
        self.R = R  # measurement noise
        self.Q = Q  # process noise
        self.N = N  # number of ensembles
        self.state_estimates = np.zeros((self.n))
        self.covariance_matrices = np.zeros((self.n, self.n))
        self.current_step = 0
        self.ensembles = np.random.multivariate_normal(x0 if x0 is not None else np.zeros(self.n), P, size=self.N).T

    def predict(self, F_results):
        if len(F_results) < self.N:
            raise ValueError("Insufficient number of simulations in F_results")
    
        for i in range(self.N):
            self.ensembles[:, i] = F_results[i]
            self.ensembles[:, i] += np.random.multivariate_normal(np.zeros(self.n), self.Q)


    def update(self, yi):
        if np.any(np.isnan(yi)):
            return np.mean(self.ensembles, axis=1)  # return current state estimate without updating

        y_ensembles = np.random.multivariate_normal(np.zeros(self.m), self.R, size=self.N).T

        ensemble_mean = np.mean(self.ensembles, axis=1)[:, np.newaxis]
        y_ensembles_mean = np.mean(y_ensembles, axis=1)[:, np.newaxis]
        Pxy = (self.ensembles - ensemble_mean).dot((y_ensembles - y_ensembles_mean).T) / (self.N - 1)
        Pyy = np.cov(y_ensembles)

        K = Pxy.dot(np.linalg.inv(Pyy + self.R))
        self.ensembles += K.dot((yi - y_ensembles_mean.T))
        return np.mean(self.ensembles, axis=1)


    def get_state_estimate(self):
        return np.mean(self.ensembles, axis=1)

    def get_covariance(self):
        return np.cov(self.ensembles)
