# EnKF.py
import numpy as np

class EnKF:
    def __init__(self, n, m, R, Q, N, x0, P):
        self.n = n  # number of states
        self.m = m  # number of measurements
        self.R = R  # measurement noise
        self.Q = Q  # process noise
        self.N = N  # number of ensembles
        self.ensembles = np.random.multivariate_normal(x0, P, size=self.N).T

    def predict(self, F_results):
        for i in range(self.N):
            self.ensembles[:, i] = F_results[i]
            self.ensembles[:, i] += np.random.multivariate_normal(np.zeros(self.n), self.Q)
        return self.ensembles

    def update(self, y, H_results):
        if np.any(np.isnan(y)):
            return np.mean(self.ensembles, axis=1)  # return current state estimate without updating

        y_ensembles = np.zeros((self.m, self.N))
        for i in range(self.N):
            y_ensembles[:, i] = H_results[i]
            y_ensembles[:, i] += np.random.multivariate_normal(np.zeros(self.m), self.R)

        ensemble_mean = np.mean(self.ensembles, axis=1)[:, np.newaxis]
        y_ensembles_mean = np.mean(y_ensembles, axis=1)[:, np.newaxis]
        Pxy = (self.ensembles - ensemble_mean).dot((y_ensembles - y_ensembles_mean).T) / (self.N - 1)
        Pyy = np.cov(y_ensembles)

        K = Pxy.dot(np.linalg.inv(Pyy + self.R))
        self.ensembles += K.dot(y[:, np.newaxis] + np.random.normal(np.zeros(self.m), self.R) - y_ensembles)
        return np.mean(self.ensembles, axis=1)

    def get_state_estimate(self):
        return np.mean(self.ensembles, axis=1)
    
    def get_covariance(self):
        return np.cov(self.ensembles)
