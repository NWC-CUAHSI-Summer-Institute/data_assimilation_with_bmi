import numpy as np

# Contact:
# Dr. Matt Bartos: mdbartos@utexas.edu

def muskingum_matrix(startnodes, endnodes, alpha, beta, chi, gamma, indegree):
    m = startnodes.size
    n = endnodes.size
    indegree_t = indegree.copy()

    A = np.zeros((n, n), dtype=np.float64)
    B = np.zeros((n, n), dtype=np.float64)
    
    # Simulate output
    for k in range(m):
        startnode = startnodes[k]
        endnode = endnodes[startnode]
        while(indegree_t[startnode] == 0):
            alpha_i = alpha[startnode]
            beta_i = beta[startnode]
            chi_i = chi[startnode]
            gamma_i = gamma[startnode]
            alpha_j = alpha[endnode]
            beta_j = beta[endnode]
            # System matrix
            A[startnode, startnode] = chi_i
            B[startnode, startnode] = gamma_i
            # Add outflow to inflow at endnode
            if startnode != endnode:
                A[endnode, startnode] += alpha_j
                A[endnode] += beta_j * A[startnode]
                B[endnode] += beta_j * B[startnode]
            indegree_t[endnode] -= 1
            startnode = endnode
            endnode = endnodes[startnode]
    return A, B