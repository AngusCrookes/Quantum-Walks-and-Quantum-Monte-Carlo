"""
Created on Sat Aug 9th 2025

@author: anguscrookes
"""

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.circuit.library import CSwapGate
import numpy as np 

# generate a circuit representing a QGB with n layers 
def generate_qgb(n):
    
    '''
    Quantum Galton Board with n layers
    
    Generates the quantum circuit for a QGB with n layers. The left, right drop
    probability of each peg is constant and such that p = 0.5, q = 0.5. 
    
    Input Parameters
    ----------
    n :  integer, representing the number of layers

    Returns 
    ----------
    qc : quantum circuit of QGB
    
    '''
    
    # number of qubits
    N = 2*n + 2 

    # number of pegs
    P = n*(n+1)/2
    
    # middle peg 
    M = N // 2
    
    # classical and quantum registers 
    q = QuantumRegister(N)
    c = ClassicalRegister(n+1) # was N
    qc = QuantumCircuit(q,c)

    # gates
    qc.x(q[M])

    # construct each layer 
    for i in range(n):

        qc.h(q[0])

        # construct each peg 
        for p in np.arange(-i,i+2,2):
            qc.append(CSwapGate(), [0, M + p, M + p - 1])
            qc.cx(q[M + p], q[0])
            qc.append(CSwapGate(), [0, M + p, M + p + 1])

            if p < i:
                qc.cx(q[M+p+1], q[0])

        if i < n:
            qc.reset(q[0])
               
    qc.barrier()

    # measure all qubits that represent buckets 
    qc.measure([q[i] for i in np.arange(M-n, M+n+2, 2)], c)
     
    return qc


# generate a circuit representing a QGB with n layers with biased pegs
def generate_custom_qgb(n, thetas):

    """ Quantum Galton Board with n layers and biased pegs
    
    Generates the quantum circuit for a QGB with n layers. The left, right drop
    probability is defined by thetas, which is a list of angles for each rotation
    matrix. Note this function is general, and reproduces the generate_qgb() function 
    when thetas = [np.pi/2, np.pi/2 ...].
    
    Input Parameters
    ----------
    n :  integer, representing the number of layers
    
    thetas: list of angles, for the parameterised rotations of each R_x gate
            Note: the left probability at each peg is such that theta_i = 2arcsin(sqrt(p))

    Returns 
    ----------
    qc : quantum circuit of QGB """
    
    # number of qubits
    N = 2*n + 2 

    # number of pegs
    P = n*(n+1)/2

    # peg index
    p_ind = 0
    
    # middle peg 
    M = N // 2
    
    # classical and quantum registers 
    q = QuantumRegister(N)
    c = ClassicalRegister(n+1) 
    qc = QuantumCircuit(q,c)

    # gates
    qc.x(q[M])

    # construct each layer 
    for i in range(n):

        qc.rx(thetas[p_ind], q[0])
        p_ind += 1
        
        # construct each peg 
        for p in np.arange(-i,i+2,2):
            qc.append(CSwapGate(), [0, M + p, M + p - 1])
            qc.cx(q[M + p], q[0])
            qc.append(CSwapGate(), [0, M + p, M + p + 1])

            if p < i:
                qc.reset(q[0])
                qc.rx(thetas[p_ind], q[0]) 
                p_ind += 1
            
        qc.barrier()

        # construct CNOT alternative
        for k in np.arange(-i+2,i+2,2):
            qc.cx(q[M+k],q[M+k-1])
            qc.reset(q[M+k])
   
        if i < n:
            qc.reset(q[0])

    qc.barrier()

    # measure all qubits that represent buckets 
    qc.measure([q[i] for i in np.arange(M-n, M+n+2, 2)], c)
     
    return qc


def generate_hadamard_qgb(n, sym=True):

    '''
    Quantum Galton Board with n layers
    
    Generates the quantum circuit for a QGB with n layers producing a quantum
    Hadamard walk distribution. 
    
    Input Parameters
    ----------
    n :  integer, representing the number of layers
    
    thetas: list of angles, for the parameterised rotations of each R_x gate
            Note: the left probability at each peg is such that theta_i = 2arcsin(sqrt(p))

    Returns 
    ----------
    qc : quantum circuit of QGB
    
    '''

    
    # number of qubits
    N = 2*n + 2 

    # number of pegs
    P = n*(n+1)/2

    print(N)
    
    # middle peg 
    M = N // 2
    
    # classical and quantum registers 
    q = QuantumRegister(N)
    c = ClassicalRegister(n+1) # was N
    qc = QuantumCircuit(q,c)

    # gates
    qc.x(q[M])

    # construct each layer 
    for i in range(n):

        if i != 0:
            qc.h(q[0])
        elif sym==True:
            qc.h(q[0])
            qc.s(q[0])
            qc.h(q[0])
        elif sym==False:
            qc.h(q[0])

        # construct each peg 
        for p in np.arange(-i,i+2,2):
            qc.append(CSwapGate(), [0, M + p, M + p - 1])
            qc.cx(q[M + p], q[0])
            qc.append(CSwapGate(), [0, M + p, M + p + 1])
             
            qc.cx(q[M+p+1], q[0])
                
    qc.barrier()

    # measure all qubits that represent buckets 
    qc.measure([q[i] for i in np.arange(M-n, M+n+2, 2)], c)
     
    return qc

def generate_thetas(n,p):

    '''
    Generates the rotation angles required for an exponential distribution

    Input parameters:

    n: Number of layers

    p: list of probabilities for each peg
    
    '''
    
    theta_p = 2*np.arcsin(np.sqrt(p))
    thetas = []
    for i in range(n):
        [thetas.append(np.pi) for i in range(i)]
        thetas.append(theta_p)
    return thetas


def normal_distribution(x, mu, std):

    '''
    Normal distribution with mean mu, and standard deviation std
    
    Input Parameters
    ----------
    x :  x data
    
    mu: mean of data

    std: standard deviation of data

    Returns 
    ----------
    normal distribution
    
    '''
    
    return (1 / np.sqrt(2*np.pi*std**2))* np.exp(-0.5*((x-mu)/std)**2)



def exponential_distribution(x, lam):
    
    '''
    Exponential distribution with mean mu, and standard deviation std
    
    Input Parameters
    ----------
    x :  x data
    
    mu: mean of data

    std: standard deviation of data

    Returns 
    ----------
    normal distribution
    '''
    
    return lam*np.exp(-lam*x)



