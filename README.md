# Quantum-Walks-and-Quantum-Monte-Carlo
The Quantum Galton Board (QGB) is capable of simulating statistical distributions, and can be used to solve differential equations. This project enables the generation of QGBs which simulate (1) normal (2) exponential (3) quantum Hadamard random walk distributions.


## Summary: 

This project enables the generation of a general QGB with arbitrary left ($p$) and right ($q$) drop probabilities for each peg, for any number of layers. It also provides functions tailored to generate exponential and Quantum Hadamard Random Walk (QHRD) distributions.

A walkthrough of all simulations is found in `qgb_simulations.ipynb`. The functions to generate each QGB are in `qgb_functions.py`. The two-page summary is in the folder `Presentations and Summary`.

### Normal distribution

The function `generate_qgb(n)` generates a quantum circuit for a QGB with $n$ layers. This specific function forces the left and right drop probabilities at each peg to be equal i.e. $p=0.5$. An example of how to use this function is shown below:

```
# number of layers 
n = 2

# create circuit
circuit = generate_qgb(n)

# draw circuit
circuit.draw(output="mpl")
```

The function `generate_custom_qgb(n, thetas)` generates a more general quantum circuit. The additional argument thetas is a list of rotation angles corresponding to the left drop probabilities on each peg.
$[\theta_1, \theta_2, \theta_3]$ corresponds $[p_1,p_2,p_3]$ where $\theta_i = 2\arcsin(\sqrt{p_i})$. Each peg in the board is numbered from top to bottom, left to right. For example, the previous circuit can be reproduced in the following way:

```
# number of layers 
n = 2

# rotation angles
angles = [np.pi/2, np.pi/2, np.pi/2]

# create circuit
circuit = generate_custom_qgb(n)

# draw circuit
circuit.draw(output="mpl")
```

This circuit will generate a normal distribution.

### Exponential Distribution

The function `generate_custom_qgb()` can also generate an exponential distribution. In particular, `generate_thetas(n,p)` forms the required rotation angles for an exponential distribution $\sim \lambda e^{-\lambda x}$ given a probabilitiy $p=1-e^{-\lambda/m}$ where $1/m$ is the size of the bucket spacing. 

```
# number of layers 
n = 6

# lambda in exponential distribution
lam = 1.5
# chooses the size of each bucket
m = 4 

# left probabilitiy
p = 1 - np.exp(-lam/m)
thetas = generate_thetas(n,p)

# Generate circuit
circuit = generate_custom_qgb(n, thetas=thetas)
```

### Hadamard Quantum Walk Distribution 

The Hadamard Quantum Walk distirbution is generated with `generate_hadamard_qgb(n, sym=True)`. $n$ is the number of layers (steps) and `sym` is the initial condition of the ancilla qubit. If set to `True` $|\psi\rangle = \frac{|0\rangle + i|1\rangle}{2}$. Otherwise if `False` $|\psi\rangle = |0\rangle$.

```
# number of layers 
n = 12

# create circuit
circuit = generate_hadamard_qgb(n, sym=True)
```

### Running in qiskit

This project uses qiskit 2.0. The circuits can be run in the following way: 

```
# Backend selection
backend = Aer.get_backend('qasm_simulator')

# Executing the scheme on the selected backend
job = backend.run(circuit, shots=20000, memory=True)

# Get the result object
result = job.result()

# Get counts
counts = result.get_counts(circuit)
```

## Other Information

Project Name: Quantum Walks and Quantum Monte Carlo 

Team Name: angusQuantum

Team Members: Angus Crookes

ID: gst-J55s2DFyZglsL4G

Presentation deck location: Quantum-Walks-and-Quantum-Monte-Carlo/Presentations and Summary/presentation_deck



