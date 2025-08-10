# Quantum-Walks-and-Quantum-Monte-Carlo
The Quantum Galton Board (QGB) is capable of simulating statistical distributions, and is used to model complex random processes, and solve differential equations. This project enables the generation of QGBs which simulate (1) normal (2) exponential (3) quantum Hadamard random walk distributions.

## Information

Project Name: Quantum Walks and Quantum Monte Carlo 

Team Name: angusQuantum

Team Members: Angus Crookes

ID: gst-J55s2DFyZglsL4G

Presentation deck location: Quantum-Walks-and-Quantum-Monte-Carlo/Presentations and Summary

## Summary: 

The Galton Board (GB) is a device capable of simulating different statistical distributions. Its quantum analogue the Quantum Galton Board (QGB) can potentially process calculations faster, offering a more efficient way to model complex random processes, and solve differential equations. This project enables the generation of a general QGB with arbitrary left ($p$) and right ($q$) drop probabilities for each peg, for any number of layers. It also provides functions tailored to generate exponential and Quantum Hadamard Random Walk (QHRD) distributions.

### QGB with n-layers

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

