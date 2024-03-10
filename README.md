# Strict hierarchy of optimal strategies for global estimations: Mapping global estimations into local ones

This repository serves as a companion to the research paper with the same name.

# Structure of this repository

The file `Data` is used to store related data. Symbolic codes are stored in the file `MATHEMATICA` which are used to provide a computer-assisted level proof for a strict hierarchy between different types of strategies. The current file stores MATLAB codes.

# Basic requirement 

- The `CVX` package should be installed. See [this website](http://cvxr.com/cvx/doc/install.html) for guidance.

- The `QETLAB`, a MATLAB toolbox for exploring quantum entanglement theory, is required. The package can be downloaded by following [this documentation](https://qetlab.com/Installation).

- The symbolic calculation using MATHEMATICA requires [`MATLink`](http://matlink.org/) to invoke MATLAB functions without leaving MATHEMATICA. One can also avoid using MATLink by directly using the data we generated in `Data`.

# List of functions
## Basic operations 

- `Max_SDP`: Primal SDP.
- `Min_SDP`: Dual SDP.
- `random_noisy_channel`: Randomly generate noisy channels to check the strict hierarchy between different types of strategies.
- `diffprior`: Study the impact of different prior distributions on the ultimate precision limit.
- `SeqNoCtrol`: Calculate the performance of sequential strategies without control.
- `Noise_diff_dist`: Study the impact of different prior distributions for noisy channels.

## Helper functions

- `myoperation`: Implement the operation $_Q\tilde{X}\coloneqq \mathrm{tr}_Q\tilde{X}\otimes \left( \mathbb{I} _Q/d_Q \right)$.
- `mysylvester`: Solve the equation $\overline{\theta C}+\\{H,\bar{C}\\}=0$.
- `myrandomChannel`: Randomly generate quantum channels $\mathcal{E}:\mathcal{L}(\mathcal{H})\rightarrow \mathcal{L}(\mathcal{H})$. Explicitly, we first generate a random unitary $U$ using [QETLAB](https://qetlab.com/RandomUnitary), uniformly according to Haar measure. Then we get Kraus operators from this unitary $U$ by Stinespring representation of $\mathcal{E}$.
- `plot_noisy`: Used to generate Fig. S2.
- `plot_unexp`: Used to generate Fig. S3.

