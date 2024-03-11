# Generating random noisy channels
- `random_noise_N2`: $\mathcal{J}_{\mathrm{max}}^{(k)}$ for different types of strategies under randomly generated noisy channels.

# Impact of prior distributions
- `GHZ_diffdelta_N2` and `GHZ_diffdelta_N3`: Data used to compare the performance of $\mathcal{J} _{\max}^{\left( i \right)}$ and $\mathcal{J} \left( \rho \right)$ when $N=2$ and $N=3$. Here, $\mathcal{J}\left( \rho \right)$ means $\mathrm{tr}\left( \bar{\rho}S^2 \right)$ when the initial state is $\rho=\left( |0\rangle ^{\otimes N}+|1\rangle ^{\otimes N} \right) /\sqrt{2}$ and the adopted strategy is of type $i$.
- `SeqCtrl_diffprior_N2` and `SeqCtrl_diffprior_N3`: Data used to compare the performance of sequential strategies with and without control when the channel $\mathcal{E}_\theta$ is unitary.
- `noise_beta_diffa_N2`, `noise_gaussian_diffdelta_N2`, `noise_gaussian_diffmean_N2` and `noise_GM_diffw_N2`: Data used to study the impact of different prior distributions on the performance of different types of strategies.

# Strict hierarchy
- `maxsdp_noise_strategyii`, `maxsdp_noise_strategyiii`, `maxsdp_noise_strategyiv`, `minsdp_noise_strategyi`, `minsdp_noise_strategyii` and `minsdp_noise_strategyiii`: Data used to check the strict hierarchy between different types of strategies when $\mathcal{E} _{\theta}=\mathcal{E} ^{(\mathrm{AD)}}\circ \mathcal{E} ^{(\mathrm{BF)}}\circ \mathcal{U} _{\theta}$.
- `maxsdp_unitary_strategyi_N2`, `maxsdp_unitary_strategyii_N2`, `maxsdp_unitary_strategyii_N3`, `minsdp_unitary_strategyi_N2`, `minsdp_unitary_strategyi_N3` and `minsdp_unitary_strategyii_N2`: Data obtained by solving primal and dual SDPs for the unitary channel $\mathcal{U}_\theta$.
- `NoisePhiandH`, `UnitaryPhiandH_N2` and `UnitaryPhiandH_N3`: Data to get $\Phi$ and $H$ for noisy channel $\mathcal{E} _{\theta}=\mathcal{E} ^{(\mathrm{AD)}}\circ \mathcal{E} ^{(\mathrm{BF)}}\circ \mathcal{U} _{\theta}$ and unitary channel when $N=2,3$.
