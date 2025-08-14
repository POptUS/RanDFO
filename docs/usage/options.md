---
title: Options
parent: Usage
nav_order: 4
---

# Solver Options

| Option | Default | Meaning |
|---|---:|---|
| `p` | 2 | Subspace dimension |
| `delta0` | 1.0 | Initial trust-region radius |
| `delta_max` | 5.0 | Maximum trust-region radius |
| `delta_end` | 1e-6 | Termination radius (if used) |
| `eta1` | 0.01 | Accept step threshold (low) |
| `eta2` | 0.9 | Accept step threshold (high) |
| `gamma` | 2.0 | Radius shrink/expand factor |
| `r` | 1.0 | Model regularization / scaling (if applicable) |
| `jlm_type` | 3 | Subspace basis type (identity matrix(n==p):2, Scaled Haar orthogonal(√(n/p) * Q):3, Haar orthogonal(Q):4) |
| `max_eval` | 1000 | Max function evaluations |
| `max_iters` | 100 | Max iterations |
| `mc_samples` | 1 | Monte-Carlo samples per estimate |
| `parallel` | False | Use joblib parallelism |
| `cpu` | 1 | Cores for parallelism |
| `adaptive_subspace` | False | Enable p-expansion after failures |
| `stochastic` | False | Treat objective as noisy |
| `seed` | 42 | RNG seed |
| **Logging Options** |||
| `deterministic_log` | False | Save deterministic evaluation log to CSV |
| `montecarlo_eval_log` | False | Save Monte-Carlo evaluation log to CSV |
| `montecarlo_log` | False | Save detailed Monte-Carlo step log |
| `progress_log` | False | Save progress per iteration to CSV |
| `objfun_log` | False | Save objective function evaluations to text file |
| **Terminal Display** |||
| `track_objfun` | False | Print each objective function evaluation to console |
| `track_progress` | False | Print iteration progress to console |


### Example
```python
import numpy as np
from core.solver import DiagHessianStarSolver, FrobStarSolver

# Noise-free deterministic function
def rosenbrock(x):
    return np.sum(100.0 * (x[1:] - x[:-1]**2)**2 + (1 - x[:-1])**2)

options = {
    "p": 2,
    "delta0": 1.0,
    "delta_max": 5.0,
    "eta1": 0.01,
    "eta2": 0.9,
    "gamma": 2.0,
    "jlm_type": 3,
    "max_eval": 1000,
    "max_iters": 100,
    "mc_samples": 1,
    "parallel": False,
    "cpu": 1,
    "adaptive_subspace": False,
    "stochastic": False,
    "seed": 42,

    # Logging options
    "deterministic_log": False,
    "montecarlo_eval_log": False,
    "montecarlo_log": False,
    "progress_log": False,
    "objfun_log": False,

    # Terminal display
    "track_objfun": False,
    "track_progress": False
}

# This is for the Diagonal Hessian Modal
diagsolver = DiagHessianStarSolver(f=rosenbrock, x0=x0, options=options)
# This is for the Frobenius Modal
frobsolver = FrobStarSolver(f=rosenbrock, x0=x0, options=options)
solution = diagsolver.solve()

```