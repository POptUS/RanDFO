---
title: PySTARS
nav_order: 1
---

# PySTARS: Python Implementation of STARS Optimization Algorithms

PySTARS is a Python-based implementation of STARS (Stochastic Trust-region Algorithm in Random Subspaces) for derivative-free optimization. It supports both deterministic and stochastic objective functions, subspace trust-region models, adaptive sampling, and comprehensive logging.

---

## Features

- Deterministic and stochastic objective function support
- Monte Carlo estimation with optional parallelism and caching
- Trust-region based subspace model
- Adaptive subspace expansion logic
- Diagonal Hessian and Frobenius models
- Logging of evaluations, progress, and objective values
- Terminal output with tracking options
- Full bounds support

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/POptUS/RanDFO.git
cd RanDFO/PySTARS
```

---

### 2. Create Virtual Environment

```bash
python -m venv venv
```

#### Activate the Virtual Environment

- **Windows (Command Prompt):**
  ```cmd
  venv\Scripts\activate
  ```

- **Windows (Git Bash / PowerShell):**
  ```bash
  source venv/Scripts/activate
  ```

- **Linux / macOS:**
  ```bash
  source venv/bin/activate
  ```

#### Install the required dependencies

```bash
pip install -r requirements.txt
```

---

## Usage Example

### 1. Define Objective Functions

```python
import numpy as np

# Noise-free deterministic function
def rosenbrock(x):
    return np.sum(100.0 * (x[1:] - x[:-1]**2)**2 + (1 - x[:-1])**2)

# Noisy version of the same function
def rosenbrock_noisy(x):
    return rosenbrock(x) * (1.0 + 0.01 * np.random.randn())
```

---

### 2. Initial Point and Bounds

```python
x0 = np.full(10, 0.0)  # 10-dimensional input
lower = np.full(10, -1.2)
upper = np.full(10, 5.0)
```

---

### 3. Solver Options (Default)

```python
options = {
    "p": 2,
    "delta0": 1.0,
    "delta_max": 5.0,
    "eta1": 0.01,
    "eta2": 0.9,
    "gamma": 2.0,
    "r": 1.0,
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

```

---

### 4. Running the Solver

#### Deterministic (Noise-Free) Mode for Diagonal Hessian and Frobenius Modal

```python
from core.solver import DiagHessianStarSolver, FrobStarSolver

diagsolver = DiagHessianStarSolver(f=rosenbrock, x0=x0)
frobsolver = FrobStarSolver(f=rosenbrock, x0=x0)
solution = diagsolver.solve()
```

#### Deterministic (Noise-Free) Mode for Diagonal Hessian and Frobenius Modal with options and bounds

```python
from core.solver import DiagHessianStarSolver, FrobStarSolver

options = {"max_eval": 100, "max_iters": 10}
diagsolver = DiagHessianStarSolver(f=rosenbrock, x0=x0, options=options, bounds=(lower, upper))
frobsolver = FrobStarSolver(f=rosenbrock, x0=x0, options=options, bounds=(lower, upper))
solution = diagsolver.solve()
```

#### Stochastic (Noisy) Mode for Diagonal Hessian and Frobenius Modal with options and bounds

```python
from core.solver import DiagHessianStarSolver, FrobStarSolver

options = {"stochastic": True}

diagsolver = DiagHessianStarSolver(f=rosenbrock_noisy, x0=x0, bounds=(lower, upper), options=options, f_true=rosenbrock)
frobsolver = DiagHessianStarSolver(f=rosenbrock_noisy, x0=x0, bounds=(lower, upper), options=options, f_true=rosenbrock)
solution = diagsolver.solve()
```

---

## Logging Output

All logs are saved to the `logs/` directory:

| File                              | Description                                      |
|-----------------------------------|--------------------------------------------------|
| `deterministic_eval_log.csv`      | Logs evaluations for deterministic runs          |
| `montecarlo_eval_log.csv`         | Logs sample-by-sample evaluations for MC runs    |
| `progress_log.csv`                | Tracks trust-region size, model updates, etc.    |
| `monte_carlo_estimates_log.csv`   | Logs mean/std of MC estimates per point          |
| `objfun_log.txt`                  | Text log of all sampled values (detailed)        |

---

## Console Display Options

| Option            | Description                                                        |
|-------------------|--------------------------------------------------------------------|
| `track_objfun`    | Prints each evaluation result to terminal                          |
| `track_progress`  | Prints each iteration’s trust region state and progress metrics    |

> Only one of these can be `True` at a time to avoid output conflicts.

---

## Available Solvers

| Solver Class             | Description                                      |
|--------------------------|--------------------------------------------------|
| `DiagHessianStarSolver`  | Uses diagonal Hessian-based model                |
| `FrobStarSolver`         | Uses Frobenius based model                       |

---

## Questions or Issues?

- Submit a GitHub [Issue](https://github.com/POptUS/RanDFO/issues)
- Start a [Discussion](https://github.com/POptUS/RanDFO/discussions)