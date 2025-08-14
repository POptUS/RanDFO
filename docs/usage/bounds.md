---
title: Bounds
parent: Usage
nav_order: 3
---

# Bounds

```python
import numpy as np
from core.solver import DiagHessianStarSolver, FrobStarSolver

x0     = np.full(10, 0.0)
lower  = np.full(10, -1.2)
upper  = np.full(10,  5.0)
options = {"max_eval": 100, "max_iters": 10}

# This is for the Diagonal Hessian Modal
diagsolver = DiagHessianStarSolver(f=rosenbrock, x0=x0, options=options, bounds=(lower, upper))
# This is for the Frobenius Modal
frobsolver = FrobStarSolver(f=rosenbrock, x0=x0, options=options, bounds=(lower, upper))
solution = diagsolver.solve()

```