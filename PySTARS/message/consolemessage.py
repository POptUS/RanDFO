import numpy as np

def print_eval_console(eval_id, x_str, fx, best_fx, from_cache):
    src = "cache" if from_cache else "compute"
    print(f"[Eval {eval_id}] f({x_str}) = {fx:.6f} | f(best so far) = {best_fx:.6f} [{src}]")

def print_eval_console_monte(eval_id, point_id, sample_id, x):
    print(f"Function eval {eval_id} at point {point_id} has f = {sample_id:.8f} at x = {x}")

def print_progress_header():
    print(f"{'Iter':>5} {'Obj':>10} {'Grad':>10} {'Delta':>10} {'rho':>8} {'Evals':>6}")

def print_progress_row(k, fks, normG, delta_k, rho_k, evals):
    print(f"{k:5d} {fks:10.2e} {normG:10.2e} {delta_k:10.2e} {rho_k:8.2e} {evals:6d}")

def print_final_summary(xk, fval, grad, hess, evals, points, header):
    print(f"\n***** {header} Model Final Summary *****")
    print(f"Solution xmin = {np.round(xk, 4)}")
    if fval is not None:
        print(f"Objective value f(xmin) = {float(fval):.9e}")
    else:
        print("Objective value not available.")
    print(f"Total evaluations = {evals} (at {points} points)")
    if grad is not None:
        print(f"Gradient = {np.round(grad, 4)}")
    else:
        print("Gradient not available.")
    print("Hessian =")
    print(np.round(hess, 4) if hess is not None else "Hessian not available.")
