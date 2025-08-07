import os
import csv
import numpy as np

def clear_log_dir(log_dir="logs"):
    if os.path.exists(log_dir):
        for filename in os.listdir(log_dir):
            file_path = os.path.join(log_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)

def ensure_dir_exists(filepath):
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

def log_eval_deterministic(filepath, eval_id, x_str, fx, best_fx):
    ensure_dir_exists(filepath)
    write_header = not os.path.exists(filepath)
    with open(filepath, "a", newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["eval", "x", "f(x)", "f(best_so_far)", "source"])
            print(f"📄 Saved file to {filepath}")
        writer.writerow([eval_id, x_str, round(fx, 6), round(best_fx, 6), "compute"])

def log_eval_montecarlo(filepath, x, new_samples, f_true, eval_counter,
                         true_fx_cache, sample_cache, best_true_fx):
    ensure_dir_exists(filepath)
    x_key_str = str(tuple(round(float(val), 4) for val in x))
    write_header = not os.path.exists(filepath)

    if x_key_str not in true_fx_cache:
        true_fx = f_true(x)
        true_fx_cache[x_key_str] = round(true_fx, 4)
    else:
        true_fx = true_fx_cache[x_key_str]

    with open(filepath, "a", newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["eval", "x", "f(x)", "f(true)", "f(best_so_far)", "source"])
            print(f"📄 Saved file to {filepath}")
        for i, sample_val in enumerate(new_samples):
            sample_key = (x_key_str, round(float(sample_val), 4))
            src = "cache" if sample_key in sample_cache else "compute"
            sample_cache.add(sample_key)
            #  # Update best using only f(true)
            # best_true_fx = min(best_true_fx, true_fx)
            writer.writerow([
                eval_counter[0] + i,
                x_key_str,
                round(sample_val, 4),
                round(true_fx, 4),
                round(best_true_fx, 4),
                src
            ])

def log_objfun_to_txt(filepath, x, fx, eval_id, point_id, sample_num):
    ensure_dir_exists(filepath)
    file_did_not_exist = not os.path.exists(filepath)
    with open(filepath, "a") as f:
        if file_did_not_exist:
            print(f"📄 Saved file to {filepath}")
        f.write(f"Function eval {eval_id} at point {point_id} sample {sample_num} has f = {fx:.8f} at x = {np.round(x, 4)}\n")

def log_mc_to_csv(filepath, x_key_str, old_samples, new_samples, combined, entry):
    ensure_dir_exists(filepath)
    write_header = not os.path.exists(filepath)
    with open(filepath, "a", newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["x_samples", "old_estimates", "new_estimates", "combined_estimates", "reuse_count", "count", "mean", "std"])
            print(f"📄 Saved file to {filepath}")
        writer.writerow([
            x_key_str,
            [round(float(v), 4) for v in old_samples],
            [round(float(v), 4) for v in new_samples],
            [round(float(v), 4) for v in combined],
            entry["reuse_count"],
            entry["count"],
            round(entry["mean"], 4),
            round(entry["std"], 4)
        ])

def log_progress_to_csv(filepath, k, xk, evals, delta_k, G, H, rho_k, success,
                        Q, interp_points, s_star, f0, fks, mks):
    ensure_dir_exists(filepath)
    write_header = not os.path.exists(filepath)
    with open(filepath, "a", newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["Iter", "xk", "evals", "delta_k", "G", "H", "rho_k", "status",
                             "Q", "interpolation_points", "sk", "f0", "fks", "mks"])
            print(f"📄 Saved file to {filepath}")
        writer.writerow([
            k, str([round(float(xi), 4) for xi in xk]), evals, round(float(delta_k), 4),
            [round(float(gi), 4) for gi in G], [[round(float(hij), 4) for hij in hi] for hi in H],
            round(float(rho_k), 4), "success" if success else "failure",
            [[round(float(qij), 4) for qij in qi] for qi in Q],
            [[round(float(ptj), 4) for ptj in pt] for pt in interp_points],
            [round(float(si), 4) for si in s_star],
            round(float(f0), 4), round(float(fks), 4), round(float(mks), 4)
        ])