import numpy as np
import multiprocessing
from models.bqmin import bqmin
from models.q_generator import generate_subspace_matrix
from models.quad import fitfroquad
from models.algorithm_6_4_poised import algorithm_6_4
from utils.common import BaseStarSolver
from message.consolemessage import (print_progress_header, print_progress_row, print_final_summary)
from utils.io import (log_progress_to_csv)

class DiagHessianStarSolver(BaseStarSolver):
    def evaluate_interpolation_points(self, xk, Q_k, delta_k, p, n_samples, track_objfun, eval_counter, point_counter):
        """
        Generates and evaluates 2p + 1 interpolation points used to construct
        a diagonal Hessian quadratic model in the current subspace.

        The first point is the current iterate xk.
        The remaining 2p points are symmetric perturbations along ±e_i directions
        in the subspace defined by Q_k.

        Parameters:
            xk (np.ndarray): Current point in R^n (base point for interpolation)
            Q_k (np.ndarray): Subspace matrix of shape (n x p)
            delta_k (float): Trust region radius
            p (int): Subspace dimension
            n_samples (int): Number of Monte Carlo samples per point
            track_objfun (bool): Whether to print each evaluation to console
            eval_counter (list[int]): Global counter for function evaluations
            point_counter (list[int]): Global counter for evaluated points

        Returns:
            y (np.ndarray): Function values at the 2p+1 interpolation points
            points (List[np.ndarray]): List of interpolation points in R^n
        """
        y = np.zeros(2 * p + 1)
        points = []
        y[0], _ = self.monte_carlo_estimate(xk, n_samples, self.montecarlo_log_file, track_objfun=track_objfun, eval_counter=eval_counter, point_counter=point_counter)
        points.append(xk.copy())
        for i in range(p):
            ei = np.zeros(p)
            ei[i] = 1
            fwd = xk + Q_k @ (delta_k * ei)
            bwd = xk - Q_k @ (delta_k * ei)
            y[1 + i], _ = self.monte_carlo_estimate(fwd, n_samples, self.montecarlo_log_file, track_objfun=track_objfun, eval_counter=eval_counter, point_counter=point_counter)
            y[1 + p + i], _ = self.monte_carlo_estimate(bwd, n_samples, self.montecarlo_log_file, track_objfun=track_objfun, eval_counter=eval_counter, point_counter=point_counter)
            points.extend([fwd, bwd])
        return y, points

    def build_C_matrix(self, delta_k, p):
        """
        Constructs matrix C used to fit a diagonal Hessian quadratic model.

        The model is built using 2p+1 function evaluations:
        - One at the center (xk)
        - p in the +delta_k direction
        - p in the -delta_k direction

        The matrix C has size (2p+1 x 2p+1), and each row corresponds to one interpolation point.
        Each row is composed of:
        - A constant term (1)
        - A linear term (±delta_k)
        - A quadratic term (0.5 * delta_k^2) on the diagonal (Hessian is assumed diagonal)

        Parameters:
            delta_k (float): Trust region radius used for interpolation spacing
            p (int): Subspace dimension (number of basis directions)

        Returns:
            C (np.ndarray): The (2p+1 x 2p+1) matrix used for model fitting
        """
        C = np.zeros((2 * p + 1, 2 * p + 1))
        C[0, 0] = 1
        for i in range(p):
            C[1 + i, 0] = 1
            C[1 + i, 1 + i] = delta_k
            C[1 + i, 1 + p + i] = 0.5 * delta_k**2
            C[1 + p + i, 0] = 1
            C[1 + p + i, 1 + i] = -delta_k
            C[1 + p + i, 1 + p + i] = 0.5 * delta_k**2
        return C
    
    def update_adaptive_model(self, xk, points_old, y_old, delta_k, p, n_samples,
                            eval_counter, point_counter, track_objfun, jlm_type):
        """
        Expands the current quadratic model by increasing the subspace dimension from `p` to `p+1`,
        reusing old interpolation points and evaluating only two new ones along the new direction.

        This adaptive model update is triggered after a failed trust-region step, allowing the 
        optimizer to improve model quality in a larger subspace without discarding prior computations.

        Parameters:
            xk (np.ndarray): Current iterate in R^n
            points_old (List[np.ndarray]): Previously evaluated 2p+1 interpolation points in R^n
            y_old (np.ndarray): Function values at those interpolation points
            delta_k (float): Current trust region radius
            p (int): Previous subspace dimension
            n_samples (int): Number of Monte Carlo samples per point
            eval_counter (list[int]): Global evaluation counter (mutable)
            point_counter (list[int]): Unique point counter (mutable)
            track_objfun (bool): Whether to print function evaluations
            jlm_type (int): Type of subspace generation (e.g., Haar = 3, Identity = 1)

        Returns:
            Q_k_new (np.ndarray): Updated subspace matrix of shape (n x (p+1))
            G_sub_new (np.ndarray): Estimated gradient vector in subspace (length p+1)
            H_sub_new (np.ndarray): Estimated diagonal Hessian matrix in subspace ((p+1) x (p+1))
            s_star_new (np.ndarray): Optimal step in subspace (length p+1)
            s_k_new (np.ndarray): Corresponding step in full space (length n)
            f0_new (float): Estimated function value at xk
            y_new (np.ndarray): Combined function values at all 2(p+1)+1 points
            points_new (List[np.ndarray]): Combined interpolation points in R^n
        """
        p_new = p + 1
        q_hat = np.sqrt(1 + 1 / p)  # Scaling for new direction

        # Step 1: Generate a new subspace matrix of size (n x (p+1))
        Q_k_new = generate_subspace_matrix(jlm_type, self.n, p_new)

        # Step 2: Reuse old interpolation points and values (2p+1 points)
        X_old = points_old
        y_old = np.array(y_old)

        # Step 3: Create two new points along ±e_{p+1} in subspace
        e_p1 = np.zeros(p_new)
        e_p1[-1] = 1
        s_plus = q_hat * delta_k * e_p1
        s_minus = -q_hat * delta_k * e_p1

        x_plus = xk + Q_k_new @ s_plus
        x_minus = xk + Q_k_new @ s_minus

        # Step 4: Evaluate function at only the two new points
        y_plus, _ = self.monte_carlo_estimate(
            x_plus, n_samples,
            track_objfun=track_objfun,
            eval_counter=eval_counter,
            point_counter=point_counter
        )
        y_minus, _ = self.monte_carlo_estimate(
            x_minus, n_samples,
            track_objfun=track_objfun,
            eval_counter=eval_counter,
            point_counter=point_counter
        )

        # Step 5: Combine old and new interpolation data
        points_new = X_old + [x_plus, x_minus]
        y_new = np.concatenate([y_old, [y_plus, y_minus]])

        # Step 6: Fit new diagonal quadratic model from updated points
        C_new = self.build_C_matrix(delta_k, p_new)
        h = np.linalg.solve(C_new, y_new)
        f0_new = h[0]
        G_sub_new = h[1:1 + p_new]
        H_diag_new = h[1 + p_new:]
        H_sub_new = np.diag(H_diag_new)

        # Step 7: Minimize the new quadratic model within the trust region
        L = -delta_k / np.sqrt(p_new) * np.ones(p_new)
        U = -L
        s_star_new, _ = bqmin(H_sub_new, G_sub_new, L, U)
        s_k_new = Q_k_new @ s_star_new

        return Q_k_new, G_sub_new, H_sub_new, s_star_new, s_k_new, f0_new, y_new, points_new

    def solve(self):
        """
        Main optimization loop for the ANASTAARS solver with optional adaptive subspace model.
        Handles trust region updates, subspace model fitting (with diagonal Hessian),
        Monte Carlo evaluation, adaptive expansion, and progress logging.
        """
        # Initialization from options
        xk = self.x0
        p = self.options.get("p", 2)
        p0 = p  # Save initial subspace dimension
        delta_k = self.options.get("delta0", 1.0)
        eta1 = self.options.get("eta1", 0.01)
        eta2 = self.options.get("eta2", 0.9)
        gamma = self.options.get("gamma", 2.0)
        max_eval = self.options.get("max_eval", 1000)
        max_iters = self.options.get("max_iters", 100)
        n_samples = self.options.get("mc_samples", 1)
        jlm_type = self.options.get("jlm_type", 3)
        r = 0 if n_samples == 1 else self.options.get("r", 1.0)

        # Flags and counters
        track_progress = self.options.get("track_progress", False)
        track_objfun = self.options.get("track_objfun", False)
        use_parallel = self.options.get("parallel", False)
        cpus = self.options.get("cpu", multiprocessing.cpu_count())
        adaptive_subspace = self.options.get("adaptive_subspace", False)

        eval_counter = [0]
        point_counter = [0]
        p_history = []

        # Bounds
        lower_bound = self.lower_bound
        upper_bound = self.upper_bound

        # Adaptive reuse cache
        reuse_model_from_adaptive = False
        cached_model = None

        # Print table header
        if track_progress and not track_objfun:
            print_progress_header()

        # Main loop
        k = 0
        while True:
            k += 1
            if max_iters is not None and k > max_iters:
                break
            if self._check_eval_limit(eval_counter, max_eval):
                break

            p_history.append(p)

            # --- Reuse adaptive model if flagged ---
            if adaptive_subspace and reuse_model_from_adaptive:
                Q_k, G_sub, H_sub, s_star, s_k, f0, y, points = cached_model
                reuse_model_from_adaptive = False
            else:
                # Step 1: Build new subspace model from scratch
                Q_k = generate_subspace_matrix(jlm_type, self.n, p)
                y, points = self.evaluate_interpolation_points(
                    xk, Q_k, delta_k, p, n_samples,
                    track_objfun, eval_counter, point_counter
                )
                if self._check_eval_limit(eval_counter, max_eval):
                    break

                # Step 2: Fit diagonal Hessian model using interpolation
                C = self.build_C_matrix(delta_k, p)
                h = np.linalg.solve(C, y)
                f0 = h[0]
                G_sub = h[1:1 + p]
                H_diag = h[1 + p:]
                H_sub = np.diag(H_diag)

                # Step 3: Solve trust-region subproblem in subspace
                L = -delta_k / np.sqrt(p) * np.ones(p)
                U = -L
                s_star, _ = bqmin(H_sub, G_sub, L, U)
                s_k = Q_k @ s_star

            # --- Bounds handling ---
            within_bounds = True
            if lower_bound is not None and upper_bound is not None:
                within_bounds = np.all((xk + s_k) >= lower_bound) and np.all((xk + s_k) <= upper_bound)

            if within_bounds:
                # Step 4: Evaluate f(xk) and f(xk + s_k)
                f_k0, std_k = self.monte_carlo_estimate(
                    xk, n_samples,
                    self.montecarlo_log_file,
                    use_parallel, cpus,
                    track_objfun, eval_counter, point_counter
                )
                if self._check_eval_limit(eval_counter, max_eval):
                    break

                f_ks, _ = self.monte_carlo_estimate(
                    xk + s_k, n_samples,
                    self.montecarlo_log_file,
                    use_parallel, cpus,
                    track_objfun, eval_counter, point_counter
                )
                if self._check_eval_limit(eval_counter, max_eval):
                    break

                # Step 5: Compute model prediction and ratio
                model_pred = f0 + G_sub @ s_star + 0.5 * s_star.T @ H_sub @ s_star
                rho_k = (f_k0 - f_ks + r * std_k) / (f0 - model_pred)
                normG = np.linalg.norm(G_sub)
                success = rho_k >= eta1 and normG >= eta2 * delta_k

                if success:
                    xk = xk + s_k
                    delta_k = min(gamma * delta_k, self.options.get("delta_max", 5.0))
                    p = p0  # Reset to base subspace dim
                else:
                    if adaptive_subspace:
                        if p >= self.n:
                            print(f"Iteration {k}: Maximum subspace dimension p = {p} reached (n = {self.n}). Shrinking delta_k instead.")
                            delta_k /= gamma
                            continue
                        Q_k, G_sub, H_sub, s_star, s_k, f0, y, points = self.update_adaptive_model(
                            xk, points, y, delta_k, p, n_samples,
                            eval_counter, point_counter,
                            track_objfun, jlm_type
                        )
                        p += 1
                        cached_model = (Q_k, G_sub, H_sub, s_star, s_k, f0, y, points)
                        reuse_model_from_adaptive = True
                    else:
                        delta_k /= gamma
            else:
                print(f"Iteration {k}: Step out of bounds. Rejecting step and shrinking trust region.")
                delta_k /= gamma
                continue

            # --- Logging: Console + CSV ---
            if track_progress and not track_objfun:
                print_progress_row(k, f_ks, normG, delta_k, rho_k, eval_counter[0])

            if self.progress_log_file:
                log_progress_to_csv(
                    self.progress_log_file, k, xk, eval_counter[0], delta_k,
                    G_sub, H_sub, rho_k, success, Q_k,
                    points, s_star, f0, f_ks, model_pred
                )

        # --- Final summary ---
        print_final_summary(xk, f0 if 'f0' in locals() else None,
                            G_sub if 'G_sub' in locals() else None,
                            H_sub if 'H_sub' in locals() else None,
                            eval_counter[0], point_counter[0], "Diagonal Hessian")

        return {
            "xmin": xk,
            "fval": float(f0) if 'f0' in locals() else None,
            "grad": G_sub if 'G_sub' in locals() else None,
            "hess": H_sub if 'H_sub' in locals() else None,
            "evals": eval_counter[0],
            "points": point_counter[0]
        }


class FrobStarSolver(BaseStarSolver):
    def update_adaptive_model_frobenius(self, xk, poised_d_old, interp_points_old, f_poised_old,
                                        delta_k, p, n_samples, eval_counter, point_counter, jlm_type):
        """
        Expands the Frobenius subspace model by increasing subspace dimension from `p` to `p+1`
        and evaluates two new points to maintain poised interpolation geometry.

        Parameters:
        - xk: Current iterate (n-dimensional vector)
        - poised_d_old: Previous poised directions in subspace (2p+1, p)
        - interp_points_old: List of previous interpolation points in R^n
        - f_poised_old: Function values at previous interpolation points (array of length 2p+1)
        - delta_k: Trust region radius
        - p: Current subspace dimension
        - n_samples: Number of Monte Carlo samples per evaluation
        - eval_counter, point_counter: Counters for function evaluations and unique points
        - jlm_type: Type of subspace matrix (e.g., Haar, identity, etc.)

        Returns:
        - Q_k_new: Updated (n x (p+1)) subspace matrix
        - poised_d_combined: Updated interpolation directions in subspace (2(p+1)+1, p+1)
        - interp_points_combined: Updated list of interpolation points in R^n
        - f_poised_combined: Function values at all interpolation points
        - G_sub, H_sub: Gradient and Hessian estimates for the quadratic model
        - s_star: Step in subspace
        - s_k: Step in ambient space
        """
        # Step 1: Increase subspace dimension
        p_new = p + 1
        q_hat = np.sqrt(1 + 1 / p)  # Scaling factor for the new direction

        # Step 2: Generate new subspace matrix Q of size (n x p+1)
        Q_k_new = generate_subspace_matrix(jlm_type, self.n, p_new)

        # Step 3: Pad old directions with zero columns to align dimensions
        poised_d_padded = np.zeros((len(poised_d_old), p_new))
        poised_d_padded[:, :p] = poised_d_old

        # Step 4: Construct new directions along ±e_{p+1}
        e_p1 = np.zeros(p_new)
        e_p1[-1] = 1  # Last dimension of the new subspace
        s_plus = q_hat * delta_k * e_p1
        s_minus = -q_hat * delta_k * e_p1

        # Step 5: Evaluate function at new interpolation points
        y_plus = xk + Q_k_new @ s_plus
        y_minus = xk + Q_k_new @ s_minus

        f_plus, _ = self.monte_carlo_estimate(y_plus, n_samples, eval_counter=eval_counter, point_counter=point_counter)
        f_minus, _ = self.monte_carlo_estimate(y_minus, n_samples, eval_counter=eval_counter, point_counter=point_counter)

        # Step 6: Combine old and new data
        poised_d_combined = np.vstack([poised_d_padded, s_plus, s_minus])  # New directions (2(p+1)+1, p+1)
        interp_points_combined = interp_points_old + [y_plus, y_minus]     # New points in R^n
        f_poised_combined = np.concatenate([f_poised_old, [f_plus, f_minus]])  # Corresponding function values

        # Step 7: Fit Frobenius quadratic model to new data
        G_sub, H_sub = fitfroquad(poised_d_combined, f_poised_combined)

        # Step 8: Solve bound-constrained trust-region subproblem in subspace
        L = -delta_k / np.sqrt(p_new) * np.ones(p_new)
        U = -L  # Symmetric bounds
        s_star, _ = bqmin(H_sub, G_sub, L, U)

        # Step 9: Map subspace step back to full space
        s_k = Q_k_new @ s_star

        # Step 10: Return all components
        return Q_k_new, poised_d_combined, interp_points_combined, f_poised_combined, G_sub, H_sub, s_star, s_k


    def solve(self):
        """
        Main optimization loop for the STAARS solver with optional adaptive subspace model.
        Handles trust region updates, model fitting using Frobenius modeling,
        Monte Carlo evaluation, and progress logging.
        """
        # Initialization from options
        xk = self.x0
        p = self.options.get("p", 2)
        p0 = p  # initial base value of subspace dimension
        delta_k = self.options.get("delta0", 1.0)
        eta1 = self.options.get("eta1", 0.01)
        eta2 = self.options.get("eta2", 0.9)
        gamma = self.options.get("gamma", 2.0)
        max_eval = self.options.get("max_eval", 1000)
        max_iters = self.options.get("max_iters", 100)
        n_samples = self.options.get("mc_samples", 1)
        jlm_type = self.options.get("jlm_type", 3)

        # Flags
        track_progress = self.options.get("track_progress", False)
        track_objfun = self.options.get("track_objfun", False)
        use_parallel = self.options.get("parallel", False)
        cpus = self.options.get("cpu", 1)
        adaptive_subspace = self.options.get("adaptive_subspace", False)

        # Counters
        eval_counter = [0]
        point_counter = [0]
        p_history = []

        # Bounds
        lower_bound = self.lower_bound
        upper_bound = self.upper_bound

        # Adaptive model reuse flag and cache
        reuse_model_from_adaptive = False
        cached_model = None

        # Progress logging header
        if track_progress and not track_objfun:
            print_progress_header()

        k = 0
        while True:
            k += 1
            if max_iters is not None and k > max_iters:
                break
            if self._check_eval_limit(eval_counter, max_eval):
                break

            p_history.append(p)

            # --- Reuse adaptive model if flagged ---
            if adaptive_subspace and reuse_model_from_adaptive:
                Q_k, poised_d, interp_points, f_poised, G_sub, H_sub, s_star, s_k = cached_model
                CurObj = f_poised[0]
                reuse_model_from_adaptive = False
            else:
                # Build interpolation set and evaluate initial point
                Q_k = generate_subspace_matrix(jlm_type, self.n, p)
                f0_scalar, _ = self.monte_carlo_estimate(
                    xk, n_samples,
                    self.montecarlo_log_file,
                    use_parallel=use_parallel,
                    cpus=cpus,
                    track_objfun=track_objfun,
                    eval_counter=eval_counter,
                    point_counter=point_counter
                )
                f0 = np.array([f0_scalar])

                if self._check_eval_limit(eval_counter, max_eval):
                    break

                poised_d, _ = algorithm_6_4(np.zeros((1, p)), delta_k, f0)
                interp_points = []
                f_poised = np.zeros(len(poised_d))
                for i in range(len(f_poised)):
                    y_i = xk + Q_k @ poised_d[i]
                    interp_points.append(y_i.copy())
                    f_val, _ = self.monte_carlo_estimate(
                        y_i, n_samples, self.montecarlo_log_file,
                        use_parallel=use_parallel, cpus=cpus,
                        track_objfun=track_objfun,
                        eval_counter=eval_counter,
                        point_counter=point_counter
                    )
                    f_poised[i] = f_val
                    if self._check_eval_limit(eval_counter, max_eval):
                        break

                G_sub, H_sub = fitfroquad(poised_d, f_poised)
                CurObj = f0.item()
                L = -delta_k / np.sqrt(p) * np.ones(p)
                U = -L
                s_star, _ = bqmin(H_sub, G_sub, L, U)
                s_k = Q_k @ s_star

            # Bound handling
            within_bounds = True
            if lower_bound is not None and upper_bound is not None:
                within_bounds = np.all((xk + s_k) >= lower_bound) and np.all((xk + s_k) <= upper_bound)

            if within_bounds:
                fks, _ = self.monte_carlo_estimate(
                    xk + s_k, n_samples,
                    self.montecarlo_log_file,
                    use_parallel=use_parallel, cpus=cpus,
                    track_objfun=track_objfun,
                    eval_counter=eval_counter,
                    point_counter=point_counter
                )

                if self._check_eval_limit(eval_counter, max_eval):
                    break

                mks = CurObj + G_sub @ s_star + 0.5 * s_star @ H_sub @ s_star
                rho_k = (CurObj - fks) / (CurObj - mks) if (CurObj - mks) != 0 else 0

                normG = np.linalg.norm(G_sub)
                success = rho_k >= eta1 and normG >= eta2 * delta_k

                if success:
                    xk += s_k
                    delta_k = min(gamma * delta_k, self.options.get("delta_max", 5.0))
                    p = p0
                else:
                    if adaptive_subspace:
                        if p >= self.n:
                            print(f"Iteration {k}: Maximum subspace dimension p = {p} reached (n = {self.n}). Shrinking delta_k instead.")
                            delta_k /= gamma
                            continue
                        Q_k, poised_d, interp_points, f_poised, G_sub, H_sub, s_star, s_k = self.update_adaptive_model_frobenius(
                            xk, poised_d, interp_points, f_poised,
                            delta_k, p, n_samples,
                            eval_counter, point_counter,
                            jlm_type
                        )
                        cached_model = (Q_k, poised_d, interp_points, f_poised, G_sub, H_sub, s_star, s_k)
                        p += 1
                        reuse_model_from_adaptive = True
                    else:
                        delta_k /= gamma
            else:
                print(f"Iteration {k}: Step out of bounds. Rejecting step and shrinking trust region.")
                delta_k /= gamma
                continue

            if track_progress and not track_objfun:
                print_progress_row(k, fks, normG, delta_k, rho_k, eval_counter[0])

            if self.progress_log_file:
                log_progress_to_csv(
                    self.progress_log_file, k, xk, eval_counter[0], delta_k,
                    G_sub, H_sub, rho_k, success, Q_k,
                    interp_points, s_star, f0, fks, mks
                )

        print_final_summary(xk, float(f0) if 'f0' in locals() else None,
                            G_sub if 'G_sub' in locals() else None,
                            H_sub if 'H_sub' in locals() else None,
                            eval_counter[0], point_counter[0], "Frobenius")

        return {
            "xmin": xk,
            "fval": float(f0) if 'f0' in locals() else None,
            "grad": G_sub if 'G_sub' in locals() else None,
            "hess": H_sub if 'H_sub' in locals() else None,
            "evals": eval_counter[0],
            "points": point_counter[0]
        }
