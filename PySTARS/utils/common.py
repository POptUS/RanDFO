import numpy as np
import multiprocessing
from collections import defaultdict
from joblib import Parallel, delayed
from utils.io import (log_eval_deterministic, log_eval_montecarlo, log_objfun_to_txt, log_mc_to_csv, clear_log_dir)
from message.consolemessage import (print_eval_console, print_eval_console_monte)
from exception.exceptionhandling import MissingTrueFunctionError, ConflictingTrackingOptionsError
clear_log_dir("logs")  # Start clean

class BaseStarSolver:
    def __init__(self, f, x0, options=None, bounds=None, f_true=None):
        """
        Initializes the Frobenius and Diagonal Hessian STAR solver.
        """
        self.f = f
        self.f_true = f_true if f_true is not None else f
        self.x0 = np.array(x0, dtype=float)
        self.options = options or {}
        self.n = len(x0)

        # Extract flags
        self.stochastic = self.options.get("stochastic", False)
        self.track_objfun = self.options.get("track_objfun", False)
        self.track_progress = self.options.get("track_progress", False)

        # Conflict checking
        if self.track_objfun and self.track_progress:
            raise ConflictingTrackingOptionsError()

        # Raise error if in stochastic mode but f_true is not provided
        if self.stochastic:
            if f_true is None:
                raise MissingTrueFunctionError()
            if f == f_true:
                raise MissingTrueFunctionError("f_true must be the *true* (noise-free) function, not the same as the noisy function.")           
        else:
            if f_true is not None:
                raise MissingTrueFunctionError("In deterministic mode, only one function (f) should be provided. Remove f_true. Set 'stochastic' to True to enable stochastic mode.")

        # Bounds
        self.lower_bound, self.upper_bound = bounds if bounds else (None, None)

        # Evaluation caches and best-tracking
        self.deterministic_cache = {}
        self.cache = defaultdict(dict)
        self.eval_display_counter = 0
        self.best_fx = float("inf")
        self.best_true_fx = float("inf")
        self.sample_cache = set()
        self.true_fx_cache = {}

        # Logging paths
        self.deterministic_log_file = "logs/deterministic_eval_log.csv" if self.options.get("deterministic_log", False) else None
        self.montecarlo_log_file = "logs/monte_carlo_estimates_log.csv" if self.options.get("montecarlo_log", False) else None
        self.montecarlo_eval_log_file = "logs/montecarlo_eval_log.csv" if self.options.get("montecarlo_eval_log", False) else None
        self.progress_log_file = "logs/progress_log.csv" if self.options.get("progress_log", False) else None
        self.track_objfun_txt_file = "logs/objfun_log.txt" if self.options.get("objfun_log", False) else None

        # Random seed
        self.seed = self.options.get("seed", 42)
        if self.seed is not None:
            np.random.seed(self.seed)
            import random
            random.seed(self.seed)

        # Force mc_samples = 1 for deterministic mode
        if not self.stochastic and self.options.get("mc_samples", 1) > 1:
            print("Warning: 'stochastic' is False but 'mc_samples' > 1. Forcing mc_samples = 1.")
            self.options["mc_samples"] = 1

    def _check_eval_limit(self, eval_counter, max_eval):
        return max_eval and eval_counter[0] >= max_eval

    def unwrap_fx(self, fx):
        """Extracts scalar from list/tuple/ndarray; returns as-is if already scalar."""
        return fx[0] if isinstance(fx, (tuple, list, np.ndarray)) else fx
    
    def deterministic_mc_eval(self, x, eval_counter, point_counter, track_objfun):
        """
        Evaluates f(x) deterministically with cache reuse and optional logging.

        Used when `mc_samples = 1` (deterministic mode). Reuses cached values if available,
        updates evaluation counters, logs to console and file if requested.

        Parameters:
            x (np.ndarray): Evaluation point.
            eval_counter (list[int]): Counter for unique evaluations.
            point_counter (list[int]): Counter for total evaluation points.
            track_objfun (bool): Whether to log to console.
        
        Returns:
            (float, float): Function value and dummy standard deviation (always 0.0).
        """
        max_eval = self.options.get("max_eval", None)
        if self._check_eval_limit(eval_counter, max_eval):
            return float("inf"), 0.0  # Stop if over eval limit

        # Prepare keys for cache and display
        x_key = tuple(round(float(val), 8) for val in x)
        x_str = tuple(round(float(i), 4) for i in x)

        # Initialize evaluation ID counter (used in logs)
        if not hasattr(self, "eval_display_counter"):
            self.eval_display_counter = 0
        eval_id = self.eval_display_counter

        # Ensure deterministic cache exists
        if not hasattr(self, "deterministic_cache"):
            self.deterministic_cache = {}

        # Try retrieving from cache
        if x_key in self.deterministic_cache:
            fx = self.deterministic_cache[x_key]
            from_cache = True
        else:
            fx = self.f(x)
            self.deterministic_cache[x_key] = fx
            from_cache = False
            eval_counter[0] += 1  # Count new evaluations only

        point_counter[0] += 1  # Every evaluated point (new or cached)

        # Track the best value seen so far
        if not hasattr(self, "best_fx"):
            self.best_fx = fx
        else:
            self.best_fx = min(self.best_fx, fx)

        # Print to terminal if requested
        if track_objfun:
            print_eval_console(eval_id, x_str, fx, self.best_fx, from_cache)

        # Log to CSV if enabled and not from cache
        if self.deterministic_log_file and not from_cache:
            log_eval_deterministic(
                filepath=self.deterministic_log_file,
                eval_id=eval_counter[0],
                x_str=x_str,
                fx=fx,
                best_fx=self.best_fx
            )

        # Update eval display ID (for terminal)
        self.eval_display_counter += 1

        return fx, 0.0  # Standard deviation is 0 for deterministic
    
    def monte_carlo_estimate(self, x, N, montecarlo_log_file=None, use_parallel=False,
                            cpus=1, track_objfun=False, eval_counter=None, point_counter=None):
        """
        Estimates f(x) using N Monte Carlo samples, with optional parallelism, caching,
        logging, and console output.

        Parameters:
            x (np.ndarray): Evaluation point
            N (int): Number of Monte Carlo samples
            montecarlo_log_file (str): Optional log file path
            use_parallel (bool): Whether to parallelize sample evaluation
            cpus (int): Number of CPUs to use if parallel
            track_objfun (bool): Whether to log to console
            eval_counter, point_counter (list[int]): Mutable counters for logging

        Returns:
            (float, float): Mean and standard deviation of the estimate
        """
        max_eval = self.options.get("max_eval", None)
        if self._check_eval_limit(eval_counter, max_eval):
            return float("inf"), 0.0  # Early stopping

        # Fall back to deterministic eval if not stochastic
        if not self.options.get("stochastic") == True and N <= 1:
            return self.deterministic_mc_eval(x, eval_counter, point_counter, track_objfun)

        x_key = tuple(round(float(val), 8) for val in x)
        x_key_str = str(tuple(round(float(val), 4) for val in x))

        # Initialize cache entry if new point
        if x_key not in self.cache:
            self.cache[x_key] = {"mean": 0.0, "std": 0.0, "count": 0, "samples": [], "reuse_count": 0}

        entry = self.cache[x_key]
        old_samples = entry["samples"][:]
        old_count = entry["count"]
        old_mean = entry["mean"]

        # Run N evaluations (parallel or serial)
        if use_parallel:
            if cpus is None:
                cpus = multiprocessing.cpu_count()
            new_samples = Parallel(n_jobs=cpus)(delayed(lambda: self.unwrap_fx(self.f(x)))() for _ in range(N))
        else:
            new_samples = [self.unwrap_fx(self.f(x)) for _ in range(N)]

        # Update estimates
        combined_samples = old_samples + new_samples
        total_count = old_count + N
        new_mean = (old_mean * old_count + sum(new_samples)) / total_count
        std_dev = np.std(combined_samples, ddof=1) if total_count > 1 else 0.0

        # Cache updated statistics
        entry.update({
            "count": total_count,
            "mean": new_mean,
            "std": std_dev,
            "samples": combined_samples,
            "reuse_count": entry["reuse_count"] + 1
        })

        # Optional convergence logging (f(x), f(true), f(best), etc.)
        if self.montecarlo_eval_log_file:
            if not hasattr(self, "best_true_fx"):
                self.best_true_fx = float("inf")
            if not hasattr(self, "sample_cache"):
                self.sample_cache = set()
            if not hasattr(self, "true_fx_cache"):
                self.true_fx_cache = {}

            log_eval_montecarlo(
                filepath=self.montecarlo_eval_log_file,
                x=x,
                new_samples=new_samples,
                f_true=self.f_true,
                eval_counter=eval_counter,
                true_fx_cache=self.true_fx_cache,
                sample_cache=self.sample_cache,
                best_true_fx=self.best_true_fx
            )
            self.best_true_fx = min(self.best_true_fx, self.true_fx_cache[str(tuple(round(float(val), 4) for val in x))])

        # track if we have already logged this point
        logged_points = set()
        for i in range(N):
            point_str = str(np.round(x, 4))  # Create a unique identifier for the point
            if point_str not in logged_points:
                if track_objfun:
                    print_eval_console_monte(eval_counter[0] + 1, point_counter[0]+1, new_samples[i], np.round(x, 4))
                logged_points.add(point_str)  # Mark the point as logged
            if self.track_objfun_txt_file:
                log_objfun_to_txt(
                    filepath=self.track_objfun_txt_file,
                    x=x,
                    fx=new_samples[i],
                    eval_id=eval_counter[0] + 1,
                    point_id=point_counter[0] + 1,
                    sample_num=i + 1
                )
            eval_counter[0] += 1
        point_counter[0] += 1

        # Optional full Monte Carlo logging
        if montecarlo_log_file:
            log_mc_to_csv(
                filepath=montecarlo_log_file,
                x_key_str=x_key_str,
                old_samples=old_samples,
                new_samples=new_samples,
                combined=combined_samples,
                entry=entry
            )

        return new_mean, std_dev