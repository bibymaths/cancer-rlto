"""
Cancer-adapted RLTO-inspired model.

This module implements a practical Python framework inspired by the
robustness-load trade-off (RLTO) model described in Choi et al.
(Science Advances, 2026), but reframed for cancer systems biology.

What is retained from the paper:
- A threshold-like fitness landscape with respect to effective protein abundance.
- Overabundance defined as o = C0 / CA, where C0 is baseline abundance and
  CA is the critical abundance threshold.
- The qualitative prediction that lower-expression essential genes tend to need
  larger overabundance to maintain robustness.

What is adapted for cancer:
- Fitness is tumor-cell proliferation instead of bacterial growth.
- Gene dosage, clone fraction, microenvironment stress, and oncogenic pressure
  are included explicitly.
- Toxicity and pathway-specific burden are modeled directly.
- The implementation is intended for simulation, sensitivity analysis,
  and optimization on bulk or single-cell omics inputs.

This is a principled engineering implementation, not a verbatim recovery of the
paper's full supplementary mathematics.
"""

from __future__ import annotations

from matplotlib import pyplot as plt

from dataclasses import dataclass, asdict, fields
from typing import Dict, Optional, Sequence, Tuple

import json
import math

import numpy as np
import pandas as pd

from scipy.optimize import minimize_scalar, minimize
from scipy.special import gammaincc

ArrayLike = Sequence[float] | np.ndarray


@dataclass(slots=True)
class Microenvironment:
    hypoxia: float = 0.0
    nutrient_limitation: float = 0.0
    immune_pressure: float = 0.0
    drug_pressure: float = 0.0
    oxidative_stress: float = 0.0

    def clipped(self) -> "Microenvironment":
        vals = {
            f.name: float(np.clip(getattr(self, f.name), 0.0, 1.0))
            for f in fields(self)
        }
        return Microenvironment(**vals)

    @property
    def stress_index(self) -> float:
        env = self.clipped()
        return float(np.mean([getattr(env, f.name) for f in fields(env)]))


@dataclass(slots=True)
class GeneConfig:
    """Per-gene parameters for the cancer RLTO-inspired model."""

    gene: str
    transcription_rate: float
    translation_efficiency: float = 1.0
    mrna_half_life: float = 6.0
    protein_half_life: float = 12.0
    copy_number: float = 2.0
    clone_fraction: float = 1.0
    essentiality_weight: float = 1.0
    threshold_abundance: Optional[float] = None
    baseline_abundance: Optional[float] = None
    burden_weight: float = 0.02
    toxicity_weight: float = 0.0
    regulation_strength: float = 0.0
    burstiness: float = 1.0
    stress_sensitivity: float = 0.2
    oncogenic_boost: float = 0.0
    pathway: str = "generic"

    def validate(self) -> None:
        if self.transcription_rate <= 0:
            raise ValueError(f"{self.gene}: transcription_rate must be > 0")
        if self.translation_efficiency <= 0:
            raise ValueError(f"{self.gene}: translation_efficiency must be > 0")
        if self.copy_number <= 0:
            raise ValueError(f"{self.gene}: copy_number must be > 0")
        if not (0 <= self.clone_fraction <= 1):
            raise ValueError(f"{self.gene}: clone_fraction must be in [0,1]")
        if self.essentiality_weight < 0:
            raise ValueError(f"{self.gene}: essentiality_weight must be >= 0")
        if self.burden_weight < 0 or self.toxicity_weight < 0:
            raise ValueError(f"{self.gene}: burden/toxicity weights must be >= 0")
        if self.burstiness <= 0:
            raise ValueError(f"{self.gene}: burstiness must be > 0")


@dataclass(slots=True)
class ModelResult:
    gene: str
    baseline_abundance: float
    threshold_abundance: float
    overabundance: float
    abundance_mean: float
    abundance_std: float
    robustness_probability: float
    burden_cost: float
    toxicity_cost: float
    net_fitness: float
    pathway: str


class CancerRLTOModel:
    """
    Cancer-adapted threshold/noise/burden model.
    """

    def __init__(
            self,
            genes: Sequence[GeneConfig],
            microenvironment: Optional[Microenvironment] = None,
            hill_n: float = 12.0,
            base_proliferation: float = 1.0,
            burden_global: float = 1.0,
            toxicity_global: float = 1.0,
            robustness_scale: float = 10.0,
            random_seed: Optional[int] = 42,
    ) -> None:
        self.genes = list(genes)
        for g in self.genes:
            g.validate()
        self.microenvironment = (microenvironment or Microenvironment()).clipped()
        self.hill_n = float(hill_n)
        self.base_proliferation = float(base_proliferation)
        self.burden_global = float(burden_global)
        self.toxicity_global = float(toxicity_global)
        self.robustness_scale = float(robustness_scale)
        self.rng = np.random.default_rng(random_seed)

    @staticmethod
    def _decay_rate_from_half_life(half_life: float) -> float:
        half_life = max(float(half_life), 1e-9)
        return math.log(2.0) / half_life

    def expected_abundance(self, gene: GeneConfig) -> float:
        if gene.baseline_abundance is not None:
            return float(max(gene.baseline_abundance, 1e-9))

        mrna_decay = self._decay_rate_from_half_life(gene.mrna_half_life)
        protein_decay = self._decay_rate_from_half_life(gene.protein_half_life)

        synthesis = (
                gene.transcription_rate
                * gene.translation_efficiency
                * gene.copy_number
                * (1.0 + gene.oncogenic_boost)
                * gene.clone_fraction
        )

        regulation_gain = 1.0 + 0.35 * gene.regulation_strength
        stress_penalty = 1.0 - gene.stress_sensitivity * self.microenvironment.stress_index
        stress_penalty = max(stress_penalty, 0.05)

        abundance = synthesis * regulation_gain * stress_penalty / max(mrna_decay * protein_decay, 1e-9)
        return float(max(abundance, 1e-9))

    def threshold_abundance(self, gene: GeneConfig, baseline: Optional[float] = None) -> float:
        baseline = self.expected_abundance(gene) if baseline is None else float(baseline)
        if gene.threshold_abundance is not None:
            return float(gene.threshold_abundance)

        stress_multiplier = 1.0 + 0.8 * gene.stress_sensitivity * self.microenvironment.stress_index
        essentiality_multiplier = 0.4 + 0.6 * gene.essentiality_weight
        regulation_discount = 1.0 / (1.0 + 0.25 * gene.regulation_strength)

        threshold = 2.5 * np.sqrt(baseline) * stress_multiplier * essentiality_multiplier * regulation_discount
        return float(max(threshold, 1e-9))

    def abundance_distribution_params(self, gene: GeneConfig, mean_abundance: Optional[float] = None) -> Tuple[
        float, float]:
        mu = self.expected_abundance(gene) if mean_abundance is None else float(mean_abundance)

        cv2 = (
                      gene.burstiness / max(gene.transcription_rate * gene.copy_number, 1e-9)
              ) / (1.0 + gene.regulation_strength)
        cv2 *= (1.0 + 0.75 * self.microenvironment.stress_index)
        cv2 = max(cv2, 1e-4)

        shape = 1.0 / cv2
        scale = mu / shape
        return float(shape), float(scale)

    def sample_abundance(self, gene: GeneConfig, n: int = 10000, mean_abundance: Optional[float] = None) -> np.ndarray:
        shape, scale = self.abundance_distribution_params(gene, mean_abundance=mean_abundance)
        return self.rng.gamma(shape=shape, scale=scale, size=int(n))

    def threshold_fitness(self, abundance: ArrayLike, threshold: float, kmax: float = 1.0) -> np.ndarray:
        abundance = np.asarray(abundance, dtype=float)
        thr = max(float(threshold), 1e-9)
        return kmax / (1.0 + (thr / np.maximum(abundance, 1e-12)) ** self.hill_n)

    def burden_cost(self, gene: GeneConfig, mean_abundance: float) -> float:
        x = mean_abundance / 1000.0
        return float(self.burden_global * gene.burden_weight * x)

    def toxicity_cost(self, gene: GeneConfig, mean_abundance: float) -> float:
        x = mean_abundance / 1000.0
        return float(
            self.toxicity_global
            * gene.toxicity_weight
            * (x ** 1.25)
            / (1.0 + gene.regulation_strength)
        )

    def robustness_probability(self, gene: GeneConfig, mean_abundance: Optional[float] = None) -> float:
        """Calculates the exact analytical probability that abundance >= threshold."""
        baseline = self.expected_abundance(gene) if mean_abundance is None else float(mean_abundance)
        threshold = self.threshold_abundance(gene, baseline=baseline)
        shape, scale = self.abundance_distribution_params(gene, mean_abundance=baseline)

        # gammaincc is the exact survival function (1 - CDF) for the Gamma distribution
        return float(gammaincc(shape, threshold / scale))

    def net_gene_fitness(self, gene: GeneConfig, mean_abundance: Optional[float] = None,
                         n_samples: int = 0) -> ModelResult:
        """
        Calculates fitness analytically, entirely eliminating Monte Carlo noise.
        (n_samples is kept in the signature for backwards compatibility but is no longer used).
        """
        baseline = self.expected_abundance(gene) if mean_abundance is None else float(mean_abundance)
        threshold = self.threshold_abundance(gene, baseline=baseline)
        shape, scale = self.abundance_distribution_params(gene, mean_abundance=baseline)

        # 1. Exact Analytical Robustness
        robustness_prob = float(gammaincc(shape, threshold / scale))

        # 2. Expected Growth
        # Because hill_n is steep (12.0), the Hill function acts essentially as a step function.
        # Thus, expected growth is perfectly approximated by the base rate times the fraction of surviving cells.
        expected_growth = self.base_proliferation * robustness_prob

        robustness_term = self.robustness_scale * expected_growth * gene.essentiality_weight

        # 3. Deterministic Costs
        burden = self.burden_cost(gene, baseline)
        toxicity = self.toxicity_cost(gene, baseline)
        net = robustness_term - burden - toxicity

        # Exact mean and standard deviation derived directly from Gamma parameters
        exact_mean = shape * scale
        exact_std = math.sqrt(shape * scale ** 2)

        return ModelResult(
            gene=gene.gene,
            baseline_abundance=baseline,
            threshold_abundance=threshold,
            overabundance=float(baseline / threshold),
            abundance_mean=exact_mean,
            abundance_std=exact_std,
            robustness_probability=robustness_prob,
            burden_cost=burden,
            toxicity_cost=toxicity,
            net_fitness=float(net),
            pathway=gene.pathway,
        )

    def evaluate(self, n_samples: int = 15000) -> pd.DataFrame:
        rows = [asdict(self.net_gene_fitness(g, n_samples=n_samples)) for g in self.genes]
        df = pd.DataFrame(rows)
        df["rank"] = df["net_fitness"].rank(ascending=False, method="dense").astype("Int64")
        return df.sort_values(["net_fitness", "overabundance"], ascending=[False, False]).reset_index(drop=True)

    # -------------------------------------------------------------------------
    # LAYER 1: Gene-Level Optimization (Replaced Grid Search)
    # -------------------------------------------------------------------------
    def optimize_gene_abundance(self, gene: GeneConfig, n_samples: int = 10000) -> Dict[str, float]:
        baseline = self.expected_abundance(gene)

        def objective(x):
            # We want to maximize net_fitness, so we minimize the negative net_fitness
            res = self.net_gene_fitness(gene, mean_abundance=float(x), n_samples=n_samples)
            return -res.net_fitness

        bounds = (1e-6, baseline * 20.0)

        res = minimize_scalar(
            objective,
            bounds=bounds,
            method="bounded"
        )

        best_abundance = float(res.x)
        best_fitness = -float(res.fun)
        thr = self.threshold_abundance(gene, baseline=best_abundance)

        return {
            "gene": gene.gene,
            "baseline_abundance": baseline,
            "optimal_abundance": best_abundance,
            "critical_threshold": thr,
            "optimal_overabundance": float(best_abundance / thr),
            "optimal_net_fitness": best_fitness,
        }

    def optimize_all(self, n_samples: int = 10000) -> pd.DataFrame:
        rows = [self.optimize_gene_abundance(g, n_samples=n_samples) for g in self.genes]
        return pd.DataFrame(rows).sort_values("optimal_net_fitness", ascending=False).reset_index(drop=True)

    # -------------------------------------------------------------------------
    # LAYER 2: Tumor-Level Global Optimization
    # -------------------------------------------------------------------------
    def optimize_global(self, n_samples: int = 0, global_resource_weight: float = 0.5) -> Tuple[np.ndarray, float]:
        """Finds the optimal abundance vector across all genes simultaneously, enforcing carrying capacity."""
        x0 = np.array([self.expected_abundance(g) for g in self.genes])
        baseline_total = float(np.sum(x0))

        def objective(x):
            # 1. Base fitness: mean of individual gene net fitnesses
            vals = [
                self.net_gene_fitness(g, mean_abundance=float(xi)).net_fitness
                for g, xi in zip(self.genes, x)
            ]
            mean_fitness = np.mean(vals)

            # 2. Global carrying capacity penalty
            # Superlinear penalty based on the total protein burden across ALL genes
            total_abundance = np.sum(x)

            # Normalize by baseline total to scale the penalty correctly
            normalized_total = total_abundance / max(baseline_total, 1.0)

            # Quadratic penalty forces a hard trade-off: you cannot maximize all genes at once
            resource_penalty = global_resource_weight * (normalized_total ** 2)

            # Minimize negative fitness
            return -(mean_fitness - resource_penalty)

        bounds = [(1e-6, xi * 20.0) for xi in x0]

        res = minimize(
            objective,
            x0,
            bounds=bounds,
            method="SLSQP"
        )
        return res.x, -res.fun

    # -------------------------------------------------------------------------
    # LAYER 3: Intervention Optimization
    # -------------------------------------------------------------------------
    def simulate_intervention(self, target_gene: str, inhibition: float, n_samples: int = 15000) -> pd.DataFrame:
        inhibition = float(np.clip(inhibition, 0.0, 1.0))
        rows = []
        for gene in self.genes:
            baseline = self.expected_abundance(gene)
            if gene.gene == target_gene:
                baseline *= (1.0 - inhibition)
            rows.append(asdict(self.net_gene_fitness(gene, mean_abundance=baseline, n_samples=n_samples)))
        return pd.DataFrame(rows).sort_values("net_fitness", ascending=False).reset_index(drop=True)

    def optimize_inhibition(self, target_gene: str, n_samples: int = 10000) -> Tuple[float, float]:
        """Finds the precise inhibition level (0 to 1) that minimizes tumor fitness."""

        def objective(u):
            df = self.simulate_intervention(target_gene, inhibition=float(u), n_samples=n_samples)
            # Minimize tumor fitness -> minimize the mean
            return df["net_fitness"].mean()

        res = minimize_scalar(
            objective,
            bounds=(0.0, 1.0),
            method="bounded"
        )
        return res.x, res.fun

    # -------------------------------------------------------------------------
    # Utilities & Setup
    # -------------------------------------------------------------------------
    def pathway_summary(self, n_samples: int = 15000) -> pd.DataFrame:
        df = self.evaluate(n_samples=n_samples)
        return (
            df.groupby("pathway", dropna=False)
            .agg(
                genes=("gene", "count"),
                mean_overabundance=("overabundance", "mean"),
                median_overabundance=("overabundance", "median"),
                mean_net_fitness=("net_fitness", "mean"),
                mean_robustness=("robustness_probability", "mean"),
            )
            .sort_values("mean_net_fitness", ascending=False)
            .reset_index()
        )

    def tumor_fitness_score(self, n_samples: int = 12000) -> float:
        df = self.evaluate(n_samples=n_samples)
        return float(np.mean(np.clip(df["net_fitness"].values, -10.0, 10.0)))

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, microenvironment: Optional[Microenvironment] = None,
                       **kwargs) -> "CancerRLTOModel":
        def f(row: pd.Series, key: str, default: float) -> float:
            val = row.get(key, default)
            return float(default if pd.isna(val) else val)

        def s(row: pd.Series, key: str, default: str) -> str:
            val = row.get(key, default)
            return str(default if pd.isna(val) else val)

        def optf(row: pd.Series, key: str) -> Optional[float]:
            val = row.get(key, np.nan)
            return None if pd.isna(val) else float(val)

        genes = []
        for _, row in df.iterrows():
            genes.append(
                GeneConfig(
                    gene=str(row["gene"]),
                    transcription_rate=float(row["transcription_rate"]),
                    translation_efficiency=f(row, "translation_efficiency", 1.0),
                    mrna_half_life=f(row, "mrna_half_life", 6.0),
                    protein_half_life=f(row, "protein_half_life", 12.0),
                    copy_number=f(row, "copy_number", 2.0),
                    clone_fraction=f(row, "clone_fraction", 1.0),
                    essentiality_weight=f(row, "essentiality_weight", 1.0),
                    threshold_abundance=optf(row, "threshold_abundance"),
                    baseline_abundance=optf(row, "baseline_abundance"),
                    burden_weight=f(row, "burden_weight", 0.02),
                    toxicity_weight=f(row, "toxicity_weight", 0.0),
                    regulation_strength=f(row, "regulation_strength", 0.0),
                    burstiness=f(row, "burstiness", 1.0),
                    stress_sensitivity=f(row, "stress_sensitivity", 0.2),
                    oncogenic_boost=f(row, "oncogenic_boost", 0.0),
                    pathway=s(row, "pathway", "generic"),
                )
            )
        return cls(genes=genes, microenvironment=microenvironment, **kwargs)

    @staticmethod
    def demo_dataset() -> pd.DataFrame:
        """
        A comprehensive synthetic dataset designed to stress-test the bounds,
        edge cases, and trade-offs of the Cancer RLTO model.
        """
        return pd.DataFrame([
            # -------------------------------------------------------------------------
            # 1. THE ANCHORS (Typical behavior, smooth curves)
            # -------------------------------------------------------------------------
            {"gene": "MYC", "transcription_rate": 8.0, "essentiality_weight": 0.8, "burden_weight": 0.015,
             "toxicity_weight": 0.020, "pathway": "proliferation"},
            {"gene": "KRAS", "transcription_rate": 5.5, "essentiality_weight": 0.9, "burden_weight": 0.012,
             "toxicity_weight": 0.015, "pathway": "MAPK"},

            # -------------------------------------------------------------------------
            # 2. THE EXTREME EDGE CASES (Testing zeros and boundaries)
            # -------------------------------------------------------------------------
            # The "Useless Passenger": Zero essentiality. Net fitness should be negative.
            {"gene": "PASSENGER", "transcription_rate": 2.0, "essentiality_weight": 0.0, "burden_weight": 0.01,
             "toxicity_weight": 0.01, "pathway": "none"},

            # The "Free Lunch": Zero toxicity and burden. Tests if the optimizer pushes to the absolute upper bound (20x baseline).
            {"gene": "MAGIC_BULLET", "transcription_rate": 3.0, "essentiality_weight": 1.0, "burden_weight": 0.0,
             "toxicity_weight": 0.0, "pathway": "synthetic"},

            # The "Micro-Subclone": Almost nonexistent clone fraction. Tests numerical stability at near-zero baseline abundances.
            {"gene": "RARE_VAR", "transcription_rate": 5.0, "clone_fraction": 0.01, "essentiality_weight": 0.8,
             "burden_weight": 0.01, "toxicity_weight": 0.01, "pathway": "subclonal"},

            # -------------------------------------------------------------------------
            # 3. NOISE AND REGULATION EXTREMES (Testing the Gamma distribution shape)
            # -------------------------------------------------------------------------
            # The "Chaos Engine": Massive burstiness, low regulation. Wide, flat Gamma distribution. Needs high overabundance to survive.
            {"gene": "CHAOS_GENE", "transcription_rate": 2.0, "burstiness": 15.0, "regulation_strength": 0.0,
             "essentiality_weight": 0.9, "burden_weight": 0.01, "toxicity_weight": 0.01, "pathway": "stress_response"},

            # The "Perfect Housekeeper": Tiny burstiness, extreme regulation. Tight spike Gamma distribution. Should hover exactly on its threshold.
            {"gene": "HOUSEKEEPER", "transcription_rate": 10.0, "burstiness": 0.1, "regulation_strength": 5.0,
             "essentiality_weight": 1.0, "burden_weight": 0.02, "toxicity_weight": 0.005, "pathway": "metabolism"},

            # -------------------------------------------------------------------------
            # 4. TOXICITY AND BURDEN EXTREMES (Testing the penalty curves)
            # -------------------------------------------------------------------------
            # The "Lethal Driver": Cranking transcription beyond what the cell can handle, coupled with massive toxicity. Optimizer should squash it down.
            {"gene": "TOXIC_DRIVER", "transcription_rate": 25.0, "oncogenic_boost": 2.0, "essentiality_weight": 0.5,
             "burden_weight": 0.05, "toxicity_weight": 0.15, "pathway": "apoptosis"},

            # The "Heavy Protein": No superlinear toxicity, but massive linear burden (e.g., massive structural protein).
            {"gene": "STRUCTURAL", "transcription_rate": 15.0, "essentiality_weight": 0.9, "burden_weight": 0.1,
             "toxicity_weight": 0.0, "pathway": "cytoskeleton"},

            # -------------------------------------------------------------------------
            # 5. MICROENVIRONMENT SENSITIVITY
            # -------------------------------------------------------------------------
            # The "Fragile Gene": Hyper-sensitive to stress. In a high-stress environment, its baseline will tank.
            {"gene": "FRAGILE", "transcription_rate": 6.0, "stress_sensitivity": 0.9, "essentiality_weight": 0.85,
             "burden_weight": 0.01, "toxicity_weight": 0.01, "pathway": "dna_repair"}
        ])

    def abundance_fitness_curve(self, gene: GeneConfig, abundance_grid: Optional[np.ndarray] = None,
                                n_samples: int = 5000) -> pd.DataFrame:
        baseline = self.expected_abundance(gene)
        if abundance_grid is None:
            abundance_grid = np.geomspace(max(baseline * 0.05, 1e-6), baseline * 20.0, 100)

        rows = []
        for x in abundance_grid:
            res = self.net_gene_fitness(gene, mean_abundance=float(x), n_samples=n_samples)
            rows.append({
                "gene": gene.gene,
                "abundance": float(x),
                "threshold": res.threshold_abundance,
                "overabundance": res.overabundance,
                "robustness_probability": res.robustness_probability,
                "burden_cost": res.burden_cost,
                "toxicity_cost": res.toxicity_cost,
                "net_fitness": res.net_fitness,
            })
        return pd.DataFrame(rows)


class ModelDiagnostics:
    """
    Diagnostic suite to verify the well-posedness and numerical stability
    of the CancerRLTOModel.
    """

    def __init__(self, model):
        self.model = model

    def _get_gene(self, gene_name: str):
        for g in self.model.genes:
            if g.gene == gene_name:
                return g
        raise ValueError(f"Gene {gene_name} not found in model.")

    def plot_component_decomposition(self, gene_name: str, n_points: int = 100, n_samples: int = 5000):
        """1. Proves well-posedness by showing the internal trade-off components."""
        gene = self._get_gene(gene_name)
        baseline = self.model.expected_abundance(gene)

        # Wide grid to capture boundaries
        x_grid = np.geomspace(max(baseline * 0.01, 1e-6), baseline * 10.0, n_points)

        net_fits, robust_terms, burdens, toxicities = [], [], [], []

        for x in x_grid:
            res = self.model.net_gene_fitness(gene, mean_abundance=float(x), n_samples=n_samples)

            net_fits.append(res.net_fitness)
            burdens.append(res.burden_cost)
            toxicities.append(res.toxicity_cost)
            # Recover the robustness term: Net = Robustness - Burden - Toxicity
            robust_terms.append(res.net_fitness + res.burden_cost + res.toxicity_cost)

        plt.figure(figsize=(8, 5))
        plt.plot(x_grid, robust_terms, label="Robustness Benefit (Plateaus)", color="tab:green", lw=2)
        plt.plot(x_grid, burdens, label="Burden Cost (Linear)", color="tab:orange", linestyle="--")
        plt.plot(x_grid, toxicities, label="Toxicity Cost (Superlinear)", color="tab:red", linestyle="-.")
        plt.plot(x_grid, net_fits, label="Net Fitness (Objective)", color="black", lw=3)

        plt.xscale("log")
        plt.xlabel("Mean Abundance (log scale)")
        plt.ylabel("Fitness Contribution")
        plt.title(f"Component Decomposition: {gene_name}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        return plt.gcf()

    def plot_noise_stability(self, gene_name: str, n_evals: int = 10, grid_points: int = 40, n_samples: int = 5000):
        """2. Checks if Monte Carlo noise overwhelms the peak curvature."""
        gene = self._get_gene(gene_name)

        # Find the approximate peak first
        opt = self.model.optimize_gene_abundance(gene, n_samples=n_samples)
        best_x = opt["optimal_abundance"]

        # Create a tight linear grid +/- 15% around the optimum
        x_grid = np.linspace(best_x * 0.85, best_x * 1.15, grid_points)

        all_x, all_y = [], []
        mean_y = []

        for x in x_grid:
            y_vals = []
            for _ in range(n_evals):
                res = self.model.net_gene_fitness(gene, mean_abundance=float(x), n_samples=n_samples)
                y_vals.append(res.net_fitness)
                all_x.append(x)
                all_y.append(res.net_fitness)
            mean_y.append(np.mean(y_vals))

        plt.figure(figsize=(8, 5))
        plt.scatter(all_x, all_y, alpha=0.4, color="tab:blue", s=10, label="Stochastic Evaluations")
        plt.plot(x_grid, mean_y, color="black", lw=2, label="Mean Landscape")
        plt.axvline(best_x, color="tab:red", linestyle="--", label="Optimizer Target")

        plt.xlabel("Mean Abundance (Linear near optimum)")
        plt.ylabel("Net Fitness")
        plt.title(f"Zoomed Peak Noise Analysis: {gene_name}\n(Band thickness vs Peak curvature)")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        return plt.gcf()

    def plot_convergence_check(self, gene_name: str, n_trials: int = 20, n_samples: int = 3000):
        """3. Verifies the optimizer finds the same peak regardless of starting guess."""
        gene = self._get_gene(gene_name)
        baseline = self.model.expected_abundance(gene)

        # Define a wide range of start points
        start_points = np.geomspace(baseline * 0.1, baseline * 5.0, n_trials)
        results = []

        def objective(x):
            val = x[0] if isinstance(x, np.ndarray) else x

            # FIX: Prevent Nelder-Mead from exploring negative abundances
            if val <= 1e-6:
                return 1e9  # Massive penalty to force the optimizer back to positive numbers

            return -self.model.net_gene_fitness(gene, mean_abundance=float(val), n_samples=n_samples).net_fitness

        for x0 in start_points:
            # Using Nelder-Mead here to strictly test start-point reliance
            res = minimize(objective, x0=[x0], method="Nelder-Mead")
            results.append(res.x[0])

        plt.figure(figsize=(8, 4))
        plt.scatter(start_points, results, color="tab:purple", alpha=0.7)
        plt.axhline(np.median(results), color="black", linestyle="--", label="Median Optimum")

        plt.xscale("log")
        plt.xlabel("Optimizer Start Point (x0)")
        plt.ylabel("Final Optimal Abundance")
        plt.title(f"Start-Point Independence Check: {gene_name}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        return plt.gcf()

    def plot_bivariate_landscape(self, gene1_name: str, gene2_name: str, grid_size: int = 40, n_samples: int = 0):
        """4. Checks for infinite ridges or compensation in the 2D global landscape."""
        g1 = self._get_gene(gene1_name)
        g2 = self._get_gene(gene2_name)

        # FIX: Find the actual peaks and center the grid around them!
        opt1 = self.model.optimize_gene_abundance(g1)["optimal_abundance"]
        opt2 = self.model.optimize_gene_abundance(g2)["optimal_abundance"]

        # Keep the background genes fixed at baseline
        bg_abundance = sum(
            self.model.expected_abundance(g)
            for g in self.model.genes
            if g.gene not in (gene1_name, gene2_name)
        )
        baseline_total = bg_abundance + self.model.expected_abundance(g1) + self.model.expected_abundance(g2)

        # Create a tight grid right around the peaks
        x1_grid = np.linspace(opt1 * 0.2, opt1 * 2.5, grid_size)
        x2_grid = np.linspace(opt2 * 0.2, opt2 * 2.5, grid_size)

        X1, X2 = np.meshgrid(x1_grid, x2_grid)
        Z = np.zeros_like(X1)

        for i in range(grid_size):
            for j in range(grid_size):
                f1 = self.model.net_gene_fitness(g1, mean_abundance=float(X1[i, j])).net_fitness
                f2 = self.model.net_gene_fitness(g2, mean_abundance=float(X2[i, j])).net_fitness

                mean_fit = (f1 + f2) / 2.0

                # Apply the global carrying capacity penalty
                total_abundance = X1[i, j] + X2[i, j] + bg_abundance
                normalized_total = total_abundance / max(baseline_total, 1.0)
                resource_penalty = 0.5 * (normalized_total ** 2)

                Z[i, j] = mean_fit - resource_penalty

        plt.figure(figsize=(7, 6))
        cp = plt.contourf(X1, X2, Z, levels=40, cmap="viridis")
        plt.colorbar(cp, label="Global Tumor Fitness (with Resource Limit)")

        plt.contour(X1, X2, Z, levels=20, colors='black', alpha=0.5, linewidths=0.5)

        plt.xlabel(f"{gene1_name} Abundance")
        plt.ylabel(f"{gene2_name} Abundance")
        plt.title(f"Bivariate Landscape: {gene1_name} vs {gene2_name}")
        plt.tight_layout()
        return plt.gcf()

    def plot_sensitivity_profile(self, gene_name: str, param_name: str = "toxicity_global", n_steps: int = 15):
        """5. Ensures smooth responses to parameter perturbations (no brittle cliffs)."""
        gene = self._get_gene(gene_name)

        # Store original value
        orig_val = getattr(self.model, param_name)
        multipliers = np.linspace(0.5, 2.0, n_steps)
        opt_abundances = []

        for mult in multipliers:
            setattr(self.model, param_name, orig_val * mult)
            opt = self.model.optimize_gene_abundance(gene, n_samples=3000)
            opt_abundances.append(opt["optimal_abundance"])

        # Restore original value
        setattr(self.model, param_name, orig_val)

        plt.figure(figsize=(7, 4))
        plt.plot(multipliers, opt_abundances, marker='o', linestyle='-', color="tab:red")

        plt.xlabel(f"{param_name} Multiplier")
        plt.ylabel("Optimal Abundance")
        plt.title(f"Sensitivity Profile: {gene_name} response to {param_name}")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        return plt.gcf()


def example_usage() -> Dict[str, object]:
    # Set up a moderately stressful environment
    env = Microenvironment(hypoxia=0.6, nutrient_limitation=0.4)
    df = CancerRLTOModel.demo_dataset()
    model = CancerRLTOModel.from_dataframe(df, microenvironment=env)

    # =========================================================
    # 1. TEST ALL EDGE CASES (Gene-Level Optimization)
    # =========================================================
    print("\n=== 1. OPTIMIZING ALL GENES (Checking Edge Cases) ===")
    opt_df = model.optimize_all()

    # Print a clean summary table so you can see the baseline vs optimal shift for every case
    summary_cols = ["gene", "baseline_abundance", "optimal_abundance", "optimal_net_fitness"]
    print(opt_df[summary_cols].to_string(index=False))

    # =========================================================
    # 2. GLOBAL & INTERVENTION TESTING
    # =========================================================
    print("\n=== 2. GLOBAL TUMOR STATE OPTIMIZATION ===")
    opt_vector, best_global_fitness = model.optimize_global()
    print(f"Max achievable global fitness (with resource limits): {best_global_fitness:.4f}")

    print("\n=== 3. INTERVENTION TESTING ===")
    # Let's see what happens if we target the TOXIC_DRIVER
    opt_inhibition, residual_fitness = model.optimize_inhibition("TOXIC_DRIVER")
    print(f"Optimal inhibition for TOXIC_DRIVER: {opt_inhibition:.4f}")

    # =========================================================
    # 4. EXHAUSTIVE MATHEMATICAL DIAGNOSTICS & PLOTTING
    # =========================================================
    print("\n=== 4. RUNNING EXHAUSTIVE DIAGNOSTICS ===")
    import os
    os.makedirs("diagnostic_plots", exist_ok=True)

    diag = ModelDiagnostics(model)

    # 1. Generate all single-gene plots
    for g in model.genes:
        name = g.gene
        print(f"Generating plots for {name}...")

        diag.plot_component_decomposition(name).savefig(f"diagnostic_plots/{name}_component.png")
        plt.close()

        diag.plot_noise_stability(name).savefig(f"diagnostic_plots/{name}_noise.png")
        plt.close()

        diag.plot_convergence_check(name).savefig(f"diagnostic_plots/{name}_convergence.png")
        plt.close()

        diag.plot_sensitivity_profile(name, param_name="toxicity_global").savefig(
            f"diagnostic_plots/{name}_sensitivity.png")
        plt.close()

    # 2. Generate specific bivariate landscapes
    print("Generating Bivariate Landscapes...")
    diag.plot_bivariate_landscape("MYC", "KRAS").savefig("diagnostic_plots/bivariate_MYC_KRAS.png")
    plt.close()

    diag.plot_bivariate_landscape("CHAOS_GENE", "HOUSEKEEPER").savefig(
        "diagnostic_plots/bivariate_CHAOS_HOUSEKEEPER.png")
    plt.close()

    diag.plot_bivariate_landscape("PASSENGER", "MAGIC_BULLET").savefig("diagnostic_plots/bivariate_PASSENGER_MAGIC.png")
    plt.close()

    diag.plot_bivariate_landscape("MAGIC_BULLET", "RARE_VAR").savefig("diagnostic_plots/bivariate_MAGIC_RARE.png")
    plt.close()

    diag.plot_bivariate_landscape("TOXIC_DRIVER", "STRUCTURAL").savefig(
        "diagnostic_plots/bivariate_TOXIC_STRUCTURAL.png")
    plt.close()

    print("All plots saved to the 'diagnostic_plots' directory!")

    # =========================================================
    # Extract the core columns from our optimization dataframe and convert to a list of dictionaries
    all_gene_results = opt_df[
        ["gene", "baseline_abundance", "optimal_abundance", "optimal_net_fitness"]
    ].to_dict(orient="records")

    return {
        "Best_global_fitness": best_global_fitness,
        "Optimal_Inhibition_Target": opt_inhibition,
        "All_Genes_Optimization": all_gene_results
    }


if __name__ == "__main__":
    summary = example_usage()
    print("\n--- Final Output ---")
    print(json.dumps(summary, indent=2))
