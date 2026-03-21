from __future__ import annotations

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

from dataclasses import dataclass, asdict, fields
from typing import Dict, Optional, Sequence, Tuple
import json
import math

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

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

    The effective abundance distribution is approximated with a gamma model,
    which is flexible and positive-valued. A gene's tumor fitness contribution
    is high when abundance remains above a critical threshold despite noise,
    but excessive abundance produces burden and toxicity costs.
    """

    def __init__(
            self,
            genes: Sequence[GeneConfig],
            microenvironment: Optional[Microenvironment] = None,
            hill_n: float = 12.0,
            base_proliferation: float = 1.0,
            burden_global: float = 1.0,
            toxicity_global: float = 1.0,
            robustness_scale: float = 5.0,
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
        """Critical threshold below which tumor-cell fitness drops sharply."""
        baseline = self.expected_abundance(gene) if baseline is None else float(baseline)
        if gene.threshold_abundance is not None:
            return float(gene.threshold_abundance)

        # Cancer reframing:
        # Higher stress and stronger essentiality raise the effective threshold.
        # Regulation partially lowers threshold through tighter control.
        stress_multiplier = 1.0 + 0.8 * gene.stress_sensitivity * self.microenvironment.stress_index
        essentiality_multiplier = 0.4 + 0.6 * gene.essentiality_weight
        regulation_discount = 1.0 / (1.0 + 0.25 * gene.regulation_strength)

        # Low-expression genes retain a lower absolute threshold but may still
        # exhibit large overabundance because baseline stays far above threshold.
        threshold = 2.5 * np.sqrt(baseline) * stress_multiplier * essentiality_multiplier * regulation_discount
        return float(max(threshold, 1e-9))

    def abundance_distribution_params(self, gene: GeneConfig, mean_abundance: Optional[float] = None) -> Tuple[
        float, float]:
        """Gamma(shape, scale) parameters for effective abundance noise."""
        mu = self.expected_abundance(gene) if mean_abundance is None else float(mean_abundance)

        # Noise decreases with higher expression and tighter regulation.
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

    def robustness_probability(
            self,
            gene: GeneConfig,
            mean_abundance: Optional[float] = None,
            n_samples: int = 20000,
    ) -> float:
        samples = self.sample_abundance(gene, n=n_samples, mean_abundance=mean_abundance)
        threshold = self.threshold_abundance(gene, baseline=mean_abundance)
        return float(np.mean(samples >= threshold))

    def threshold_fitness(
            self,
            abundance: ArrayLike,
            threshold: float,
            kmax: float = 1.0,
    ) -> np.ndarray:
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

    def net_gene_fitness(
            self,
            gene: GeneConfig,
            mean_abundance: Optional[float] = None,
            n_samples: int = 15000,
    ) -> ModelResult:
        baseline = self.expected_abundance(gene) if mean_abundance is None else float(mean_abundance)
        threshold = self.threshold_abundance(gene, baseline=baseline)
        samples = self.sample_abundance(gene, n=n_samples, mean_abundance=baseline)
        growth = self.threshold_fitness(samples, threshold, kmax=self.base_proliferation)
        robustness_prob = float(np.mean(samples >= threshold))
        robustness_term = self.robustness_scale * float(np.mean(growth)) * gene.essentiality_weight
        burden = self.burden_cost(gene, baseline)
        toxicity = self.toxicity_cost(gene, baseline)
        net = robustness_term - burden - toxicity

        return ModelResult(
            gene=gene.gene,
            baseline_abundance=baseline,
            threshold_abundance=threshold,
            overabundance=float(baseline / threshold),
            abundance_mean=float(np.mean(samples)),
            abundance_std=float(np.std(samples, ddof=1)),
            robustness_probability=robustness_prob,
            burden_cost=burden,
            toxicity_cost=toxicity,
            net_fitness=float(net),
            pathway=gene.pathway,
        )

    def evaluate(self, n_samples: int = 15000) -> pd.DataFrame:
        rows = [asdict(self.net_gene_fitness(g, n_samples=n_samples)) for g in self.genes]
        df = pd.DataFrame(rows)

        bad = ~np.isfinite(df["net_fitness"])
        if bad.any():
            print("Non-finite net_fitness rows:")
            print(df.loc[bad, ["gene", "baseline_abundance", "threshold_abundance",
                               "overabundance", "abundance_mean", "abundance_std",
                               "robustness_probability", "burden_cost",
                               "toxicity_cost", "net_fitness"]])
            raise ValueError("Non-finite values detected in net_fitness")

        df["rank"] = df["net_fitness"].rank(ascending=False, method="dense").astype("Int64")
        return df.sort_values(["net_fitness", "overabundance"], ascending=[False, False]).reset_index(drop=True)

    def optimize_gene_abundance(
            self,
            gene: GeneConfig,
            abundance_grid: Optional[np.ndarray] = None,
            n_samples: int = 10000,
    ) -> Dict[str, float]:
        baseline = self.expected_abundance(gene)
        if abundance_grid is None:
            abundance_grid = np.geomspace(max(baseline * 0.05, 1e-6), baseline * 20.0, 200)

        fitness = []
        for x in abundance_grid:
            res = self.net_gene_fitness(gene, mean_abundance=float(x), n_samples=n_samples)
            fitness.append(res.net_fitness)
        idx = int(np.argmax(fitness))
        best = float(abundance_grid[idx])
        thr = self.threshold_abundance(gene, baseline=best)
        return {
            "gene": gene.gene,
            "baseline_abundance": baseline,
            "optimal_abundance": best,
            "critical_threshold": thr,
            "optimal_overabundance": float(best / thr),
            "optimal_net_fitness": float(fitness[idx]),
        }

    def optimize_all(self, n_samples: int = 10000) -> pd.DataFrame:
        rows = [self.optimize_gene_abundance(g, n_samples=n_samples) for g in self.genes]
        return pd.DataFrame(rows).sort_values("optimal_net_fitness", ascending=False).reset_index(drop=True)

    def pathway_summary(self, n_samples: int = 15000) -> pd.DataFrame:
        df = self.evaluate(n_samples=n_samples)
        out = (
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
        return out

    def simulate_intervention(
            self,
            target_gene: str,
            inhibition: float,
            n_samples: int = 15000,
    ) -> pd.DataFrame:
        inhibition = float(np.clip(inhibition, 0.0, 1.0))
        rows = []
        for gene in self.genes:
            baseline = self.expected_abundance(gene)
            if gene.gene == target_gene:
                baseline *= (1.0 - inhibition)
            rows.append(asdict(self.net_gene_fitness(gene, mean_abundance=baseline, n_samples=n_samples)))
        return pd.DataFrame(rows).sort_values("net_fitness", ascending=False).reset_index(drop=True)

    def tumor_fitness_score(self, n_samples: int = 12000) -> float:
        df = self.evaluate(n_samples=n_samples)
        # Mean positive fitness contribution across genes, clipped to avoid negatives dominating numerically.
        return float(np.mean(np.clip(df["net_fitness"].values, -10.0, 10.0)))

    @classmethod
    def from_dataframe(
            cls,
            df: pd.DataFrame,
            microenvironment: Optional[Microenvironment] = None,
            **kwargs,
    ) -> "CancerRLTOModel":
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
        return pd.DataFrame(
            [
                {
                    "gene": "MYC",
                    "transcription_rate": 8.0,
                    "translation_efficiency": 1.8,
                    "copy_number": 6.0,
                    "essentiality_weight": 0.8,
                    "burden_weight": 0.015,
                    "toxicity_weight": 0.020,
                    "regulation_strength": 0.3,
                    "burstiness": 2.2,
                    "stress_sensitivity": 0.3,
                    "oncogenic_boost": 0.7,
                    "pathway": "proliferation",
                },
                {
                    "gene": "KRAS",
                    "transcription_rate": 5.5,
                    "translation_efficiency": 1.4,
                    "copy_number": 3.0,
                    "essentiality_weight": 0.9,
                    "burden_weight": 0.012,
                    "toxicity_weight": 0.015,
                    "regulation_strength": 0.2,
                    "burstiness": 1.6,
                    "stress_sensitivity": 0.25,
                    "oncogenic_boost": 0.5,
                    "pathway": "MAPK",
                },
                {
                    "gene": "RRM2",
                    "transcription_rate": 3.2,
                    "translation_efficiency": 1.5,
                    "copy_number": 2.5,
                    "essentiality_weight": 1.0,
                    "burden_weight": 0.018,
                    "toxicity_weight": 0.004,
                    "regulation_strength": 0.1,
                    "burstiness": 1.8,
                    "stress_sensitivity": 0.45,
                    "pathway": "replication",
                },
                {
                    "gene": "POLR2A",
                    "transcription_rate": 2.4,
                    "translation_efficiency": 1.1,
                    "copy_number": 2.0,
                    "essentiality_weight": 1.0,
                    "burden_weight": 0.025,
                    "toxicity_weight": 0.002,
                    "regulation_strength": 0.4,
                    "burstiness": 1.2,
                    "stress_sensitivity": 0.35,
                    "pathway": "transcription",
                },
                {
                    "gene": "MCL1",
                    "transcription_rate": 2.0,
                    "translation_efficiency": 1.2,
                    "copy_number": 4.0,
                    "essentiality_weight": 0.95,
                    "burden_weight": 0.010,
                    "toxicity_weight": 0.010,
                    "regulation_strength": 0.2,
                    "burstiness": 2.8,
                    "stress_sensitivity": 0.5,
                    "oncogenic_boost": 0.4,
                    "pathway": "apoptosis_evasion",
                },
                {
                    "gene": "DHODH",
                    "transcription_rate": 1.4,
                    "translation_efficiency": 1.0,
                    "copy_number": 2.0,
                    "essentiality_weight": 0.9,
                    "burden_weight": 0.015,
                    "toxicity_weight": 0.003,
                    "regulation_strength": 0.05,
                    "burstiness": 3.0,
                    "stress_sensitivity": 0.55,
                    "pathway": "metabolism",
                },
                {
                    "gene": "PSMA1",
                    "transcription_rate": 1.1,
                    "translation_efficiency": 1.0,
                    "copy_number": 2.0,
                    "essentiality_weight": 1.0,
                    "burden_weight": 0.011,
                    "toxicity_weight": 0.001,
                    "regulation_strength": 0.15,
                    "burstiness": 2.4,
                    "stress_sensitivity": 0.4,
                    "pathway": "proteostasis",
                },
            ]
        )

    def save_results(self, output_prefix: str, n_samples: int = 15000) -> None:
        eval_df = self.evaluate(n_samples=n_samples)
        opt_df = self.optimize_all(n_samples=max(n_samples // 2, 2000))
        path_df = self.pathway_summary(n_samples=n_samples)
        eval_df.to_csv(f"{output_prefix}_gene_fitness.csv", index=False)
        opt_df.to_csv(f"{output_prefix}_optimal_abundance.csv", index=False)
        path_df.to_csv(f"{output_prefix}_pathway_summary.csv", index=False)

    def plot_expression_vs_overabundance(self, n_samples: int = 12000):
        if plt is None:
            raise RuntimeError("matplotlib is not available")
        df = self.evaluate(n_samples=n_samples)
        src = []
        for g in self.genes:
            src.append({"gene": g.gene, "transcription_rate": g.transcription_rate})
        src_df = pd.DataFrame(src)
        merged = df.merge(src_df, on="gene", how="left")
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(merged["transcription_rate"], merged["overabundance"])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Transcription rate")
        ax.set_ylabel("Overabundance")
        ax.set_title("Cancer RLTO: expression vs overabundance")
        fig.tight_layout()
        return fig, ax

    def abundance_fitness_curve(
            self,
            gene: GeneConfig,
            abundance_grid: Optional[np.ndarray] = None,
            n_samples: int = 5000,
    ) -> pd.DataFrame:
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


def example_usage() -> Dict[str, object]:
    """Run a demonstration and return serializable summaries."""
    env = Microenvironment(
        hypoxia=0.6,
        nutrient_limitation=0.4,
        immune_pressure=0.2,
        drug_pressure=0.1,
        oxidative_stress=0.5,
    )
    df = CancerRLTOModel.demo_dataset()
    model = CancerRLTOModel.from_dataframe(df, microenvironment=env, hill_n=10.0)
    for g in model.genes:
        b = model.expected_abundance(g)
        print(g.gene, "expected_abundance =", b)
    evaluation = model.evaluate(n_samples=8000)
    optimal = model.optimize_all(n_samples=5000)
    score = model.tumor_fitness_score(n_samples=6000)

    plot_fig, plot_ax = model.plot_expression_vs_overabundance(n_samples=12000)
    plt.show()

    abundance_df = model.abundance_fitness_curve(model.genes[0], n_samples=10000)
    plt.figure()
    plt.plot(abundance_df["abundance"], abundance_df["net_fitness"])
    plt.xlabel("Abundance")
    plt.ylabel("Net fitness")
    plt.title("Abundance vs net fitness")
    plt.show()

    print(evaluation[[
        "gene", "baseline_abundance", "threshold_abundance", "overabundance",
        "robustness_probability", "burden_cost", "toxicity_cost", "net_fitness", "rank"
    ]].to_string(index=False))

    return {
        "tumor_fitness_score": score,
        "top_gene_fitness": evaluation.head(5).to_dict(orient="records"),
        "top_optimal_abundance": optimal.head(5).to_dict(orient="records"),
    }


if __name__ == "__main__":
    summary = example_usage()
    print(json.dumps(summary, indent=2))
