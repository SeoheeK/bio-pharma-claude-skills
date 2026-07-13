"""PopPK case-study database and similarity recommender.

A compact, structured subset of the ``pk-case-study-db`` skill: real published
population-PK models keyed by drug class, modality, and route. Given a new
drug's attributes, the recommender scores and ranks similar cases and proposes
initial parameter estimates from their median -- the starting point a modeller
would use to seed an estimation run.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Sequence


@dataclass
class PKCase:
    drug: str
    drug_class: str
    modality: str            # "mAb" | "small_molecule"
    route: str               # "iv" | "oral" | "sc"
    structure: str           # e.g. "2-compartment linear"
    CL: float                # L/h
    V1: float                # L
    Q: float | None = None
    V2: float | None = None
    iiv_cl_cv: float | None = None   # %
    key_covariates: list[str] = field(default_factory=list)
    reference: str = ""
    nct: str = ""


# A curated subset (units normalised to L/h and L).
CASES: list[PKCase] = [
    PKCase("Nivolumab", "immune_checkpoint", "mAb", "iv", "2-compartment linear",
           CL=0.0095, V1=3.5, Q=0.0046, V2=5.6, iiv_cl_cv=35,
           key_covariates=["body_weight", "ECOG", "albumin", "ADA"],
           reference="Bajaj G et al. J Clin Pharmacol 2017", nct="NCT01667809"),
    PKCase("Pembrolizumab", "immune_checkpoint", "mAb", "iv", "2-compartment time-varying CL",
           CL=0.209, V1=3.45, Q=0.555, V2=4.10, iiv_cl_cv=24,
           key_covariates=["body_weight", "tumor_burden", "ECOG"],
           reference="Ahamadi M et al. CPT:PSP 2017", nct="NCT01866319"),
    PKCase("Atezolizumab", "immune_checkpoint", "mAb", "iv", "2-compartment linear",
           CL=0.200 / 24, V1=3.28, Q=0.485 / 24, V2=3.40, iiv_cl_cv=33,
           key_covariates=["body_weight", "albumin", "ADA"],
           reference="FDA clinical pharmacology review", nct="NCT01375842"),
    PKCase("Bevacizumab", "anti_vegf", "mAb", "iv", "2-compartment linear",
           CL=0.231 / 24, V1=2.73, Q=0.462 / 24, V2=1.69, iiv_cl_cv=32,
           key_covariates=["body_weight", "sex", "ADA"],
           reference="Avastin label", nct=""),
    PKCase("Trastuzumab", "anti_her2", "mAb", "iv", "2-compartment (TMDD approx)",
           CL=0.225 / 24, V1=2.95, Q=0.614 / 24, V2=2.84, iiv_cl_cv=44,
           key_covariates=["sHER2", "ALT", "body_weight"],
           reference="Herceptin label", nct=""),
    PKCase("Adalimumab", "anti_tnf", "mAb", "sc", "1-compartment SC",
           CL=0.476 / 24, V1=8.17, iiv_cl_cv=36,
           key_covariates=["body_weight", "methotrexate", "ADA"],
           reference="Humira label", nct=""),
    PKCase("Imatinib", "tki", "small_molecule", "oral", "1-compartment first-order abs",
           CL=14.3, V1=347, iiv_cl_cv=40,
           key_covariates=["body_weight", "AGP", "age"],
           reference="Widmer N et al. Clin Pharmacokinet 2008", nct=""),
    PKCase("Ibrutinib", "btk_inhibitor", "small_molecule", "oral", "1-compartment Weibull abs",
           CL=62.0, V1=683, iiv_cl_cv=66,
           key_covariates=["food", "CYP3A4_inhibitor"],
           reference="Imbruvica PCYC-1102", nct="NCT01105716"),
    PKCase("Venetoclax", "bcl2_inhibitor", "small_molecule", "oral", "2-compartment first-order abs",
           CL=389.0, V1=256, Q=85.7, V2=897, iiv_cl_cv=52,
           key_covariates=["food", "CYP3A4_inhibitor"],
           reference="Venclexta M14-032", nct="NCT02141282"),
    PKCase("Busulfan", "alkylator", "small_molecule", "iv", "1-compartment",
           CL=7.82, V1=48.3, iiv_cl_cv=28,
           key_covariates=["body_weight", "GSTA1", "azole", "bilirubin"],
           reference="HSCT conditioning case", nct="NCT01052883"),
    PKCase("Cyclosporine", "calcineurin_inhibitor", "small_molecule", "oral", "2-compartment oral",
           CL=35.2, V1=989, Q=28.3, V2=1280, iiv_cl_cv=45,
           key_covariates=["CYP3A5", "MDR1", "voriconazole", "hematocrit"],
           reference="HSCT/transplant case", nct="NCT00823940"),
    PKCase("Tacrolimus", "calcineurin_inhibitor", "small_molecule", "oral", "1-compartment CYP3A5 mixture",
           CL=17.8, V1=590, iiv_cl_cv=60,
           key_covariates=["CYP3A5", "hematocrit", "azole"],
           reference="Transplant case", nct="NCT01094119"),
    PKCase("Vancomycin", "glycopeptide", "small_molecule", "iv", "2-compartment infusion",
           CL=3.68, V1=28.6, Q=2.11, V2=51.4, iiv_cl_cv=48,
           key_covariates=["creatinine_clearance", "body_weight", "sepsis", "CRRT"],
           reference="ICU PopPK", nct="NCT02535949"),
    PKCase("Voriconazole", "azole_antifungal", "small_molecule", "oral", "1-compartment nonlinear (MM)",
           CL=2.8, V1=4.6 * 70, iiv_cl_cv=65,
           key_covariates=["CYP2C19", "ABCB1"],
           reference="Vfend label", nct=""),
    PKCase("Fluconazole", "azole_antifungal", "small_molecule", "oral", "1-compartment first-order abs",
           CL=1.25, V1=41.7, iiv_cl_cv=28,
           key_covariates=["creatinine_clearance", "age"],
           reference="Diflucan label", nct=""),
]


@dataclass
class Recommendation:
    matches: list[dict]
    suggested_CL: float
    suggested_V1: float
    suggested_structure: str
    suggested_covariates: list[str]

    def summary(self) -> str:
        lines = ["Top matching PopPK cases:"]
        for m in self.matches:
            lines.append(f"  {m['drug']:14s} (score {m['score']}) "
                         f"CL={m['CL']:.4g} L/h, V1={m['V1']:.4g} L, {m['structure']}")
        lines.append(f"\nSuggested initial estimates: "
                     f"CL={self.suggested_CL:.4g} L/h, V1={self.suggested_V1:.4g} L")
        lines.append(f"Suggested structure: {self.suggested_structure}")
        lines.append(f"Consider covariates: {', '.join(self.suggested_covariates)}")
        return "\n".join(lines)


def recommend_cases(
    drug_class: str | None = None,
    modality: str | None = None,
    route: str | None = None,
    top_n: int = 5,
) -> Recommendation:
    """Score cases by attribute overlap and return the closest matches with
    median-based initial parameter suggestions."""
    scored = []
    for c in CASES:
        score = 0
        if drug_class and c.drug_class == drug_class:
            score += 40
        if modality and c.modality == modality:
            score += 30
        if route and c.route == route:
            score += 20
        scored.append((score, c))
    scored.sort(key=lambda x: x[0], reverse=True)
    top = [c for s, c in scored[:top_n]]

    def med(vals):
        vals = sorted(v for v in vals if v is not None)
        if not vals:
            return 0.0
        n = len(vals)
        return vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])

    covs: list[str] = []
    for c in top:
        for cov in c.key_covariates:
            if cov not in covs:
                covs.append(cov)

    structures = [c.structure for c in top]
    struct = max(set(structures), key=structures.count) if structures else "1-compartment"

    matches = [{"drug": c.drug, "score": scored[i][0] if i < len(scored) else 0,
                "CL": c.CL, "V1": c.V1, "structure": c.structure}
               for i, c in enumerate(top)]

    return Recommendation(
        matches=matches,
        suggested_CL=med([c.CL for c in top]),
        suggested_V1=med([c.V1 for c in top]),
        suggested_structure=struct,
        suggested_covariates=covs[:6],
    )
