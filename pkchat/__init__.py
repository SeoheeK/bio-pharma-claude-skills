"""pkchat -- a chatbot-driven pharmacokinetic modelling engine.

Layers:
  models / regimen / simulate  -- structural PK models and simulation
  nca                          -- non-compartmental analysis
  estimation                   -- individual + two-stage population PK fitting
  ml.allometry                 -- animal-to-human allometric scaling ensemble
  casedb                       -- PopPK case-study database + recommender
  chat                         -- natural-language layer (Claude tool-use +
                                  dependency-free rule-based fallback)

The numeric core is pure standard library, so it runs with no external
packages; the chat layer uses the Anthropic SDK when available and an API key
is set, otherwise a local intent router drives the same tools.
"""
from .models import PKParameters, CompartmentModel, analytical_concentration
from .regimen import DosingRegimen, DoseEvent
from .simulate import simulate, simulate_population, SimulationResult, PopulationResult
from .nca import nca, NCAResult
from .estimation import fit_individual, fit_population, FitResult, PopPKResult
from .ml.allometry import allometry_ensemble, SpeciesPK, AllometryResult
from .casedb import recommend_cases, CASES, PKCase

__version__ = "0.1.0"

__all__ = [
    "PKParameters", "CompartmentModel", "analytical_concentration",
    "DosingRegimen", "DoseEvent",
    "simulate", "simulate_population", "SimulationResult", "PopulationResult",
    "nca", "NCAResult",
    "fit_individual", "fit_population", "FitResult", "PopPKResult",
    "allometry_ensemble", "SpeciesPK", "AllometryResult",
    "recommend_cases", "CASES", "PKCase",
]
