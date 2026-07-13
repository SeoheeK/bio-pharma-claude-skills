"""Test suite for the pkchat engine (pure stdlib, run with `python -m unittest`)."""
import math
import os
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pkchat import (
    PKParameters, DosingRegimen, analytical_concentration, simulate,
    simulate_population, nca, fit_individual, fit_population,
    allometry_ensemble, SpeciesPK, recommend_cases,
)
from pkchat.chat import dispatch, tool_schemas
from pkchat.chat.agent import LocalRouter


class TestModels(unittest.TestCase):
    def test_1cmt_iv_ode_matches_analytical(self):
        p = PKParameters(CL=5.0, V1=30.0)
        reg = DosingRegimen.single(100, "iv_bolus")
        t = [0, 0.5, 1, 2, 4, 8, 12, 24]
        ana = analytical_concentration(p, t, 100, "iv_bolus")
        ode = simulate(p, reg, times=t, dt=0.01).concentrations
        for a, o in zip(ana, ode):
            self.assertAlmostEqual(a, o, places=4)

    def test_2cmt_iv_ode_matches_analytical(self):
        p = PKParameters(CL=4.0, V1=8.0, Q=3.0, V2=10.0)
        reg = DosingRegimen.single(100, "iv_bolus")
        t = [0.25, 0.5, 1, 2, 4, 8, 12, 24]
        ana = analytical_concentration(p, t, 100, "iv_bolus")
        ode = simulate(p, reg, times=t, dt=0.005).concentrations
        for a, o in zip(ana, ode):
            self.assertAlmostEqual(a, o, delta=1e-2)

    def test_oral_ode_matches_analytical(self):
        p = PKParameters(CL=5.0, V1=30.0, ka=1.2, F=0.8)
        reg = DosingRegimen.single(200, "oral")
        t = [0.25, 0.5, 1, 2, 4, 8, 12, 24]
        ana = analytical_concentration(p, t, 200, "oral")
        ode = simulate(p, reg, times=t, dt=0.005).concentrations
        for a, o in zip(ana, ode):
            self.assertAlmostEqual(a, o, delta=1e-2)

    def test_infusion_steady_state(self):
        # 1-cmt infusion: plateau approaches R0/CL.
        p = PKParameters(CL=5.0, V1=30.0)
        reg = DosingRegimen.single(500, "infusion", tinf=100.0)
        res = simulate(p, reg, times=[99.9], t_end=100, dt=0.01)
        plateau = (500 / 100) / 5.0
        self.assertAlmostEqual(res.concentrations[0], plateau, delta=0.05)

    def test_multiple_dose_accumulation(self):
        p = PKParameters(CL=5.0, V1=30.0)
        single = DosingRegimen.single(100, "iv_bolus")
        multi = DosingRegimen.multiple(100, interval=4, n_doses=5, route="iv_bolus")
        c_single = simulate(p, single, times=[16.0], t_end=20).concentrations[0]
        c_multi = simulate(p, multi, times=[16.0], t_end=20).concentrations[0]
        self.assertGreater(c_multi, c_single)  # accumulation

    def test_nonlinear_michaelis_menten_runs(self):
        p = PKParameters(CL=1.0, V1=30.0, Vmax=50.0, Km=2.0)
        reg = DosingRegimen.single(200, "iv_bolus")
        res = simulate(p, reg, t_end=24, dt=0.02)
        self.assertGreater(res.cmax(), 0)
        self.assertTrue(all(c >= 0 for c in res.concentrations))


class TestNCA(unittest.TestCase):
    def test_nca_recovers_cl_and_thalf(self):
        p = PKParameters(CL=5.0, V1=30.0)
        t = [0.5, 1, 2, 4, 8, 12, 24]
        c = analytical_concentration(p, t, 100, "iv_bolus")
        r = nca(t, c, dose=100, route="iv_bolus")
        self.assertAlmostEqual(r.clearance, 5.0, delta=0.6)
        self.assertAlmostEqual(r.t_half, math.log(2) / (5 / 30), delta=0.2)


class TestEstimation(unittest.TestCase):
    def test_fit_recovers_truth(self):
        truth = PKParameters(CL=5.0, V1=30.0)
        t = [0.5, 1, 2, 4, 8, 12, 24]
        obs = analytical_concentration(truth, t, 100, "iv_bolus")
        fit = fit_individual(t, obs, DosingRegimen.single(100, "iv_bolus"),
                             estimate=("CL", "V1"))
        self.assertAlmostEqual(fit.estimates["CL"], 5.0, delta=0.1)
        self.assertAlmostEqual(fit.estimates["V1"], 30.0, delta=1.0)

    def test_two_stage_population(self):
        subjects = []
        for i, cl in enumerate([4.0, 5.0, 6.0, 5.5, 4.5]):
            p = PKParameters(CL=cl, V1=30.0)
            t = [0.5, 1, 2, 4, 8, 12, 24]
            c = analytical_concentration(p, t, 100, "iv_bolus")
            subjects.append({"times": t, "concentrations": c,
                             "regimen": DosingRegimen.single(100, "iv_bolus")})
        pop = fit_population(subjects, estimate=("CL", "V1"))
        self.assertAlmostEqual(pop.typical["CL"], 4.9, delta=0.5)
        self.assertGreater(pop.iiv_cv_pct["CL"], 0)


class TestAllometry(unittest.TestCase):
    def test_ensemble_reasonable(self):
        sp = [SpeciesPK("mouse", 0.02, 0.8), SpeciesPK("rat", 0.25, 3.5),
              SpeciesPK("dog", 10, 45), SpeciesPK("monkey", 5, 28)]
        res = allometry_ensemble(sp)
        self.assertGreater(res.ensemble, 0)
        self.assertTrue(0.4 < res.exponent < 1.2)
        self.assertGreater(res.r2, 0.8)

    def test_needs_two_species(self):
        with self.assertRaises(ValueError):
            allometry_ensemble([SpeciesPK("rat", 0.25, 3.5)])


class TestCaseDB(unittest.TestCase):
    def test_recommend_ranks_matching_class(self):
        rec = recommend_cases(drug_class="immune_checkpoint", modality="mAb", route="iv")
        self.assertEqual(rec.matches[0]["score"], 90)
        self.assertGreater(rec.suggested_CL, 0)


class TestTools(unittest.TestCase):
    def test_all_tools_have_schemas(self):
        schemas = tool_schemas()
        self.assertEqual(len(schemas), 6)
        for s in schemas:
            self.assertIn("name", s)
            self.assertIn("input_schema", s)

    def test_dispatch_simulate(self):
        out = dispatch("simulate_pk", {"CL": 5, "V1": 30, "dose": 100, "route": "iv_bolus"})
        self.assertIn("Cmax", out)

    def test_dispatch_unknown(self):
        self.assertIn("Unknown tool", dispatch("nope", {}))

    def test_dispatch_error_is_caught(self):
        # Missing required arg surfaces as readable text, not an exception.
        out = dispatch("simulate_pk", {})
        self.assertIn("failed", out.lower())


class TestLocalRouter(unittest.TestCase):
    def setUp(self):
        self.r = LocalRouter()

    def test_simulate_intent(self):
        out = self.r.handle("simulate CL=5 V1=30 dose=100 iv")
        self.assertIn("Cmax", out)

    def test_nca_intent(self):
        out = self.r.handle("run NCA times=[0.5,1,2,4,8] conc=[8,6,4,2,0.5] dose=100 iv")
        self.assertIn("AUC", out)

    def test_allometry_intent(self):
        out = self.r.handle("scale to human: rat 0.25 3.5, dog 10 45, monkey 5 28")
        self.assertIn("ENSEMBLE", out)

    def test_recommend_intent(self):
        out = self.r.handle("recommend a PopPK model for a PD-1 mAb given IV")
        self.assertIn("Nivolumab", out)

    def test_population_intent(self):
        out = self.r.handle("population CL=5 V1=30 dose=100 iiv=30% n=50")
        self.assertIn("prediction interval", out)

    def test_help_fallback(self):
        out = self.r.handle("hello there")
        self.assertIn("I can", out)

    def test_oral_population_peaks_after_zero(self):
        # Regression: oral ka must be forwarded, so the peak is not at t=0.
        out = self.r.handle("population CL=6 V1=40 dose=200 oral ka=1.0 iiv=45% n=50")
        m = re.search(r"at t = ([0-9.]+) h", out)
        self.assertIsNotNone(m)
        self.assertGreater(float(m.group(1)), 0.5)


class TestRouteGuards(unittest.TestCase):
    def test_oral_without_ka_rejected(self):
        out = dispatch("simulate_pk", {"CL": 5, "V1": 30, "dose": 100, "route": "oral"})
        self.assertIn("ka", out)
        self.assertIn("failed", out.lower())

    def test_infusion_without_tinf_rejected(self):
        out = dispatch("simulate_pk", {"CL": 5, "V1": 30, "dose": 100, "route": "infusion"})
        self.assertIn("tinf", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
