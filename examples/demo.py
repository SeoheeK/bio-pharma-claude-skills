"""End-to-end walkthrough of the pkchat engine (Python API + chatbot).

Run:  python examples/demo.py
Requires nothing beyond the standard library.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pkchat import (
    PKParameters, DosingRegimen, simulate, simulate_population, nca,
    fit_individual, fit_population, allometry_ensemble, SpeciesPK,
    recommend_cases, analytical_concentration,
)
from pkchat.chat.agent import PKChatAgent


def section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main():
    section("1. Simulate a 2-compartment IV infusion, multiple doses")
    p = PKParameters(CL=4.0, V1=8.0, Q=3.0, V2=10.0)
    reg = DosingRegimen.multiple(500, interval=12, n_doses=4,
                                 route="infusion", tinf=1.0)
    res = simulate(p, reg, t_end=60)
    print(f"Cmax over regimen = {res.cmax():.3g} mg/L at t = {res.tmax():.3g} h")

    section("2. Non-compartmental analysis of a single IV-bolus profile")
    truth = PKParameters(CL=5.0, V1=30.0)
    t = [0.5, 1, 2, 4, 8, 12, 24]
    obs = analytical_concentration(truth, t, 100, "iv_bolus")
    r = nca(t, obs, dose=100, route="iv_bolus")
    print(f"AUC(0-inf) = {r.auc_inf:.3g}, CL = {r.clearance:.3g} L/h, "
          f"t1/2 = {r.t_half:.3g} h")

    section("3. Fit CL and V1 back from the observed data")
    fit = fit_individual(t, obs, DosingRegimen.single(100, "iv_bolus"),
                         estimate=("CL", "V1"))
    print(fit.summary())

    section("4. Two-stage population PK across 5 subjects")
    subjects = []
    for cl in [4.0, 5.0, 6.0, 5.5, 4.5]:
        c = analytical_concentration(PKParameters(CL=cl, V1=30.0), t, 100, "iv_bolus")
        subjects.append({"times": t, "concentrations": c,
                         "regimen": DosingRegimen.single(100, "iv_bolus")})
    pop = fit_population(subjects, estimate=("CL", "V1"))
    print(pop.summary())

    section("5. Allometric scaling: animal -> human clearance")
    sp = [SpeciesPK("mouse", 0.02, 0.8), SpeciesPK("rat", 0.25, 3.5),
          SpeciesPK("dog", 10, 45), SpeciesPK("monkey", 5, 28)]
    al = allometry_ensemble(sp)
    print(al.summary())

    section("6. Recommend a starting PopPK model for a new PD-1 mAb")
    rec = recommend_cases(drug_class="immune_checkpoint", modality="mAb", route="iv")
    print(rec.summary())

    section("7. Population simulation with 30% IIV on CL")
    popsim = simulate_population(truth, DosingRegimen.single(100, "iv_bolus"),
                                 n_subjects=200, iiv_cv={"CL": 0.3, "V1": 0.2},
                                 residual_cv=0.1)
    med = popsim.median()
    imax = max(range(len(med)), key=lambda i: med[i])
    print(f"median Cmax = {med[imax]:.3g} mg/L; "
          f"90% PI = [{popsim.percentile(5)[imax]:.3g}, "
          f"{popsim.percentile(95)[imax]:.3g}]")

    section("8. The chatbot layer (local backend)")
    agent = PKChatAgent(use_claude=False)
    for q in [
        "simulate CL=5 V1=30 dose=100 iv",
        "scale to human: mouse 0.02 0.8, rat 0.25 3.5, dog 10 45, monkey 5 28",
        "recommend a model for a TKI given orally",
    ]:
        print(f"\n> {q}")
        print(agent.chat(q))


if __name__ == "__main__":
    main()
