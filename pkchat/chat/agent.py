"""The chatbot layer.

``PKChatAgent`` turns natural-language requests into calls against the PK tool
surface. Two backends, selected automatically:

* **Claude** -- when the ``anthropic`` SDK is installed and credentials are
  available, the agent runs a tool-use loop: Claude reads the request, calls
  the registered PK tools with structured arguments, and explains the results.
* **Local router** -- a dependency-free intent parser that recognises common PK
  requests and dispatches to the same tools. This keeps the system fully
  functional with zero external packages.

Both backends share the tool registry in ``tools.py``, so the numeric behaviour
is identical; only the language understanding differs.
"""
from __future__ import annotations

import re

from .tools import TOOLS, dispatch, tool_schemas

SYSTEM_PROMPT = (
    "You are a pharmacokinetics (PK) modelling assistant. You help pharmacometricians "
    "simulate concentration-time profiles, run non-compartmental analysis, estimate "
    "compartmental model parameters, scale animal data to humans allometrically, run "
    "population simulations, and find similar published PopPK cases. "
    "Use the provided tools to do the actual computation -- never invent numeric "
    "results. Report units, state assumptions, and remind the user that predictions "
    "are model-based and require clinical/expert review before any decision."
)

DEFAULT_MODEL = "claude-opus-4-8"


class PKChatAgent:
    def __init__(self, model: str = DEFAULT_MODEL, use_claude: bool | None = None,
                 max_steps: int = 6):
        self.model = model
        self.max_steps = max_steps
        self._client = None
        self.backend = "local"
        if use_claude is not False:
            self._client = _try_claude_client()
            if self._client is not None:
                self.backend = "claude"

    # -- public API ---------------------------------------------------------

    def chat(self, message: str) -> str:
        if self.backend == "claude" and self._client is not None:
            try:
                return self._chat_claude(message)
            except Exception as e:
                return (f"[Claude backend error: {e}]\n"
                        + LocalRouter().handle(message))
        return LocalRouter().handle(message)

    # -- Claude backend -----------------------------------------------------

    def _chat_claude(self, message: str) -> str:
        client = self._client
        tools = tool_schemas()
        messages = [{"role": "user", "content": message}]
        transcript: list[str] = []

        for _ in range(self.max_steps):
            resp = client.messages.create(
                model=self.model,
                max_tokens=16000,
                system=SYSTEM_PROMPT,
                tools=tools,
                messages=messages,
            )
            if resp.stop_reason == "refusal":
                return "[Request declined by safety policy.]"

            # Collect text and any tool calls.
            tool_results = []
            for block in resp.content:
                if block.type == "text":
                    transcript.append(block.text)
                elif block.type == "tool_use":
                    result = dispatch(block.name, dict(block.input))
                    tool_results.append({
                        "type": "tool_result",
                        "tool_use_id": block.id,
                        "content": result,
                    })

            if resp.stop_reason != "tool_use":
                break

            messages.append({"role": "assistant", "content": resp.content})
            messages.append({"role": "user", "content": tool_results})

        return "\n".join(t for t in transcript if t).strip() or "(no response)"


def _try_claude_client():
    try:
        import anthropic  # noqa
    except Exception:
        return None
    try:
        client = anthropic.Anthropic()
        # Cheap credential probe deferred to first real call; assume ok if
        # constructor succeeds and a key source is present.
        import os
        if not (os.environ.get("ANTHROPIC_API_KEY")
                or os.environ.get("ANTHROPIC_AUTH_TOKEN")):
            return None
        return client
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Local rule-based router (no external dependencies).
# ---------------------------------------------------------------------------

_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


class LocalRouter:
    """Parses common PK requests and dispatches to the tool surface.

    This is intentionally simple keyword/regex intent detection -- it covers the
    documented request patterns, not arbitrary free text. When it cannot map a
    request it returns guidance listing what it understands.
    """

    def handle(self, message: str) -> str:
        text = message.lower()
        intent = self._intent(text)
        if intent == "nca":
            return self._nca(message)
        if intent == "fit":
            return self._fit(message)
        if intent == "allometry":
            return self._allometry(message)
        if intent == "recommend":
            return self._recommend(message)
        if intent == "population":
            return self._population(message)
        if intent == "simulate":
            return self._simulate(message)
        return self._help()

    # -- intent detection ---------------------------------------------------

    def _intent(self, text: str) -> str:
        species_hits = sum(1 for a in ("mouse", "rat", "rabbit", "monkey", "dog")
                           if a in text)
        if ("allometr" in text or "scale to human" in text or "scaling" in text
                or "animal to human" in text or species_hits >= 2):
            return "allometry"
        if any(k in text for k in ["nca", "non-compartmental", "noncompartmental", "auc from", "observed"]):
            return "nca"
        if any(k in text for k in ["fit", "estimate", "regress"]):
            return "fit"
        if any(k in text for k in ["population", "virtual", "iiv", "variability", "prediction interval"]):
            return "population"
        if any(k in text for k in ["recommend", "similar", "case", "initial estimate", "which model"]):
            return "recommend"
        if any(k in text for k in ["simulate", "profile", "concentration", "cmax", "dose", "predict conc"]):
            return "simulate"
        return "help"

    # -- extractors ---------------------------------------------------------

    def _num(self, text, *labels, default=None):
        for lab in labels:
            m = re.search(rf"{lab}\s*[=:]?\s*({_NUM})", text, re.I)
            if m:
                return float(m.group(1))
        return default

    def _route(self, text):
        t = text.lower()
        if "infus" in t:
            return "infusion"
        if "oral" in t or "po " in t or " po" in t or "ka" in t:
            return "oral"
        return "iv_bolus"

    def _arrays(self, text):
        """Extract two numeric arrays (times, concentrations) from bracketed
        lists or 'times=... conc=...' phrasing."""
        arrs = re.findall(r"\[([^\]]+)\]", text)
        parsed = []
        for a in arrs:
            nums = [float(x) for x in re.findall(_NUM, a)]
            if nums:
                parsed.append(nums)
        if len(parsed) >= 2:
            return parsed[0], parsed[1]
        return None, None

    # -- handlers -----------------------------------------------------------

    def _simulate(self, text):
        args = {
            "CL": self._num(text, "cl", "clearance", default=5.0),
            "V1": self._num(text, "v1", "v", "volume", default=30.0),
            "dose": self._num(text, "dose", "mg", default=100.0),
            "route": self._route(text),
        }
        q = self._num(text, "q")
        v2 = self._num(text, "v2")
        if q and v2:
            args["Q"], args["V2"] = q, v2
        ka = self._num(text, "ka")
        if ka:
            args["ka"] = ka
            args["route"] = "oral"
        tinf = self._num(text, "tinf", "infusion")
        if tinf and args["route"] == "infusion":
            args["tinf"] = tinf
        interval = self._num(text, "interval", "q", "every")
        n = self._num(text, "n_doses", "doses")
        if interval and n and n > 1:
            args["interval"] = interval
            args["n_doses"] = int(n)
        return dispatch("simulate_pk", args)

    def _nca(self, text):
        times, concs = self._arrays(text)
        if not times or not concs:
            return ("For NCA give two bracketed lists, e.g.\n"
                    "  'run NCA times=[0.5,1,2,4,8] conc=[8,6,4,2,0.5] dose=100 iv'")
        args = {"times": times, "concentrations": concs,
                "route": self._route(text)}
        d = self._num(text, "dose", "mg")
        if d:
            args["dose"] = d
        return dispatch("run_nca", args)

    def _fit(self, text):
        times, concs = self._arrays(text)
        if not times or not concs:
            return ("For fitting give times and concentrations as bracketed lists, "
                    "plus the dose, e.g.\n"
                    "  'fit CL,V1 to times=[0.5,1,2,4,8,12,24] "
                    "conc=[3,2.8,2.4,1.7,0.9,0.45,0.06] dose=100 iv'")
        args = {"times": times, "concentrations": concs,
                "dose": self._num(text, "dose", "mg", default=100.0),
                "route": self._route(text)}
        est = re.search(r"(?:fit|estimate)\s+([a-z0-9, ]+?)\s+(?:to|from|for)", text, re.I)
        if est:
            params = [p.strip().upper().replace("VD", "V1").replace("V", "V1") if p.strip().lower() in ("v", "vd") else p.strip().capitalize()
                      for p in re.split(r"[ ,]+", est.group(1)) if p.strip()]
            # Normalise common names.
            norm = []
            for p in params:
                pl = p.lower()
                if pl in ("cl", "clearance"):
                    norm.append("CL")
                elif pl in ("v", "v1", "vd", "volume"):
                    norm.append("V1")
                elif pl in ("q",):
                    norm.append("Q")
                elif pl in ("v2",):
                    norm.append("V2")
                elif pl in ("ka",):
                    norm.append("ka")
            if norm:
                args["estimate"] = norm
        ka = self._num(text, "ka")
        if ka:
            args["ka"] = ka
        return dispatch("fit_pk_model", args)

    def _allometry(self, text):
        # Parse "species weight clearance" triples like "rat 0.25 3.5".
        species = []
        for name in ["mouse", "rat", "rabbit", "monkey", "dog"]:
            m = re.search(rf"{name}\s*[:=]?\s*({_NUM})\s*(?:kg)?\s*[, ]\s*({_NUM})", text, re.I)
            if m:
                species.append({"species": name,
                                "body_weight_kg": float(m.group(1)),
                                "clearance": float(m.group(2))})
        if len(species) < 2:
            return ("For allometric scaling list >= 2 species as "
                    "'name weight clearance', e.g.\n"
                    "  'scale to human: mouse 0.02 0.8, rat 0.25 3.5, dog 10 45, monkey 5 28'")
        return dispatch("allometric_scaling", {"species": species})

    def _recommend(self, text):
        text = text.lower()
        args = {}
        classes = {
            "checkpoint": "immune_checkpoint", "pd-1": "immune_checkpoint",
            "pd-l1": "immune_checkpoint", "tki": "tki", "kinase": "tki",
            "calcineurin": "calcineurin_inhibitor", "tacrolimus": "calcineurin_inhibitor",
            "cyclosporine": "calcineurin_inhibitor", "vancomycin": "glycopeptide",
            "antifungal": "azole_antifungal", "azole": "azole_antifungal",
            "alkylat": "alkylator", "busulfan": "alkylator", "her2": "anti_her2",
            "vegf": "anti_vegf", "tnf": "anti_tnf",
        }
        for k, v in classes.items():
            if k in text:
                args["drug_class"] = v
                break
        if "mab" in text or "antibod" in text or "biologic" in text:
            args["modality"] = "mAb"
        elif "small molecule" in text or "small-molecule" in text:
            args["modality"] = "small_molecule"
        # Word-boundary route match so "given" does not read as "iv".
        if re.search(r"\biv\b|intravenous", text):
            args["route"] = "iv"
        elif re.search(r"\boral\b|orally|\bpo\b", text):
            args["route"] = "oral"
        elif re.search(r"\bsc\b|subcut", text):
            args["route"] = "sc"
        return dispatch("recommend_pk_case", args)

    def _population(self, text):
        args = {
            "CL": self._num(text, "cl", "clearance", default=5.0),
            "V1": self._num(text, "v1", "v", "volume", default=30.0),
            "dose": self._num(text, "dose", "mg", default=100.0),
            "route": self._route(text),
            "n_subjects": int(self._num(text, "subjects", "n", default=200)),
        }
        cv_cl = self._num(text, "iiv", "cv")
        iiv = {"CL": (cv_cl / 100 if cv_cl and cv_cl > 1 else cv_cl) if cv_cl else 0.3,
               "V1": 0.2}
        args["iiv_cv"] = iiv
        return dispatch("simulate_population", args)

    def _help(self):
        return (
            "PK chat -- I can:\n"
            "  * simulate a profile: 'simulate CL=5 V1=30 dose=100 iv'\n"
            "  * run NCA: 'NCA times=[0.5,1,2,4,8] conc=[8,6,4,2,0.5] dose=100 iv'\n"
            "  * fit a model: 'fit CL,V1 to times=[...] conc=[...] dose=100'\n"
            "  * allometric scaling: 'scale to human: rat 0.25 3.5, dog 10 45, monkey 5 28'\n"
            "  * population sim: 'population CL=5 V1=30 dose=100 iiv=30% n=200'\n"
            "  * find cases: 'recommend a PopPK model for a PD-1 mAb given IV'\n"
        )
