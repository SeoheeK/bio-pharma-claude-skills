"""Dosing regimens: a sequence of dose events applied during simulation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

DoseRoute = Literal["iv_bolus", "infusion", "oral"]


@dataclass
class DoseEvent:
    time: float               # dose time [h]
    amount: float             # dose amount [mg]
    route: DoseRoute = "iv_bolus"
    tinf: float = 0.0         # infusion duration [h]; only for route="infusion"


@dataclass
class DosingRegimen:
    events: list[DoseEvent]

    @classmethod
    def single(cls, amount: float, route: DoseRoute = "iv_bolus",
               time: float = 0.0, tinf: float = 0.0) -> "DosingRegimen":
        return cls([DoseEvent(time, amount, route, tinf)])

    @classmethod
    def multiple(cls, amount: float, interval: float, n_doses: int,
                 route: DoseRoute = "iv_bolus", start: float = 0.0,
                 tinf: float = 0.0) -> "DosingRegimen":
        """Repeated fixed dose every ``interval`` hours."""
        return cls([
            DoseEvent(start + i * interval, amount, route, tinf)
            for i in range(n_doses)
        ])

    def duration(self) -> float:
        """Time of the last dosing-related event."""
        last = 0.0
        for e in self.events:
            last = max(last, e.time + (e.tinf if e.route == "infusion" else 0.0))
        return last
