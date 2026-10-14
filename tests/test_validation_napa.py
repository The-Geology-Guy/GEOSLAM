"""Regression test on real terrain: 2014 South Napa (Worrell 2016 inputs).

The fault-plane solution (NP1, strike 345) must place most of the mapped West
Napa fault inside its seismo-lineament; the auxiliary plane (NP2) must not.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))


def test_napa_np1_contains_west_napa_fault():
    import matplotlib
    matplotlib.use("Agg")
    import napa_2014_validation as v

    rows = {(t, m): (a, f) for t, m, a, f, _ in v.main()}
    assert rows[("worrell_np1", "envelope")][1] > 0.7
    assert rows[("worrell_np2", "envelope")][1] < 0.2
    assert rows[("comcat_np1", "envelope")][1] > 0.6
    # the envelope always contains Cronin's 7-slab union
    assert rows[("worrell_np1", "envelope")][0] > rows[("worrell_np1", "cronin")][0]
