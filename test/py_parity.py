"""Runs the ACTUAL Python functions from code/ on the shared test dataset and
prints JSON, so the qFlow browser core can be compared against them number for
number. Imports are safe: both scripts guard their work behind __main__."""
import json, sys, os

CODE_DIR = os.path.join(os.path.dirname(__file__), "..", "code")
sys.path.insert(0, CODE_DIR)

import receptor_quantification as rq
import flow_cytometry_validation as fcv

with open(os.path.join(os.path.dirname(__file__), "test_data.json")) as f:
    data = json.load(f)

beads = [(p[0], p[1]) for p in data["beads"]]
lin_beads = [(p[0], p[1]) for p in data["linear_beads"]]

curve = rq.fit_calibration_curve(beads)          # auto model selection
lin_curve = rq.fit_calibration_curve(lin_beads)

def quant(c):
    out = []
    for m in data["measurements"]:
        at = rq.apply_curve(c, m["mfi"])
        ai = rq.apply_curve(c, m["isotype"])
        out.append({"abc_target": at, "abc_isotype": ai, "net": at - ai})
    return out

result = {
    "curve": curve,
    "linear_curve": lin_curve,
    "quant": quant(curve),
    "normalize": {
        "MDA-MB-453": fcv.normalize_cell_line("MDA-MB-453"),
        "mda mb 453": fcv.normalize_cell_line("mda mb 453"),
        "HCC70": fcv.normalize_cell_line("HCC70"),
    },
    "subtype_lookup": {k: v for k, v in fcv.CELL_LINE_SUBTYPE.items()},
}
print(json.dumps(result))
