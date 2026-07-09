# qFlow test suite

Automated tests for the qFlow browser tool (`website/qflow.html` / `code/qflow_tool.html`).

The point of this suite is that qFlow's numbers are trustworthy: the receptor-quantification math in the browser is the *same* math as the Python scripts in `code/`, and every other module (unit conversion, gating, indices, compensation, FCS parsing) is checked against known-correct values.

## What it does

`qflow_core.test.js`:

1. Extracts the pure compute core straight out of `website/qflow.html` (the block between the `QFLOW-CORE-START` / `QFLOW-CORE-END` markers) and loads it in Node. There is no second copy of the math — the test runs the exact code the page ships, so the two can't drift.
2. Runs `py_parity.py`, which calls the **real** functions in `code/receptor_quantification.py` and `code/flow_cytometry_validation.py` on a shared dataset, and compares the browser output to the Python output number for number (calibration curve, ABC conversion, net receptors, cell-line normalization, subtype lookup).
3. Exercises the rest: calibration edge cases, unit-conversion round trips, descriptive statistics, stain/separation index, compensation (true → measured → recovered round trip, matrix inversion, singular-matrix handling), and the FCS parser (builds a synthetic FCS 3.1 file and parses it back).

## Run it

```
cd ~/Desktop/tnbc-project
node test/qflow_core.test.js
```

Requires `node` and `python3` on PATH. Expected output: `96 passed, 0 failed`.

To test a specific HTML file (e.g. the repo dev copy) instead of the deployed one:

```
QFLOW_HTML=code/qflow_tool.html node test/qflow_core.test.js
```

## Files

- `qflow_core.test.js` — the harness and all assertions.
- `py_parity.py` — imports the real Python functions and prints their output as JSON.
- `test_data.json` — the shared calibration/measurement dataset used by both sides.
