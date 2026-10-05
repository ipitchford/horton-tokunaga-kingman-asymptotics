# Agent entry point

Version 1.1.0-candidate; creator Anonymous; DOI 10.5281/zenodo.23160583; unrefereed analytic candidate.

## Claims and assumptions
Read [manuscript.tex](manuscript.tex), Theorem 2.1, Proposition 3.4, the Section 4 corollary and Theorem 5.1. Critical Tokunaga c>1, gamma>0, leaves of order one; Kingman deterministic limiting densities, with size limit before order limit. [STATUS.json](STATUS.json) records separate claim and diagnostic states.

## Evidence and replay
[verify_exact.py](verify_exact.py) checks symbolic identities and finite rational probabilities, not universal analytic proofs. [verify_numeric.py](verify_numeric.py) performs noncertified finite ODE diagnostics. [verify_stability.py](verify_stability.py) checks step/order consistency, not rigorous accuracy. [negative_controls.py](negative_controls.py) requires six corruptions to be rejected. Commands and environment are in [README.md](README.md) and [requirements.txt](requirements.txt).

## Boundaries and objections
[ASSURANCE.md](ASSURANCE.md) separates internal checking, unavailable independent assurance and formal verification. [REVIEW_RESPONSE.md](REVIEW_RESPONSE.md) records all supplied review dispositions, including the original numerical drift. [SOURCES.md](SOURCES.md) describes the bounded primary-source comparison; [PROVENANCE.md](PROVENANCE.md) binds supplied archives. The c=2 classical case and existing Kingman exponent are not new claims. No empirical application is established.

## Reuse
[LICENSES.md](LICENSES.md) maps original prose/data CC0 and code MIT; third-party review and literature are not relicensed. [MANIFEST.sha256](MANIFEST.sha256) covers exact shipped files. New work should focus on independent analytic scrutiny, certified constants or joint limits, not infer these from current diagnostics.
