# Leaf-conditioned order in critical Tokunaga trees and strong Horton asymptotics for Kingman's coalescent

Anonymous. Version 1.1.0-candidate. Unrefereed candidate. DOI: 10.5281/zenodo.23160583.

Start with [AI_INDEX.md](AI_INDEX.md), [manuscript.tex](manuscript.tex) and [manuscript.pdf](manuscript.pdf).

Two analytic results: an exact-leaf-conditioned Tokunaga order law for every c>1, including phase-dependent fluctuations and fixed moments; and a positive geometric asymptotic for deterministic Kingman branch densities. The latter takes tree size to infinity first at fixed order. No empirical river model, joint finite-size/order limit or certified decimal constant is claimed.

## Replay

Python 3.13; install requirements.txt in an isolated environment.

```sh
python verify_exact.py
python -O verify_exact.py
python verify_numeric.py --orders 40 --output numeric40.json
python -O verify_numeric.py --orders 60 --step .04 --output numeric60.json
python verify_stability.py numeric40.json numeric60.json
python negative_controls.py
python -O negative_controls.py
```

Generated JSON diagnostics may be overwritten by replay; retain the immutable release ZIP separately. Finite computation corroborates, but does not establish, the infinite analytical claims. See [REVIEW_RESPONSE.md](REVIEW_RESPONSE.md), [SOURCES.md](SOURCES.md), [ASSURANCE.md](ASSURANCE.md) and [LICENSES.md](LICENSES.md).
