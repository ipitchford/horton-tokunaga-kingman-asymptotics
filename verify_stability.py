#!/usr/bin/env python3
"""Compare retained 40/60-order diagnostics; not a certified error enclosure."""
import argparse
import json
import math
from pathlib import Path

def assess(records):
    if len(records)<2:
        raise ValueError('at least two runs required')
    for r in records:
        if any(not math.isfinite(r['estimates'][key]) for key in ['rho','B_star','C_star','H']):
            raise ValueError('nonfinite constant estimate')
        if r.get('warnings') or not r.get('all_finite_checks_passed'):
            raise ValueError('warnings or failed finite checks')
        if not r.get('quadrature') or any(not math.isfinite(q['estimated_absolute_error']) or q['estimated_absolute_error']>2e-10 for q in r['quadrature']):
            raise ValueError('missing or failed quadrature diagnostics')
    if {r['parameters']['orders'] for r in records}!={40,60}:
        raise ValueError('require both 40 and 60 orders')
    if len({r['parameters']['max_step'] for r in records})<2:
        raise ValueError('require different integration steps')
    errors={key:max(abs(r['estimates'][key]-records[0]['estimates'][key]) for r in records) for key in ['rho','B_star','C_star','H']}
    if any(not math.isfinite(x) or x>2e-7 for x in errors.values()):
        raise ValueError('constant estimates failed cross-run stability: '+str(errors))
    return {'cross_run_stability_passed':True,'absolute_tolerance':2e-7,'maximum_differences':errors,
            'scope':'Empirical step/order sensitivity only; not interval certification or a bound on asymptotic bias.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records',nargs='+',type=Path)
    args=parser.parse_args()
    print(json.dumps(assess([json.loads(p.read_text()) for p in args.records]),indent=2))
