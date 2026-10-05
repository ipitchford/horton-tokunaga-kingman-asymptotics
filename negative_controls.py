"""Deliberate semantic corruption tests; active under Python -O."""
import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from verify_stability import assess

root=Path(__file__).resolve().parent
original=(root/'verify_exact.py').read_text()
for before,after in [('2*p*q-p*p','2*p*q+p*p'),('n*2**(2*n-1)','n*2**(2*n)')]:
    if before not in original:
        raise ValueError('mutation anchor missing')
    with tempfile.TemporaryDirectory() as d:
        target=Path(d)/'mutated.py'
        target.write_text(original.replace(before,after,1))
        p=subprocess.run([sys.executable,*(['-O'] if not __debug__ else []),str(target)],capture_output=True,text=True)
        if p.returncode==0:
            raise ValueError('semantic mutation incorrectly accepted: '+before)
records=[json.loads((root/n).read_text()) for n in ['numeric40.json','numeric60.json']]
assess(records)
for kind in ['warning','drift','missing-error','nan']:
    bad=copy.deepcopy(records)
    if kind=='warning':bad[1]['warnings']=[{'category':'IntegrationWarning','message':'injected'}]
    elif kind=='drift':bad[1]['estimates']['H']*=1.00245
    elif kind=='missing-error':bad[1]['quadrature']=[]
    else:bad[1]['estimates']['H']=float('nan')
    try:assess(bad)
    except ValueError:pass
    else:raise ValueError('accuracy mutation incorrectly accepted: '+kind)
print(json.dumps({'negative_controls':'passed','mutations':6,'optimized':not __debug__}))
