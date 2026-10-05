#!/usr/bin/env python3
"""Floating-point diagnostics for the Kingman proof.

No output of this program is used as a proof of an infinite limit.
Dependencies: numpy, scipy, mpmath. Output is written next to this script.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import math
import warnings
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
import mpmath as mp

def require(condition, message):
    if not condition:
        raise ValueError(message)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orders',type=int,default=40)
    parser.add_argument('--step',type=float,default=0.08)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('numeric_checks.json'))
    parser.add_argument('--window',type=float,default=25.)
    args=parser.parse_args()
    if not 15 <= args.window <= 30:
        parser.error('window must lie in [15,30].')
    if not 8 <= args.orders <= 80 or not 0 < args.step <= .25:
        parser.error('orders must lie in [8,80], step in (0,0.25].')
    K=args.orders
    S=max(65.,2.0*K+15.)
    def fun(s:float,z:np.ndarray)->np.ndarray:
        v=np.exp(z-s)
        ret=np.zeros(K)
        ret[1:]=2*v[:-1]-np.exp(2*z[:-1]-z[1:]-s)
        return ret
    sol=solve_ivp(fun,[0.,S],np.zeros(K),method='DOP853',rtol=2e-12,
                  atol=2e-13,dense_output=True,max_step=args.step)
    if not sol.success:
        raise RuntimeError(sol.message)
    avec=np.exp(sol.y[:,-1]); ratios=avec[1:]/avec[:-1]
    rho=float(np.median(ratios[-min(8,len(ratios)):]))
    B=[]; quadrature=[]; warning_records=[]
    for k in range(K):
        center=float(np.log(avec[k]))
        def integrand(r:float)->float:
            ss=r+center
            if ss<0:return 0.
            defect=-float(np.expm1(sol.sol(ss)[k]-ss))
            return float(np.exp(-r)*defect**2)
        lo=max(-center,-args.window); hi=min(S-center,args.window)
        points=[lo]+[p for p in [-15.,-8.,0.,8.,15.] if lo<p<hi]+[hi]
        val=0.; err=0.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for left,right in zip(points,points[1:]):
                piece,piece_err=quad(integrand,left,right,epsabs=2e-12,epsrel=2e-12,limit=300)
                val+=piece; err+=piece_err
        warning_records.extend({'order':k+1,'category':w.category.__name__,'message':str(w.message)} for w in caught)
        # The exact profile envelope bounds only omitted integration tails.
        tail_bound=(1/(1+math.exp(-lo)) if lo>-center else 0.)+1/(1+math.exp(hi))
        quadrature.append({'order':k+1,'window':[lo,hi],'estimated_absolute_error':err,
                           'analytic_omitted_tail_bound':tail_bound})
        B.append(val)
    B=np.asarray(B); C=avec/rho**np.arange(1,K+1)
    scaledN=rho*B/C
    # Independent high-precision quadrature for the first nontrivial endpoint.
    mp.mp.dps=60
    a2exact=(1+mp.exp(2))/2
    def h2(x):return (1+mp.exp(2*x))/2
    def primitive(x):return x+(mp.exp(2*x)-1)/2
    a3exact=mp.exp(primitive(1))*(1-mp.quad(lambda x:h2(x)**2*mp.exp(-primitive(x)),[0,1]))
    N3exact=mp.exp(4)/128-mp.exp(2)/8+mp.mpf(233)/384
    a2err=abs(avec[1]-float(a2exact));a3err=abs(avec[2]-float(a3exact))
    N3err=abs(B[1]/avec[1]-float(N3exact))
    require(a2err<2e-10 and a3err<2e-9 and N3err<2e-10,'low-order comparison failed')
    require(abs(B[0]-1/3)<2e-10,'first integral comparison failed')
    # Full finite endpoint supermultiplicativity grid.
    super_gaps=[]
    for l in range(1,K):
        for j in range(1,K-l+1):
            lhs=avec[l+j-1]+1;rhs=(avec[l-1]+1)*(avec[j-1]+1)
            super_gaps.append((lhs-rhs)/max(lhs,rhs))
    require(min(super_gaps)>-2e-10,'endpoint supermultiplicativity failed')
    # Endpoint-normalized profiles, including their reciprocal extensions.
    us=np.r_[0.,np.geomspace(1e-6,1e6,2500)]
    def profile(k:int)->np.ndarray:
        aa=avec[k-1]
        ans=np.ones(len(us))
        inside=(us>0)&(us<=aa)
        ss=np.log(aa)-np.log(us[inside])
        ans[inside]=np.exp(sol.sol(ss)[k-1]-np.log(aa))
        outside=us>aa
        ans[outside]=1/us[outside]
        return ans
    qs=np.array([profile(k) for k in range(1,min(K,24)+1)])
    qgap=float(np.min(qs[:-1]-qs[1:]))
    lowergap=float(np.min(qs-1/(1+us)))
    upper=np.ones_like(us);upper[us>1]=1/us[us>1]
    uppergap=float(np.min(upper-qs))
    require(min(qgap,lowergap,uppergap)>-3e-10,'profile comparison failed')
    # Fractional-linear inequality throughout a finite x grid.
    xs=np.r_[np.linspace(0,.99,250),1.]
    def hval(k:int,xx:np.ndarray)->np.ndarray:
        out=np.empty_like(xx)
        mask=xx<1
        out[mask]=np.exp(sol.sol(-np.log1p(-xx[mask]))[k-1])
        out[~mask]=avec[k-1]
        return out
    mobiusgaps=[]
    for l in range(1,min(8,K//2)+1):
        b=avec[l-1]/(1+avec[l-1]);a=1/(1+avec[l-1])
        den=1-b*xs;tt=a*xs/den
        tt[-1]=1. # exact endpoint identity; avoid roundoff in a/(1-b).
        for j in range(1,min(8,K-l)+1):
            lhs=hval(l+j,xs)
            rhs=b/den+a/den**2*hval(j,tt)
            mobiusgaps.extend(((lhs-rhs)/np.maximum(1,np.maximum(lhs,rhs))).tolist())
    mobiusgap=float(min(mobiusgaps))
    require(mobiusgap>-5e-10,'Mobius grid comparison failed')
    out={
      'status':'floating-point diagnostics, not interval-certified',
      'parameters':{'orders':K,'s_endpoint':S,'max_step':args.step,
                    'method':'DOP853','rtol':2e-12,'atol':2e-13},
      'software':{'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__},
      'diagnostics':{'a2_absolute_error':a2err,'a3_independent_quadrature_absolute_error':a3err,
       'N3_absolute_error':N3err,'minimum_relative_endpoint_supermultiplicativity_gap':float(min(super_gaps)),
       'minimum_profile_monotonicity_gap':qgap,'minimum_profile_lower_bound_gap':lowergap,
       'minimum_profile_upper_bound_gap':uppergap,'minimum_relative_mobius_grid_gap':mobiusgap},
      'estimates':{'rho':rho,'B_star':float(B[-1]),'C_star':float(C[-1]),'H':float(scaledN[-1])},
      'a':avec.tolist(),'endpoint_ratios':ratios.tolist(),'B':B.tolist(),'C_using_estimated_rho':C.tolist(),
      'scaled_N_using_estimated_rho':scaledN.tolist(),
      'quadrature':quadrature,'warnings':warning_records,
      'constant_estimate_stability':'not assessed by a single run; use verify_stability.py',
      'error_scope':'Quadrature estimates and analytic tail bounds do not bound ODE or floating-point error.',
      'all_finite_checks_passed':True}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','parameters','diagnostics','estimates','all_finite_checks_passed']},indent=2))
    require(not warning_records,'numerical warnings recorded: accuracy status failed')
    require(all(math.isfinite(q['estimated_absolute_error']) and q['estimated_absolute_error']<2e-10 for q in quadrature),'quadrature accuracy status failed')

if __name__=='__main__':
    main()
