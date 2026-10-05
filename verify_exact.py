#!/usr/bin/env python3
"""Exact algebra and finite-distribution checks. These do not replace the proofs."""
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path
import sympy as s

def require(condition, message='exact check failed'):
    if not condition:
        raise ValueError(message)

x,a,b,p,q,pprime,qprime=s.symbols('x a b p q pprime qprime')
D=1-b*x
A=a/D**2
B=b/D
# q' = 2pq-p^2 after the independent-variable substitution t'=A.
transformed_derivative=s.diff(B,x)+s.diff(A,x)*q+A*A*(2*p*q-p*p)
transformed_rhs=2*(B+A*p)*(B+A*q)-(B+A*p)**2
require(s.factor(transformed_derivative-transformed_rhs)==0)
# Cooperative comparison identity for two ordered hierarchies.
u,v=s.symbols('u v')
require(s.expand(2*p*u-p*p-2*q*v+q*q - ((p+q)*(u-v)+(p-q)*(u+v-p-q)))==0)
z,c=s.symbols('z c')
F2=z*z/(c-(c-1)*z)
require(s.simplify(s.diff(F2,z).subs(z,1)-(c+1))==0)
# First-moment recurrence, with R=2c.
k=s.symbols('k',integer=True,positive=True)
R=s.symbols('R',positive=True)
mean=(R**k+R-2)/(2*(R-1))
require(s.simplify(mean.subs(k,k+1)-(R*mean-(R-2)/2))==0)

N=64
out={"symbolic_mobius_identity":True,"symbolic_comparison_identity":True,
     "symbolic_mean_recurrence":True,"maximum_leaf_count":N,"cases":[]}
for c0 in [Q(5,4),Q(3,2),Q(2),Q(7,2)]:
    pp=1/c0; qq=1-pp
    cur=[Q(0) for _ in range(N+1)];cur[1]=Q(1)
    levels=[cur]
    for lev in range(2,N.bit_length()+1):
        nxt=[Q(0) for _ in range(N+1)]
        for n in range(2,N+1):
            nxt[n]=sum((cur[m]*comb(n-2,2*m-2)*pp**(2*m-1)*qq**(n-2*m)
                         for m in range(1,n//2+1) if cur[m]),Q(0))
        require(all(t>=0 for t in nxt))
        require(all(nxt[n]==0 for n in range(min(2**(lev-1),N+1))))
        levels.append(nxt);cur=nxt
    marginal=[sum((row[n]/2**(lev+1) for lev,row in enumerate(levels)),Q(0))
              for n in range(N+1)]
    if c0==2:
        for n in range(1,N+1):
            expected=Q(comb(2*n-2,n-1),n*2**(2*n-1))
            require(marginal[n]==expected,str((n,marginal[n],expected)))
    out['cases'].append({"c":str(c0),"finite_probabilities_nonnegative":True,
                          "minimum_order_support_correct":True,
                          "c2_catalan_identity":c0==2,
                          "P_N_2":str(marginal[2]),"P_N_3":str(marginal[3])})
p=Path(__file__).with_name('exact_checks.json')
p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
