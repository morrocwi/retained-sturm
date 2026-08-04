import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ.setdefault(v,"1")
import sys, statistics, time
import os as _os
sys.path.insert(0, _os.environ.get("IDM_REPO", _os.getcwd()))
import numpy as np, scipy.linalg as sla
from retained_spectral.engine import (raw_benchmark_targets, retained_raw_input_readout,
    native_eigvals_from_tridiagonal, warm_native_kernel)
from retained_spectral.competition.executor_audit import retained_tridiagonal
warm_native_kernel()

def hot(fn, reps=15):
    fn()
    s=[]
    for _ in range(reps):
        t=time.perf_counter(); r=fn(); s.append(time.perf_counter()-t)
    return statistics.median(s), r

INT=768
rows=[]
for tgt in raw_benchmark_targets():
    p=tgt.problem
    win=retained_raw_input_readout(p).window
    d,e,_=retained_tridiagonal(p,win,INT)
    k=p.modes
    eps=p.tolerance
    w=max(eps*0.01,2e-12)
    norm1=float(np.max(np.abs(d))+2*np.max(np.abs(e)))
    res={}
    # native at the two candidate widths
    for lbl,tw in (("nat@eps",eps),("nat@0.01eps",w)):
        res[lbl]=hot(lambda tw=tw: native_eigvals_from_tridiagonal(d,e,k,tw))[0]
    # scipy: default, matched-eps, matched-0.01eps
    for lbl,kw in (("sci@default",{}),("sci@eps",{"tol":eps}),("sci@0.01eps",{"tol":w})):
        res[lbl]=hot(lambda kw=kw: sla.eigh_tridiagonal(d,e,select="i",select_range=(0,k-1),
                       eigvals_only=True,check_finite=False,lapack_driver="stebz",**kw))[0]
    rows.append((p.name,eps,norm1,res))
    print(f"{p.name:32s} eps={eps:.1e} |T|1={norm1:.3e} epsmach*|T|1={2.22e-16*norm1:.2e}")
    print(f"   AS-SHIPPED  nat@eps/sci@default = {res['sci@default']/res['nat@eps']:.3f}")
    print(f"   MATCHED eps      = {res['sci@eps']/res['nat@eps']:.3f}")
    print(f"   MATCHED 0.01eps  = {res['sci@0.01eps']/res['nat@0.01eps']:.3f}")
    sys.stdout.flush()

import math
def gm(x): return math.exp(sum(math.log(v) for v in x)/len(x))
print()
print("GEOMEAN as-shipped:", round(gm([r[3]['sci@default']/r[3]['nat@eps'] for r in rows]),4))
print("GEOMEAN matched eps:", round(gm([r[3]['sci@eps']/r[3]['nat@eps'] for r in rows]),4))
print("GEOMEAN matched 0.01eps:", round(gm([r[3]['sci@0.01eps']/r[3]['nat@0.01eps'] for r in rows]),4))
