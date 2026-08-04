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

def hot(fn,reps):
    fn(); s=[]
    for _ in range(reps):
        t=time.perf_counter(); fn(); s.append(time.perf_counter()-t)
    return statistics.median(s)

tgt=[t for t in raw_benchmark_targets() if t.problem.name=="harmonic_low4"][0]
p=tgt.problem; win=retained_raw_input_readout(p).window
w=max(p.tolerance*0.01,2e-12); k=p.modes
print(f"harmonic_low4  matched width w={w:.1e}  k={k}")
print(f"{'n':>9} {'RMS ms':>10} {'DSTEBZ ms':>10} {'ratio':>7} {'reps':>5}")
for M in (768, 3072, 12288, 49152, 196608, 786432):
    d,e,_=retained_tridiagonal(p,win,M)
    reps = 11 if M<=12288 else (5 if M<=196608 else 3)
    tn=hot(lambda: native_eigvals_from_tridiagonal(d,e,k,w),reps)
    ts=hot(lambda: sla.eigh_tridiagonal(d,e,select="i",select_range=(0,k-1),
             eigvals_only=True,check_finite=False,tol=w,lapack_driver="stebz"),reps)
    print(f"{M-1:>9} {tn*1e3:>10.3f} {ts*1e3:>10.3f} {ts/tn:>7.3f} {reps:>5}")
    sys.stdout.flush()
