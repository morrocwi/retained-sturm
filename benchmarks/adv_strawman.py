import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ.setdefault(v,"1")
import sys,statistics,time
import os as _os
sys.path.insert(0, _os.environ.get("IDM_REPO", _os.getcwd()))
import numpy as np, scipy.linalg as sla
from retained_spectral.engine import raw_benchmark_targets, retained_raw_input_readout
from retained_spectral.competition.executor_audit import retained_tridiagonal
def hot(fn,r=9):
    fn(); s=[]
    for _ in range(r):
        t=time.perf_counter(); fn(); s.append(time.perf_counter()-t)
    return statistics.median(s)
tgt=[t for t in raw_benchmark_targets() if t.problem.name=="harmonic_low4"][0]
p=tgt.problem; k=p.modes; w=max(p.tolerance*0.01,2e-12)
def cost(win,M):
    d,e,_=retained_tridiagonal(p,win,M)
    return hot(lambda: sla.eigh_tridiagonal(d,e,select="i",select_range=(0,k-1),
        eigvals_only=True,check_finite=False,tol=w,lapack_driver="stebz"))
L1=[(( -8.0,8.0),M) for M in (128,256,512,1024,2048,4096)]
L2=[((-10.0,10.0),M) for M in (160,320,640,1280,2560,5120)]
c1=[cost(*x) for x in L1]; c2=[cost(*x) for x in L2]
tot_shipped=sum(c1)+sum(c2)
tot_retained=sum(c1)+c2[-1]
print("first ladder solves (ms):", [round(x*1e3,3) for x in c1])
print("second ladder solves (ms):", [round(x*1e3,3) for x in c2])
print(f"comparator solve time, as shipped (12 solves): {tot_shipped*1e3:.3f} ms")
print(f"comparator solve time, if it retained (7 solves): {tot_retained*1e3:.3f} ms")
print(f"fraction removed by retention: {(1-tot_retained/tot_shipped)*100:.1f}%")
print(f"measured e2e: RMS 4.439 ms, SciPy 20.262 ms, ratio 4.565")
sci_ret = 20.262 - (tot_shipped-tot_retained)*1e3
print(f"SciPy e2e if it retained (subtracting removed solves): {sci_ret:.3f} ms -> ratio {sci_ret/4.439:.3f}")
