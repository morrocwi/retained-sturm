import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ.setdefault(v,"1")
import sys,statistics,time,random,json
from pathlib import Path

_THIS_DIR = Path(__file__).resolve().parent
_DEFAULT_IDM_REPO = _THIS_DIR.parent.parent / "information-discrete-math"
sys.path.insert(0, os.environ.get("IDM_REPO_PATH", str(_DEFAULT_IDM_REPO)))
import numpy as np, scipy.linalg as sla
from retained_spectral.engine import native_eigvals_from_tridiagonal, warm_native_kernel
sys.path.insert(0, str(_THIS_DIR)); from adv_kernel import wilkinson, glued_wilkinson, toeplitz_121, scaled_extremes
warm_native_kernel(); rng=random.Random(1)
def paired(f,g,reps):
    f();g();a=[];b=[]
    for _ in range(reps):
        for w in ([0,1] if rng.random()<.5 else [1,0]):
            t=time.perf_counter(); (f if w==0 else g)(); (a if w==0 else b).append(time.perf_counter()-t)
    return statistics.median(a)*1e3, statistics.median(b)*1e3
cases=[]
d,e=glued_wilkinson(8); cases += [("glued_wilkinson k=64",d,e,64,2e-10),("glued_wilkinson k=168 (all)",d,e,168,2e-10)]
d,e=toeplitz_121(4000); cases += [("toeplitz_121 k=256",d,e,256,2e-10),("toeplitz_121 k=1024",d,e,1024,2e-10)]
d,e=scaled_extremes(2000); nrm=float(np.max(np.abs(d))+2*np.max(np.abs(e)))
cases += [("scaled_extremes k=4 REL tol",d,e,4,2e-10*nrm)]
d,e=wilkinson(200); cases += [("wilkinson_w401 k=64",d,e,64,2e-10)]
print(f"{'case':30s} {'n':>6} {'k':>5} {'nat err':>11} {'ste err':>11} {'nat ms':>9} {'ste ms':>9} {'ratio':>6}")
out={}
for name,d,e,k,tol in cases:
    n=len(d); dense=np.diag(d)+np.diag(e,1)+np.diag(e,-1)
    ref=np.sort(np.linalg.eigvalsh(dense))[:k]
    try:
        nat=np.sort(np.asarray(native_eigvals_from_tridiagonal(d,e,k,tol))); ne=float(np.max(np.abs(nat-ref)))
    except Exception as ex: ne=float('inf'); print("  native raised",type(ex).__name__,ex)
    try:
        ste=np.sort(np.asarray(sla.eigh_tridiagonal(d,e,select="i",select_range=(0,k-1),eigvals_only=True,check_finite=False,tol=tol,lapack_driver="stebz"))); se=float(np.max(np.abs(ste-ref)))
    except Exception as ex: se=float('inf'); print("  dstebz raised",type(ex).__name__,ex)
    mn,ms=paired(lambda: native_eigvals_from_tridiagonal(d,e,k,tol),
        lambda: sla.eigh_tridiagonal(d,e,select="i",select_range=(0,k-1),eigvals_only=True,check_finite=False,tol=tol,lapack_driver="stebz"),5)
    print(f"{name:30s} {n:>6} {k:>5} {ne:>11.3e} {se:>11.3e} {mn:>9.3f} {ms:>9.3f} {ms/mn:>6.3f}")
    out[name]={"n":n,"k":k,"native_err":ne,"dstebz_err":se,"native_ms":mn,"dstebz_ms":ms,"ratio":ms/mn}
    sys.stdout.flush()
_out_path = _THIS_DIR / "results" / "adversarial_kernel2.json"
_out_path.parent.mkdir(parents=True, exist_ok=True)
json.dump(out,open(_out_path,"w"),indent=2)
