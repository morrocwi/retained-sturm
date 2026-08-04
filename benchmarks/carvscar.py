import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ.setdefault(v,"1")
import sys, time, random, statistics, warnings, json
warnings.filterwarnings("ignore")
import os as _os
sys.path.insert(0, _os.environ.get("IDM_REPO", _os.getcwd()))
import numpy as np
from pyslise import Pyslise
from retained_spectral.engine import (RawSpectralProblem, raw_benchmark_targets,
    retained_raw_input_readout, potential_values, warm_native_kernel)
from scipy.special import ai_zeros
warm_native_kernel(); DIR=np.array([0.0,1.0])

def dvr(V,L,N):
    x=np.linspace(-L,L,N+2)[1:-1]; h=x[1]-x[0]; i=np.arange(N); d=i[:,None]-i[None,:]
    T=np.where(d==0,0.0,2.0*(-1.0)**d/np.maximum(d**2,1))/(2*h**2); np.fill_diagonal(T,np.pi**2/(6*h**2))
    return np.linalg.eigvalsh(T+np.diag(V(x)))
a,ap,_,_=ai_zeros(8)
extra=[(RawSpectralProblem('abs_linear_low4','abs_linear',(('g',1.0),),4,2e-8),
        np.sort(np.concatenate([-a/2**(1/3),-ap/2**(1/3)]))[:4], None),
       (RawSpectralProblem('double_well_a2_9','symmetric_double_well',(('lam',1.0),('a2',9.0)),4,2e-8),
        dvr(lambda x:(x*x-9.0)**2,8.0,6000)[:4], (-8.0,8.0))]
items=[(t.problem,np.asarray(t.reference),None) for t in raw_benchmark_targets()]+extra

def matslise_values(p, win):
    V2=lambda x: 2.0*float(potential_values(p,np.array([x]))[0])
    s=Pyslise(V2,win[0],win[1],tolerance=p.tolerance)
    got=s.eigenvaluesByIndex(0,p.modes,DIR)
    if len(got)<p.modes:                       # fall back to an energy window
        lo=min(V2(x) for x in np.linspace(win[0],win[1],2001))
        got=s.eigenvalues(lo-1.0, lo+400.0, DIR)
    return np.array(sorted(e for i,e in got)[:p.modes])/2.0

REPEATS=21
def run_launch(seed):
    rng=random.Random(seed); out={}
    for p,ref,ow in items:
        r=retained_raw_input_readout(p); win=ow or r.window
        try:
            mv=matslise_values(p,win); merr=float(np.max(np.abs(mv-ref)))
        except Exception as ex:
            out[p.name]={"matslise_error":f"{type(ex).__name__}: {str(ex)[:60]}"}; continue
        rerr=float(np.max(np.abs(np.asarray(r.values)-ref)))
        arms={"rms": lambda: retained_raw_input_readout(p),
              "matslise": lambda: matslise_values(p,win)}
        for f in arms.values(): f()
        s={k:[] for k in arms}; names=list(arms)
        for _ in range(REPEATS):
            rng.shuffle(names)
            for n in names:
                t=time.perf_counter(); arms[n](); s[n].append(time.perf_counter()-t)
        out[p.name]={"rms_ms":statistics.median(s["rms"])*1e3,
                     "mat_ms":statistics.median(s["matslise"])*1e3,
                     "rms_err":rerr,"mat_err":merr,
                     "rms_samples":[x*1e3 for x in s["rms"]],
                     "mat_samples":[x*1e3 for x in s["matslise"]],
                     "window":list(win),"hand_window":ow is not None}
    return out

if __name__=="__main__":
    L=int(sys.argv[1]); json.dump(run_launch(20260804+L), open(f"carvscar{L}.json","w"), indent=2)
    print("launch",L,"done")
