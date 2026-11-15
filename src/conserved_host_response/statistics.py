"""Statistical kernels used by current score analyses; no legacy Wilcoxon default."""
from __future__ import annotations
import itertools,math
import numpy as np
from scipy import stats

def adjust(p,method="BH",*,family_size=None):
    p=np.asarray(p,dtype=float);n=len(p) if family_size is None else family_size
    if n<len(p) or method not in ("BH","holm"):raise ValueError("Declare the original BH/Holm family size,including unavailable positions")
    valid=np.isfinite(p);out=np.full(len(p),np.nan)
    ix=np.flatnonzero(valid);ix=ix[np.argsort(p[ix],kind="stable")]
    if not len(ix):return out
    r=np.arange(1,len(ix)+1)
    q=np.maximum.accumulate((n-r+1)*p[ix]) if method=="holm" else np.minimum.accumulate((p[ix]*n/r)[::-1])[::-1]
    out[ix]=np.minimum(q,1);return out

def welch(case,control):
    case=np.asarray(case,float);control=np.asarray(control,float)
    if min(len(case),len(control))<2 or not np.isfinite(np.r_[case,control]).all():raise ValueError("Finite independent biological-unit values required")
    eff=case.mean()-control.mean();a=case.var(ddof=1)/len(case);b=control.var(ddof=1)/len(control);se=np.sqrt(a+b)
    if se<1e-15:return dict(effect=float(eff),CI95_low=float(eff),CI95_high=float(eff),p_value=1. if abs(eff)<1e-15 else 0.,df=len(case)+len(control)-2)
    df=(a+b)**2/(a*a/(len(case)-1)+b*b/(len(control)-1));q=stats.t.ppf(.975,df)
    return dict(effect=eff,CI95_low=eff-q*se,CI95_high=eff+q*se,p_value=2*stats.t.sf(abs(eff/se),df),df=df)

def paired(changes):
    v=np.asarray(changes,float);n=len(v)
    if n<2 or not np.isfinite(v).all():raise ValueError("Finite paired biological-unit changes required")
    eff=v.mean();se=v.std(ddof=1)/np.sqrt(n);q=stats.t.ppf(.975,n-1)
    return dict(effect=eff,CI95_low=eff-q*se,CI95_high=eff+q*se,p_value=2*stats.t.sf(abs(eff/se),n-1),df=n-1)

def exact_signflip(changes):
    v=np.asarray(changes,float)
    if not len(v) or not np.isfinite(v).all():raise ValueError('Finite nonempty paired changes required')
    if len(v)>20:raise ValueError("Exact enumeration is restricted to20pairs; use the declared Monte Carlo workflow for larger analyses")
    signs=np.array(list(itertools.product([-1,1],repeat=len(v))))
    return float(np.mean(abs(signs@v/len(v))>=abs(v.mean())-1e-12))

def hc3(y,X,contrast):
    y=np.asarray(y,float);X=np.asarray(X,float);c=np.asarray(contrast,float)
    if np.linalg.matrix_rank(X)!=X.shape[1] or len(y)<=X.shape[1]:raise ValueError("Full rank and positive residual degrees of freedom required")
    pinv=np.linalg.pinv(X);res=y-X@(pinv@y);h=np.einsum('ij,ji->i',X,pinv);w=c@pinv
    effect=float(w@y);se=float(np.sqrt(np.sum((w*res/(1-h))**2)));df=len(y)-X.shape[1];q=stats.t.ppf(.975,df)
    return dict(effect=effect,CI95_low=effect-q*se,CI95_high=effect+q*se,p_value=2*stats.t.sf(abs(effect/se),df),df=df)
