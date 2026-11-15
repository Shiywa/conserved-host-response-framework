"""Executed statistical functions only;current clinical entry is the262-profile workflow."""
import itertools,math
import numpy as np,pandas as pd
from scipy import stats
rng=np.random.default_rng(20261120)
def truth(x):return pd.Series(x).astype(str).str.lower().eq('true').to_numpy()

def padjust(p,method,n=None):
 p=np.asarray(p,float);good=np.flatnonzero(np.isfinite(p));out=np.full(len(p),np.nan);n=max(n or len(good),len(good))
 if not len(good):return out
 order=good[np.argsort(p[good],kind='stable')];rank=np.arange(1,len(order)+1);v=np.maximum.accumulate((n-rank+1)*p[order]) if method=='holm' else np.minimum.accumulate((p[order]*n/rank)[::-1])[::-1];out[order]=np.minimum(v,1);return out

def welch(y,g):
 y=np.asarray(y,float);g=np.asarray(g,bool);a=y[g];b=y[~g];empty=dict(effect=np.nan,se=np.nan,df=np.nan,CI95_low=np.nan,CI95_high=np.nan,p_value=np.nan)
 if min(len(a),len(b))<2 or not np.isfinite(y).all():return empty
 v1=np.var(a,ddof=1)/len(a);v0=np.var(b,ddof=1)/len(b);se=np.sqrt(v1+v0);eff=float(a.mean()-b.mean())
 if se<1e-15:return dict(effect=eff,se=0.,df=float(len(y)-2),CI95_low=eff,CI95_high=eff,p_value=1. if abs(eff)<1e-15 else 0.)
 df=(v1+v0)**2/(v1*v1/(len(a)-1)+v0*v0/(len(b)-1));q=stats.t.ppf(.975,df);return dict(effect=eff,se=float(se),df=float(df),CI95_low=float(eff-q*se),CI95_high=float(eff+q*se),p_value=float(2*stats.t.sf(abs(eff/se),df)))

def paired(v):
 v=np.asarray(v,float);n=len(v)
 if n<2 or not np.isfinite(v).all():return dict(effect=np.nan,se=np.nan,df=np.nan,CI95_low=np.nan,CI95_high=np.nan,p_value=np.nan)
 mean=float(v.mean());se=float(np.std(v,ddof=1)/np.sqrt(n));q=stats.t.ppf(.975,n-1);return dict(effect=mean,se=se,df=float(n-1),CI95_low=float(mean-q*se),CI95_high=float(mean+q*se),p_value=float(2*stats.t.sf(abs(mean/se),n-1)))

def ols_hc3(y,X,c):
 y=np.asarray(y,float);X=np.asarray(X,float);c=np.asarray(c,float);pinv=np.linalg.pinv(X);res=y-X@(pinv@y);h=np.einsum('ij,ji->i',X,pinv);w=c@pinv;se=float(np.sqrt(np.sum((w*res/(1-h))**2)));eff=float(w@y);df=len(y)-X.shape[1];q=stats.t.ppf(.975,df);return dict(effect=eff,se=se,df=float(df),CI95_low=float(eff-q*se),CI95_high=float(eff+q*se),p_value=float(2*stats.t.sf(abs(eff/se),df)))

def independent_permutation(y,g):
 y=np.asarray(y,float);g=np.asarray(g,bool);n=len(y);nt=int(g.sum());obs=abs(y[g].mean()-y[~g].mean());ncomb=math.comb(n,nt)
 if ncomb<=100000:
  ix=np.asarray(list(itertools.combinations(range(n),nt)));sums=y[ix].sum(axis=1);null=sums/nt-(y.sum()-sums)/(n-nt);return float(np.mean(np.abs(null)>=obs-1e-12)),len(null),'Exact group-label allocations'
 count=0
 for start in range(0,19999,1000):
  size=min(1000,19999-start);ix=np.argsort(rng.random((size,n)),axis=1)[:,:nt];sums=y[ix].sum(axis=1);null=sums/nt-(y.sum()-sums)/(n-nt);count+=int(np.sum(np.abs(null)>=obs-1e-12))
 return (1+count)/20000,19999,'Monte Carlo group-label permutation'

def signflip(v):
 v=np.asarray(v,float);n=len(v);obs=abs(v.mean());sign=np.asarray(list(itertools.product([-1,1],repeat=n)));null=sign@v/n;return float(np.mean(np.abs(null)>=obs-1e-12))
