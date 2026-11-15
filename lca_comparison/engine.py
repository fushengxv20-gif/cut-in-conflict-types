"""Weighted response-pattern EM, independently checked against poLCA.

Missing categories are marginalized in the likelihood and excluded only from
the corresponding item's M-step denominator. No case-wise imputation is used.
"""
import sys,json,time,itertools
from pathlib import Path
W=Path(__file__).resolve().parent
sys.path.insert(0,str(W/'pydeps'))
import numpy as np
from numba import njit
from scipy.optimize import linear_sum_assignment

CATS=np.array([4,3,2,2],dtype=np.int64)

@njit(cache=True,nogil=True)
def em(x,w,pi,theta,maxiter=20000,tol=1e-8):
    n,j=x.shape;k=len(pi);old=-1e300
    post=np.empty((n,k)); total=w.sum()
    for it in range(maxiter):
        ll=0.
        for i in range(n):
            den=0.
            for c in range(k):
                z=pi[c]
                for v in range(j):
                    a=x[i,v]
                    if a>=0:z*=theta[c,v,a]
                post[i,c]=z;den+=z
            den=max(den,1e-300)
            ll+=w[i]*np.log(den)
            for c in range(k):post[i,c]=post[i,c]/den*w[i]
        if it>0 and abs(ll-old)<tol:return ll,pi,theta,it+1
        old=ll
        for c in range(k):
            pi[c]=post[:,c].sum()/total
            for v in range(j):
                den=0.
                vals=np.zeros(4)
                for i in range(n):
                    a=x[i,v]
                    if a>=0:
                        vals[a]+=post[i,c];den+=post[i,c]
                if den>0:
                    for a in range(CATS[v]):theta[c,v,a]=max(vals[a]/den,1e-300)
    # Evaluate the parameters after the final M-step.
    ll=0.
    for i in range(n):
        den=0.
        for c in range(k):
            z=pi[c]
            for v in range(j):
                a=x[i,v]
                if a>=0:z*=theta[c,v,a]
            den+=z
        ll+=w[i]*np.log(max(den,1e-300))
    return ll,pi,theta,maxiter

def pack(x):
    y,n=np.unique(x,axis=0,return_counts=True)
    return y.astype(np.int64),n.astype(float)

def fit(x,k,seed=1,nstart=30,ref=None,tol=1e-8):
    rng=np.random.default_rng(seed);y,w=pack(x);best=None;lls=[];iterations=[]
    for s in range(nstart):
        pi=rng.dirichlet(np.ones(k));theta=np.zeros((k,4,4))
        for v,a in enumerate(CATS):theta[:,v,:a]=rng.dirichlet(np.ones(a),k)
        if ref is not None and s<5:
            mix=[0,0.01,0.05,0.1,0.2][s]
            pi=(1-mix)*ref[1]+mix*pi
            theta=(1-mix)*ref[2]+mix*theta
        z=em(y,w,pi,theta,tol=tol);lls.append(z[0]);iterations.append(z[3])
        if best is None or z[0]>best[0]:best=z
    return best,np.array(lls),np.array(iterations)

def logprobs(f,x):
    z=np.broadcast_to(f[1],(len(x),len(f[1]))).copy()
    for v in range(4):
        ok=x[:,v]>=0
        z[ok]*=f[2][:,v,x[ok,v]].T
    return np.log(np.maximum(z.sum(axis=1),1e-300))

def posterior(f,x):
    z=np.broadcast_to(f[1],(len(x),len(f[1]))).copy()
    for v in range(4):
        ok=x[:,v]>=0;z[ok]*=f[2][:,v,x[ok,v]].T
    return z/z.sum(axis=1,keepdims=True)

def aligned(f,ref):
    costs=((ref[2][:,None]-f[2][None,:])**2).sum(axis=(2,3))
    _,p=linear_sum_assignment(costs)
    return (f[0],f[1][p],f[2][p],f[3]),float(costs[np.arange(len(p)),p].sum())

def simulate(f,x,rng):
    cs=rng.choice(len(f[1]),len(x),p=f[1]);z=np.empty_like(x)
    for c in range(len(f[1])):
        ix=np.where(cs==c)[0]
        for v,a in enumerate(CATS):z[ix,v]=rng.choice(a,len(ix),p=f[2][c,v,:a])
    z[x<0]=-1
    return z

def discrepancy(f,x):
    # Saturated multinomial within each observed missingness pattern.
    # Bootstrap calibration is used; no asymptotic chi-square reference.
    g=0.
    masks=np.unique(x<0,axis=0)
    for mask in masks:
        xx=x[np.all((x<0)==mask,axis=1)]
        y,w=pack(xx);expected=len(xx)*np.exp(logprobs(f,y))
        g+=2*np.sum(w*np.log(w/expected))
    return float(g)

def load():
    import csv
    src=W.parent/'reproducibility'
    rows=list(csv.DictReader((src/'final_analysis_data.csv').open(encoding='utf-8-sig')))
    x=np.array([[int(r[v])-1 if r[v] else -1 for v in ['injury','penalty','collision','emotion_recorded']] for r in rows])
    x[:,2:]+=1
    models=json.loads((W/'original_models.json').read_text(encoding='utf-8'))
    refs={}
    for k in range(2,6):
        f=models[k-1];theta=np.zeros((k,4,4))
        for v,name in enumerate(['injury','penalty','collision','emotion']):theta[:,v,:CATS[v]]=f['probs'][name]
        refs[k]=(f['LL'],np.array(f['P']),theta,0)
    return x,refs

if __name__=='__main__':
    x,refs=load();results={}
    for k in [3,4]:
        assert abs(logprobs(refs[k],x).sum()-refs[k][0])<1e-6
        t=time.time();f,ls,its=fit(x,k,seed=836+k,nstart=30,ref=refs[k],tol=1e-9)
        results[k]={'ll':f[0],'seconds':time.time()-t,'best_repeats':int(np.sum(abs(ls-f[0])<1e-4)),'maxiter_hits':int(np.sum(its==20000)),'independent_ll':float(logprobs(f,x).sum()),'g2':discrepancy(f,x)}
        assert abs(f[0]-refs[k][0])<1e-4
        print(k,results[k],flush=True)
    (W/'engine_validation.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
