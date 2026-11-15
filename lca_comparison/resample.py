from engine import *
from concurrent.futures import ThreadPoolExecutor,as_completed
import csv,hashlib
x,refs=load();N=len(x)

def robust_fit(z,k,seed,nstart=30,use_ref=True):
    f,ls,its=fit(z,k,seed,nstart,refs[k] if use_ref else None)
    repeats=int(np.sum(abs(ls-f[0])<1e-4))
    if repeats<2:
        ff,ll2,it2=fit(z,k,seed+10000000,30,f)
        if ff[0]>f[0]:f=ff
        ls=np.r_[ls,ll2];its=np.r_[its,it2]
        repeats=int(np.sum(abs(ls-f[0])<1e-4))
    if f[3]>=20000:
        y,w=pack(z);f=em(y,w,f[1].copy(),f[2].copy(),maxiter=100000,tol=1e-9)
    return f,{'best_repeats':repeats,'starts':len(ls),'best_iterations':int(f[3]),'start_iteration_limits':int(np.sum(its>=20000))}

def param_task(i):
    rng=np.random.default_rng(31106000+i);z=simulate(refs[3],x,rng)
    a,da=robust_fit(z,3,41106000+i)
    b,db=robust_fit(z,4,51106000+i)
    assert b[0]>=a[0]-1e-4
    return {'replicate':i,'lr':2*(b[0]-a[0]),'g3':discrepancy(a,z),'diagnostics3':da,'diagnostics4':db}

def gof4_task(i):
    rng=np.random.default_rng(61106000+i);z=simulate(refs[4],x,rng)
    a,da=robust_fit(z,4,71106000+i)
    return {'replicate':i,'g4':discrepancy(a,z),'diagnostics':da}

def nonparam_task(i):
    rng=np.random.default_rng(81106000+i);z=x[rng.integers(0,N,N)]
    a,da=robust_fit(z,3,91106000+i);b,db=robust_fit(z,4,101106000+i)
    a,ad=aligned(a,refs[3]);b,bd=aligned(b,refs[4])
    return {'replicate':i,'lr':2*(b[0]-a[0]),'bic4_minus3':-2*(b[0]-a[0])+8*np.log(N),'a3':a[1].tolist(),'a4':b[1].tolist(),'theta3':a[2].tolist(),'theta4':b[2].tolist(),'distance3':ad,'distance4':bd,'diagnostics3':da,'diagnostics4':db}

def cv_task(job):
    repeat,fold,indices=job;test=indices[fold::10];train=np.ones(N,bool);train[test]=False
    a,da=robust_fit(x[train],3,121106000+repeat*10+fold,use_ref=False)
    b,db=robust_fit(x[train],4,131106000+repeat*10+fold,use_ref=False)
    la=logprobs(a,x[test]);lb=logprobs(b,x[test])
    return {'repeat':repeat,'fold':fold,'test_n':len(test),'ll3':float(la.sum()),'ll4':float(lb.sum()),'delta4_minus3':float((lb-la).sum()),'test_indices':test.tolist(),'pointwise_delta':(lb-la).tolist(),'diagnostics3':da,'diagnostics4':db}

def run(name,jobs,fn):
    p=W/(name+'.json');done=json.loads(p.read_text()) if p.exists() else []
    seen={(r['repeat'],r['fold']) if name=='cv' else r['replicate'] for r in done}
    jobs=[j for j in jobs if ((j[0],j[1]) if name=='cv' else j) not in seen]
    print(name,'remaining',len(jobs),flush=True);t=time.time()
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(fn,j) for j in jobs]
        for fut in as_completed(futures):
            done.append(fut.result())
            if len(done)%25==0:
                p.write_text(json.dumps(done),encoding='utf-8')
                print(name,len(done),'elapsed',round(time.time()-t,1),flush=True)
    done.sort(key=lambda r:(r['repeat'],r['fold']) if name=='cv' else r['replicate'])
    p.write_text(json.dumps(done),encoding='utf-8')

if __name__=='__main__':
    # Compile before using worker threads.
    fit(x,3,1,1)
    run('parametric_bootstrap',range(499),param_task)
    run('nonparametric_bootstrap',range(300),nonparam_task)
    jobs=[]
    for r in range(5):
        inds=np.random.default_rng(141106000+r).permutation(N)
        jobs.extend((r,f,inds) for f in range(10))
    run('cv',jobs,cv_task)
    run('parametric_gof4',range(199),gof4_task)
    print('ALL RESAMPLING FINISHED',flush=True)
