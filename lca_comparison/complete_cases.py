from engine import *
from concurrent.futures import ThreadPoolExecutor,as_completed
x,old=load();x=x[np.all(x>=0,axis=1)];refs={};fitinfo={}
for k in [3,4]:
    f,ls,its=fit(x,k,151106000+k,100,old[k],tol=1e-9)
    refs[k]=f;fitinfo[k]={'N':len(x),'LL':f[0],'AIC':-2*f[0]+2*(8*k-1),'BIC':-2*f[0]+(8*k-1)*np.log(len(x)),'g2':discrepancy(f,x),'P':f[1].tolist(),'theta':f[2].tolist(),'best_repeats':int(np.sum(abs(ls-f[0])<1e-4))}
print('COMPLETE CASE FITS',fitinfo,flush=True)
(W/'complete_case_fits.json').write_text(json.dumps(fitinfo),encoding='utf-8')
def work(job):
    k,i=job;rng=np.random.default_rng(161106000+k*1000+i);z=simulate(refs[k],x,rng)
    a,_,_=fit(z,k,171106000+k*1000+i,30,refs[k]);res={'null_k':k,'replicate':i,'g':discrepancy(a,z)}
    if k==3:
        b,_,_=fit(z,4,181106000+i,30,refs[4]);res['lr']=2*(b[0]-a[0])
        assert res['lr']>=-1e-4
    return res
p=W/'complete_case_bootstrap.json';done=json.loads(p.read_text()) if p.exists() else [];seen={(r['null_k'],r['replicate']) for r in done}
jobs=[(k,i) for k in [3,4] for i in range(199) if (k,i) not in seen]
with ThreadPoolExecutor(max_workers=3) as pool:
    for f in as_completed([pool.submit(work,j) for j in jobs]):
        done.append(f.result())
        if len(done)%50==0:p.write_text(json.dumps(done));print('CC BOOT',len(done),flush=True)
p.write_text(json.dumps(done));print('COMPLETE CASE DONE',flush=True)
