from engine import *
x,refs=load()
rows=[]
for i in [0,17,44]:
    z=simulate(refs[3],x,np.random.default_rng(31106000+i))
    np.savetxt(W/f'simulation_{i}.csv',np.where(z<0,np.nan,z+1),delimiter=',',header='injury,penalty,collision,emotion',comments='')
    for k in [3,4]:
        f,_,_=fit(z,k,202611090+i+k,30,refs[k])
        rows.append({'replicate':i,'K':k,'LL':f[0]})
(W/'simulation_engine_crosschecks.json').write_text(json.dumps(rows),encoding='utf-8')
