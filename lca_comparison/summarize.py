from engine import *
import csv,hashlib,platform
x,refs=load();result={'N':len(x),'missing_counts':(x<0).sum(axis=0).tolist()}
def read(n):return json.loads((W/(n+'.json')).read_text())
def pboot(vals,obs):
    a=np.array(vals);exc=int(np.sum(a>=obs));B=len(a)
    return {'observed':float(obs),'B':B,'exceedances':exc,'p_plus_one':(exc+1)/(B+1),'null_quantiles':np.quantile(a,[.025,.5,.95,.975]).tolist()}
a=read('parametric_bootstrap');g=read('parametric_gof4');b=read('nonparametric_bootstrap');cv=read('cv')
result['BLRT_3_vs_4']=pboot([r['lr'] for r in a],2*(refs[4][0]-refs[3][0]))
result['missing_pattern_gof3']=pboot([r['g3'] for r in a],discrepancy(refs[3],x))
result['missing_pattern_gof4']=pboot([r['g4'] for r in g],discrepancy(refs[4],x))
result['missing_pattern_gof_note']='Conditional simulation fixes observed missingness masks and assumes simulated indicators independent of those masks. These checks can reflect the missingness mechanism as well as model misfit. Complete-case refits are separately calibrated.'
profiles=[]
for k in [3,4]:
    f=refs[k];order=np.lexsort((f[2][:,2,1],f[2][:,1,1]));pi=np.array([r[f'a{k}'] for r in b]);th=np.array([r[f'theta{k}'] for r in b]);post=posterior(f,x);labels=post.argmax(axis=1)
    result[f'bootstrap{k}']={'B':len(b),'profile_squared_distance_quantiles':np.quantile([r[f'distance{k}'] for r in b],[.025,.5,.975]).tolist(),'classes':[]}
    for c,j in enumerate(order):
        rec={'class':f'C{c+1}','original_model_index':int(j),'share':float(f[1][j]),'share_bootstrap_percentiles':np.quantile(pi[:,j],[.025,.5,.975]).tolist(),'modal_n':int(np.sum(labels==j)),'mean_max_posterior':float(post[labels==j,j].mean()),'injury':f[2][j,0,:4].tolist(),'penalty':f[2][j,1,:3].tolist(),'collision':float(f[2][j,2,1]),'emotion':float(f[2][j,3,1]),'no_injury_bootstrap_percentiles':np.quantile(th[:,j,0,0],[.025,.5,.975]).tolist(),'no_injury_over_90_bootstrap_fraction':float(np.mean(th[:,j,0,0]>.9))}
        result[f'bootstrap{k}']['classes'].append(rec)
        for v,name in enumerate(['injury','penalty','collision','emotion']):
            for cat in range(CATS[v]):profiles.append({'K':k,'class':f'C{c+1}','indicator':name,'category':cat+1,'estimate':f[2][j,v,cat],'bootstrap_p025':np.quantile(th[:,j,v,cat],.025),'bootstrap_p50':np.quantile(th[:,j,v,cat],.5),'bootstrap_p975':np.quantile(th[:,j,v,cat],.975)})
result['bootstrap_BIC4_selected']=sum(r['bic4_minus3']<0 for r in b)
result['bootstrap_BIC_difference_quantiles']=np.quantile([r['bic4_minus3'] for r in b],[.025,.5,.975]).tolist()
result['cv']={'repeats':5,'folds':10,'delta_definition':'held-out LL4 minus held-out LL3; positive favors four','repeat_totals':[sum(r['delta4_minus3'] for r in cv if r['repeat']==i) for i in range(5)],'positive_folds':sum(r['delta4_minus3']>0 for r in cv),'note':'Folds and repeats share observations; no independent-sample t-test is applied.'}
result['cv']['average_delta_per_case']=sum(result['cv']['repeat_totals'])/(5*len(x))
cc=read('complete_case_fits');cboot=read('complete_case_bootstrap')
result['complete_case']={'N':cc['3']['N'],'BIC3':cc['3']['BIC'],'BIC4':cc['4']['BIC'],'AIC3':cc['3']['AIC'],'AIC4':cc['4']['AIC'],'BLRT':pboot([r['lr'] for r in cboot if r['null_k']==3],2*(cc['4']['LL']-cc['3']['LL'])),'gof3':pboot([r['g'] for r in cboot if r['null_k']==3],cc['3']['g2']),'gof4':pboot([r['g'] for r in cboot if r['null_k']==4],cc['4']['g2'])}
result['diagnostics']={}
for name,rows in [('parametric',a),('nonparametric',b),('cv',cv),('gof4',g)]:
    ds=[v for r in rows for key,v in r.items() if key.startswith('diagnostics')]
    result['diagnostics'][name]={'fits':len(ds),'best_at_100000_iteration_limit':sum(d['best_iterations']>=100000 for d in ds),'best_replicated_fewer_than_twice':sum(d['best_repeats']<2 for d in ds),'max_best_iterations':max(d['best_iterations'] for d in ds)}
en=read('simulation_engine_crosschecks');rr=list(csv.DictReader((W/'simulation_poLCA_crosschecks.csv').open()))
errors=[abs(e['LL']-float(next(r['LL'] for r in rr if int(r['replicate'])==e['replicate'] and int(r['K'])==e['K']))) for e in en]
result['simulation_engine_poLCA_max_abs_LL_difference']=max(errors)
assert max(errors)<1e-4,errors
with (W/'profiles_and_bootstrap_ranges.csv').open('w',newline='',encoding='utf-8-sig') as f:
    wr=csv.DictWriter(f,fieldnames=list(profiles[0]));wr.writeheader();wr.writerows(profiles)
(W/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
