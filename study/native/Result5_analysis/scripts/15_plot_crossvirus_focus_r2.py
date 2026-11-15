"""Six-panel main figure and four-module supplement from unchanged Result5 outputs."""
from pathlib import Path
import ast,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,Patch
from matplotlib.text import Text
from matplotlib.colors import ListedColormap

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result5_analysis');F=OUT/'figures';TAG='crossvirus_focus_r2';R4=(OUT.parent/'Result4') if (OUT.parent/'Result4').exists() else OUT.parent/'Result4_analysis'
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
code=(R4/'scripts/20_plot_cellular_focus_r1.py').read_text();tree=ast.parse(code)
for name in ['new','ax','letter','clean','zero','heatmap','save']:
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name);exec(ast.get_source_segment(code,node))
original_save=save
def save(fig,stem):
    extra={a.get_gid():dict(xlim=list(a.get_xlim()),ylim=list(a.get_ylim()),scatter_offsets=[np.asarray(c.get_offsets()).tolist() for c in a.collections if type(c).__name__=='PathCollection'],scatter_facecolors=[np.asarray(c.get_facecolors()).tolist() for c in a.collections if type(c).__name__=='PathCollection'],annotations=[dict(text=t.get_text(),xy=list(t.xy),text_position=list(t.xyann)) for t in a.texts if hasattr(t,'xyann')]) for a in fig.axes}
    original_save(fig,stem);p=F/(stem+'_layout.json');d=json.loads(p.read_text());d['actual_artist_extra']=extra;dump(p,d)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'legend.title_fontsize':7,'axes.titlesize':7,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
TEAL='#287D8E';ORANGE='#D7813D';BLUE='#486F9A';PURPLE='#81649A';GRAY='#777777';RED='#BB4F4A'
E=pd.read_csv(OUT/'results/virus_comparison_effects.csv');S=pd.read_csv(OUT/'results/profile_scores.csv');M=pd.read_csv(OUT/'results/sensitivity_effects.csv');R=pd.read_csv(OUT/'results/locked_gene_replication.csv');R=R[R.replicated].sort_values('adjusted_P_discovery');D=pd.read_csv(OUT/'results/DC_time_matched_donor_values.csv');U=pd.read_csv(OUT/'results/participant_donor_values.csv');RN=pd.read_csv(OUT/'results/rank_induction_and_replication_effects.csv')
branches=['IFN-alpha','IFN-gamma','TNF / NF-kB','TLR signaling','Neutrophil degranulation','RTK signaling','Microtubule'];bn=['IFN-α','IFN-γ','TNF / NF-κB','TLR','Degranulation-annotated','RTK (down-directed)','Microtubule (down-directed)'];ct=['B_cell','T_cell','NK_cell','Monocyte'];cs=['B','T','NK','Mono.'];donorcolors=[TEAL,ORANGE,PURPLE]
payload={'version':TAG,'statistics_recomputed':False,'main':{},'supplement':{},'core_members':388,'main_panel_count':6,'supplement_modules':4}
def get(key,target='CRRM388'):return E[E.key.eq(key)&E.score.eq(target)].iloc[0]
def rec(panel,rows,supp=False):payload['supplement' if supp else 'main'][panel]=rows
def ef(a,r,y,color=TEAL,marker='o',gid=None,qcol='adjusted_P',ms=3):
    if not np.isfinite(r.estimate):return
    q=r.get(qcol,np.nan);item=a.errorbar(r.estimate,y,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt=marker,color=color,mfc=color if q<.05 else 'white',mec=color,mew=.75,ms=ms,elinewidth=.7,capsize=1.3)
    item.lines[0].set_gid(gid or str(r.get('key',''))+'|'+str(r.get('score','')))

fig=new(197,'Fig5');letter(fig,'A',4,10)
for source,x,w,groups,label in [('GSE38900_GPL6884',16,49,['Healthy','RSV','Influenza','Rhinovirus'],'Dallas (GSE38900)'),('GSE42026',85,29,['Healthy','RSV','Influenza'],'UK (GSE42026)')]:
    a=ax(fig,x,18,w,36,'5A_'+source);rows=[];random=np.random.default_rng(100)
    for i,group in enumerate(groups):
        g=S[S.source_label.eq(source)&S.analysis_group.eq(group)];v=g.CRRM388.to_numpy();c={'Healthy':GRAY,'RSV':TEAL,'Influenza':ORANGE,'Rhinovirus':PURPLE}[group]
        a.scatter(i+random.uniform(-.18,.18,len(v)),v,s=3.5,facecolors=c,edgecolors='none',alpha=.65,zorder=2)
        a.boxplot(v,positions=[i],widths=.56,showfliers=False,patch_artist=True,boxprops=dict(facecolor='none',edgecolor=c,lw=.65),whiskerprops=dict(color=c,lw=.6),capprops=dict(color=c,lw=.6),medianprops=dict(color=c,lw=1.))
        rows.append(dict(source=source,group=group,n=len(v),values=v.tolist(),median=float(np.median(v))))
    a.set(xticks=range(len(groups)),xticklabels=[{'Healthy':'HC','RSV':'RSV','Influenza':'Flu','Rhinovirus':'RV'}[g]+'\n'+str(len(S[S.source_label.eq(source)&S.analysis_group.eq(g)])) for g in groups],xlabel=label,ylim=(-.62,1.54),yticks=[-.5,0,.5,1,1.5]);a.axhline(0,color='#aaaaaa',lw=.5,zorder=0);clean(a)
    if source.startswith('GSE38900'):a.set_ylabel('CRRM score\n(within-study reference)')
    rec('5A_'+source,rows)

letter(fig,'B',119,10);a=ax(fig,160,18,42,36,'5B');rr=[]
for i,(key,color,marker) in enumerate([('Dallas_RSV_Influenza',BLUE,'o'),('UK_RSV_Influenza',TEAL,'s'),('Dallas_RSV_Rhinovirus',BLUE,'o')]):
    r=get(key);ef(a,r,i,color,marker);rr.append(r.to_dict())
zero(a);clean(a);a.set(yticks=range(3),yticklabels=['Dallas: RSV−flu','UK: RSV−flu','Dallas: RSV−RV'],ylim=(2.5,-.5),xlabel='Within-study CRRM difference',xlim=(-.65,.15),xticks=[-.6,-.3,0]);a.legend([Line2D([],[],ls='',marker=m,color=c,ms=3) for c,m in [(BLUE,'o'),(TEAL,'s')]],['Dallas','UK'],frameon=False,ncol=2,handlelength=.7,handletextpad=.3,columnspacing=.8,loc='lower right',bbox_to_anchor=(1.03,1.05),borderpad=0);rec('5B',rr)

letter(fig,'C',4,72);a=ax(fig,53,80,39,30,'5C');keys=['Dallas_RSV_Influenza','UK_RSV_Influenza','Dallas_RSV_Rhinovirus'];rr=[[get(k,t) for k in keys] for t in branches];mesh=heatmap(a,[[r.estimate for r in row] for row in rr],'RdBu_r',-1.1,1.1)
for i,row in enumerate(rr):
    for j,r in enumerate(row):
        if r.adjusted_P<.05:a.text(j,i,'*',ha='center',va='center',color='white' if abs(r.estimate)>.65 else 'black',fontsize=7)
a.add_patch(Rectangle((-.5,1.5),2,1,fill=False,lw=.7,edgecolor='#333333'))
a.set(yticks=range(7),yticklabels=bn,xticks=range(3),xticklabels=['Dallas\nvs flu','UK\nvs flu','Dallas\nvs RV']);a.get_yticklabels()[2].set_fontweight('bold');a.tick_params(length=0);cax=ax(fig,96,80,2.4,30,'5C_colorbar');cb=fig.colorbar(mesh,cax=cax,ticks=[-1,0,1]);cb.set_label('RSV−comparator',labelpad=2);rec('5C',dict(rows=branches,columns=keys,values=[r.to_dict() for row in rr for r in row]))

letter(fig,'D',119,72);a=ax(fig,155,80,47,30,'5D');rr=[]
for i,cell in enumerate(ct):
    for off,target,color,marker,size in [(-.20,'CRRM388','black','D',3),(0,'IFN-list overlap',PURPLE,'o',3),(.20,'Outside IFN lists',TEAL,'s',3.7)]:
        r=get('PBMC_'+cell+'_RSV_COVID',target);ef(a,r,i+off,color,marker,ms=size);rr.append(r.to_dict())
zero(a);clean(a);a.set(yticks=range(4),yticklabels=['B (19/17)','T (19/17)','NK (19/17)','Mono. (19/15)'],ylim=(3.6,-.6),xlabel='PBMC RSV−COVID',xticks=[-.6,-.3,0]);a.legend([Line2D([],[],ls='',marker=m,color=c,ms=3) for c,m in [(TEAL,'s'),('black','D'),(PURPLE,'o')]],['Outside IFN lists','Total','IFN-list overlap'],frameon=False,ncol=3,handlelength=.7,handletextpad=.25,columnspacing=.5,loc='lower right',bbox_to_anchor=(1.03,1.06),borderpad=0);rec('5D',rr)

letter(fig,'E',4,128);a=ax(fig,18,144,84,36,'5E');gg=pd.read_csv(OUT/'results/gene_induction_and_difference.csv');xx=gg[gg.key.eq('Dallas_Influenza_Healthy')].set_index('human_ensembl');yy=gg[gg.key.eq('Dallas_RSV_Healthy')].set_index('human_ensembl');diff=pd.read_csv(OUT/'results/discovery_gene_differences.csv');diff=diff[diff.key.eq('Dallas_RSV_Influenza')].set_index('human_ensembl');ids=xx.index[xx.measured&yy.measured];v=diff.loc[ids];cols=np.where(v.adjusted_P.lt(.05),BLUE,'#b7b7b7');a.scatter(xx.loc[ids,'aligned_gene_effect'],yy.loc[ids,'aligned_gene_effect'],s=7,c=cols,edgecolors='none',alpha=.85)
rid=R.human_ensembl;a.scatter(xx.loc[rid,'aligned_gene_effect'],yy.loc[rid,'aligned_gene_effect'],s=18,facecolors='none',edgecolors=RED,lw=.85)
extent=np.r_[xx.loc[ids,'aligned_gene_effect'],yy.loc[ids,'aligned_gene_effect']];lims=(min(-1.3,float(np.nanmin(extent))-.2),max(2.,float(np.nanmax(extent))+.2));a.plot(lims,lims,color=GRAY,lw=.5,ls='--');a.axhline(0,color=GRAY,lw=.4);a.axvline(0,color=GRAY,lw=.4);clean(a);a.set(xlim=lims,ylim=lims,xlabel='Influenza−HC signed gene change',ylabel='RSV−HC signed gene change');a.locator_params(axis='both',nbins=4)
positions={'RRBP1':(1.7,-.92),'CFLAR':(-.2,-1.15),'ZCCHC2':(1.30,1.65),'LYN':(.30,1.05)};annotations=[]
for r in R.itertuples():
    xy=(float(xx.loc[r.human_ensembl,'aligned_gene_effect']),float(yy.loc[r.human_ensembl,'aligned_gene_effect']));textpos=positions[r.display_symbol];a.annotate(r.display_symbol,xy=xy,xytext=textpos,fontsize=7,va='center',ha='left',color='#222222',bbox=dict(boxstyle='square,pad=.13',facecolor='white',edgecolor='none',alpha=.94),arrowprops=dict(arrowstyle='-',lw=.55,color='#555555',shrinkA=1.5,shrinkB=3.0));annotations.append(dict(gene=r.display_symbol,point=xy,label=textpos))
a.legend([Line2D([],[],ls='',marker='o',color=c,ms=3,mfc='none' if c==RED else c) for c in ['#b7b7b7',BLUE,RED]],['Other members','Dallas discoveries (97)','UK-replicated (4)'],frameon=False,loc='lower left',bbox_to_anchor=(0,1.03),ncol=2,handlelength=.7,handletextpad=.25,columnspacing=.6,labelspacing=.15,borderpad=.1);rec('5E',dict(genes=ids.tolist(),x=xx.loc[ids,'aligned_gene_effect'].tolist(),y=yy.loc[ids,'aligned_gene_effect'].tolist(),direct_adjusted_P=v.adjusted_P.tolist(),annotations=annotations,gene_family=776,replication_family=388))

letter(fig,'F',119,128);a=ax(fig,150,144,52,36,'5F');rr=[]
for i,r in enumerate(R.itertuples()):
    for off,source,color,marker in [(-.13,'discovery',BLUE,'o'),(.13,'replication',TEAL,'s')]:
        if source=='discovery':z=diff.loc[r.human_ensembl];q=pd.Series(dict(estimate=r.aligned_gene_effect_discovery,lower=z.lower,upper=z.upper,adjusted_P=r.adjusted_P_discovery))
        else:q=pd.Series(dict(estimate=r.aligned_gene_effect_replication,lower=r.lower,upper=r.upper,adjusted_P=r.adjusted_P_replication))
        ef(a,q,i+off,color,marker,gid=r.display_symbol+'|'+source);rr.append(dict(gene=r.display_symbol,source=source,**q.to_dict()))
zero(a);clean(a);a.set(yticks=range(4),yticklabels=R.display_symbol.tolist(),ylim=(3.5,-.6),xlabel='RSV−flu signed gene difference\n(within-study reference)',xlim=(-2.8,.15),xticks=[-2,-1,0]);a.legend([Line2D([],[],ls='',marker=m,color=c,ms=3) for c,m in [(BLUE,'o'),(TEAL,'s')]],['Dallas','UK'],frameon=False,ncol=2,loc='lower right',bbox_to_anchor=(1.02,1.06),handlelength=.7,handletextpad=.3,borderpad=0);rec('5F',rr);save(fig,'Figure5_'+TAG)

# Supplement: full-width sensitivity, paired PBMC modules, then one DC module.
fig=new(204,'FigS6');letter(fig,'A',4,10);a=ax(fig,57,18,145,42,'S6A')
keys=['Dallas_RSV_Influenza','Dallas_RSV_Rhinovirus']+['PBMC_'+cell+'_RSV_COVID' for cell in ct]+['DC_RSV_Flu','UK_RSV_Influenza'];methods=['Primary','Signed transcriptome rank','Age/sex HC3','Common measured members','Minimum15cells','Minimum50cells','Detection-supported members'];vals=np.full((8,7),np.nan);rr=[]
for i,key in enumerate(keys):
    for j,method in enumerate(methods):
        if key=='UK_RSV_Influenza':h=E[E.key.eq(key)&E.score.eq('CRRM388')] if method=='Primary' else RN[RN.key.eq(key)&RN.score.eq('CRRM388')] if method=='Signed transcriptome rank' else pd.DataFrame()
        else:h=E[E.key.eq(key)&E.score.eq('CRRM388')] if method=='Primary' else M[M.key.eq(key)&M.score.eq('CRRM388')&M.method.eq(method)]
        if len(h) and np.isfinite(h.iloc[0].estimate):
            r=h.iloc[0];q=r.get('adjusted_P',r.get('adjusted_P_BH7',np.nan));vals[i,j]=np.sign(r.estimate)*(1 if q<.05 else .35);rr.append(dict(key=key,method=method,estimate=r.estimate,adjusted_P=q,scope='UK replication BH3' if key.startswith('UK_') else 'Primary direct family BH7'))
mesh=heatmap(a,vals,ListedColormap([BLUE,'#c7d6e5','white','#e7cfb7',ORANGE]),-1,1)
for i in range(8):
    for j in range(7):
        if not np.isfinite(vals[i,j]):a.text(j,i,'—',ha='center',va='center',fontsize=7,color='#777777')
        elif abs(vals[i,j])==1:a.text(j,i,'*',ha='center',va='center',fontsize=7,color='white')
a.axhline(6.5,color='#333333',lw=.85);a.set(yticks=range(8),yticklabels=['Dallas: RSV−flu','Dallas: RSV−RV','PBMC: B','PBMC: T','PBMC: NK','PBMC: Mono.','DC: RSV−flu','UK: RSV−flu (replication)'],xticks=range(7),xticklabels=['Primary','Rank','Age/sex','313 common\nmembers','≥15 cells','≥50 cells','Detection-\nsupported']);a.tick_params(length=0)
a.legend([Patch(color=BLUE),Patch(color=ORANGE),Line2D([],[],ls='',marker='*',color='#444444',ms=4),Line2D([],[],ls='-',color=GRAY,lw=.7)],['Lower in RSV','Higher in RSV','Adjusted P<0.05','Unavailable'],ncol=4,frameon=False,loc='lower right',bbox_to_anchor=(1.01,1.05),handlelength=.8,handletextpad=.3,columnspacing=.9,borderpad=0);rec('S6A',dict(rows=rr,matrix=vals.tolist(),keys=keys,methods=methods),True)

letter(fig,'B',4,74);rr=[]
for j,(method,x) in enumerate([('Primary',23),('Signed transcriptome rank',71)]):
    a=ax(fig,x,82,33,32,'S6B_'+str(j))
    for i,cell in enumerate(ct):
        for off,target,color,marker in [(-.13,'Outside IFN lists',TEAL,'o'),(.13,'Neutrophil degranulation',ORANGE,'s')]:
            key='PBMC_'+cell+'_RSV_COVID';r=get(key,target) if method=='Primary' else M[M.key.eq(key)&M.score.eq(target)&M.method.eq(method)].iloc[0];qcol='adjusted_P' if method=='Primary' else 'adjusted_P_BH7';ef(a,r,i+off,color,marker,qcol=qcol);rr.append({**r.to_dict(),'display_method':method})
    zero(a);clean(a);a.set(yticks=range(4),yticklabels=cs if j==0 else [],ylim=(3.6,-.6),xlabel=('Gene-z' if j==0 else 'Rank')+'\nRSV−COVID');a.locator_params(axis='x',nbins=3)
a.legend([Line2D([],[],ls='',marker=m,color=c,ms=3) for c,m in [(TEAL,'o'),(ORANGE,'s')]],['Outside IFN lists','Degranulation annotation'],frameon=False,ncol=2,loc='upper right',bbox_to_anchor=(1.02,-.50),borderpad=0,handlelength=.7,handletextpad=.3,columnspacing=.6);rec('S6B',rr,True)

letter(fig,'C',119,74);comp=pd.read_csv(OUT/'results/PBMC_composition_effects.csv');rankcomp=pd.read_csv(OUT/'results/PBMC_rank_composition_effects.csv');rr=[]
for j,(table,x) in enumerate([(comp,154),(rankcomp,184)]):
    a=ax(fig,x,82,19,32,'S6C_'+str(j))
    for i,target in enumerate(['CRRM388','IFN-list overlap','Outside IFN lists']):
        for off,method,color,marker in [(-.15,'Observed composition',GRAY,'o'),(.15,'Fixed healthy composition',TEAL,'s')]:
            r=table[table.score.eq(target)&table.method.eq(method)].iloc[0];r['CI_exclusion_flag']=0 if r.lower*r.upper>0 else 1;ef(a,r,i+off,color,marker,qcol='CI_exclusion_flag');rr.append({**r.to_dict(),'scoring':'Gene-z' if j==0 else 'Rank'})
    zero(a);clean(a);a.set(yticks=range(3),yticklabels=['Total','IFN overlap','Outside IFN'] if j==0 else [],ylim=(2.6,-.6),xlabel=('Gene-z' if j==0 else 'Rank')+'\nRSV−COVID');a.locator_params(axis='x',nbins=2)
a.legend([Line2D([],[],ls='',marker=m,color=c,ms=3) for c,m in [(GRAY,'o'),(TEAL,'s')]],['Observed composition','Fixed healthy composition'],frameon=False,ncol=1,loc='upper right',bbox_to_anchor=(1.04,-.50),borderpad=0,handlelength=.7,handletextpad=.3,labelspacing=.1);rec('S6C',rr,True)

letter(fig,'D',4,143);rr=[]
for j,target in enumerate(['CRRM388','IFN-list overlap','Outside IFN lists']):
    a=ax(fig,16+35*j,155,25,32,'S6D_time_'+str(j));g=D[D.cell.eq('Dendritic_cell')&D.case.eq('RSV')&D.control.eq('Flu')&D.score.eq(target)]
    for k,(donor,h) in enumerate(g.groupby('donor')):
        h=h.set_index('time').loc[['4h','8h','18h']];line=a.plot(range(3),h.difference,color=donorcolors[k],lw=.75,marker='o',ms=2.3)[0];line.set_gid(target+'|'+donor);rr+=h.reset_index().to_dict('records')
    a.axhline(0,color=GRAY,lw=.5);clean(a);a.set(xticks=range(3),xticklabels=['4','8','18'],xlabel={'CRRM388':'Total CRRM\nTime (h)','IFN-list overlap':'IFN-list overlap\nTime (h)','Outside IFN lists':'Outside IFN lists\nTime (h)'}[target]);a.locator_params(axis='y',nbins=3)
    if j==0:
        a.set_ylabel('RSV−flu difference');a.legend([Line2D([],[],color=c,marker='o',ms=2.3,lw=.75) for c in donorcolors],['Donor 1','Donor 2','Donor 3'],frameon=False,ncol=3,loc='lower left',bbox_to_anchor=(-.04,1.08),borderpad=0,handlelength=.8,handletextpad=.25,columnspacing=.7)
rec('S6D_time',rr,True)
a=ax(fig,150,155,52,32,'S6D_mean');rr=[]
for i,key in enumerate(['DC_RSV_Flu','DC_RSV_HK']):
    ca=U[U.key.eq(key)&U.score.eq('CRRM388')&U.role.eq('case')].set_index('analysis_unit').value;co=U[U.key.eq(key)&U.score.eq('CRRM388')&U.role.eq('control')].set_index('analysis_unit').value;delta=ca-co
    for j,(u,v) in enumerate(delta.items()):a.plot(v,i+[-.15,0,.15][j],'o',color=donorcolors[j],ms=3)
    r=get(key);ef(a,r,i+.30,'black','D');rr.append(dict(key=key,donors=delta.to_dict(),model=r.to_dict()))
zero(a);clean(a);a.set(yticks=[0,1],yticklabels=['RSV−flu','RSV−HK'],ylim=(1.8,-.5),xlabel='Time-averaged total CRRM\nRSV−comparator',xlim=(-.6,.4),xticks=[-.5,0,.3]);a.legend([Line2D([],[],color='black',marker='D',ls='-',lw=.7,ms=3)],['Mean (95% CI)'],frameon=False,loc='lower right',bbox_to_anchor=(1.02,1.08),handlelength=1.2,borderpad=0);rec('S6D_mean',rr,True);save(fig,'FigureS6_'+TAG)
dump(OUT/'audit'/f'figure_payload_{TAG}.json',payload)
