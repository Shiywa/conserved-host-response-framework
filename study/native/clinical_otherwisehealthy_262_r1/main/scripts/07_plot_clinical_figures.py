from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1/main')
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.55,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':0,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
UP='#BA603A';DOWN='#387C98';INK='#244B60';LIGHT='#A6B2B9';COL={'Healthy':'#8D969E','Mild':'#558F75','Moderate':'#D99B3B','Severe':'#B74742'};DS=['GSE246622','GSE69606','GSE105450'];branches=['IFN-alpha','IFN-gamma','TNF / NF-kB','TLR signaling','Neutrophil degranulation','RTK signaling','Microtubule'];blabel=['IFN-α','IFN-γ','TNF/NF-κB','TLR signaling','Neutrophil\ndegranulation','RTK signaling','Microtubule'];rng=np.random.default_rng(20261120);audit=[]
read=lambda f:pd.read_csv(OUT/'tables'/f)
scores=read('07_sample_scores.csv');effects=read('09_total_clinical_effects.csv').set_index('dataset').loc[DS];fun=read('10_functional_clinical_effects.csv');comp=read('11_directional_contributions.csv');adj=read('17_covariate_adjusted_effects.csv');rec=read('18_paired_recovery_values.csv');recfx=read('19_paired_recovery_effects.csv');gene=read('27_infection_severity_gene_map.csv');relation=json.loads((OUT/'gene_relationship_summary.json').read_text())
def axmm(f,x,y,w,h):
 W,H=np.asarray(f.get_size_inches())*25.4;return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def textmm(f,x,y,s,**kw):
 W,H=np.asarray(f.get_size_inches())*25.4;return f.text(x/W,1-y/H,s,**kw)
def letter(f,x,y,s):textmm(f,x,y,s,fontsize=12,fontweight='bold',va='top')
def clean(ax,zero=None):
 ax.spines[['top','right']].set_visible(False)
 if zero=='x':ax.axvline(0,color='#ADB4B8',lw=.6,zorder=0)
 if zero=='y':ax.axhline(0,color='#ADB4B8',lw=.6,zorder=0)
def err(ax,r,y,color=INK,offset=0,marker='o',mfc=None):
 if not hasattr(ax,'_effect_intervals'):ax._effect_intervals=[]
 ax._effect_intervals.append((float(r.CI95_low),float(r.CI95_high)))
 ax.errorbar(r.effect,y+offset,xerr=[[r.effect-r.CI95_low],[r.CI95_high-r.effect]],fmt=marker,color=color,mfc=color if mfc is None else mfc,ms=3.3,lw=.85,capsize=1.5)
def cb(f,im,x,y,h,label,ticks):
 c=f.colorbar(im,cax=axmm(f,x,y,2,h));c.set_ticks(ticks);c.ax.tick_params(labelsize=7,length=1.5,pad=1);c.outline.set_linewidth(.4);c.set_label(label,labelpad=3)
def finish(f,stem,tag):
 # Extend an existing range only if an actual confidence interval would be clipped.
 for a in f.axes:
  if hasattr(a,'_effect_intervals'):
   oldlo,oldhi=a.get_xlim();lo=min(v[0] for v in a._effect_intervals);hi=max(v[1] for v in a._effect_intervals)
   if lo<oldlo or hi>oldhi:
    span=max(oldhi-oldlo,hi-lo);a.set_xlim(min(oldlo,lo-.04*span),max(oldhi,hi+.04*span))
    from matplotlib.ticker import MaxNLocator
    a.xaxis.set_major_locator(MaxNLocator(nbins=4))
 textmm(f,207,3,tag,fontsize=12,ha='right',va='top');f.canvas.draw();rend=f.canvas.get_renderer();W,H=f.bbox.width,f.bbox.height;outside=[];texts=[t for t in f.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
 for t in texts:
  b=t.get_window_extent(rend)
  if not np.isfinite(b.extents).all() or b.x0<-.2 or b.y0<-.2 or b.x1>W+.2 or b.y1>H+.2:outside.append(dict(text=t.get_text(),bbox_px=[float(v) for v in b.extents]))
 for ext in ['pdf','svg','png']:f.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=600 if ext=='png' else 300)
 audit.append(dict(figure=stem,width_mm=f.get_figwidth()*25.4,height_mm=f.get_figheight()*25.4,text_sizes_pt=sorted(set(t.get_fontsize() for t in texts)),all_text_inside_canvas=not outside,outside_text=outside,identifier=tag,identifier_pt=12,ordinary_pt=7,panel_letter_pt=12,single_compact_layout=True));plt.close(f)

f=plt.figure(figsize=(210/25.4,191/25.4))
letter(f,2.5,12,'A');ax=axmm(f,13,19,82,40);m=scores[scores.dataset.eq('GSE246622')];levels=['Healthy','Mild','Moderate','Severe'];vals=[m.loc[m.group.eq(g),'score_Total'].to_numpy() for g in levels]
b=ax.boxplot(vals,positions=np.arange(4),widths=.5,patch_artist=True,showfliers=False,medianprops=dict(color=INK,lw=.9),boxprops=dict(lw=.7),whiskerprops=dict(lw=.7),capprops=dict(lw=.7))
for i,g in enumerate(levels):
 b['boxes'][i].set_facecolor(COL[g]);b['boxes'][i].set_alpha(.3);ax.scatter(i+rng.uniform(-.19,.19,len(vals[i])),vals[i],s=3.5,c=COL[g],alpha=.65,edgecolors='none',zorder=2)
ax.set(xticks=np.arange(4),xticklabels=[f'{g}\nn = {len(v)}' for g,v in zip(levels,vals)],ylabel='CRRM score',xlabel='Clinical group');clean(ax,'y');textmm(f,15,16,'GSE246622; acute RSV and healthy reference',va='bottom')

letter(f,105,12,'B');ax=axmm(f,153,19,48,40)
for i,ds in enumerate(DS):err(ax,effects.loc[ds],i,color=INK if i<2 else LIGHT,marker='o' if i<2 else 's')
ax.axhline(1.5,color=LIGHT,lw=.6,ls=(0,(2,2)));ax.set(yticks=np.arange(3),yticklabels=['Whole blood: severe − mild\nGSE246622','PBMC: severe − mild\nGSE69606','Whole blood: inpatient − outpatient\nGSE105450'],ylim=(2.55,-.55),xlim=(-.25,.43),xticks=[-.2,0,.2,.4],xlabel='Mean score change (95% CI)');ax.tick_params(axis='y',pad=4);clean(ax,'x')

letter(f,2.5,76,'C');ax=axmm(f,36,82,59,35)
for metric,color,off in [('Up contribution',UP,-.14),('Down contribution',DOWN,.14)]:
 for i,ds in enumerate(DS):r=comp[comp.dataset.eq(ds)&comp.component.eq(metric)].iloc[0];err(ax,r,i,color=color,offset=off)
ax.axhline(1.5,color=LIGHT,lw=.6,ls=(0,(2,2)));ax.set(yticks=np.arange(3),yticklabels=DS,ylim=(2.5,-.5),xlim=(-.18,.46),xticks=[-.1,0,.2,.4],xlabel='Contribution to score change');clean(ax,'x');ax.legend(handles=[Line2D([],[],color=UP,marker='o',ms=3,label='Up'),Line2D([],[],color=DOWN,marker='o',ms=3,label='Down')],loc='lower right',bbox_to_anchor=(1,1.03),ncol=2,frameon=False,handlelength=.8,columnspacing=1)

letter(f,105,76,'D');ax=axmm(f,138,82,56,35);mat=fun.pivot(index='branch',columns='dataset',values='effect').reindex(index=branches,columns=DS);q=fun.pivot(index='branch',columns='dataset',values='BH_adjusted_P').reindex(index=branches,columns=DS);assert mat.notna().all().all();im=ax.imshow(mat,aspect='auto',interpolation='nearest',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1));ax.axvline(1.5,color='white',lw=1.2)
for i in range(7):
 for j in range(3):
  if q.iloc[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=7,color='white' if abs(mat.iloc[i,j])>.6 else 'black')
ax.set(yticks=np.arange(7),yticklabels=blabel,xticks=np.arange(3),xticklabels=['Whole blood\n246622','PBMC\n69606','Hospitalization\n105450']);ax.tick_params(axis='x',length=0,pad=3);ax.spines[:].set_visible(False);cb(f,im,197,82,35,'Branch score change',[-1,0,1]);textmm(f,135,129,'* BH-adjusted P < 0.05',va='top')

letter(f,2.5,138,'E');ax=axmm(f,13,145,54,35)
for direction,color in [('up',UP),('down',DOWN)]:
 v=gene[gene.frozen_direction.eq(direction)];ax.scatter(v.Mild_minus_Healthy_adjusted_log_effect,v.Severe_minus_Mild_adjusted_log_effect,s=4.5,c=color,alpha=.7,edgecolors='none')
ax.axvline(0,color=LIGHT,lw=.6,zorder=0);ax.axhline(0,color=LIGHT,lw=.6,zorder=0);ax.set(xlabel='Mild − healthy (log effect)',ylabel='Severe − mild (log effect)');clean(ax);ax.text(.03,.96,f'ρ = {relation["Spearman_rho"]:.2f}',transform=ax.transAxes,va='top')
label_positions={'ICAM1':(.49,.72),'TNFAIP3':(.50,.40),'CXCL10':(.50,.20),'ISG15':(.73,.56)}
for name,pos in label_positions.items():
 r=gene[gene.display_symbol.eq(name)].iloc[0];ax.annotate(name,(r.Mild_minus_Healthy_adjusted_log_effect,r.Severe_minus_Mild_adjusted_log_effect),xytext=pos,textcoords='axes fraction',fontsize=7,arrowprops=dict(arrowstyle='-',lw=.4,color=LIGHT))

letter(f,74,138,'F');ax=axmm(f,98,145,40,35)
for model,color,off in [('unadjusted',LIGHT,-.12),('adjusted',INK,.12)]:
 for i,ds in enumerate(['GSE246622','GSE105450']):r=adj[adj.dataset.eq(ds)&adj.model.eq(model)].iloc[0];err(ax,r,i,color=color,offset=off)
ax.set(yticks=[0,1],yticklabels=['GSE246622','GSE105450'],ylim=(1.5,-.5),xlim=(-.2,.27),xticks=[-.2,0,.2],xlabel='Mean score change');clean(ax,'x');ax.legend(handles=[Line2D([],[],color=LIGHT,marker='o',ms=3,label='Unadj.'),Line2D([],[],color=INK,marker='o',ms=3,label='Adjusted')],loc='lower right',bbox_to_anchor=(1.02,1.02),ncol=2,frameon=False,handlelength=.7,columnspacing=.6)

letter(f,145,138,'G');ax=axmm(f,158,145,45,35)
for _,r in rec.iterrows():ax.plot([0,1],[r.acute_Total,r.recovery_Total],'-o',color=COL[r.group],lw=.7,ms=2.2,alpha=.7)
ax.plot([0,1],[rec.acute_Total.mean(),rec.recovery_Total.mean()],'-D',color=INK,lw=1.4,ms=3.8);ax.set(xticks=[0,1],xticklabels=['Acute','Recovery'],xlim=(-.2,1.2),ylabel='CRRM score',xlabel='Paired visit');clean(ax,'y');textmm(f,160,141,'GSE69606; 15 pairs',va='bottom')
finish(f,'Figure3_CRRM_clinical_associations','Liu et al. Fig3')

f=plt.figure(figsize=(210/25.4,186/25.4))
letter(f,2.5,12,'A');ax=axmm(f,14,19,47,35);sub=scores[scores.dataset.eq('GSE246622')];vals=[sub.loc[sub.group.eq(g),'age_days'].dropna().to_numpy()/30.4375 for g in levels];b=ax.boxplot(vals,positions=np.arange(4),widths=.55,patch_artist=True,showfliers=False,medianprops=dict(color=INK,lw=.7))
for patch,g in zip(b['boxes'],levels):patch.set_facecolor(COL[g]);patch.set_alpha(.45)
labels=[]
for short,g in zip(['HC','Mild','Mod.','Sev.'],levels):
 v=sub[sub.group.eq(g)];labels.append(f'{short}\nF {v.sex.eq("Female").mean()*100:.0f}%')
ax.set(xticks=np.arange(4),xticklabels=labels,ylabel='Age (months)');clean(ax);textmm(f,15,16,'GSE246622',va='bottom');ax=axmm(f,74,19,24,35);sub=scores[scores.dataset.eq('GSE105450')];b=ax.boxplot([sub.loc[sub.group.eq(g),'age_months'].to_numpy() for g in ['Outpatient','Inpatient']],widths=.5,patch_artist=True,showfliers=False,medianprops=dict(color=INK,lw=.7))
for patch in b['boxes']:patch.set_facecolor(LIGHT);patch.set_alpha(.5)
ax.set(xticks=[1,2],xticklabels=[f'Out\nF {sub.loc[sub.group.eq("Outpatient"),"sex"].eq("Female").mean()*100:.0f}%',f'In\nF {sub.loc[sub.group.eq("Inpatient"),"sex"].eq("Female").mean()*100:.0f}%']);clean(ax);textmm(f,74,16,'GSE105450',va='bottom')

letter(f,105,12,'B');ax=axmm(f,139,19,55,35);coverage=read('06_coverage_summary.csv').set_index('dataset');fc=read('08_functional_branch_coverage.csv');colDS=['GSE246622','GSE246622_RMA','GSE69606','GSE105450'];matrix=[]
for metric,denom in [('n_scored',388),('n_up',280),('n_down',108)]:matrix.append([coverage.loc[ds,metric]/denom for ds in colDS])
for branch in branches:matrix.append([fc.loc[fc.dataset.eq(ds)&fc.branch.eq(branch),'fraction_measured'].iloc[0] for ds in colDS])
im=ax.imshow(matrix,aspect='auto',interpolation='nearest',cmap='YlGnBu',vmin=.5,vmax=1);ax.set(yticks=np.arange(10),yticklabels=['Total','Up','Down','IFN-α','IFN-γ','TNF/NF-κB','TLR signaling','Degranulation','RTK signaling','Microtubule'],xticks=np.arange(4),xticklabels=['246622','RMA','69606','105450']);ax.tick_params(axis='x',length=0,pad=3);ax.spines[:].set_visible(False);cb(f,im,197,19,35,'Measured fraction',[.5,.75,1]);textmm(f,137,62,' / '.join(str(int(coverage.loc[d,"n_scored"])) for d in colDS)+' scored members',va='top')

letter(f,2.5,73,'C');ax=axmm(f,41,79,55,35);common=read('21_common_member_sensitivity.csv').set_index('dataset').loc[DS]
for i,ds in enumerate(DS):err(ax,effects.loc[ds],i,color=INK,offset=-.12);err(ax,common.loc[ds],i,color=UP,offset=.12)
ax.set(yticks=np.arange(3),yticklabels=DS,ylim=(2.5,-.5),xlim=(-.25,.43),xticks=[-.2,0,.2,.4],xlabel='Mean score change (95% CI)');clean(ax,'x');ax.legend(handles=[Line2D([],[],color=INK,marker='o',ms=3,label='Full'),Line2D([],[],color=UP,marker='o',ms=3,label='Common '+str(int(common.n_common.iloc[0])))],loc='lower right',bbox_to_anchor=(1.03,1.02),ncol=2,frameon=False,handlelength=.8,columnspacing=.8)

letter(f,105,73,'D');ax=axmm(f,145,79,57,35);st=read('24_hospitalization_study_strata.csv')
err(ax,effects.loc['GSE105450'],0,color=INK)
for i,(_,r) in enumerate(st.iterrows(),1):err(ax,r,i,color=LIGHT)
ax.set(yticks=[0,1,2],yticklabels=['Overall','Discovery','Validation'],ylim=(2.5,-.5),xlim=(-.31,.2),xticks=[-.3,0,.2],xlabel='Inpatient − outpatient');clean(ax,'x');textmm(f,146,76,'GSE105450; original study strata',va='bottom')

letter(f,2.5,134,'E');ax=axmm(f,34,141,30,32);loo=read('13_patient_deletion_influence.csv').groupby('contrast_key').effect.agg(['min','max'])
for i,ds in enumerate(DS):
 k=effects.loc[ds,'contrast_key'];ax.plot([loo.loc[k,'min'],loo.loc[k,'max']],[i,i],lw=2.2,color=LIGHT,solid_capstyle='round');ax.plot(effects.loc[ds,'effect'],i,'o',ms=3.2,color=INK)
ax.set(yticks=np.arange(3),yticklabels=DS,ylim=(2.5,-.5),xlabel='Deletion effect');clean(ax,'x');ax=axmm(f,73,141,25,32);boot=read('28_patient_bootstrap_gene_map.csv');split=read('29_disjoint_Mild_split_sensitivity.csv');bins=np.linspace(-1,1,31);ax.hist(boot.Spearman_rho,bins=bins,density=True,color=INK,alpha=.45,label='Bootstrap');ax.hist(split.mean_swapped_rho,bins=bins,density=True,color=UP,alpha=.6,label='Split mild');ax.axvline(relation['Spearman_rho'],color=INK,lw=.8,ls='--');ax.set(xlim=(-.35,.7),xticks=[0,.5],yticks=[],xlabel='Gene-effect ρ');clean(ax);ax.legend(loc='upper left',frameon=True,facecolor='white',edgecolor='none',framealpha=.95,handlelength=.7,handletextpad=.4,borderpad=.2,labelspacing=.2)

letter(f,105,134,'F');ax=axmm(f,141,141,60,32)
for i,b in enumerate(branches):
 r=recfx[recfx.metric.eq(b)].iloc[0];err(ax,r,i,color=INK if r.branch_BH7_P<.05 else LIGHT)
ax.set(yticks=np.arange(7),yticklabels=blabel,ylim=(6.5,-.5),xlim=(-1.75,.5),xticks=[-1.5,-1,-.5,0,.5],xlabel='Recovery − acute branch score');clean(ax,'x')
finish(f,'FigureS3_CRRM_clinical_sensitivity','Liu et al. FigS3')
(OUT/'logs/figure_layout_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
for r in audit:print(r)
