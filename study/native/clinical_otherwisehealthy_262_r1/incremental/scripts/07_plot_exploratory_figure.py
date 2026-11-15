from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1/incremental')
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.55,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':0,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#387C98';ORANGE='#BA603A';INK='#244B60';GREY='#A6B2B9';colors={'GSE246622':BLUE,'GSE69606':ORANGE,'GSE105450':GREY};DS=['GSE246622','GSE69606'];read=lambda f:pd.read_csv(OUT/'tables'/f)
ab=read('07_LCN2_ablation_effects.csv');models=read('08_conditional_clinical_models.csv');infl=read('10_every_gene_omission_influence.csv');cm=read('17_measured_composition_models.csv');perf=read('20_descriptive_out_of_fold_comparison.csv');rec=read('14_paired_recovery_ablation.csv')
f=plt.figure(figsize=(210/25.4,187/25.4))
def axmm(x,y,w,h):return f.add_axes([x/210,1-(y+h)/187,w/210,h/187])
def textmm(x,y,text,**kw):return f.text(x/210,1-y/187,text,**kw)
def letter(x,y,text):textmm(x,y,text,fontsize=12,fontweight='bold',va='top')
def clean(ax):
    ax.spines[['top','right']].set_visible(False);ax.axvline(0,color=GREY,lw=.6,zorder=0)
def err(ax,r,y,color,marker='o',mfc=None):
    if not hasattr(ax,'_effect_intervals'):ax._effect_intervals=[]
    ax._effect_intervals.append((float(r.CI95_low),float(r.CI95_high)))
    ax.errorbar(r.effect,y,xerr=[[r.effect-r.CI95_low],[r.CI95_high-r.effect]],fmt=marker,ms=3.3,lw=.85,color=color,mfc=color if mfc is None else mfc,capsize=1.5)

letter(2.5,13,'A');ax=axmm(43,23,52,36)
for k,ds in enumerate(DS):
    for j,target in enumerate(['CRRM25','CRRM24','CRRMDisjoint']):err(ax,ab[ab.dataset.eq(ds)&ab.target.eq(target)].iloc[0],j+(-.12 if k==0 else .12),colors[ds])
ax.set(yticks=range(3),yticklabels=['Full subset','LCN2 omitted','LCN2 and IFN\nmembers omitted'],ylim=(2.48,-.48),xlim=(-.05,1.52),xticks=[0,.5,1,1.5],xlabel='Severe − mild (95% CI)');clean(ax)
handles=[Line2D([],[],color=colors[ds],marker='o',ms=3,lw=.8,label=ds) for ds in DS]
ax.legend(handles=handles,loc='lower right',bbox_to_anchor=(1,1.04),ncol=2,frameon=False,handlelength=.6,columnspacing=.7,borderaxespad=0)

letter(105,13,'B');ax=axmm(154,23,49,36)
for k,ds in enumerate(DS):
    for j,model in enumerate(['Clinical','Clinical+IFNUnion','Clinical+Zhang10','Clinical+Zhang10+IFNUnion']):err(ax,models[models.dataset.eq(ds)&models.target.eq('CRRM24')&models.model.eq(model)].iloc[0],j+(-.12 if k==0 else .12),colors[ds])
ax.set(yticks=range(4),yticklabels=['Clinical factors','+ IFN reference','+ Published panel','+ Both references'],ylim=(3.5,-.5),xlim=(-1.2,1.85),xticks=[-1,0,1],xlabel='Conditional severity difference');clean(ax)

letter(2.5,75,'C');ax=axmm(43,86,52,30);rng=np.random.default_rng(20261122)
for j,ds in enumerate(DS):
    d=infl[infl.dataset.eq(ds)&infl.model.eq('Clinical')];r=models[models.dataset.eq(ds)&models.target.eq('CRRM24')&models.model.eq('Clinical')].iloc[0]
    ax.plot([d.effect.min(),d.effect.max()],[j,j],color=colors[ds],lw=1,zorder=1);ax.scatter(d.effect,j+rng.uniform(-.12,.12,len(d)),color=colors[ds],s=7,alpha=.65,edgecolors='none');ax.plot(r.effect,j,marker='D',color=INK,ms=4,ls='none',zorder=3)
ax.set(yticks=[0,1],yticklabels=['Whole blood\nGSE246622','PBMC\nGSE69606'],ylim=(1.48,-.48),xlim=(.18,.91),xticks=[.2,.4,.6,.8],xlabel='Severity effect after member omission');clean(ax)
ax.legend(handles=[Line2D([],[],marker='D',color=INK,ls='none',ms=3,label='All remaining'),Line2D([],[],marker='o',color=BLUE,ls='none',ms=3,label='One more omitted')],loc='lower right',bbox_to_anchor=(1,1.05),ncol=2,frameon=False,handlelength=.5,columnspacing=.7,borderaxespad=0)

letter(105,75,'D');ax=axmm(154,86,49,30)
for j,typ in enumerate(['Measured_myeloid','Measured_CLR3']):
    for k,state in enumerate(['before','after']):
        r=cm[cm.model.eq(typ+'__Same_cases_'+state+'_composition')].iloc[0];err(ax,r,j+(-.12 if k==0 else .12),GREY if k==0 else BLUE,marker='s' if k==0 else 'o',mfc='white' if k==0 else BLUE)
ax.set(yticks=[0,1],yticklabels=['Myeloid fractions\n'+str(int(cm.loc[cm.model.eq("Measured_myeloid__Same_cases_after_composition"),"n"].iloc[0]))+' matched cases','All-cluster CLR PCs\n'+str(int(cm.loc[cm.model.eq("Measured_CLR3__Same_cases_after_composition"),"n"].iloc[0]))+' complete cases'],ylim=(1.5,-.5),xlim=(-.4,.82),xticks=[-.4,0,.4,.8],xlabel='Conditional severity difference');clean(ax)
ax.legend(handles=[Line2D([],[],marker='s',color=GREY,mfc='white',ms=3,lw=.8,label='Before'),Line2D([],[],marker='o',color=BLUE,ms=3,lw=.8,label='After')],loc='lower right',bbox_to_anchor=(1,1.05),ncol=2,frameon=False,handlelength=.6,columnspacing=1,borderaxespad=0)

letter(2.5,132,'E');ax=axmm(43,143,52,30)
for j,ds in enumerate(DS+['GSE105450']):
    r=perf[perf.dataset.eq(ds)].iloc[0]
    for k,(metric,mark) in enumerate([('logloss','s'),('Brier','o')]):
        gain=100*r[metric+'_improvement']/r['Baseline_'+metric];ax.plot(gain,j+(-.1 if k==0 else .1),marker=mark,color=colors[ds],ms=4,ls='none',mfc='white' if k==0 else colors[ds])
ax.axhline(1.5,color=GREY,lw=.5,ls=(0,(2,2)));ax.set(yticks=range(3),yticklabels=['GSE246622','GSE69606','Hospitalization\nGSE105450'],ylim=(2.5,-.5),xlim=(-2.2,7.0),xticks=[-2,0,2,4,6],xlabel='Held-out loss reduction (%)');clean(ax)
ax.legend(handles=[Line2D([],[],marker='s',color=INK,mfc='white',ls='none',ms=3,label='Log loss'),Line2D([],[],marker='o',color=INK,ls='none',ms=3,label='Brier')],loc='lower right',bbox_to_anchor=(1,1.05),ncol=2,frameon=False,handlelength=.5,columnspacing=1,borderaxespad=0)

letter(105,132,'F');ax=axmm(154,143,49,30)
metrics=['CRRM25','CRRM24','CRRMDisjoint','Zhang10','IFNUnion'];labels=['Full subset','LCN2 omitted','LCN2 and IFN\nmembers omitted','Published panel','IFN reference']
for j,metric in enumerate(metrics):err(ax,rec[rec.metric.eq(metric)].iloc[0],j,INK if metric.startswith('CRRM') else GREY)
ax.set(yticks=range(5),yticklabels=labels,ylim=(4.5,-.5),xlim=(-1.85,.1),xticks=[-1.5,-1,-.5,0],xlabel='Recovery − acute (95% CI)');clean(ax)

textmm(207,3,'Liu et al. FigS4',fontsize=12,ha='right',va='top')
for a in f.axes:
    if hasattr(a,'_effect_intervals'):
        oldlo,oldhi=a.get_xlim();lo=min(v[0] for v in a._effect_intervals);hi=max(v[1] for v in a._effect_intervals)
        if lo<oldlo or hi>oldhi:
            span=max(oldhi-oldlo,hi-lo);a.set_xlim(min(oldlo,lo-.04*span),max(oldhi,hi+.04*span))
            from matplotlib.ticker import MaxNLocator
            a.xaxis.set_major_locator(MaxNLocator(nbins=4))
f.canvas.draw();rend=f.canvas.get_renderer();outside=[];texts=[t for t in f.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
for t in texts:
    b=t.get_window_extent(rend)
    if not np.isfinite(b.extents).all() or b.x0<-.2 or b.y0<-.2 or b.x1>f.bbox.width+.2 or b.y1>f.bbox.height+.2:outside.append(dict(text=t.get_text(),bbox_px=b.extents.tolist()))
for ext in ['pdf','svg','png']:f.savefig(OUT/'figures'/f'FigureS4_CRRM_incremental_information.{ext}',dpi=600 if ext=='png' else 300)
f.savefig(OUT/'logs/FigureS4_preview.png',dpi=150)
audit=dict(width_mm=210,height_mm=187,ordinary_pt=7,panel_letters_pt=12,identifier_pt=12,text_sizes_pt=sorted(set(t.get_fontsize() for t in texts)),all_text_inside_canvas=not outside,outside_text=outside,single_compact_layout=True)
(OUT/'logs/figure_layout_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n');print(json.dumps(audit,ensure_ascii=False))
