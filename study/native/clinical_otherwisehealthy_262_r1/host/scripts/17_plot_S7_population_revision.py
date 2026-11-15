"""Redraw S7 with the updated otherwise-healthy clinical vector; all other panels reuse original results."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
from matplotlib.text import Text

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1/host');F=OUT/'figures';TAG='otherwisehealthy_262_r1';FS=8.2;WIDTH=210
plt.rcParams.update({'font.family':'Arial','font.size':FS,'axes.labelsize':FS,'xtick.labelsize':FS,'ytick.labelsize':FS,'legend.fontsize':FS,'legend.title_fontsize':FS,'axes.titlesize':FS,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
TEAL='#287D8E';ORANGE='#D7813D';BLUE='#486F9A';PURPLE='#81649A';GRAY='#777777'
E=pd.read_csv(OUT/'results/score_effects.csv');S=pd.read_csv(OUT/'results/profile_scores.csv');C=pd.read_csv(OUT/'results/mouse_identity_composition_effects.csv');A=pd.read_csv(OUT/'results/clinical_gene_alignment.csv');G=pd.read_csv(OUT/'results/clinical_gene_alignment_points.csv');O=pd.read_csv(OUT/'results/every_gene_omission.csv');TH=pd.read_csv(OUT/'results/threshold_sensitivities.csv')
branches=['IFN-alpha','IFN-gamma','TNF / NF-kB','TLR signaling','Neutrophil degranulation','RTK signaling','Microtubule']
bn=['IFN-α','IFN-γ','TNF/NF-κB','TLR','Degranulation annotation','RTK (down-directed)','Microtubule (down-directed)']
lineages=['Epithelial_cell','Endothelial_cell','Fibroblast','Mural_cell','Monocyte','Macrophage','Dendritic_cell','Neutrophil','B_cell','T_cell','NK_cell']
ln=['Epithelium','Endothelium','Fibroblasts','Mural cells','Monocytes','Macrophages','DCs','Neutrophils','B cells','T cells','NK cells']
relative=['TNF / NF-kB minus IFN-list overlap','Neutrophil degranulation minus IFN-list overlap','TNF / NF-kB outside IFN lists minus IFN-list overlap','Neutrophil degranulation outside IFN lists minus IFN-list overlap']
payload={'revision':TAG,'statistics_recomputed':True,'recomputed_scope':'S7E whole-blood clinical vector and shared BH3 only','main':{},'supplement':{},'panel_mapping':{'6A':'old6A','6B':'old6B','6C':'old6C','6D':'old6D_epithelium + oldS7C','S7A':'oldS7A','S7B':'old6D_all11lineages','S7C':'oldS7B_rank_only','S7D':'oldS7D','S7E':'old6E','removed':'old6F manuscript-summary boxes'}}
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda a:a.item() if isinstance(a,np.generic) else str(a))+'\n')
def new(height,label):
    fig=plt.figure(figsize=(WIDTH/25.4,height/25.4));fig._mm_height=height;fig.text(204/WIDTH,1-2/height,'Liu et al. '+label,ha='right',va='top',fontsize=12);return fig
def ax(fig,x,y,w,h,gid):
    a=fig.add_axes([x/WIDTH,1-(y+h)/fig._mm_height,w/WIDTH,h/fig._mm_height]);a.set_gid(gid);return a
def letter(fig,label,x,y):fig.text(x/WIDTH,1-y/fig._mm_height,label,fontsize=12,fontweight='bold',va='top')
def clean(a):a.spines[['top','right']].set_visible(False)
def zero(a):a.axvline(0,color=GRAY,lw=.55,zorder=0)
def get(key,score='CRRM388',method='Primary'):
    x=E[E.key.eq(key)&E.score.eq(score)&E.method.eq(method)&E.scope.eq('Intervention')];assert len(x)==1;return x.iloc[0]
def ef(a,r,y,color=TEAL,marker='o',q='adjusted_P'):
    if not np.isfinite(r.estimate):return
    h=a.errorbar(r.estimate,y,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt=marker,ms=3.3,color=color,mfc=color if r[q]<.05 else 'white',mec=color,mew=.75,elinewidth=.75,capsize=1.4)
    h.lines[0].set_gid(str(r.get('test_id',r.get('method','')))+'|'+str(r.get('score','')))
def rec(key,value,supp=False):payload['supplement' if supp else 'main'][key]=value
def score_legend(a,anchor=(1.02,1.04)):
    a.legend([Line2D([],[],marker=m,color=c,ls='',ms=3.3) for c,m in [('black','D'),(PURPLE,'o'),(TEAL,'s')]],['Total','IFN overlap','Outside IFN lists'],frameon=False,ncol=3,loc='lower right',bbox_to_anchor=anchor,borderpad=0,handlelength=.7,handletextpad=.25,columnspacing=.6)
def component_groups(a):
    """One heading per two-row block; use existing row order and colors."""
    a.set(yticks=[0,1,3,4],
          yticklabels=['TNF/NF-κB','Degranulation annotation']*2,
          ylim=(4.65,-1.15))
    a.axhline(1.65,color='#bbbbbb',lw=.5)
    # Headings sit in otherwise empty inter-row bands inside the widened plot.
    for y,label in [(-.8,'Original subsets'),(2.2,'Shared IFN members removed')]:
        a.text(.015,y,label,transform=a.get_yaxis_transform(),
               ha='left',va='center',fontweight='bold',zorder=4,
               bbox=dict(facecolor='white',edgecolor='none',pad=.6))
def save(fig,stem):
    fig.canvas.draw();ren=fig.canvas.get_renderer();off=set();boxes=[];clipped=[];artists={}
    for a in fig.axes:
        for axis in [a.xaxis,a.yaxis]:
            lo,hi=sorted(axis.get_view_interval())
            for t in axis.get_major_ticks():
                if not lo-1e-8<=t.get_loc()<=hi+1e-8:off.update([id(t.label1),id(t.label2)])
        artists[a.get_gid()]=dict(xlim=list(a.get_xlim()),ylim=list(a.get_ylim()),lines=[dict(gid=l.get_gid(),x=np.asarray(l.get_xdata()).tolist(),y=np.asarray(l.get_ydata()).tolist(),marker=l.get_marker(),mfc=str(l.get_markerfacecolor())) for l in a.lines],intervals=[[np.asarray(v).tolist() for v in c.get_segments()] for c in a.collections if hasattr(c,'get_segments')],scatters=[np.asarray(c.get_offsets()).tolist() for c in a.collections if type(c).__name__=='PathCollection'],meshes=[np.asarray(c.get_array()).tolist() for c in a.collections if type(c).__name__=='QuadMesh'])
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text() or id(t) in off:continue
        b=t.get_window_extent(ren);boxes.append(dict(text=t.get_text(),size_pt=t.get_fontsize(),bbox_px=[b.x0,b.y0,b.x1,b.y1]))
        assert t.get_fontsize()*180/WIDTH>=6.99,('Small at180mm',t.get_text())
        if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:clipped.append(t.get_text())
    assert not clipped,(stem,clipped)
    for ext in ['pdf','svg','png']:fig.savefig(F/(stem+'.'+ext),dpi=600)
    dump(F/(stem+'_layout.json'),dict(width_mm=WIDTH,height_mm=fig._mm_height,ordinary_pt=FS,effective_pt_at180mm=FS*180/WIDTH,panel_pt=12,identifier_pt=12,text_boxes=boxes,clipped=clipped,artists=artists,axes=[dict(id=a.get_gid(),position=a.get_position().bounds) for a in fig.axes]))
    plt.close(fig);print('Exported',stem,flush=True)

# One supplement with complete checks, full lineages, rank-only component tests,
# omission ranges, and all three original clinical gene-effect scatterplots.
fig=new(258,'FigS7');letter(fig,'A',4,10);a=ax(fig,51,21,151,43,'S7A');keys=['BRD4_Airway']+['RelA_'+ct for ct in lineages];methods=['Primary','Rank','Exact test','Minimum15cells','Minimum50cells','Detection-supported','Abundance rescaling'];matrix=np.full((12,7),np.nan);rr=[]
for i,key in enumerate(keys):
    for j,method in enumerate(methods):
        if method.startswith('Minimum'):
            if key.startswith('BRD4'):continue
            r=TH[TH.key.eq(key)&TH.score.eq('CRRM388')&TH.method.eq(method)].iloc[0]
        else:r=get(key,method='Primary' if method=='Exact test' else method)
        if not np.isfinite(r.estimate):continue
        q=r.exact_adjusted_P if method=='Exact test' else r.adjusted_P
        matrix[i,j]=np.sign(r.estimate)*(1 if q<.05 else .35);rr.append(dict(key=key,method=method,estimate=r.estimate,adjusted_P=q))
a.pcolormesh(np.arange(8)-.5,np.arange(13)-.5,np.ma.masked_invalid(matrix),cmap=ListedColormap([BLUE,'#c7d6e5','white','#e7cfb7',ORANGE]),vmin=-1,vmax=1,shading='flat')
for i in range(12):
    for j in range(7):
        if not np.isfinite(matrix[i,j]):a.text(j,i,'—',ha='center',va='center',color=GRAY)
        elif abs(matrix[i,j])==1:a.text(j,i,'*',ha='center',va='center',color='white')
a.set(xlim=(-.5,6.5),ylim=(11.5,-.5),yticks=range(12),yticklabels=['Human BRD4i']+['Mouse: '+n for n in ln],xticks=range(7),xticklabels=['Gene-z','Rank','Exact\npermutation','≥15 cells','≥50 cells','Detection-\nsupported','TPM\nrescaling']);a.tick_params(length=0)
for y in [.5,1.5]:a.axhline(y,color='#555555',lw=.7 if y==.5 else .5)
a.legend([Patch(color=BLUE),Patch(color=ORANGE),Line2D([],[],marker='*',color=GRAY,ls='',ms=4),Line2D([],[],color=GRAY,lw=.7)],['Lower','Higher','BH < 0.05','Unavailable / inapplicable'],frameon=False,ncol=4,loc='lower right',bbox_to_anchor=(1.01,1.055),borderpad=0,handlelength=.7,handletextpad=.3,columnspacing=1);rec('S7A',dict(matrix=matrix.tolist(),rows=rr,keys=keys,methods=methods),True)

letter(fig,'B',4,82);a=ax(fig,41,94,61,68,'S7B');rr=[]
for i,ct in enumerate(lineages):
    for off,target,col,marker in [(-.20,'CRRM388','black','D'),(0,'IFN-list overlap',PURPLE,'o'),(.20,'Outside IFN lists',TEAL,'s')]:
        r=get('RelA_'+ct,target);ef(a,r,i+off,col,marker);rr.append(r.to_dict())
zero(a);clean(a);a.set(yticks=range(11),yticklabels=ln,ylim=(10.7,-.7),xlabel='RelA CKO − RSV vehicle, day 21\n(gene-z; same 4 vs 3 lung samples)',xlim=(-.9,.6),xticks=[-.8,-.4,0,.4]);score_legend(a);rec('S7B',rr,True)

letter(fig,'C',113,82);a=ax(fig,148,94,54,38,'S7C');rr=[]
for i,(y,t) in enumerate(zip([0,1,3,4],relative)):
    r=get('BRD4_Airway',t,'Rank');ef(a,r,y,BLUE if i<2 else TEAL);rr.append(r.to_dict())
component_groups(a);zero(a);clean(a)
# Preserve the exact r2 auto-range; add a readable positive tick within it.
rank_xlim=json.loads((F/'FigureS7_component_focus_r2_layout.json').read_text())['artists']['S7C']['xlim']
a.set(xlabel='Effect difference vs IFN overlap\n(percentile-rank units)',
      xlim=rank_xlim,xticks=[-.05,-.025,0,.025]);rec('S7C',rr,True)

letter(fig,'D',113,145);rr=[]
for j,(key,x,w,col,label) in enumerate([('BRD4_Airway',136,26,TEAL,'Human BRD4i'),('RelA_Epithelial_cell',178,24,ORANGE,'Mouse epithelium')]):
    a=ax(fig,x,157,w,18,'S7D_'+str(j))
    for i,target in enumerate(['CRRM388','Neutrophil degranulation']):
        g=O[O.key.eq(key)&O.score.eq(target)];full=float(g.full_estimate.iloc[0]);lo=g.estimate_without_gene.min();hi=g.estimate_without_gene.max();a.plot([lo,hi],[i,i],color=col,lw=2);a.plot(full,i,'D',color='black',ms=3.3);rr.append(dict(key=key,score=target,full=full,minimum=lo,maximum=hi,gene_omissions=len(g)))
    clean(a);a.set(yticks=[0,1],yticklabels=['Total','Degranulation'] if j==0 else [],ylim=(1.5,-.5),xlabel='Gene-z effect\nInfluence range');a.locator_params(axis='x',nbins=3);a.text(.5,1.04,label,transform=a.transAxes,ha='center',va='bottom')
rec('S7D',rr,True)

letter(fig,'E',4,195);rr=[]
for j,(key,x,label) in enumerate([('GSE246622_Severe_Mild',18,'Whole blood: Severe−Mild'),('GSE69606_Severe_Mild',85,'PBMC: Severe−Mild'),('GSE69606_Recovery_Acute',152,'PBMC: Recovery−Acute')]):
    a=ax(fig,x,207,50,32,'S7E_'+str(j));g=G[G.clinical_contrast.eq(key)];r=A[A.clinical_contrast.eq(key)].iloc[0]
    a.scatter(g.clinical_log_effect,g.intervention_log_effect,c=GRAY,s=10,lw=0,alpha=.85);a.axhline(0,color='#aaaaaa',lw=.5);a.axvline(0,color='#aaaaaa',lw=.5);clean(a);a.set(xlabel=label+'\nLog-expression difference',ylim=(-6.7,1.5));a.locator_params(axis='x',nbins=3);a.set_yticks([-6,-3,0])
    if j==0:a.set_ylabel('BRD4i − RSV\nLog-expression difference')
    a.text(.02,1.04,'n = '+str(r.members)+' genes\nρ = '+f'{r.rho:.2f}'+'; BH P = '+f'{r.adjusted_P_BH3:.2f}',transform=a.transAxes,ha='left',va='bottom');rr.append(dict(statistic=r.to_dict(),points=g.to_dict('records')))
rec('S7E',rr,True);save(fig,'FigureS7_'+TAG)
dump(OUT/'audit'/f'figure_payload_{TAG}.json',payload)
