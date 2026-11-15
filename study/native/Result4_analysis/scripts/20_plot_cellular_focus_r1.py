"""Reorganize existing CRRM results; no score fitting or statistical tests."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,zipfile
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle
from matplotlib.text import Text

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis');ROOT=Path(str(Path(__import__('os').environ['CRRM_DATA_ROOT']) / ''));REV=Path(__import__('os').environ['CRRM_WORKSPACE']);CURRENT=ROOT/'当前稿件';TAG='cellular_focus_r1';F=OUT/'figures'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=True,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
# Archival of local author documents is outside the public figure replay.
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'legend.title_fontsize':7,'axes.titlesize':7,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
P=pd.read_csv(OUT/'results/primary_cellular_effects.csv');M=pd.read_csv(OUT/'results/measurement_method_sensitivities.csv');V=pd.read_csv(OUT/'results/minimum_cell_sensitivities.csv');G=pd.read_csv(OUT/'results/fixed_gene_contributions.csv')
ids=['PBMC_B_cell','PBMC_T_cell','PBMC_NK_cell','PBMC_Monocyte','Berg_Epithelial_cell','DC_Dendritic_cell','DC_Monocyte']
short=['B','T','NK','Mono.','Airway','DC','Mono. (DC)'];heat_labels=['B','T','NK','Mono.','Airway','DC','Mono.\n(DC)']
source_order=['Berg_ALI','CD1c_DC','Down_syndrome','HLCA','Infant_PBMC','Infant_adult_ALI','RSV_history_wheeze','Sunshine','Yoshida_airway','iAT2_macrophage','iPSC_organoid','nasal_submerged_model','pathogen_specific_TRM']
source_names=['Berg ALI','CD1c DC','Down syndrome','HLCA','Infant PBMC','Infant/adult ALI','RSV history/wheeze','Sunshine','Yoshida airway','iAT2/macrophage','iPSC organoid','Nasal submerged','Pathogen-specific TRM']
TEAL='#287D8E';ORANGE='#D7813D';BLUE='#486F9A';PURPLE='#81649A';GRAY='#666666';colors=[TEAL]*4+['#C76730',PURPLE,PURPLE]
branches=['IFN-alpha','IFN-gamma','TNF / NF-kB','TLR signaling','Neutrophil degranulation','RTK signaling','Microtubule'];branch_labels=['IFN-α','IFN-γ','TNF / NF-κB','TLR','Degranulation-annotated','RTK (down-directed)','Microtubule (down-directed)']
payload={'version':TAG,'statistics_recomputed':False,'panel_mapping':{'4A':'old4B','4B':'old4G','4C':'old4C','4D':'old4D','4E':'old4E','4F':'old4F','S5A':'old4A','S5B':'oldS5B','S5C':'oldS5C','S5D':'oldS5D','S5E':'oldS5E','S5F':'oldS5F'},'panels':{},'omitted_primary_ciliated':P[P.comparison_id.eq('Berg_Ciliated')&P.score.eq('CRRM388')].to_dict('records')}
def new(height,label):
    fig=plt.figure(figsize=(210/25.4,height/25.4));fig._mm_height=height;fig.text(204/210,1-2/height,'Liu et al. '+label,ha='right',va='top',fontsize=12);return fig
def ax(fig,x,y,w,h,gid):
    a=fig.add_axes([x/210,1-(y+h)/fig._mm_height,w/210,h/fig._mm_height]);a.set_gid(gid);return a
def letter(fig,label,x,y):fig.text(x/210,1-y/fig._mm_height,label,fontsize=12,fontweight='bold',va='top')
def clean(a):a.spines[['top','right']].set_visible(False)
def zero(a):a.axvline(0,color='#777777',lw=.6,zorder=0)
def heatmap(a,values,cmap,vmin,vmax):
    values=np.asarray(values);nr,nc=values.shape
    mesh=a.pcolormesh(np.arange(nc+1)-.5,np.arange(nr+1)-.5,np.ma.masked_invalid(values),cmap=cmap,vmin=vmin,vmax=vmax,shading='flat',edgecolors='none',rasterized=False)
    a.set(xlim=(-.5,nc-.5),ylim=(nr-.5,-.5));return mesh
def effect(a,r,y,color,marker='o',filled=True):
    result=a.errorbar(r.estimate,y,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt=marker,ms=3,color=color,mfc=color if filled else 'white',mec=color,mew=.7,elinewidth=.7,capsize=1.4)
    result.lines[0].set_gid(str(r.comparison_id) if hasattr(r,'comparison_id') else str(r.score)+'|'+str(getattr(r,'model','')));return result
def save(fig,stem):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();boxes=[];clipped=[];tiny=[];off_axis_labels=set();artists={}
    for a in fig.axes:
        for c in a.collections:c.set_rasterized(False)
        for axis in [a.xaxis,a.yaxis]:
            lo,hi=sorted(axis.get_view_interval())
            for tick in list(axis.get_major_ticks())+list(axis.get_minor_ticks()):
                if not lo-1e-9<=tick.get_loc()<=hi+1e-9:off_axis_labels.update([id(tick.label1),id(tick.label2)])
        artists[a.get_gid()]={'lines':[{'gid':l.get_gid(),'x':np.asarray(l.get_xdata()).tolist(),'y':np.asarray(l.get_ydata()).tolist(),'marker':l.get_marker(),'facecolor':str(l.get_markerfacecolor())} for l in a.lines], 'interval_segments':[[np.asarray(v).tolist() for v in c.get_segments()] for c in a.collections if hasattr(c,'get_segments')], 'meshes':[np.asarray(c.get_array()).tolist() for c in a.collections if type(c).__name__=='QuadMesh'], 'bars':[{'x':p.get_x(),'y':p.get_y(),'width':p.get_width(),'height':p.get_height()} for p in a.patches if isinstance(p,Rectangle)]}
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text() or id(t) in off_axis_labels:continue
        b=t.get_window_extent(renderer)
        if t.get_fontsize()<6.99:tiny.append(t.get_text())
        if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:clipped.append(t.get_text())
        boxes.append({'text':t.get_text(),'size_pt':t.get_fontsize(),'bbox_px':[b.x0,b.y0,b.x1,b.y1]})
    assert not tiny,('Text below7pt',tiny)
    assert not clipped,('Clipped text',stem,clipped)
    for ext in ['pdf','svg','png']:fig.savefig(F/(stem+'.'+ext),dpi=600)
    dump(F/(stem+'_layout.json'),dict(width_mm=210,height_mm=fig._mm_height,ordinary_pt=7,panel_pt=12,identifier_pt=12,clipped=clipped,text_boxes=boxes,artists=artists,axes=[dict(id=a.get_gid(),position=a.get_position().bounds) for a in fig.axes]))
    plt.close(fig);print('Exported',stem,flush=True)

# Figure4: the three inferential studies and their biological findings.
fig=new(184,'Fig4');letter(fig,'A',4,8);a=ax(fig,49,14,56,38,'4A')
labels=['B cells (19 RSV, 17 HC)','T cells (19 RSV, 17 HC)','NK cells (19 RSV, 17 HC)','Monocytes (19 RSV, 17 HC)','Airway epithelium (6 donors)','DCs (3 donors)','Monocytes, DC culture (3 donors)']
totals=P[P.score.eq('CRRM388')].set_index('comparison_id').loc[ids].reset_index();payload['panels']['4A']=totals.to_dict('records')
for y,r in enumerate(totals.itertuples(index=False)):effect(a,r,y,colors[y],'o' if y<4 else 's',r.adjusted_P<.05)
zero(a);a.set(yticks=range(7),yticklabels=labels,ylim=(6.65,-.65),xlim=(-.08,1.12),xticks=[0,.5,1],xlabel='Within-study CRRM change');clean(a)
a.legend([Line2D([],[],marker=m,color=c,ls='',ms=3) for c,m in [(TEAL,'o'),('#C76730','s'),(PURPLE,'s')]],['Infant PBMC','Airway culture','DC culture'],ncol=3,frameon=False,loc='lower left',bbox_to_anchor=(-.02,1.02),borderpad=0,handlelength=.8,handletextpad=.25,columnspacing=.6)

letter(fig,'B',113,8);a=ax(fig,137,14,65,38,'4B');contributions=[]
for y,c in enumerate(ids):
    up=P[P.comparison_id.eq(c)&P.score.eq('Up contribution')].iloc[0];down=P[P.comparison_id.eq(c)&P.score.eq('Down contribution')].iloc[0];total=totals.iloc[y]
    assert abs(up.estimate+down.estimate-total.estimate)<1e-12
    a.barh(y-.16,up.estimate,height=.29,color=ORANGE);a.barh(y+.16,down.estimate,height=.29,color=BLUE);a.plot(up.estimate+down.estimate,y,'D',color='black',ms=2.3)
    contributions.append(dict(comparison_id=c,up=up.estimate,down=down.estimate,total=total.estimate,case_units=int(total.case_units),control_units=int(total.control_units)))
payload['panels']['4B']=contributions
zero(a);a.set(yticks=range(7),yticklabels=short,ylim=(6.65,-.65),xlim=(-.06,1.0),xticks=[0,.4,.8],xlabel='Contribution to total CRRM change');clean(a)
a.legend([Patch(color=ORANGE),Patch(color=BLUE),Line2D([],[],marker='D',color='black',ls='',ms=2.5)],['Up-directed','Down-directed','Total'],frameon=False,ncol=3,loc='lower right',bbox_to_anchor=(1.01,1.02),borderpad=0,handlelength=.8,handletextpad=.3,columnspacing=.7)

letter(fig,'C',4,65);a=ax(fig,43,71,64,31,'4C')
z=P[P.comparison_id.isin(ids)&P.score.isin(branches)].pivot(index='score',columns='comparison_id',values='estimate').reindex(index=branches,columns=ids)
q=P[P.comparison_id.isin(ids)&P.score.isin(branches)].pivot(index='score',columns='comparison_id',values='adjusted_P').reindex(index=branches,columns=ids)
im=heatmap(a,z,'RdBu_r',-2,2);a.set(xticks=range(7),xticklabels=heat_labels,yticks=range(7),yticklabels=branch_labels);a.tick_params(length=0)
for y,x in np.argwhere(q.to_numpy()<.05):a.plot(x,y,'o',ms=1.5,color='black')
payload['panels']['4C']={'columns':ids,'rows':branches,'estimates':z.to_numpy().tolist(),'q':q.to_numpy().tolist(),'original_correction_family_tests':56}
cb=fig.colorbar(im,cax=ax(fig,43,112,64,2,'4C:colorbar'),orientation='horizontal',ticks=[-2,0,2]);cb.set_label('Signed branch change',labelpad=0);cb.ax.tick_params(pad=1,length=1)

letter(fig,'D',113,65);a=ax(fig,135,71,67,31,'4D');rows=[]
for j,target in enumerate(['IFN-list overlap','Outside IFN lists']):
    for y,c in enumerate(ids):
        r=P[P.comparison_id.eq(c)&P.score.eq(target)].iloc[0];effect(a,r,y+(j-.5)*.25,TEAL if j==0 else ORANGE,'o' if j==0 else 's');rows.append(r.to_dict())
payload['panels']['4D']=rows
zero(a);a.set(yticks=range(7),yticklabels=short,ylim=(6.6,-.6),xlim=(-.1,2.35),xticks=[0,1,2],xlabel='Within-study component change');clean(a)
a.legend([Line2D([],[],color=TEAL,marker='o',ls='',ms=3),Line2D([],[],color=ORANGE,marker='s',ls='',ms=3)],['IFN-list overlap','Outside IFN lists'],frameon=False,ncol=2,loc='lower right',bbox_to_anchor=(1.01,1.02),borderpad=0,handletextpad=.25,columnspacing=.6)

letter(fig,'E',4,128);a=ax(fig,29,135,53,33,'4E');cc=pd.read_csv(OUT/'results/PBMC_composition_standardization.csv');cc=cc[cc.endpoint.eq('RSV_vs_Healthy')];rows=[]
for y,target in enumerate(['CRRM388','Outside IFN lists']):
    for j,model in enumerate(['Observed composition','Fixed healthy composition']):
        r=cc[cc.score.eq(target)&cc.model.eq(model)].iloc[0];effect(a,r,y+(j-.5)*.25,ORANGE if j==0 else TEAL,'o' if j==0 else 's');rows.append(r.to_dict())
payload['panels']['4E']=rows
zero(a);a.set(yticks=[0,1],yticklabels=['CRRM','Outside IFN lists'],ylim=(1.6,-.7),xlim=(0,.28),xticks=[0,.1,.2],xlabel='RSV − healthy, four-lineage PBMC');clean(a)
a.legend([Line2D([],[],color=ORANGE,marker='o',ls='',ms=3),Line2D([],[],color=TEAL,marker='s',ls='',ms=3)],['Observed composition','Fixed healthy composition'],frameon=False,ncol=1,loc='lower left',bbox_to_anchor=(-.02,1.04),borderpad=0,handletextpad=.25,labelspacing=.3)

letter(fig,'F',92,128);balanced=pd.read_csv(OUT/'results/Berg_three_state_matched_profiles.csv');state=['ctrl','bystdr','inf'];state_names=['Control','Bystander','Infected'];assert balanced.analysis_unit.nunique()==6
payload['panels']['4F']={'profiles':balanced.to_dict('records'),'independent_donors':6,'matched_donor_time_strata':22,'aggregation':'equal matched-time mean within donor'}
for k,(target,x) in enumerate(zip(['CRRM388','IFN-alpha','Outside IFN lists'],[105,140,175])):
    a=ax(fig,x,135,27,33,'4F:'+target);pivot=balanced.pivot(index='analysis_unit',columns='infection_status',values=target).reindex(columns=state)
    for donor,r in pivot.iterrows():a.plot(range(3),r.to_numpy(),color='#b8b8b8',lw=.6,marker='o',ms=1.8,gid='donor:'+donor)
    a.plot(range(3),pivot.mean().to_numpy(),color=TEAL,lw=1.2,marker='o',ms=3,gid='mean')
    a.set(xticks=range(3),xticklabels=state_names,xlim=(-.25,2.25),xlabel=['CRRM','IFN-α','Outside IFN\nlists'][k]);a.xaxis.labelpad=1;clean(a)
    if k==0:
        a.set_ylabel('Relative score',labelpad=1)
        a.legend([Line2D([],[],color='#b8b8b8',lw=.6,marker='o',ms=2),Line2D([],[],color=TEAL,lw=1.2,marker='o',ms=3)],['Donor','Mean'],frameon=False,ncol=2,loc='lower left',bbox_to_anchor=(1.45,1.04),borderpad=0,handletextpad=.3,columnspacing=.7)
save(fig,'Figure4_'+TAG)

# FigureS5: descriptive background and all existing substantive checks.
fig=new(214,'FigS5');letter(fig,'A',4,8);a=ax(fig,35,16,62,41,'S5A')
d=pd.read_csv(OUT/'results/atlas_human_lineage_footprint.csv');local=pd.read_csv(OUT/'results/atlas_source_relative_localization.csv');local=local.set_index(['source_label','display_lineage']).relative_CRRM388
lineages=['Epithelial','Myeloid','B','T / NK','Stromal','Other'];norm=plt.Normalize(-1.5,1.5)
for iy,source in enumerate(source_order):
    for ix,lin in enumerate(lineages):
        z=d[d.source_label.eq(source)&d.display_lineage.eq(lin)]
        if len(z) and int(z.cells.iloc[0]):
            n=int(z.cells.iloc[0]);s=3+7*np.log10(max(n,1))**1.35;value=local.get((source,lin),np.nan)
            a.scatter(ix,iy,s=s,c=[value] if np.isfinite(value) else ['#dddddd'],cmap='RdBu_r' if np.isfinite(value) else None,norm=norm if np.isfinite(value) else None,lw=0)
a.set(xticks=range(6),xticklabels=['Epi.','Myeloid','B','T/NK','Stroma','Other'],yticks=range(13),yticklabels=source_names,xlim=(-.65,5.65),ylim=(12.65,-.65),xlabel='Descriptive cellular background');a.tick_params(length=0)
for spine in a.spines.values():spine.set_visible(False)
handles=[a.scatter([],[],s=3+7*np.log10(n)**1.35,c=TEAL,lw=0,label=f'{n:,}') for n in [1000,10000,100000]]
a.legend(handles=handles,title='Indexed cells',loc='lower center',bbox_to_anchor=(.5,1.02),ncol=3,frameon=False,handletextpad=.25,columnspacing=.7,borderpad=0,labelspacing=.1)
sm=plt.cm.ScalarMappable(norm=norm,cmap='RdBu_r');cb=fig.colorbar(sm,cax=ax(fig,100,16,2,41,'S5A:colorbar'),ticks=[-1.5,0,1.5]);cb.set_label('Within-source CRRM (descriptive)',labelpad=1);cb.ax.tick_params(pad=1,length=1)
payload['panels']['S5A']={'cells':d.to_dict('records'),'localization':local.reset_index().to_dict('records'),'role':'descriptive only','source_objects':13,'human_cells':int(d.cells.sum())}

letter(fig,'B',116,8);a=ax(fig,140,16,62,41,'S5B');order=ids[:5]+['Berg_Ciliated']+ids[5:];rows=[]
for j,threshold in enumerate([15,30,50]):
    frame=P if threshold==30 else V[V.minimum_cells.eq(threshold)]
    for y,c in enumerate(order):
        r=frame[frame.comparison_id.eq(c)&frame.score.eq('CRRM388')].iloc[0];rows.append(r.to_dict())
        if np.isfinite(r.estimate):effect(a,r,y+(j-1)*.2,[ORANGE,TEAL,BLUE][j],['^','o','s'][j])
payload['panels']['S5B']=rows
zero(a);a.set(yticks=range(8),yticklabels=['B','T','NK','Monocytes','Airway','Ciliated','DC','Mono. (DC)'],ylim=(7.7,-.7),xlim=(-.1,1.7),xticks=[0,.5,1,1.5],xlabel='Within-study CRRM change');clean(a)
a.legend([Line2D([],[],color=c,marker=m,ms=3,ls='') for c,m in zip([ORANGE,TEAL,BLUE],['^','o','s'])],['15 cells','30 cells','50 cells'],frameon=False,ncol=3,loc='lower right',bbox_to_anchor=(1.01,1.02),borderpad=0,handletextpad=.1,columnspacing=.5)

letter(fig,'C',4,70);a=ax(fig,36,74,65,29,'S5C');matrix=[];rlab=[];records=[]
for target,tname in [('CRRM388','Total'),('Outside IFN lists','Outside'),('Neutrophil degranulation','Degran.')]:
    for method,mname in [('Primary','Z'),('Common expression-supported genes','Supported'),('Signed transcriptome rank','Rank')]:
        frame=P if method=='Primary' else M[M.measurement_method.eq(method)];qcol='adjusted_P' if method=='Primary' else 'descriptive_BH8_P';values=[]
        for c in ids:
            r=frame[frame.comparison_id.eq(c)&frame.score.eq(target)].iloc[0];positive=r.estimate>=0;support=r[qcol]<.05;values.append(3 if positive and support else 2 if positive else 0 if support else 1);records.append(dict(comparison_id=c,score=target,method=method,estimate=r.estimate,q=r[qcol],color_category=values[-1]))
        matrix.append(values);rlab.append(f'{tname} / {mname}')
im=heatmap(a,matrix,ListedColormap(['#b34f26','#f3d3bb','#c8dfe4',TEAL]),-.5,3.5);a.set(xticks=range(7),xticklabels=heat_labels,yticks=range(9),yticklabels=rlab);a.tick_params(length=0)
for y in [2.5,5.5]:a.axhline(y,color='white',lw=1)
fig.text(36/210,1-114/214,'Dark: BH < 0.05; light: BH ≥ 0.05',fontsize=7);fig.text(36/210,1-118/214,'Teal: positive; orange: negative',fontsize=7)
payload['panels']['S5C']={'rows':rlab,'columns':ids,'categories':matrix,'statistics':records}

letter(fig,'D',109,70);a=ax(fig,131,74,71,29,'S5D');adjusted=pd.read_csv(OUT/'results/PBMC_covariate_adjusted_effects.csv');author=pd.read_csv(OUT/'results/PBMC_author_scope_effects.csv');rows=[]
for j,kind in enumerate(['Primary','Age / sex','Author subset']):
    for y,c in enumerate(ids[:4]):
        if j==0:r=P[P.comparison_id.eq(c)&P.score.eq('CRRM388')].iloc[0]
        elif j==1:r=adjusted[adjusted.endpoint.eq('RSV_vs_Healthy')&adjusted.cell_identity.eq(c[5:])&adjusted.score.eq('CRRM388')].iloc[0]
        else:r=author[author.comparison_id.eq(c)&author.score.eq('CRRM388')].iloc[0]
        effect(a,r,y+(j-1)*.2,[TEAL,BLUE,ORANGE][j],['o','s','^'][j]);rows.append({**r.to_dict(),'comparison_id':c,'method':kind})
payload['panels']['S5D']=rows
zero(a);a.set(yticks=range(4),yticklabels=['B','T','NK','Monocytes'],ylim=(3.6,-.6),xlim=(0,.55),xticks=[0,.25,.5],xlabel='RSV − healthy CRRM change');clean(a)
a.legend([Line2D([],[],color=c,marker=m,ms=3,ls='') for c,m in zip([TEAL,BLUE,ORANGE],['o','s','^'])],['Primary','Age / sex','Author subset'],frameon=False,loc='lower right',bbox_to_anchor=(1.01,1.03),ncol=3,borderpad=0,handletextpad=.1,columnspacing=.7)

degran=G[G.in_degranulation].display_symbol.drop_duplicates().sort_values().tolist();gd=G[G.in_degranulation].copy();gd['mean_detection']=(gd.case_detection_fraction+gd.control_detection_fraction)/2;assert len(degran)==25
for side,x,target,scale,cmap_name,limits,clabel in [('E',23,'mean_detection',100,'Blues',(0,100),'Cells with detected RNA (%)'),('F',124,'degranulation_score_contribution',1,'RdBu_r',(-.1,.1),'Gene contribution to branch change')]:
    letter(fig,side,4 if side=='E' else 109,123);a=ax(fig,x,126,77,65,'S5'+side)
    z=gd.pivot(index='display_symbol',columns='comparison_id',values=target).reindex(index=degran,columns=ids)*scale;cm=plt.get_cmap(cmap_name).copy();cm.set_bad('#bdbdbd')
    im=heatmap(a,z,cm,limits[0],limits[1]);a.set(yticks=range(25),yticklabels=degran,xticks=range(7),xticklabels=heat_labels);a.tick_params(length=0)
    y=degran.index('PIGR');a.add_patch(Rectangle((-.5,y-.5),1,1,fill=False,edgecolor='black',lw=.65));a.get_yticklabels()[y].set_fontweight('bold')
    cb=fig.colorbar(im,cax=ax(fig,x,200,77,2,'S5'+side+':colorbar'),orientation='horizontal',ticks=[0,50,100] if side=='E' else [-.1,0,.1]);cb.set_label(clabel,labelpad=0);cb.ax.tick_params(pad=1,length=1)
    payload['panels']['S5'+side]={'genes':degran,'columns':ids,'values':z.to_numpy().tolist(),'source_column':target,'PIGR_B_cell_value':float(z.loc['PIGR','PBMC_B_cell'])}
save(fig,'FigureS5_'+TAG)

# Coverage is already a compact resource field in unchanged TableS5. Add one
# source-data summary including reference-stratum variable-member ranges.
coverage=pd.read_csv(OUT/'results/score_measurement_coverage.csv');rows=[]
for source in source_order:
    g=coverage[coverage.source_label.eq(source)&coverage.score.eq('CRRM388')];assert g.measured_members.nunique()==1
    rows.append(dict(source_label=source,nominal_members=388,common_measured_members=int(g.measured_members.iloc[0]),scored_variable_members_min=int(g.scored_variable_members.min()),scored_variable_members_max=int(g.scored_variable_members.max()),reference_strata=int(g.study_id.nunique()),role='Descriptive measurement coverage; not independent validation count',table_location='TableS5: Resource rows, Measured genes'))
pd.DataFrame(rows).to_csv(OUT/'results'/f'CRRM_measurement_coverage_{TAG}.csv',index=False)
dump(OUT/'audit'/f'plot_payload_{TAG}.json',payload)
print('No statistics recomputed. Seven finite main comparisons retain the original BH8/BH56/BH16 families.',flush=True)
