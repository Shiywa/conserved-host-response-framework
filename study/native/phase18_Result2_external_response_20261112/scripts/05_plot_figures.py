from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase18_Result2_external_response_20261112')
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.55,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':0,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
UP='#BA603A';DOWN='#387C98';INK='#244B60';LIGHT='#A6B2B9';BG='#E6E6E6'
keys=['Airway24h','Mouse_multistrain_D1','Mouse119849','Mouse140226_D5','Mouse18170','Hamster_lung','Human_D3']
short=['Airway, 24 h','Mouse, day 1','Mouse, GSE119849','Mouse, day 5','Mouse, GSE18170','Hamster lung','Human, day 3']
columns=['Airway','M1','M2','M3','M4','Ham.','Human']
branches=['IFN-alpha','IFN-gamma','TNF / NF-kB','TLR signaling','Neutrophil degranulation','RTK signaling','Microtubule']
blabels=['IFN-α','IFN-γ','TNF/NF-κB','TLR signaling','Neutrophil\ndegranulation','RTK signaling','Microtubule']
primary=pd.read_csv(OUT/'tables/07_primary_response_effects.csv').set_index('contrast_key').loc[keys]
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');functions=pd.read_csv(OUT/'tables/08_primary_functional_effects.csv');components=pd.read_csv(OUT/'tables/09_directional_component_effects.csv');genes=pd.read_csv(OUT/'tables/10_external_gene_effects.csv')
time=pd.read_csv(OUT/'tables/15_airway_time_effects.csv');donors=pd.read_csv(OUT/'tables/16_airway_donor_changes.csv');human=pd.read_csv(OUT/'tables/19_human_challenge_effects.csv');hpairs=pd.read_csv(OUT/'tables/20_human_challenge_donor_changes.csv')
audit=[]
def axmm(f,x,y,w,h):
 W,H=np.asarray(f.get_size_inches())*25.4
 return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def textmm(f,x,y,s,**kw):
 W,H=np.asarray(f.get_size_inches())*25.4
 return f.text(x/W,1-y/H,s,**kw)
def letter(f,x,y,s):textmm(f,x,y,s,fontsize=12,fontweight='bold',va='top')
def clean(ax,zero='x'):
 ax.spines[['top','right']].set_visible(False)
 if zero=='x':ax.axvline(0,color='#ADB4B8',lw=.6,zorder=0)
 if zero=='y':ax.axhline(0,color='#ADB4B8',lw=.6,zorder=0)
def errors(ax,df,y,color=INK,offset=0,markersize=3.3,**kw):
 for i,(_,r) in enumerate(df.iterrows()):
  ax.errorbar(r.effect,y[i]+offset,xerr=[[r.effect-r.CI95_low],[r.CI95_high-r.effect]],fmt='o',color=color,ms=markersize,lw=.8,capsize=1.4,**kw)
def vertical_colorbar(f,im,x,y,h,label):
 cb=f.colorbar(im,cax=axmm(f,x,y,2,h));cb.set_ticks([-2,0,2] if im.norm.vmax==2 else [-1.5,0,1.5]);cb.ax.tick_params(labelsize=7,length=1.5,pad=1);cb.outline.set_linewidth(.4);cb.set_label(label,labelpad=3)
def finish(f,stem,tag):
 textmm(f,207,3,tag,fontsize=12,ha='right',va='top')
 f.canvas.draw();rend=f.canvas.get_renderer();W,H=f.bbox.width,f.bbox.height
 text_objects=[t for t in f.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
 outside=[]
 for t in text_objects:
  b=t.get_window_extent(rend)
  if b.x0<-.2 or b.y0<-.2 or b.x1>W+.2 or b.y1>H+.2:outside.append({'text':t.get_text(),'bbox_px':list(b.extents)})
 sizes=sorted(set(round(t.get_fontsize(),4) for t in text_objects))
 for ext in ['pdf','svg','png']:f.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=600 if ext=='png' else 300)
 audit.append({'figure':stem,'width_mm':f.get_figwidth()*25.4,'height_mm':f.get_figheight()*25.4,'text_sizes_pt':sizes,'all_text_inside_canvas':not outside,'outside_text':outside,'identifier':tag,'identifier_pt':12,'ordinary_pt':7,'panel_letter_pt':12,'single_compact_layout':True})
 plt.close(f)

f=plt.figure(figsize=(210/25.4,193/25.4))
letter(f,2.5,12,'A');ax=axmm(f,48,17,49,44);y=np.arange(7)
for i,(_,r) in enumerate(primary.iterrows()):
 ax.errorbar(r.effect,i,xerr=[[r.effect-r.CI95_low],[r.CI95_high-r.effect]],fmt='o',color=INK,mfc=INK if r.primary_Holm_P<.05 else 'white',ms=3.8,lw=.9,capsize=1.5)
labels=['Human airway, 24 h\nGSE247298; 4 donor pairs','Mouse strains, day 1\nGSE112956; 15 RSV / 14 vehicle','Mouse lung\nGSE119849; 10 / 10','Mouse lung, day 5\nGSE140226; 3 / 3','Mouse lung\nGSE18170; 3 / 3','Hamster lung, day 5\nPRJNA1271255; 4 / 4','Human challenge, day 3\n22 infected / 16 uninfected']
ax.set(yticks=y,yticklabels=labels,ylim=(6.6,-.6),xlim=(-.25,2.12),xticks=[0,.5,1,1.5,2],xlabel='Mean score change (95% CI)');ax.tick_params(axis='y',pad=4);clean(ax)
textmm(f,5,72,'Filled marker: primary Holm P < 0.05',va='top')

letter(f,105,12,'B');ax=axmm(f,126,17,68,44)
order=pd.concat([core[core.frozen_direction.eq(d)].sort_values('display_symbol') for d in ['up','down']]).human_ensembl.tolist()
matrix=genes.pivot(index='human_ensembl',columns='contrast_key',values='aligned_effect').reindex(index=order,columns=keys)
cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad(BG)
im=ax.imshow(np.ma.masked_invalid(matrix.to_numpy()),aspect='auto',interpolation='nearest',cmap=cmap,vmin=-2,vmax=2,rasterized=True)
ax.axhline(279.5,color='white',lw=.7);ax.set(xticks=np.arange(7),xticklabels=columns,yticks=[139.5,333.5],yticklabels=['Up\n280','Down\n108']);ax.tick_params(axis='x',pad=3,length=0);ax.tick_params(axis='y',pad=4);ax.spines[:].set_visible(False)
vertical_colorbar(f,im,197,17,44,'Sign-aligned log effect')

letter(f,2.5,81,'C');ax=axmm(f,30,86,67,35)
v=components.pivot(index='contrast_key',columns='component',values='effect').loc[keys]
ax.barh(y,v['Up contribution'],height=.61,color=UP,label='Up members');ax.barh(y,v['Down contribution'],left=v['Up contribution'],height=.61,color=DOWN,label='Down members')
ax.set(yticks=y,yticklabels=short,ylim=(6.6,-.6),xlim=(0,1.55),xticks=[0,.5,1,1.5],xlabel='Contribution to total score change');clean(ax,zero=None)
ax.legend(handles=[Patch(color=UP,label='Up'),Patch(color=DOWN,label='Down')],loc='lower right',bbox_to_anchor=(1.02,1.03),ncol=2,frameon=False,handlelength=.8,columnspacing=1)

letter(f,105,81,'D');ax=axmm(f,134,86,60,35)
mat=functions.pivot(index='branch',columns='contrast_key',values='effect').reindex(index=branches,columns=keys)
q=functions.pivot(index='branch',columns='contrast_key',values='FDR_BH49').reindex(index=branches,columns=keys)
assert mat.notna().all().all() and q.notna().all().all()
im=ax.imshow(mat,aspect='auto',interpolation='nearest',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-1.8,vcenter=0,vmax=1.8))
for i in range(7):
 for j in range(7):
  if q.iloc[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',color='white' if mat.iloc[i,j]>1.1 else 'black',fontsize=7)
ax.set(yticks=np.arange(7),yticklabels=blabels,xticks=np.arange(7),xticklabels=columns);ax.tick_params(axis='x',pad=3,length=0);ax.tick_params(axis='y',pad=3);ax.spines[:].set_visible(False)
vertical_colorbar(f,im,197,86,35,'Branch score change')
textmm(f,133,134,'* BH-adjusted P < 0.05 (49 tests)',va='top')

letter(f,2.5,142,'E');ax=axmm(f,13,148,49,32)
for d,g in donors[donors.background.eq('Non-CF')].groupby('donor_id'):
 ax.plot(np.r_[0,g.time_h],np.r_[0,g.delta],lw=.55,color=LIGHT,alpha=.8)
t=time[time.background.eq('Non-CF')];ax.errorbar(t.time_h,t.effect,yerr=[t.effect-t.CI95_low,t.CI95_high-t.effect],fmt='o-',ms=3,color=INK,lw=1,capsize=1.6)
ax.set(xticks=[0,24,48,72,120],xlim=(-3,125),xlabel='Hours after inoculation',ylabel='Change from baseline',ylim=(-.8,1.7),yticks=[-.5,0,.5,1,1.5]);clean(ax,'y');textmm(f,15,144,'Non-CF airway; 4 donors',va='bottom')

letter(f,70,142,'F');ax=axmm(f,91,148,46,32);st=pd.read_csv(OUT/'tables/18_mouse_strain_time_effects.csv');strains=['AJ','AKRJ','Balb','C57BL6J','DBA2J'];display=['A/J','AKR/J','BALB/cJ','C57BL/6J','DBA/2J']
assert set(strains)==set(st.strain)
for day,col,off in [(1,UP,-.13),(5,DOWN,.13)]:
 ss=st[st.time_d.eq(day)].set_index('strain').loc[strains];errors(ax,ss,np.arange(5),color=col,offset=off,markersize=2.8)
ax.set(yticks=np.arange(5),yticklabels=display,ylim=(4.6,-.6),xlim=(-1.3,2.9),xticks=[-1,0,1,2],xlabel='Mean score change');clean(ax)
ax.legend(handles=[Line2D([],[],color=UP,marker='o',ms=3,label='Day 1'),Line2D([],[],color=DOWN,marker='o',ms=3,label='Day 5')],loc='lower right',bbox_to_anchor=(1.04,1.02),ncol=2,frameon=False,handlelength=.6,columnspacing=.7)

letter(f,145,142,'G');ax=axmm(f,158,148,46,32);rng=np.random.default_rng(20261112)
for i,group in enumerate(['Uninfected','Infected']):
 z=hpairs[hpairs.time_d.eq(3)&hpairs.eventual_outcome.eq(group)];ax.scatter(i-.12+rng.uniform(-.14,.14,len(z)),z.delta,s=7,c=LIGHT,alpha=.9,edgecolors='none');r=human[human.time_d.eq(3)&human.target.eq(group.lower())].iloc[0]
 ax.errorbar(i+.23,r.effect,yerr=[[r.effect-r.CI95_low],[r.CI95_high-r.effect]],fmt='D',ms=3.7,color=INK,lw=.9,capsize=2)
ax.set(xticks=[0,1],xticklabels=['Uninfected\nn = 16','Infected\nn = 22'],xlim=(-.6,1.65),xlabel='Infection outcome',ylabel='Day 3 − baseline',ylim=(-1.4,1.25),yticks=[-1,0,1]);clean(ax,'y');textmm(f,159,144,'Human challenge',va='bottom')
finish(f,'Figure2_CRRM_external_response','Liu et al. Fig2')

# Supplementary figure: coverage and prespecified sensitivity/context analyses.
f=plt.figure(figsize=(210/25.4,179/25.4))
datasets=[primary.loc[k,'dataset'] for k in keys]+['PRJNA1271255_Blood'];slabels=short+['Hamster blood'];scols=['Airway','M1','M2','M3','M4','Ham.\nlung','Human','Ham.\nblood'];coverage=pd.read_csv(OUT/'tables/04_coverage_summary.csv').set_index('dataset').loc[datasets]
letter(f,2.5,12,'A');ax=axmm(f,31,18,66,38);ys=np.arange(8)
ax.barh(ys,coverage.n_up,height=.65,color=UP);ax.barh(ys,coverage.n_down,left=coverage.n_up,height=.65,color=DOWN)
for i,(_,r) in enumerate(coverage.iterrows()):ax.text(r.n_scored+5,i,str(int(r.n_scored)),va='center')
ax.set(yticks=ys,yticklabels=slabels,ylim=(7.6,-.6),xlim=(0,440),xticks=[0,100,200,300,388],xlabel='Scored members (of 388)');clean(ax,zero=None)
ax.legend(handles=[Patch(color=UP,label='Up'),Patch(color=DOWN,label='Down')],loc='lower right',bbox_to_anchor=(1,1.02),ncol=2,frameon=False,handlelength=.8,columnspacing=1)

letter(f,106,12,'B');ax=axmm(f,141,18,62,38);common=pd.read_csv(OUT/'tables/13_common_member_sensitivity.csv').set_index('contrast_key').loc[keys]
errors(ax,primary,np.arange(7),color=INK,offset=-.12,markersize=3);errors(ax,common,np.arange(7),color=UP,offset=.12,markersize=3)
ax.set(yticks=np.arange(7),yticklabels=short,ylim=(6.6,-.6),xlim=(-.25,2.12),xticks=[0,1,2],xlabel='Mean score change (95% CI)');clean(ax)
ax.legend(handles=[Line2D([],[],color=INK,marker='o',ms=3,label='Full'),Line2D([],[],color=UP,marker='o',ms=3,label='Common 273')],loc='lower right',bbox_to_anchor=(1.02,1.02),ncol=2,frameon=False,handlelength=.8,columnspacing=.8)

letter(f,2.5,69,'C');ax=axmm(f,31,75,62,36);cov=pd.read_csv(OUT/'tables/06_functional_branch_coverage.csv');mat=cov.pivot(index='branch',columns='dataset',values='fraction_measured').reindex(index=branches,columns=datasets)
assert mat.notna().all().all()
im=ax.imshow(mat,aspect='auto',interpolation='nearest',cmap='YlGnBu',vmin=.5,vmax=1)
ax.set(yticks=np.arange(7),yticklabels=blabels,xticks=np.arange(8),xticklabels=scols);ax.tick_params(axis='x',length=0,pad=3);ax.spines[:].set_visible(False)
cb=f.colorbar(im,cax=axmm(f,96,75,2,36));cb.set_ticks([.5,.75,1]);cb.ax.tick_params(labelsize=7,length=1.5,pad=1);cb.outline.set_linewidth(.4)

letter(f,111,69,'D');ax=axmm(f,141,75,62,36)
for i,group in enumerate(['Uninfected','Infected']):
 z=hpairs[hpairs.time_d.eq(7)&hpairs.eventual_outcome.eq(group)];ax.scatter(i-.12+rng.uniform(-.13,.13,len(z)),z.delta,s=10,c=LIGHT,edgecolors='none');r=human[human.time_d.eq(7)&human.target.eq(group.lower())].iloc[0];ax.errorbar(i+.23,r.effect,yerr=[[r.effect-r.CI95_low],[r.CI95_high-r.effect]],fmt='D',ms=3.7,color=INK,lw=.9,capsize=2)
ax.set(xticks=[0,1],xticklabels=['Uninfected\nn = 4','Infected\nn = 9'],xlim=(-.6,1.65),xlabel='Infection outcome',ylabel='Day 7 − baseline',yticks=[-1,0,1,2]);clean(ax,'y')

letter(f,2.5,131,'E');ax=axmm(f,14,137,83,30)
for bg,col in [('Non-CF',INK),('CF',UP)]:
 for d,g in donors[donors.background.eq(bg)].groupby('donor_id'):ax.plot(np.r_[0,g.time_h],np.r_[0,g.delta],lw=.5,color=col,alpha=.2)
 t=time[time.background.eq(bg)];ax.errorbar(t.time_h,t.effect,yerr=[t.effect-t.CI95_low,t.CI95_high-t.effect],fmt='o-',ms=3,color=col,lw=.9,capsize=1.5,label=bg)
ax.set(xticks=[0,12,24,48,72,120],xlim=(-3,125),xlabel='Hours after inoculation',ylabel='Change from baseline',ylim=(-1,3),yticks=[-1,0,1,2,3]);clean(ax,'y');ax.legend(loc='upper left',frameon=False,ncol=2,handlelength=1,columnspacing=1)

letter(f,106,131,'F');ax=axmm(f,141,137,62,30);loo=pd.read_csv(OUT/'tables/11_unit_deletion_influence.csv').groupby('contrast_key').effect.agg(['min','max']).loc[keys]
for i,k in enumerate(keys):ax.plot([loo.loc[k,'min'],loo.loc[k,'max']],[i,i],lw=2.1,color=LIGHT,solid_capstyle='round');ax.plot(primary.loc[k,'effect'],i,'o',ms=3.3,color=INK)
ax.set(yticks=np.arange(7),yticklabels=short,ylim=(6.6,-.6),xlim=(0,1.7),xticks=[0,.5,1,1.5],xlabel='Effect after deleting one unit');clean(ax)
finish(f,'FigureS2_CRRM_external_sensitivity','Liu et al. FigS2')
(OUT/'logs/figure_layout_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
for record in audit:print(record)
