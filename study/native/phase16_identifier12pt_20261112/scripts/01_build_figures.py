from pathlib import Path
import json,os,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/crrm_phase16_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors,font_manager,patches
import numpy as np,pandas as pd

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase16_identifier12pt_20261112')
REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
P7=REV/'phase07_A4_7pt_figures';P6=REV/'phase06_external_validation_figure1'
font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'legend.title_fontsize':7,'figure.titlesize':7,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.formatter.use_mathtext':False,'savefig.facecolor':'white'})
UP='#B54C48';DOWN='#347CA6';TEAL='#398979';GOLD='#C99744';GREY='#8B969B';CMAP=colors.LinearSegmentedColormap.from_list('fc',[DOWN,'#FFFFFF',UP])
W=210;H=174
def read(name):return pd.read_csv(OUT/'tables'/name)
atlas=read('01_all388_gene_function_atlas.csv');bp=read('17_Figure1_GO_terms.csv');pr=read('33_Figure1_interpretive_pathways.csv');diag=read('34_Figure1_reason_examples.csv');decomp=read('15_IFN_component_summary.csv');cov=read('05_annotation_coverage.csv');cc=read('19_FigureS1_MF_CC_terms.csv');matrix=read('20_FigureS1_all388_function_matrix.csv')
keys=[('GSE231666','Y0'),('GSE231666','Y1'),('GSE270463','All'),('GSE309353','All')]
build=[]
def newfig(height):
    global H
    H=height
    return plt.figure(figsize=(W/25.4,H/25.4))
def axmm(fig,x,top,w,h,label):return fig.add_axes([x/W,(H-top-h)/H,w/W,h/H],label=label)
def txt(fig,x,top,s,**kw):return fig.text(x/W,(H-top)/H,s,fontsize=7,va='top',**kw)
def tag(fig,l,x,top):return fig.text(x/W,(H-top)/H,l,fontsize=12,fontweight='bold',va='top')
def clean(ax):ax.spines[['top','right']].set_visible(False);ax.tick_params(pad=2)
def wr(s,width=36):return '\n'.join(textwrap.wrap(s,width=width,break_long_words=False))
def save(fig,stem):
    config=json.loads((OUT/'parameters/figure_identifier.json').read_text())
    number='S1' if stem.startswith('FigureS1') else '1'
    identifier=config['format'].format(surname=config['surname'],number=number)
    fig.text(config['right_mm']/W,(H-config['top_mm'])/H,identifier,ha='right',va='top',fontfamily=config['font'],fontsize=config['font_pt'])
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    axes=[dict(label=ax.get_label(),position=list(ax.get_position().bounds),tight_bbox_mm=[v/fig.dpi*25.4 for v in ax.get_tightbbox(renderer).bounds]) for ax in fig.axes]
    for ext in ['pdf','svg','png']:fig.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=600)
    build.append(dict(figure=stem,artwork_mm=[210,H],ordinary_pt=7,panel_letter_pt=12,identifier=identifier,identifier_font_pt=config['font_pt'],single_layout=True,axes=axes))
    plt.close(fig)

fig=newfig(172);tag(fig,'A',8,3)
ax=axmm(fig,8,9,128,26,'A');ax.set(xlim=(-1,101),ylim=(-1,51));ax.axis('off')
labels=['Human A549\nGSE231788; 48 h\n3 Mock / 3 RSV','Mouse lung\nGSE197868; 12 h\n3 control / 3 RSV','Cotton-rat lung\nGSE297903\nPBS / D1; n = 24']
for x,label in zip([0,35,70],labels):
    ax.add_patch(patches.FancyBboxPatch((x,32),30,17,boxstyle='round,pad=.35',facecolor='#EDF3F6',edgecolor='#7895A6',lw=.5));ax.text(x+15,40.5,label,ha='center',va='center',linespacing=1.08)
    ax.annotate('',xy=(50,28),xytext=(x+15,31),arrowprops=dict(arrowstyle='->',lw=.5,color='#7895A6'))
ax.add_patch(patches.FancyBboxPatch((4,15),92,13,boxstyle='round,pad=.35',facecolor='#F3F4F2',edgecolor=GREY,lw=.5));ax.text(50,21.5,'Native DE; unique human mapping\nFDR < 0.05; |log2FC| > 0.58; same direction',ha='center',va='center',linespacing=1.1)
ax.annotate('',xy=(50,10.5),xytext=(50,14.5),arrowprops=dict(arrowstyle='->',lw=.5,color='#7895A6'))
ax.add_patch(patches.FancyBboxPatch((6,0),88,10,boxstyle='round,pad=.35',facecolor='#E2EEF4',edgecolor=DOWN,lw=.7));ax.text(50,5,'3 / 3 models: 388 genes; 280 up / 108 down',ha='center',va='center',linespacing=1.04)
txt(fig,8,36.5,'Broad pool: 1,980 genes supported by at least 2 models')

tag(fig,'B',150,3);ax=axmm(fig,167,10,35,24,'B')
am=atlas.sort_values(['frozen_direction','human_ensembl'],ascending=[False,True]);hm=am[['A549_log2FC','mouse_log2FC','cotton_rat_log2FC']].to_numpy();ax.imshow(hm,aspect='auto',interpolation='nearest',cmap=CMAP,vmin=-4,vmax=4);ax.axhline(279.5,color='white',lw=.6);ax.set_xticks([0,1,2],['A549','Mouse','Cotton rat']);ax.xaxis.tick_top();ax.set_yticks([139.5,333.5],['280 up','108 down']);ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)
cax=axmm(fig,167,37,35,1.5,'B_colorbar');cb=fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(-4,4),cmap=CMAP),cax=cax,orientation='horizontal');cb.set_ticks([-4,0,4]);cb.outline.set_linewidth(.4);txt(fig,149,36,'log2FC')

tag(fig,'C',8,43);ax=axmm(fig,64,50,38,30,'C');yy=np.arange(len(bp))[::-1]
ax.scatter(bp.fold_enrichment,yy,s=10+bp.k_overlap.to_numpy()*.55,c=[UP if v=='up' else DOWN for v in bp.direction],edgecolors='white',linewidths=.3)
ax.set_yticks(yy,[wr(t,36) for t in bp.term]);ax.set_ylim(-.8,len(bp)-.2);ax.set_xlim(0,max(bp.fold_enrichment)*1.24);ax.set_xlabel('Fold enrichment (GO BP)',labelpad=1.5);clean(ax);ax.tick_params(axis='y',length=0)
ax.plot([],[],marker='o',linestyle='',color=UP,label='Up');ax.plot([],[],marker='o',linestyle='',color=DOWN,label='Down');ax.legend(loc='lower right',bbox_to_anchor=(1.02,1.035),ncol=2,frameon=False,handletextpad=.3,columnspacing=.7,borderpad=0)

tag(fig,'D',111,43);ax=axmm(fig,165,50,28,30,'D');norm=colors.Normalize(0,50);cm=plt.get_cmap('YlOrRd');yy=np.arange(len(pr))[::-1]
for y,(_,r) in zip(yy,pr.iterrows()):
    sig=r.FDR_global_BH<.05;ax.scatter(r.fold_enrichment,y,s=12+r.k_overlap*.7,c=[cm(norm(min(50,-np.log10(r.FDR_global_BH))))] if sig else ['white'],edgecolors='#898989' if not sig else 'white',linewidths=.6 if not sig else .3,zorder=3)
ax.set_yticks(yy,[wr(s,33) for s in pr.display_label]);ax.set_ylim(-.8,len(pr)-.2);ax.set_xlim(0,max(pr.fold_enrichment)*1.18);ax.set_xticks([0,5,10,15]);ax.set_xlabel('Fold enrichment',labelpad=1.5);ax.tick_params(axis='y',length=0);clean(ax)
handles=[ax.scatter([],[],s=12+n*.7,facecolor='#B0B0B0',edgecolors='white',linewidths=.3,label=str(n)) for n in [5,20,60]]
ax.legend(handles=handles,title='Genes',loc='lower right',bbox_to_anchor=(1.02,1.025),ncol=3,frameon=False,handletextpad=.2,columnspacing=.5,labelspacing=.1,borderpad=0)
cax=axmm(fig,198,50.5,1.5,29,'D_colorbar');cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),cax=cax);cb.set_ticks([0,25,50]);cb.ax.yaxis.set_label_position('left');cb.set_label('−log10 FDR',labelpad=2);cb.outline.set_linewidth(.4)

# The eleven-gene diagnostic matrix is taller;F/G now share the right column.
tag(fig,'E',8,91)
codes=['UP','DOWN','EFFECT','FDR','FILTER','MAP','NO_ANNOT'];palette=[TEAL,DOWN,'#D8A24E','#D8DDE1','#EDF1F4','#B49ACB','#FFFFFF'];glyph={'UP':'↑','DOWN':'↓','EFFECT':'Δ','FDR':'q','FILTER':'F','MAP':'M','NO_ANNOT':'?'}
data=np.array([[codes.index(r[fam+'_evidence_code']) for fam in ['A549','mouse','cotton_rat']] for _,r in diag.iterrows()]);ax=axmm(fig,55,102,43,43,'E');ax.imshow(data,aspect='auto',interpolation='nearest',cmap=colors.ListedColormap(palette),vmin=-.5,vmax=len(codes)-.5)
ax.set_yticks(range(len(diag)),diag.gene);ax.set_xticks([0,1,2],['A549','Mouse','Cotton rat']);ax.xaxis.tick_top();ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)
for y,r in diag.iterrows():
    for x,fam in enumerate(['A549','mouse','cotton_rat']):
        code=r[fam+'_evidence_code'];ax.text(x,y,glyph[code],ha='center',va='center',color='white' if code in ['UP','DOWN'] else '#202020',fontsize=7)
for sep in [1.5,4.5,5.5,8.5]:ax.axhline(sep,color='white',lw=.8)
for y,label in [(0.5,'Shared'),(3,'Mapping /\ncoverage'),(5,'Effect gate'),(7,'Expression\nfilter'),(9.5,'Opposite')]:ax.text(-2.55,y,label,ha='center',va='center',clip_on=False)
ma=axmm(fig,104,102,4.5,43,'E_members');ma.imshow(diag.CRRM_member.astype(int).to_numpy()[:,None],cmap=colors.ListedColormap(['#EEEEEE',TEAL]),vmin=0,vmax=1,aspect='auto');ma.set_xticks([0],['CRRM'],rotation=45,ha='left');ma.xaxis.tick_top();ma.set_yticks([]);ma.tick_params(length=0,pad=1);ma.spines[:].set_visible(False)
for x,label,color in [(8,'↑ Up',TEAL),(28,'↓ Down',DOWN),(54,'Δ Effect','#9F7027'),(82,'q FDR','#5F6B72')]:txt(fig,x,148,label,color=color)
for x,label,color in [(8,'F Filtered','#5F6B72'),(42,'M Mapping','#806397'),(78,'? Coverage','#5F6B72')]:txt(fig,x,153,label,color=color)

tag(fig,'F',119,91);txt(fig,130,92,'79 of 388 in IFN reference; 309 outside')
ax=axmm(fig,137,100,65,10,'F');ax.barh([1,0],[78,1],color=GOLD,height=.5,label='IFN overlap');ax.barh([1,0],[202,107],left=[78,1],color=TEAL,height=.5,label='Outside IFN reference')
ax.text(39,1,'78',ha='center',va='center',color='white');ax.text(179,1,'202',ha='center',va='center',color='white');ax.text(54,0,'107',ha='center',va='center',color='white');ax.text(111,0,'+ 1',va='center',ha='left');ax.set_yticks([1,0],['Up','Down']);ax.set_xlim(0,295);ax.set_xticks([0,100,200,280]);ax.set_xlabel('Genes',labelpad=1);clean(ax)
fh,fl=ax.get_legend_handles_labels();fig.legend(fh,fl,loc='upper right',bbox_to_anchor=(202/W,1-117/H),ncol=2,frameon=False,handlelength=1,columnspacing=.7,borderpad=0,handletextpad=.4)

tag(fig,'G',119,123);ax=axmm(fig,160,133,42,18,'G')
labels=['Nasal 96 h, Y0; n = 11','Nasal 96 h, Y1; n = 11','Lung ALI 72 h; n = 3','Nasal ALI 20 h; n = 3']
for j,(component,color,label) in enumerate([('CRRM388',DOWN,'CRRM 388'),('Outside_IFN_reference',TEAL,'Outside IFN reference')]):
    for y,(cohort,stratum) in zip(np.arange(4)[::-1],keys):
        z=decomp[(decomp.cohort==cohort)&(decomp.stratum==stratum)&(decomp.component==component)].iloc[0];ax.errorbar(z.mean_delta,y+(j-.5)*.22,xerr=[[z.mean_delta-z.CI95_low],[z.CI95_high-z.mean_delta]],fmt='o',markersize=3,color=color,lw=.8,capsize=1.5,label=label if y==3 else None)
ax.axvline(0,ls='--',lw=.5,color=GREY);ax.set_yticks(np.arange(4)[::-1],labels);ax.set_ylim(-.5,3.5);ax.set_xlim(-.08,1.7);ax.set_xticks([0,.5,1,1.5]);ax.set_xlabel('Mean paired signed-score change\n95% CI',labelpad=1.5);clean(ax)
gh,gl=ax.get_legend_handles_labels();fig.legend(gh,gl,loc='upper right',bbox_to_anchor=(202/W,1-125/H),ncol=2,frameon=False,handlelength=1,columnspacing=1,borderpad=0,handletextpad=.4)
save(fig,'Figure1_CRRM_7pt')

fig=newfig(172);tag(fig,'A',8,3);pp=pd.read_csv(P7/'tables/11_Figure1_external_paired_scores.csv');levels=pp.panel.drop_duplicates().tolist()
for i,l in enumerate(levels):
    ax=axmm(fig,20+i*48,14,37,21,'A'+str(i));x=pp[pp.panel==l]
    for unit,s in x.groupby('unit_id'):
        v=s.set_index('condition').score;ax.plot([0,1],[v['mock'],v['RSV']],color='#ADB5BA',lw=.5,zorder=1);ax.scatter([0,1],[v['mock'],v['RSV']],c=['#697782',DOWN],s=8,zorder=2)
    ax.set_xticks([0,1],['Mock','RSV']);ax.set_xlim(-.25,1.25);ax.set_ylabel('Signed score' if i==0 else '');clean(ax);ax.set_title(['Nasal 96 h, Y0\nn = 11 pairs','Nasal 96 h, Y1\nn = 11 pairs','Lung ALI 72 h\nn = 3 donors','Nasal ALI 20 h\nn = 3 donors'][i],pad=3)

tag(fig,'B',8,42);ax=axmm(fig,38,49,61,28,'B');yy=np.arange(len(cov))[::-1];ax.barh(yy,cov.members_annotated,color=TEAL,height=.55);ax.set_yticks(yy,['GO BP','GO MF','GO CC','Reactome 98','Hallmark']);ax.set_xlim(0,435);ax.set_xticks([0,100,200,300,388]);ax.set_xlabel('Annotated CRRM members (of 388)',labelpad=1.5);clean(ax)
for y,n in zip(yy,cov.members_annotated):ax.text(n+5,y,str(n),va='center')
tag(fig,'C',111,42);ax=axmm(fig,173,49,29,28,'C');yy=np.arange(len(cc))[::-1];ax.scatter(cc.fold_enrichment,yy,s=14+cc.k_overlap*.6,c=[UP if d=='up' else DOWN for d in cc.direction],edgecolor='white',lw=.3);ax.set_yticks(yy,cc.display_label);ax.set_ylim(-.7,len(cc)-.3);ax.set_xlim(0,cc.fold_enrichment.max()*1.2);ax.set_xticks([0,5,10]);ax.set_xlabel('Fold enrichment',labelpad=1.5);clean(ax);ax.tick_params(axis='y',length=0)

tag(fig,'D',8,90);ax=axmm(fig,28,100,74,37,'D');ax.imshow(matrix.iloc[:,2:].to_numpy(),aspect='auto',interpolation='nearest',cmap=colors.ListedColormap(['white',TEAL]),vmin=0,vmax=1);ax.axhline(279.5,color=UP,lw=.7);ax.set_yticks([139.5,333.5],['280 up','108 down']);ax.set_xticks(range(matrix.shape[1]-2),matrix.columns[2:],rotation=43,ha='right');ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)

tag(fig,'E',119,90);ee=pd.read_csv(P7/'tables/12_all_external_effects.csv');ax=axmm(fig,160,101,42,20,'E')
for j,(geneset,color,label) in enumerate([('priority388',DOWN,'CRRM 388'),('broad1980',GOLD,'Broad 1,980')]):
    for y,(cohort,stratum) in zip(np.arange(4)[::-1],keys):
        r=ee[(ee.cohort==cohort)&(ee.subset==stratum)&(ee['set']==geneset)].iloc[0];ax.errorbar(r.mean_delta,y+(j-.5)*.22,xerr=[[r.mean_delta-r.CI95_low],[r.CI95_high-r.mean_delta]],fmt='o',markersize=3,color=color,lw=.8,capsize=1.4,label=label if y==3 else None)
ax.axvline(0,color=GREY,ls='--',lw=.5);ax.set_yticks(np.arange(4)[::-1],['Nasal 96 h, Y0','Nasal 96 h, Y1','Lung ALI 72 h','Nasal ALI 20 h']);ax.set_ylim(-.5,3.5);ax.set_xlim(0,1.7);ax.set_xlabel('Paired score change; 95% CI',labelpad=1.5);clean(ax);ax.legend(loc='lower right',bbox_to_anchor=(1,1.025),ncol=1,frameon=False,handlelength=1,columnspacing=.8,borderpad=0,labelspacing=.2)

tag(fig,'F',119,130);nest=pd.read_csv(P6/'tables/05_external_donor_pairs.csv');nest=nest[nest.cohort=='GSE309353']
for i,geneset in enumerate(['priority388','broad1980']):
    ax=axmm(fig,130+i*43,139,29,18,'F'+str(i));z=nest[nest['set']==geneset]
    for unit,s in z.groupby('donor_id'):
        v=s.set_index('method').delta;ax.plot([0,1],[v['primary'],v['culture_profile_mean']],color='#ADB5BA',lw=.6);ax.scatter([0,1],[v['primary'],v['culture_profile_mean']],color=DOWN,s=11)
    ax.axhline(0,color=GREY,ls='--',lw=.5);ax.set_xticks([0,1],['Sum counts\nthen score','Score profiles\nthen average']);ax.set_xlim(-.25,1.25);ax.set_ylabel('Paired change' if i==0 else '',labelpad=1.5);ax.set_title('CRRM 388' if i==0 else 'Broad 1,980',pad=3);clean(ax)
save(fig,'FigureS1_CRRM_7pt')
(OUT/'logs/figure_build.json').write_text(json.dumps(dict(layout_only=True,figures=build,ordinary_pt=7,panel_letter_pt=12),ensure_ascii=False,indent=2)+'\n')
print('Single compact figure layout with top-right author/number identifiers built;source values unchanged.')
