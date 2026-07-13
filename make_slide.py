import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

NAVY="#0f2c4d"; BLUE="#2166ac"; RED="#b2182b"; GREEN="#1a7a3a"; GREY="#4d4d4d"; LGREY="#eceff3"

bs = pd.read_csv("/Users/kherlambang/.claude-science/orgs/985db041-bc74-48ae-9fad-4462234be58a/artifacts/proj_566d88b73186/97e6c5e5-119d-4e07-8204-5281443b1860/v17f884e8_Pancreatic_Beta_Cell_catalog_strict.csv")
es = pd.read_csv("/Users/kherlambang/.claude-science/orgs/985db041-bc74-48ae-9fad-4462234be58a/artifacts/proj_566d88b73186/6dd838bb-90b4-423f-945a-9df0d87c9246/vb71ead3d_Small_Intestinal_Enterocyte_catalog_strict.csv")

fig=plt.figure(figsize=(13.33,7.5)); fig.patch.set_facecolor("white")
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")
def box(x,y,w,h,fc,ec="none",lw=0,r=0.02):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0,rounding_size={r}",fc=fc,ec=ec,lw=lw,transform=ax.transData))

box(0,90.5,100,9.5,NAVY)
ax.text(4,96.7,"Can we computationally find \u2014 and pre-screen for cross-species translatability \u2014",fontsize=14,color="white",fontweight="bold",va="center")
ax.text(4,93.9,"tissue-specific AAV \u201con-switches\u201d for the extrahepatic tissues where AAV actually matters?",fontsize=14,color="white",fontweight="bold",va="center")
ax.text(4,91.6,"RESEARCH TRACK  \u00b7  a scientific question, answered with a purpose-built, reusable method",fontsize=8.3,color="#cdd8e6",va="center",style="italic")

box(3,71,44,16.5,LGREY)
ax.text(5,85.4,"WHY EXTRAHEPATIC, WHY AAV",fontsize=10.5,color=NAVY,fontweight="bold")
for i,t in enumerate(["• LNP already delivers to liver \u2014 AAV's real niche is",
    "   heart, CNS, muscle, retina, pancreas (LNP-hard tissues)",
    "• Cross-species failure is the field's costliest failure mode:",
    "   elements tuned in mouse routinely fail in human/NHP",
    "• We flag likely non-translating parts computationally, up front"]):
    ax.text(5,82.6-i*2.55,t,fontsize=7.8,color=GREY,va="center")

box(3,50.5,44,18.5,LGREY)
ax.text(5,66.7,"ONE METHOD, RUN ACROSS 7 CELL TYPES",fontsize=10.5,color=NAVY,fontweight="bold")
ax.text(5,63.3,"$ python run.py --auto \"heart muscle\"",fontsize=8.4,color=RED,family="monospace")
ax.text(5,60.4,"$ python run.py --auto \"insulin-producing cells\"",fontsize=8.4,color=BLUE,family="monospace")
ax.text(5,56.6,"Code lines changed between any two runs:",fontsize=8.4,color=GREY)
ax.text(41,56.6,"0",fontsize=15,color=GREEN,fontweight="bold",va="center")
ax.text(5,53.2,"Plain-English cell type \u2192 atlas name (Claude); no per-tissue code.",fontsize=7.3,color="#777",style="italic")

box(3,38,44,10.5,LGREY)
ax.text(5,46.2,"OPEN DATA",fontsize=10.5,color=NAVY,fontweight="bold")
ax.text(5,43.1,"CATLAS scATAC (1.14M cCRE × 222 cell types) · JASPAR · GENCODE",fontsize=8,color=GREY)
ax.text(5,40.4,"hg38 · CELLxGENE Census · Ensembl orthologs   |   code: MIT",fontsize=8,color=GREY)

box(50,38,47,49.5,"white",ec="#c9d3df",lw=1.2,r=0.015)
ax.text(52.5,84.5,"THE FUNNEL  (identical thresholds, every cell type)",fontsize=10.5,color=NAVY,fontweight="bold")
stages=[("Union cCREs",1143424),("Accessible",181897),("Specific (\u03c4\u22650.9)",82038),("Size \u2264800 bp",65101),("Ranked shortlist",25)]
x0=53; bmax=26; vmax=1143424; ytop=79
shortlist_w=max(1.6,bmax*(np.log10(25)/np.log10(vmax)))
for i,(lab,val) in enumerate(stages):
    y=ytop-i*5.0; w=max(1.6,bmax*(np.log10(val)/np.log10(vmax))); c=NAVY if i==4 else BLUE
    box(x0,y-1.5,w,3.0,c,r=0.008)
    ax.text(x0-0.5,y,lab,fontsize=8,color=GREY,ha="right",va="center")
    ax.text(x0+bmax+2,y,f"{val:,}",fontsize=9.2,color=c,fontweight="bold",va="center")
y=ytop-25.0
seg=[("conserved",34,GREEN),("risk-flagged",5,"#d73027"),("unknown",11,"#9aa5b1")]
xseg=x0
for lbl,n,c in seg:
    w=shortlist_w*(n/50.0)
    box(xseg,y-1.5,w,3.0,c,r=0.006); xseg+=w
ax.text(x0-0.5,y,"Cross-species labelled",fontsize=8,color=GREY,ha="right",va="center")
ax.text(x0+bmax+2,y,"34 conserved / 5 risk / 11 unknown  (n=50)",fontsize=7.4,color=GREY,va="center")
ax.text(52.5,40.2,"~1.1M candidate elements  \u2192  a ranked, translation-aware shortlist per tissue",fontsize=8,color="#777",style="italic")

cards=[
 ("POSITIVE CONTROL \u2713 (extrahepatic)",GREEN,["Beta \u2192 INS rank 1  \u00b7  Alpha \u2192 GCG  \u00b7  Heart \u2192 RYR2","Gut \u2192 OLFM4/APOA4/PIGR  \u00b7  each recovers identity",
   "Housekeeping (ACTB\u2026) \u2192 \u03c4\u22480.02, filtered out","Alpha config was AI-proposed, held out from curation"]),
 ("MOTIF ENRICHMENT \u2713",BLUE,["Beta shortlist 4.7\u00d7 enriched","for PDX1/NKX6-1/PAX6/NEUROD1",
   "p = 1.6\u00d710\u207b\u00b9\u00b3 (72% vs 12% PDX1)","Identity-TF PWMs from JASPAR, per tissue"]),
 ("CROSS-SPECIES \u2605 headline",RED,["Every part flagged for translatability","Enterocyte PHGR1 \u2192 no mouse ortholog:",
   "human-specific; a mouse-only run misses it","MALRD1 \u2192 divergent mouse expression"]),
]
cw=30; gap=1.5; x=3; ctop=32.5; ch=25
for title,col,lines in cards:
    box(x,ctop-ch,cw,ch,"white",ec=col,lw=1.6,r=0.02)
    box(x,ctop-3.6,cw,3.6,col,r=0.02)
    ax.text(x+1.5,ctop-1.8,title,fontsize=9,color="white",fontweight="bold",va="center")
    for j,ln in enumerate(lines):
        ax.text(x+1.5,ctop-6.5-j*3.9,ln,fontsize=7.9,color=GREY,va="center")
    x+=cw+gap

fig.savefig("summary_slide_strict.png",dpi=200,facecolor="white",bbox_inches="tight")
fig.savefig("summary_slide_strict.pdf",facecolor="white",bbox_inches="tight")