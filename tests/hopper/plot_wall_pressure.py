"""Wall pressure of the Couto silo runs against Couto et al. (2013) and
Zheng and Yu (2015) Fig. 11 (WP H.7). WORKDIR is a test_hopper_pressure.sh
work folder holding couto_mc_aligned and/or couto_dp_aligned.

    python3 plot_wall_pressure.py WORKDIR OUT.png "title"
"""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ["CEL_CELL"]="0.025"
import hopper_pressure_checks as hp, hopper_case as hc
c=dict(hc.COUTO); rows=hp.table()
root=sys.argv[1]; out=sys.argv[2]; title=sys.argv[3]
runs=[("MC", f"{root}/couto_mc_aligned", "#c0392b"), ("DP plane strain", f"{root}/couto_dp_aligned", "#1f6fb4")]
runs=[r for r in runs if os.path.isdir(r[1])]
def prof(p, name):
    np.savetxt(f"{p}/cells.dat", hp.read_cells(p, name))
    cyl, hop = hc.wall_pressures(p, c, 0.025)
    # hopper: average per height band
    zs=np.unique(np.round(hop[:,0],3)); hp_=np.array([[z, hop[np.abs(hop[:,0]-z)<1e-3,1].mean()] for z in zs])
    return np.vstack([hp_, cyl])
cells={n:z for n,z in hp.GA.items()}; cells["GA0"]=0.41
fem=np.loadtxt(os.path.join(hp.ROOT,'docs','spec','hopper-data','zheng2015_fig11_fem.csv'),delimiter=',',comments=('#','z_m'))
fig,axs=plt.subplots(1,2,figsize=(12,7.5),dpi=150,sharey=True)
for ax,phase,label in ((axs[0],"static","At rest (after filling, t = 0.5 s)"),(axs[1],"discharge","Discharge (mean over 0.25-1.0 s after opening)")):
    for k,(n,z) in enumerate(cells.items()):
        lo,hi=hp.exp_range(rows,phase,n)
        ax.plot([lo,hi],[z,z],color='k',lw=5,alpha=0.35,solid_capstyle='butt',label='Couto et al. 2013 (range of 7 tests)' if k==0 else None)
        ax.text(hi+0.4,z,n,va='center',fontsize=8,color='#444')
    for name,p,col in runs:
        if phase=="static":
            pr=prof(p,"cells_0.5"); ax.plot(pr[:,1],pr[:,0],color=col,lw=1.6,label=f"celFoam {name}")
        else:
            ws=[w for w in hp.WRITES if hp.AVG[0]<=w<=hp.AVG[1]]
            P=[prof(p,f"cells_{w:g}") for w in ws]
            m=min(len(q) for q in P); A=np.array([q[:m,1] for q in P]); z=P[0][:m,0]
            ax.fill_betweenx(z,A.min(0),A.max(0),color=col,alpha=0.12,lw=0)
            ax.plot(A.mean(0),z,color=col,lw=1.6,label=f"celFoam {name} (band: min-max)")
    if phase=="discharge":
        ax.plot(fem[:,1],fem[:,0],'--',color='#2a9d55',lw=1.4,label='Zheng & Yu 2015 Abaqus CEL, 0.5 s')
    ax.axhline(0.4764,color='#999',lw=0.6,ls=':'); ax.text(0.2,0.50,'transition',fontsize=7,color='#777')
    ax.set_title(label,fontsize=11); ax.set_xlabel('normal wall pressure (kPa)'); ax.grid(alpha=0.25)
    ax.legend(fontsize=8,loc='upper right')
axs[0].set_ylabel('height above outlet z (m)'); axs[0].set_ylim(0,2.5); axs[0].set_xlim(0,16); axs[1].set_xlim(0,35)
fig.suptitle(title,fontsize=12); fig.tight_layout(); fig.savefig(out); print("ok")
