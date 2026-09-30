"""Produces results/calibration_curve.json - the compute-response of the repaired
operator on the calibration cell. Written because an earlier version of that file
was hand-assembled with no generating script, which CLAUDE.md rule 7 forbids."""

import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from search.edit4 import Edit4Graph, all_subgroups_s4, word_transforms
from search.edit4_orbit_ls import build_instance, ils

g=Edit4Graph(6)
triv=next(s for s in all_subgroups_s4() if len(s)==1)
usable,adj,w,_=build_instance(g, word_transforms(triv,6,False))
print("calibration cell (trivial = unrestricted problem), REPAIRED operator")
print("budget(core-min)  best   seeds")
out=[]
for seconds,nseeds in [(120.0,3),(300.0,3),(300.0,4),(900.0,4)]:
    best=0
    for s in range(1,nseeds+1):
        _,wt,_=ils(w,adj,seed=s,time_limit=seconds)
        best=max(best,wt)
    cm=seconds*nseeds/60.0
    print(f"{cm:>14.0f}  {best:>5}   {nseeds}x{seconds:.0f}s", flush=True)
    out.append({"core_minutes":cm,"best":best,"seeds":nseeds,"seconds":seconds})
json.dump({"cell":"trivial (unrestricted)","operator":"repaired","points":out,
           "broken_operator_reference":{"core_minutes":6.0,"best":99},
           "repo_ils_reference":{"core_minutes":150.0,"best":111},
           "kamis_reference":{"best":113}},
          open("results/calibration_curve.json","w"),indent=2)
