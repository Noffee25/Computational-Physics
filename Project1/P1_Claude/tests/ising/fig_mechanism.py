"""Mechanism figure: conditional acceptance probability of a Wolff bond test
given the outcomes of the tests made p and q draws earlier (generator property),
small multiples per generator; plus how often those earlier draws fall in the
same cluster (algorithm property)."""
import sys, os, re
import numpy as np
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import plotstyle as ps
import matplotlib.pyplot as plt
ps.apply()
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
txt = open(os.path.join(HERE, "scratch/check_triples.txt")).read()
blocks = txt.split("== ")[1:]
data = {}
for b in blocks:
    cfg = b.splitlines()[0].strip()
    m = re.search(r"kappa3\(I\) =\s*([-\dE.+]+) \+-\s*([\dE.+-]+).*=\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", b)
    data[cfg] = dict(k3=float(m.group(1)), sk3=float(m.group(2)), rr=float(m.group(3)), aa=float(m.group(4)), mix=float(m.group(5)))
padd = 2 - np.sqrt(2)
panels = [("MT19937", "src=7 lags=103,250,0,0", ps.COLORS["MT19937"], "lags (103, 250)"),
          ("R250", "src=5 lags=103,250,0,0", ps.COLORS["R250"], "lags (103, 250)"),
          ("R250 + R250 (mod $2^{32}$)", "src=23 lags=103,250,0,0", ps.SINGLE, "lags (103, 250)"),
          ("additive LFG", "src=24 lags=103,250,0,0", ps.SINGLE, "lags (103, 250)"),
          ("RANDU", "src=1 lags=1,2,0,0", ps.COLORS["RANDU"], "lags (1, 2)")]
fig = plt.figure(figsize=(10, 6.6))
gs = fig.add_gridspec(2, 5, height_ratios=[1, 0.95], hspace=0.6, wspace=0.55)
cats = ["R, R", "R, A", "A, A"]
for i, (title, key, col, lagtxt) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i])
    d = data[key]
    vals = [d["rr"], d["mix"], d["aa"]]
    ax.bar(range(3), vals, color=col, width=0.62)
    ax.axhline(padd, color=ps.EXPECT, ls="--", lw=1.0)
    ax.set_xticks(range(3)); ax.set_xticklabels(cats, fontsize=8.5)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title(title + "\n" + lagtxt, fontsize=9)
    ax.text(1.0, 1.1, f"$\\kappa_3$ = {d['k3']:+.4f}", ha="center", va="center", fontsize=8.5, color=ps.INK2)
    ax.set_xlabel("earlier outcomes", fontsize=8.5)
    if i == 0:
        ax.set_ylabel("P(accept | earlier two tests)")
    else:
        ax.set_yticklabels([])
# bottom: consumption geometry of the Wolff algorithm
def load_diag(fn):
    a = np.loadtxt(os.path.join(HERE, fn), comments="#")
    return a[:, 0], a[:, 1], a[:, 2]
ax = fig.add_subplot(gs[1, :3])
lag, pt, psame = load_diag("runs/diag_lifo.dat")
ax.plot(lag, psame, color=ps.SINGLE, label="Wolff (stack): same cluster")
lag2, pt2, psame2 = load_diag("runs/diag_metr.dat")
ax.axvline(103, color=ps.MUTED, lw=0.8, ls=":"); ax.axvline(250, color=ps.MUTED, lw=0.8, ls=":")
ax.text(103, 1.02, "103", ha="center", fontsize=8, color=ps.INK2); ax.text(250, 1.02, "250", ha="center", fontsize=8, color=ps.INK2)
ax.set_xlabel(r"lag $\ell$ (number of random numbers back)")
ax.set_ylabel(r"P(draw $n-\ell$ in same cluster)")
ax.set_ylim(0, 1.08); ax.set_xlim(0, 400)
ax.set_title("Wolff: draws that share a cluster update (L = 16)", fontsize=10)
# histogram of draws per cluster
ax2 = fig.add_subplot(gs[1, 3:])
H = {}
for line in open(os.path.join(HERE, "runs/diag_lifo.dat")):
    if line.startswith("#H"):
        a = line.split(); H[int(a[1])] = int(a[2])
dd = np.array(sorted(H)); cc = np.array([H[k] for k in dd], float)
w = dd * cc / (dd * cc).sum()
bins = np.arange(0, 520, 10)
ax2.hist(dd, bins=bins, weights=w, color=ps.SINGLE)
ax2.axvline(250, color=ps.MUTED, lw=0.8, ls=":")
ax2.set_xlabel("draws used by one cluster")
ax2.set_ylabel("fraction of all draws")
ax2.set_title("Draws per cluster (draw weighted)", fontsize=10)
fig.savefig  # noqa
print(ps.save(fig, "ising_mechanism"))
