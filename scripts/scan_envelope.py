from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

SCAN = f'{Path(__file__).resolve().parents[1]}/03-scan/'  # scan directory, repo-relative
EH = 2625.4996  # Eh -> kJ/mol


def read_allxyz(name):
    # all frames of an ORCA relaxed scan: list of (energy, atom rows)
    frames = []
    for block in open(f'{SCAN}{name}/{name}.allxyz').read().split('>'):  # frames are separated by '>'
        lines = block.strip().split('\n')
        n = int(lines[0])  # number of atoms
        e = float(lines[1].split()[-1])  # energy is the last word of the comment line
        atoms = [l.split() for l in lines[2:2 + n]]
        frames.append((e, atoms))
    return frames


def xyz(atoms):
    # coordinates as an (n, 3) array
    return np.array([[float(v) for v in a[1:4]] for a in atoms])


def angle(x, a, b, c):
    # bond angle a-b-c in degrees
    u, v = x[a] - x[b], x[c] - x[b]
    return np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))


def dihedral(x, a, b, c, d):
    # dihedral a-b-c-d in degrees (-180, 180]
    b0, b1, b2 = x[a] - x[b], x[c] - x[b], x[d] - x[c]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1  # projections perpendicular to the central bond
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))


def branch(name, coord):
    # (coordinate measured from the geometry, energy, run name, atoms) for every scan point
    return [(round(coord(xyz(a)), 1), e, name, a) for e, a in read_allxyz(name)]


def envelope(points):
    # keep the lowest energy for each coordinate value
    best = {}
    for p in points:
        if p[0] not in best or p[1] < best[p[0]][1]:
            best[p[0]] = p
    return [best[q] for q in sorted(best)]


def write(env, tag, ref):
    # envelope energies (.dat) and geometries (.xyz)
    with open(f'{SCAN}{tag}_envelope.dat', 'w') as f:
        f.write('# coordinate/deg  E-E(trans)/kJmol-1  E/Eh  source\n')
        for q, e, src, _ in env:
            f.write(f'{q:8.1f} {(e - ref) * EH:10.2f} {e:16.9f}  {src}\n')
    with open(f'{SCAN}{tag}_envelope.xyz', 'w') as f:
        for q, e, src, a in env:
            f.write(f'{len(a)}\n{tag} {q:.1f} E {e:.9f} {src}\n' + '\n'.join(' '.join(r) for r in a) + '\n')


torsion = lambda x: abs(dihedral(x, 2, 0, 1, 13))  # CNNC dihedral, C2-N0=N1-C13
cnn_trans = lambda x: angle(x, 2, 0, 1)  # CNN angle C2-N0=N1
cnn_cis = lambda x: 360 - angle(x, 2, 0, 1)  # cis side mapped to 180-235 so the path runs trans -> linear -> cis

# all branches of both scans
tors = {n: branch(n, torsion) for n in ['torsion', 'torsion_from75']}
inv = {n: branch(n, cnn_trans) for n in ['inversion_trans', 'inversion_trans_from160']}
inv['inversion_cis'] = branch('inversion_cis', cnn_cis)

ref = tors['torsion'][0][1]  # trans energy (torsion 180°) is the zero
env_t = envelope(sum(tors.values(), []))[::-1]  # torsion from 180° down to 0°
env_i = envelope(sum(inv.values(), []))
write(env_t, 'torsion', ref)
write(env_i, 'inversion', ref)

plt.rcParams['font.size'] = 11  # larger font for the report
styles = [('#0072B2', 'o', '-'), ('#D55E00', 's', '--'), ('#009E73', '^', ':')]  # Okabe-Ito colour + marker + line
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, branches, env in [(ax1, tors, env_t), (ax2, inv, env_i)]:
    # grey band = lower envelope, lines = individual branches
    ax.plot([p[0] for p in env], [(p[1] - ref) * EH for p in env], color='0.85', lw=8, label='lower envelope')
    for (c, m, ls), (name, pts) in zip(styles, branches.items()):
        pts = sorted(pts)
        ax.plot([p[0] for p in pts], [(p[1] - ref) * EH for p in pts], color=c, marker=m, ls=ls, ms=4, lw=1.5,
                label=name.replace('_', ' '))
    ax.legend(frameon=False, fontsize=9)
ax1.set_xlim(180, 0)
ax1.set_xlabel('CNNC dihedral / °')
ax1.set_ylabel('$E - E$(trans) / kJ mol$^{-1}$')
ax2.set_xticks([120, 140, 160, 180, 200, 220], ['120', '140', '160', '180', '160', '140'])  # cis side labelled as real CNN
ax2.set_xlabel('CNN angle / °  (trans side → linear → cis side)')
fig.tight_layout()
fig.savefig(f'{SCAN}scan_envelope.png', dpi=200)
print(len(env_t), len(env_i))  # number of envelope points (37, 25)
