import glob
import re
from pathlib import Path
import matplotlib.pyplot as plt

R = f'{Path(__file__).resolve().parents[1]}/'  # repo root
EV = 27.211386  # Eh -> eV


def parse(out):
    # first two CIS singlets: (dE/eV, dominant i->a, weight) and oscillator strengths
    txt = open(out).read()
    assert 'ORCA TERMINATED NORMALLY' in txt
    states = []
    for block in re.split(r'\nSTATE\s+\d+:', txt.split('CIS-EXCITED STATES (SINGLETS)')[1])[1:3]:
        de = float(re.search(r'E=\s+([\d.]+) au', block).group(1)) * EV  # excitation energy
        exc = re.findall(r'(\d+)a ->\s+(\d+)a\s+:\s+([\d.]+)', block.split('\n\n')[0])  # i -> a : weight
        i, a, w = max(exc, key=lambda t: float(t[2]))  # dominant excitation
        states.append((de, f'{i}->{a}', float(w)))
    spec = txt.split('ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS')[1]
    f = [float(m) for m in re.findall(r'0-1A\s+->\s+\d-1A\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)', spec)[:2]]
    return states, f


def zindo_s0(out):
    # ZINDO/S energy of the ground state (SCF "Total Energy"), Eh
    return float(re.search(r'\nTotal Energy\s+:\s+(-?[\d.]+) Eh', open(out).read()).group(1))


data = {}
for scan in ['torsion', 'inversion']:
    rows = []
    outs = sorted(glob.glob(f'{R}05-zindo/{scan}/p*.out'))
    ref = zindo_s0(outs[0])  # reference: první bod skenu (trans)
    for out in outs:
        q = float(open(out.replace('.out', '.xyz')).read().split('\n')[1].split()[1])  # coordinate from the xyz comment
        st, f = parse(out)
        rows.append((q, (zindo_s0(out) - ref) * EV, st, f))
    data[scan] = rows
    with open(f'{R}05-zindo/{scan}_zindo.dat', 'w') as fh:
        fh.write('# coord/deg  S0(ZINDO/S)  S1  S2 (eV, rel. to trans; S_n = S0(ZINDO/S) + dE_n(ZINDO/S))  '
                 'dE1  dE2 (eV)  f1  f2  main excitation S1, S2 (weight)\n')
        for q, e0, st, (f1, f2) in rows:
            fh.write(f'{q:7.1f} {e0:8.3f} {e0 + st[0][0]:8.3f} {e0 + st[1][0]:8.3f} {st[0][0]:7.3f} {st[1][0]:7.3f} '
                     f'{f1:7.4f} {f2:7.4f}  {st[0][1]} ({st[0][2]:.2f})  {st[1][1]} ({st[1][2]:.2f})\n')

plt.rcParams['font.size'] = 11  # larger font for the report
styles = [('#0072B2', 'o', '-'), ('#D55E00', 's', '--'), ('#009E73', '^', ':')]  # S0, S1, S2: colour + marker + line
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex='col', sharey='row')
for col, scan in enumerate(['torsion', 'inversion']):
    rows = data[scan]
    q = [r[0] for r in rows]
    top, bottom = axes[0, col], axes[1, col]
    c, m, ls = styles[0]
    top.plot(q, [r[1] for r in rows], color=c, marker=m, ls=ls, ms=4, lw=1.5, label='S$_0$ (ZINDO/S)')
    for k in range(2):
        c, m, ls = styles[k + 1]
        # top: state energies on the S0 scale, bottom: vertical excitation energies alone
        top.plot(q, [r[1] + r[2][k][0] for r in rows], color=c, marker=m, ls=ls, ms=4, lw=1.5,
                 label=f'S$_{k + 1}$ = S$_0$ + Δ$E_{k + 1}$(ZINDO/S)')
        bottom.plot(q, [r[2][k][0] for r in rows], color=c, marker=m, ls=ls, ms=4, lw=1.5,
                    label=f'Δ$E_{k + 1}$ (ZINDO/S)')
    top.legend(frameon=False, fontsize=9)
    bottom.legend(frameon=False, fontsize=9)
axes[0, 0].set_ylabel('$E - E_{S_0}$(trans) / eV')
axes[1, 0].set_ylabel('vertical excitation energy / eV')
axes[1, 0].set_xlim(180, 0)
axes[1, 0].set_xlabel('CNNC dihedral / °')
axes[1, 1].set_xticks([120, 140, 160, 180, 200, 220], ['120', '140', '160', '180', '160', '140'])  # cis side as real CNN
axes[1, 1].set_xlabel('CNN angle / °  (trans side → linear → cis side)')
fig.tight_layout()
fig.savefig(f'{R}05-zindo/zindo_states.png', dpi=200)

# Composite figure: S0 from the ωB97X-D3/def2-TZVP scan envelope, S_n = S0(DFT) + ΔE_n(ZINDO/S).
# The ZINDO/S ground-state energy is not usable along the path (cis 157 eV above trans), see 05-zindo/NOTES.md.
KJ2EV = 1 / 96.485332  # kJ/mol -> eV
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for col, scan in enumerate(['torsion', 'inversion']):
    env = [l.split() for l in open(f'{R}03-scan/{scan}_envelope.dat') if not l.startswith('#')]
    rows = data[scan]
    assert len(env) == len(rows)
    for e, r in zip(env, rows):
        assert abs(float(e[0]) - r[0]) < 0.11 or abs(360 - float(e[0]) - r[0]) < 0.11, (e[0], r[0])
    q = [float(e[0]) for e in env]
    s0 = [float(e[1]) * KJ2EV for e in env]
    ax = axes[col]
    c, m, ls = styles[0]
    ax.plot(q, s0, color=c, marker=m, ls=ls, ms=4, lw=1.5, label='S$_0$ (ωB97X-D3)')
    for k in range(2):
        c, m, ls = styles[k + 1]
        ax.plot(q, [s + r[2][k][0] for s, r in zip(s0, rows)], color=c, marker=m, ls=ls, ms=4, lw=1.5,
                label=f'S$_{k + 1}$ = S$_0$ + Δ$E_{k + 1}$(ZINDO/S)')
    ax.legend(frameon=False, fontsize=9)
axes[0].set_ylabel('$E - E_{S_0}$(trans) / eV')
axes[0].set_xlim(180, 0)
axes[0].set_xlabel('CNNC dihedral / °')
axes[1].set_xticks([120, 140, 160, 180, 200, 220], ['120', '140', '160', '180', '160', '140'])
axes[1].set_xlabel('CNN angle / °  (trans side → linear → cis side)')
fig.tight_layout()
fig.savefig(f'{R}05-zindo/states_dft_zindo.png', dpi=200)
