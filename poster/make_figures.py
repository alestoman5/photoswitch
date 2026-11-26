"""Figures for the A3 poster (vector SVG + PNG crops). Run from the repo root: python poster/make_figures.py"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import qrcode
import qrcode.image.svg

R = Path(__file__).resolve().parents[1]
OUT = R / 'poster'
BLUE, ORANGE, INK, GREY = '#0072B2', '#D55E00', '#1d2733', '#8a949e'
ORANGE_TXT, GREY_TXT = '#A84300', '#5f6b76'  # darker shades for text (contrast >= 4.5:1 on white)
KJ2EV = 1 / 96.485332
plt.rcParams.update({'font.family': 'Inter', 'font.size': 9, 'axes.edgecolor': INK, 'axes.labelcolor': INK,
                     'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False,
                     'axes.spines.right': False, 'svg.fonttype': 'none', 'axes.linewidth': 0.8})


def rows(path):
    return [l.split() for l in open(path) if not l.startswith('#')]


# 1) main figure: S0 (DFT) and S1, S2 = S0 + dE(ZINDO/S) along torsion and inversion
fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.75), sharey=True, gridspec_kw={'wspace': 0.08})
for ax, scan in zip(axes, ['torsion', 'inversion']):
    env, zin = rows(R / f'03-scan/{scan}_envelope.dat'), rows(R / f'05-zindo/{scan}_zindo.dat')
    q = [float(e[0]) for e in env]
    s0 = [float(e[1]) * KJ2EV for e in env]
    s1 = [a + float(z[4]) for a, z in zip(s0, zin)]
    s2 = [a + float(z[5]) for a, z in zip(s0, zin)]
    for y, c, lab in [(s2, BLUE, 'S$_2$'), (s1, ORANGE, 'S$_1$'), (s0, INK, 'S$_0$')]:
        ax.plot(q, y, color=c, lw=2.2)
        x_lab = q[-1]
        ax.text(x_lab + (-2 if scan == 'torsion' else 2), y[-1], ' ' + lab, color=ORANGE_TXT if c == ORANGE else c, va='center',
                ha='left', fontsize=9, fontweight='semibold')
    ax.set_ylim(-0.2, 6.3)
    ax.grid(axis='y', color='#e3e7eb', lw=0.6)
    ax.set_axisbelow(True)
axes[0].set_xlim(180, 0)
axes[0].set_xticks([180, 135, 90, 45, 0])
axes[0].set_xlabel('torze: úhel CNNC / °')
axes[0].set_title('torze kolem N=N', fontsize=9.5, color=INK, loc='left', fontweight='semibold')
axes[0].set_ylabel('E − E(S$_0$, trans) / eV')
i = [float(e[0]) for e in rows(R / '03-scan/torsion_envelope.dat')].index(90.0)
axes[0].annotate('ΔE$_1$ = 0,08 eV', xy=(90, 2.36), xytext=(178, 1.55), fontsize=8.5, color=INK,
                 arrowprops=dict(arrowstyle='-', color=INK, lw=0.7))
axes[1].set_xlim(112, 250)
axes[1].set_xticks([120, 150, 180, 210, 235], ['120', '150', '180', '150', '125'])
axes[1].set_xlabel('inverze: úhel CNN / °  (trans → lineární → cis)')
axes[1].set_title('inverze na dusíku', fontsize=9.5, color=INK, loc='left', fontweight='semibold')
axes[1].text(180, 0.3, 'mezera S$_1$–S$_0$ ≥ 0,86 eV', ha='center', fontsize=8.5, color=INK)
for a in axes:
    a.text(0.0 if a is axes[0] else 1.0, -0.02, '', transform=a.transAxes)
axes[0].text(178, 0.22, 'trans', fontsize=8, color=GREY_TXT)
axes[0].text(2, 0.78, 'cis', fontsize=8, color=GREY_TXT, ha='right')
fig.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.18)
fig.savefig(OUT / 'fig_states.svg')

# 2) isomerization enthalpy vs experiment
fig, ax = plt.subplots(figsize=(3.4, 2.1))
meth = [('PM6', -16.5, GREY), ('B3LYP/6-31G*', 61.9, GREY), ('ωB97X-D3/def2-TZVP', 50.0, BLUE)]
ax.axvspan(47.4 - 3.2, 47.4 + 3.2, color=ORANGE, alpha=0.18, lw=0)
ax.axvline(47.4, color=ORANGE, lw=1.6)
ax.text(47.4, 2.75, 'experiment\n47,4 ± 3,2', color=ORANGE_TXT, ha='center', va='bottom', linespacing=1.1, bbox=dict(fc='white', ec='none', pad=1.5), fontsize=8.5, fontweight='semibold')
for k, (m, v, c) in enumerate(meth):
    ax.barh(k, v, color=c, height=0.55)
    ax.text(v + 2 if v > 0 else 2, k, f'{v:.1f}'.replace('.', ',').replace('-', '−'),
            va='center', ha='left', fontsize=9, color=INK, fontweight='semibold')
ax.set_yticks(range(3), [m for m, _, _ in meth])
ax.axvline(0, color=INK, lw=0.8)
ax.set_xlim(-32, 75)
ax.set_ylim(-0.5, 3.6)
ax.set_xlabel('ΔH(cis − trans) / kJ mol$^{-1}$')
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
fig.tight_layout()
fig.savefig(OUT / 'fig_isomerization.svg')

# 3) ground-state pathway levels (electronic energy relative to trans)
fig, ax = plt.subplots(figsize=(3.4, 2.5))
levels = [(0, 0.0, 'trans', INK), (5.6, 50.9, 'cis', INK), (2, 166.9, 'TS inverze', BLUE),
          (2, 222.6, 'torze 90°, uzavřené slupky', GREY), (2, 139.6, 'torze 90°, BS singlet', ORANGE),
          (2, 124.2, 'torze 90°, triplet', ORANGE)]
for x, e, lab, c in levels:
    ax.plot([x - 0.55, x + 0.55], [e, e], color=c, lw=2.6, solid_capstyle='butt')
    side = x + 0.65 if x != 4 else x - 0.65
    ax.text(x + 0.7 if x == 2 else x, e + (0 if x == 2 else (-9 if x == 0 else 9)),
            f'{lab}  {e:.0f}'.replace('.', ','), ha='left' if x == 2 else 'center',
            va='center' if x == 2 else ('top' if x == 0 else 'bottom'), fontsize=8, color={GREY: GREY_TXT, ORANGE: ORANGE_TXT}.get(c, c))
for x0, e0, x1, e1, c in [(0.55, 0, 1.45, 166.9, BLUE), (0.55, 0, 1.45, 124.2, ORANGE)]:
    ax.plot([x0, x1], [e0, e1], color=c, lw=0.8, ls=(0, (3, 2)))
ax.set_xlim(-0.8, 6.4)
ax.set_ylim(-25, 245)
ax.set_xticks([])
ax.spines['bottom'].set_visible(False)
ax.set_ylabel('ΔE vůči trans / kJ mol$^{-1}$')
fig.tight_layout()
fig.savefig(OUT / 'fig_pathways.svg')

# 4) orbitals of trans (n, pi, pi*) cropped from 06-orbitals/orbitals.png
im = Image.open(R / '06-orbitals/orbitals.png').convert('RGB')
for k, mo in enumerate([29, 33, 34]):
    im.crop((80, 500 * k + 120, 610, 500 * k + 360)).save(OUT / f'mo{mo}.png')

# 5) QR code to the repository
qrcode.make('https://github.com/alestoman5/photoswitch', image_factory=qrcode.image.svg.SvgPathImage,
            border=0).save(OUT / 'qr.svg')
print('ok')
