import sys
import numpy as np

# usage: displace_mopac.py <FORCE output> <input .mop> <amplitude/Å> <new .mop>
out, mop, amp, new = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]

# lowest normal mode = first column of the normal coordinate analysis
txt = open(out).read()
lines = txt.split('NORMAL COORDINATE ANALYSIS')[1].split('\n')
rows = [l.split() for l in lines[8:] if l.strip()]
mode = []
for r in rows:
    if not r[0].isdigit():  # end of the first column block
        break
    mode.append(float(r[1]))
mode = np.array(mode).reshape(-1, 3)

# keywords/title lines and Cartesian coordinates of the input (x 1 y 1 z 1 format)
head = open(mop).read().split('\n')[:3]
atoms = [l.split() for l in open(mop).read().split('\n')[3:] if l.strip()]
x = np.array([[float(a[1]), float(a[3]), float(a[5])] for a in atoms])
assert len(x) == len(mode)

# MOPAC reorients the molecule for FORCE; rotate the mode back into the input frame
orient = txt.split('ORIENTATION OF MOLECULE IN FORCE CALCULATION')[1].split('\n')[4:4 + len(x)]
y = np.array([[float(v) for v in l.split()[2:5]] for l in orient])
p, q = y - y.mean(0), x - x.mean(0)  # centred FORCE and input geometries
u, _, vt = np.linalg.svd(p.T @ q)  # Kabsch rotation p -> q
rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt  # proper rotation (no reflection)
print('max alignment error / A:', np.abs(p @ rot - q).max())
x += amp * (mode @ rot) / np.abs(mode).max()  # amp = largest atomic displacement
with open(new, 'w') as f:
    f.write('\n'.join(head) + '\n')
    for a, r in zip(atoms, x):
        f.write(f'{a[0]:2s} {r[0]:14.8f} 1 {r[1]:14.8f} 1 {r[2]:14.8f} 1\n')
