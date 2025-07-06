# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2025 Vissarion Fisikopoulos
# Copyright (c) 2018-2025 Apostolos Chalkis
# Copyright (c) 2025-2025 Iva Janković

# Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import matplotlib.pyplot as plt

shape = input("Name of the polytope ('cube', 'simplex', 'birkhoff', or anything else): ").strip()
dim   = int(input("Number of dimensions: ").strip())

# Loading the data
filename = f"build/sb_{shape}_{dim}_run.txt"
data     = np.loadtxt(filename)
n, d     = data.shape
eps = 1e-7

# Finding the facet
if shape == "cube":
    facet_col = np.random.randint(0, d)
    facet_val = np.random.choice([+1.0, -1.0])
    desc      = f"x_{facet_col}≈{facet_val}"
    mask      = np.abs(data[:, facet_col] - facet_val) < eps

elif shape in ("birkhoff", "simplex"):
    facet_col = np.random.randint(0, d)
    facet_val = 0.0
    desc      = f"x_{facet_col}≈0"
    mask      = np.abs(data[:, facet_col] - facet_val) < eps

else:
    facet_col = None
    desc      = "full cloud"
    mask      = np.ones(n, dtype=bool)

facet_pts = data[mask]

# "Center" of the facet
centroid = facet_pts.mean(axis=0)

# Projection on two random axes
axes = [i for i in range(d) if i != facet_col] or list(range(d))
i, j = np.random.choice(axes, size=2, replace=False)

# Scatter of points in the facet

plt.figure(figsize=(6,6))
plt.scatter(facet_pts[:, i], facet_pts[:, j], s=12, alpha=0.7, label="Samples")
plt.scatter(centroid[i], centroid[j], c='red', marker='*', s=150, label='Center')
plt.xlabel(f"x_{i}")
plt.ylabel(f"x_{j}")
plt.title(f"{shape.capitalize()} {dim}D — facet {desc} with {len(facet_pts)} points ")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Scaling ratio testing
coverage = {}
filename = f"build/sb_{shape}_{dim}_coverage.txt"
with open(filename) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] != 'Facet':
            continue
        facet = int(parts[1])
        tail = ' '.join(parts[4:])
        pairs = [p.strip().rstrip(',') for p in tail.split(',') if p.strip()]
        xs, covs = [], []
        for p in pairs:
            x_str, cov_str = p.split(':', 1)
            x = float(x_str)
            cov = float(cov_str)
            xs.append(x)
            covs.append(cov)
        coverage[facet] = (np.array(xs), np.array(covs))

plt.figure(figsize=(8,6))
plt.grid(True) 
for facet, (xs, covs) in coverage.items():
    plt.plot(xs, covs, label=f'Facet {facet}')

plt.xlabel(r'$x^{\,d}$')
plt.ylabel('Coverage ratio')
plt.grid(True) 
plt.title('Scaling coverage per facet (vs $x^{d}$)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()