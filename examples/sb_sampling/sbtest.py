import numpy as np
import matplotlib.pyplot as plt

shape = input("Name of the polytope ('cube', 'simplex', 'birkhoff', or anything else): ")\
            .strip()
print(shape)
dim   = int(input("Number of dimensions: ").strip())

filename = f"build/sb_{shape}_{dim}_run.txt"
data     = np.loadtxt(filename)
n, d     = data.shape

eps = 1e-1

if shape == "cube":
    facet_col = np.random.randint(0, d)
    facet_val = np.random.choice([+1.0, -1.0])
    desc      = f"x_{facet_col}≈{facet_val}"
    mask      = np.abs(data[:, facet_col] - facet_val) < eps

elif shape == "simplex":
    choice = np.random.randint(0, d+1)
    if choice < d:
        facet_col = choice
        facet_val = 0.0
        desc      = f"x_{facet_col}≈0"
        mask      = np.abs(data[:, facet_col] - 0.0) < eps
    else:
        facet_col = None
        facet_val = None
        desc      = "Σx≈1"
        sums      = data.sum(axis=1)
        mask      = np.abs(sums - 1.0) < eps

elif shape == "birkhoff":
    facet_col = np.random.randint(0, d)
    facet_val = 0.0
    desc      = f"x_{facet_col}≈0"
    mask      = np.abs(data[:, facet_col] - 0.0) < eps

else:
    facet_col = None
    facet_val = None
    desc      = "full point‐cloud"
    mask      = np.ones(n, dtype=bool)

facet_pts = data[mask]

print(f"Loaded {n} points from '{filename}'")
print(f"Selected slice: {desc}, found {len(facet_pts)} points")

available = list(range(d))
if facet_col is not None and facet_col in available:
    available.remove(facet_col)
i, j = np.random.choice(available, size=2, replace=False)
print(f"Projecting on coords x_{i} vs x_{j}")

plt.figure(figsize=(6,6))
plt.scatter(facet_pts[:, i], facet_pts[:, j], s=12, alpha=0.7)
plt.xlabel(f"x_{i}")
plt.ylabel(f"x_{j}")
plt.title(f"{shape.capitalize()} {dim}D — {desc}")
plt.grid(True)
plt.tight_layout()
plt.show()
