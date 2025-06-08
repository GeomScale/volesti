import numpy as np
import matplotlib.pyplot as plt

shape = input("Name of the polytope ('cube' or 'simplex'): ").strip().lower()
dim   = int(input("Number of dimensions: ").strip())

# ime fajla za Running
filename = f"build/sb_{shape}_{dim}_run.txt"
data     = np.loadtxt(filename)
d        = data.shape[1]

eps = 1e-1

# === Izbor nasumične fasete ===
if shape == "cube":
    # za svaku koordinatu postoje 2 facete: x_i = +1 ili x_i = -1
    facet_col = np.random.randint(0, d)
    facet_val = np.random.choice([+1.0, -1.0])
    desc      = f"x_{facet_col}≈{facet_val}"
elif shape == "simplex":
    # simpleks ima d “osi” x_i=0 i jednu “sumu” Σx=1
    choice = np.random.randint(0, d+1)
    if choice < d:
        facet_col = choice
        facet_val = 0.0
        desc      = f"x_{facet_col}≈0"
    else:
        facet_col = None   # sumarna faseta
        facet_val = None
        desc      = "Σx≈1"
else:
    raise ValueError("Polytope mora biti 'cube' ili 'simplex'")

# Filtriranje tačaka na izabranoj faseti
if shape == "simplex" and facet_col is None:
    sums = data.sum(axis=1)
    mask = np.abs(sums - 1.0) < eps
else:
    mask = np.abs(data[:, facet_col] - facet_val) < eps

facet_pts = data[mask]

print(f"Loaded {len(data)} points from '{filename}'")
print(f"Random facet: {desc}, found {len(facet_pts)} points")

# Izbor dve ose za projekciju koje nisu konstantne na toj faseti
available = list(range(d))
if facet_col is not None:
    available.remove(facet_col)
axes = np.random.choice(available, size=2, replace=False)
print(f"Projection on axes x_{axes[0]} and x_{axes[1]}")

# Crtanje
plt.figure(figsize=(6,6))
plt.scatter(facet_pts[:, axes[0]], facet_pts[:, axes[1]], s=12, alpha=0.7)
plt.xlabel(f"x_{axes[0]}")
plt.ylabel(f"x_{axes[1]}")
plt.title(f"{shape.capitalize()} {dim}D — facet {desc}")
plt.grid(True)
plt.show()
