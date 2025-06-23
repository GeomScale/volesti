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

# Uniformity test for the facet 
axis = np.random.choice([k for k in range(d) if k != facet_col])
vals = facet_pts[:, axis]
cent = centroid[axis]
R = np.max(np.abs(vals - cent))


Ss = np.linspace(0, 1, 1000)
emp = [(np.abs(vals - cent) <= s * R).mean() for s in Ss]

if shape == "cube":
    F_theo=Ss
elif shape== "simplex":
    m = d-1
    F_theo = 1 - (1 - Ss)**m 

plt.figure(figsize=(6,6))
plt.plot(Ss, emp,  label="Empirical")
plt.plot(Ss, F_theo,   '--', label="Uniform")
plt.xlabel("x")
plt.ylabel("F(x)")
plt.title(f"Uniformity test for the {shape}-facet {desc}, axes x_{axis}")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# Uniformity test for the polytope 
filename = f"build/sb_{shape}_{dim}_tvals.txt"
t = np.loadtxt(filename)   
n = len(t)

t_sorted = np.sort(t)
F_emp = np.arange(1, n+1)/n


plt.figure(figsize=(6,4))
plt.step(t_sorted, F_emp, where='post', label="Boundary sampling")
plt.plot([0,1],[0,1],'--', label="Interior sampling")
plt.xlabel("x$")
plt.ylabel("F(x)")
plt.title("Uniformity test for the polytope ")
plt.legend()
plt.grid(True)
plt.show()