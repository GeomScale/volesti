# Count Linear Extensions using OrderPolytope

Count the number of linear extensions of a partial order through volume
approximation of the corresponding order polytope.

## Key Identity

```
#LinearExtensions(P) = n! × Vol(OrderPolytope(P))
```
*(Stanley, 1986)*

## Usage

```bash
mkdir build && cd build
cmake ..
make

# Chain of 4 elements (exact = 1 linear extension)
./count_le_orderpolytope ../instances/chain_4.txt cb

# Antichain of 4 elements (exact = 24 linear extensions)
./count_le_orderpolytope ../instances/antichain_4.txt sob 20

# Diamond poset (exact = 3 linear extensions)
./count_le_orderpolytope ../instances/diamond_4.txt cg
```

## Input Format

The poset file has the format:
```
n
i1 j1
i2 j2
...
```
Where `n` is the number of elements and each pair `i j` means `a_i ≤ a_j`.

## Algorithms

- `sob` — Sequence of Balls
- `cg`  — Cooling Gaussians
- `cb`  — Cooling Balls (recommended)
