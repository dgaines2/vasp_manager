# DFT+U Configuration

When `"hubbards": "wang"` is set, we apply hubbard corrections to certain
transition metals in the presence of oxygen.

!!! Warning
    Only PBE (`"gga": "PE"`) is supported with DFT+U in `vasp_manager`

We use [`LDAUTYPE=2`](https://www.vasp.at/wiki/index.php/Category:DFT%2BU): The
simplified (rotationally invariant) approach to DFT+U, introduced by Dudarev et
al [^1] and U-J coefficients from Wang et al [^2]. If no elements needing DFT+U
are found in the structure, DFT+U is automatically not used.

In the generated INCAR, this looks like:
```
# Example for FeO
#= DFT+U =#
LDAU = .TRUE.
LDAUPRINT = 1
LDAUU = 4.0 0.0
LDAUJ = 0.0 0.0
LDAUL = 2 -1
LMAXMIX = 4
```

The DFT+U elements and their settings can be found in
[`hubbards.json`](https://github.com/dgaines2/vasp_manager/blob/main/vasp_manager/static_files/hubbards.json).

```title="hubbards.json"
--8<-- "vasp_manager/static_files/hubbards.json"
```

[^1]:
	https://doi.org/10.1103/PhysRevB.57.1505
[^2]:
	https://doi.org/10.1103/PhysRevB.73.195107
