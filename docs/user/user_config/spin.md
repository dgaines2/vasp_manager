# Spin Configuration

Spin polarization is controlled by the `"ispin"` tag in `calc_config.json`.

To include spin polarization, set `"ispin": "auto"` in `calc_config.json`;
otherwise set `"ispin": 1`. With this setting, all elements with valence *d* or
*f* electrons will start with initial magnetic moments of 5 and 7 $\mu_B$,
respectively.

See more about `ISPIN` in VASP: [ISPIN](https://www.vasp.at/wiki/index.php/ISPIN).

The elements considered to be $d$ and $f$ block can be found in
[`d_f_block.json`](https://github.com/dgaines2/vasp_manager/blob/main/vasp_manager/static_files/d_f_block.json).

??? info "d- and f-block elements"

    ```title="d_f_block.json"
    --8<-- "vasp_manager/static_files/d_f_block.json"
    ```

[`VaspManager`][vasp_manager.VaspManager] also accepts an additional argument
`magmom_per_atom_cutoff` which defaults to `0.0`, although I often use a value
around `0.1` to avoid numerical noise and only use spin polarization in further
calculations if they really need it. If this argument is passed, `rlx`
calculations that finish with a magmom per atom less than this value with be
re-run without spin polarization. This argument only affects `rlx`
calculations, and the spin setting for following `static`, `bulkmod`, or
`elastic` calculations is inferred from the final `rlx` calculation.
