# Calculation Modes

We include calculation modes `"rlx-coarse"`, `"rlx"`, `"static"`, `"bulkmod"`,
and `"elastic"`.  The desired modes to calculate are specified when
initializing a `VaspManager` object.

* `rlx-coarse`: lower precision energy-based relaxation
* `rlx`: tighter force-based relaxation
* `static`: high accuracy static SCF calculation
* `bulkmod`: bulk modulus calculation using an Equation of State (EOS) fit to an
energy-volume curve
* `elastic`: Determination of elastic constants using the strain/deformation
method built into `VASP`

I generally recommend starting from `rlx-coarse`, although the functionality is
there to start a `rlx` calculation from the initially provided POSCAR.

Most users' workflows follow `rlx-coarse` &#8594; `rlx` &#8594; `static`. The
modes `static`, `bulkmod`, and `elastic` can all be run independently of each
other. For example, workflows might look like `rlx-coarse` &#8594; `rlx`
&#8594; `static` &#8594; `bulkmod`, or `rlx` &#8594; `elastic`, or simply
`static` or `bulkmod`.  The `elastic` mode requires at least `rlx` preceding it
in order to guarantee converged lattice parameters and atomic positions.

Example workflows are shown below:

=== "Basic"

    ```mermaid
    graph LR
    rlx-coarse --> rlx --> static
    ```

=== "Full"

    ```mermaid
    graph LR
    rlx-coarse --> rlx --> static & bulkmod & elastic
    ```

=== "Elastic"

    ```mermaid
    graph LR
    rlx --> elastic
    ```
