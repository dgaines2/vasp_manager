# Installation

## Install python

1. Install python version $\geq$ 3.10.
    - Preferably, you should use some sort of environment manager like
      [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2)
      and create a new environment.

## Install vasp_manager

The simplest way is to use pip and install directly from
[PyPi](https://pypi.org/project/vasp-manager/#description).
To install, run:
```bash
pip install vasp-manager
```

Alternatively, start by cloning the github repository:
```bash
git clone https://github.com/dgaines2/vasp_manager.git
```

And then cd into the `vasp_manager/` directory and install locally with pip:

=== "Local"
    ```bash
    pip install .
    ```

=== "Editable"
    ```bash
    pip install -e .
    ```

=== "Editable for Development"
    ```bash
    pip install -e '.[dev]'
    ```
