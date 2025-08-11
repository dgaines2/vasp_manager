# Pseudopotential Configuration

We use VASP's [recommended
pseudopotentials](https://www.vasp.at/wiki/index.php/Available_pseudopotentials)
for each element, specifically those from potpaw.54.

This is stored in the
[`pot_dict.json`](https://github.com/dgaines2/vasp_manager/blob/main/vasp_manager/static_files/pot_dict.json)
file.

??? info "Recommended Pseudopotentials"
    ```title="pot_dict.json"
    --8<-- "vasp_manager/static_files/pot_dict.json"
    ```

!!! warning "Warning"
    This is technically an implementation detail, but `pot_dict.json` is a
    symbolic link to `pot_dict_54.json` in the same folder. We also provide
    `pot_dict_52.json` which contains the recommended pseudopotentials for
    potpaw.52.

    To change the symbolic link for potpaw.52, run the following commands in the
    [`static_files`](https://github.com/dgaines2/vasp_manager/blob/main/vasp_manager/static_files/)
    directory:
    ```bash
    rm pot_dict.json
    ln -s pot_dict_52.json pot_dict.json
    ```

    We also include a legacy version of pseudopotential choices for potpaw.52
    to match the OQMD in `pot_dict_oqmd.json`.
