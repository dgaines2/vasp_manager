"""Generate the code reference pages and navigation."""

from __future__ import annotations

from pathlib import Path

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()

for path in sorted(Path("vasp_manager").rglob("*.py")):
    module_path = path.with_suffix("")
    doc_path = path.with_suffix(".md")
    full_doc_path = Path("reference", doc_path)

    parts = tuple(module_path.parts)

    ignore = ["__init__", "__main__", "tests", "static_files"]
    skip = any(p in ignore for p in parts)
    if skip:
        continue

    nav[parts] = doc_path.as_posix()

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        ident = ".".join(parts)
        fd.write(f"::: {ident}")

    mkdocs_gen_files.set_edit_path(full_doc_path, Path("..") / path)

with mkdocs_gen_files.open(Path("reference", "SUMMARY.md"), "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
