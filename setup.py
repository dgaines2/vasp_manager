from setuptools import setup

setup(
    name="vasp_manager",
    version="0.1.0",
    description="A simple package to run and analyze vasp calculations",
    packages=["vasp_manager"],
    url="https://github.com/dgaines2/vasp_manager",
    author="Dale Gaines II",
    author_email="dalegainesii@gmail.com",
    license="MIT",
    include_package_data=True,
    install_requires=[
        "numpy>=1.22",
        "pandas>=1.4",
        "pymatgen>=2022.5",
    ],
    extras_require={"dev": ["black", "isort", "pre-commit"]},
)
