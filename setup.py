from setuptools import find_packages, setup

setup(
    name="vasp_manager",
    version="1.0.0",
    description="A simple package to run and analyze vasp calculations",
    packages=["vasp_manager"],
    url="https://github.com/dgaines2/vasp_manager",
    author="Dale Gaines II",
    author_email="dalegainesii@gmail.com",
    license="MIT",
    install_requires=[
        "numpy",
        "pandas",
        "pymatgen",
    ],
)
