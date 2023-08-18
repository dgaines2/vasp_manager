from setuptools import setup

setup(
    name="vasp_manager",
    version="1.0.3",
    description="A simple package to run and analyze VASP calculations",
    packages=["vasp_manager"],
    url="https://github.com/dgaines2/vasp_manager",
    author="Dale Gaines II",
    author_email="dalegainesii@gmail.com",
    license="MIT",
    include_package_data=True,
    python_requires=">=3.10.0",
    install_requires=[
        "numpy>=1.22",
        "pandas>=1.4",
        "pymatgen>=2022.5",
    ],
    extras_require={"dev": ["black", "isort", "pre-commit", "pytest"]},
)
