from setuptools import setup

with open("README.md") as file:
    long_description = file.read()

setup(
    name="vasp_manager",
    version="1.1.0",
    description="A simple package to run and analyze VASP calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
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
        "importlib_resources",
    ],
    extras_require={"dev": ["black", "coverage", "isort", "pre-commit", "pytest"]},
)
