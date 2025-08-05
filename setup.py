from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="simspace",
    version="0.2.4",
    description="SimSpace: a comprehensive in-silico spatial omics data simulation framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Tianxiao Zhao",
    packages=find_packages(include=["simspace", "simspace.*"]),
    include_package_data=True,
    package_data={
        "simspace": [
            "R/*.R",
            "R/.Rprofile",
            "R/renv.lock",
            "R/renv/activate.R",
        ],
    },
    install_requires=[
        'numpy >=2.2.6',
        'scipy >=1.15.2',
        'pandas >=2.3.1',
        'matplotlib >=3.10.3',
        'seaborn >=0.13.2',
        'scikit-learn >=1.7.0',
        'colorcet >=3.1.0',
        'esda >=2.7.1',
        'libpysal >=4.12.1'
    ],
    extras_require={
        "dev": ["pytest"]
    },
    python_requires=">=3.10",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)