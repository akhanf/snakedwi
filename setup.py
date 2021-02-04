import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snakedwi", 
    version="0.1.0",
    author="Ali Khan",
    author_email="alik@robarts.ca",
    description="Snakemake BIDS app for dwi pre-processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/akhanf/snakedwi",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={'console_scripts': [
        'snakedwi=run:main'
    ]},
    install_requires=[
        "snakebids==0.2.1",
        "snakemake==5.28.0",
        "pandas",
        "nibabel",
        "numpy"
    ],
    python_requires='>=3.7'
)
