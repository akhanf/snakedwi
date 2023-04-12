
## Contributing to snakedwi

Snakedwi python package dependencies are managed with Poetry 
which you will need installed on your machine. You can find installation 
instructions on the [Poetry website](https://python-poetry.org/).

Snakedwi also has a few dependencies outside of Python, including popular 
neuroimaging packages like ANTs, Freesurfer, FSL, MRtrix3, and others. We 
strongly recommend using snakedwi with the `--use-singularity` flag, 
which will pull and use the required containers, unless you are comfortable 
using installing and using all of these tools yourself.

## Setup the development environment

Clone the repository and install all dependencies (including dev) with poetry:

```bash
git clone https://github.com/khanlab/snakedwi.git 
cd snakedwi
poetry install --with dev 
```

Poetry will automatically create a virtual environment. To customize where these 
virtual environments are stored, see the poetry docs here

Then, you can run snakedwi with:

```bash
poetry run snakedwi
```

or you can activate a virtual environment shell and run snakedwi directly:

```bash
poetry shell
snakedwi
```

You can exit the poetry shell with `exit`

## Running and fixing code format quality

Snakedwi uses [poethepoet](https://github.com/nat-n/poethepoet) as a task 
runner. You can see what commands are available by running:

```bash
poetry run poe 
```

We use a few tools, including `black`, `flake8`, `isort`, `snakefmt`, and 
`yamlfix` to ensure formatting and style of our codebase is consistent. There 
are two task runners you can use to check and fix your code, which can be 
invoked with:

```bash
poetry run poe quality_check
poetry run poe quality_fix
```

Note: If you are in a poetry shell, you do not need to prepend poetry run to the
command.

## Dry-run / testing your workflow

Using Snakemake’s dry-run option (`--dry-run`/`-n`) is an easy way to verify any 
changes made to the workflow are working correctly. The test_data folder 
contains example BIDS datasets that are useful for verifying different workflow
conditions. These dry-run tests are part of the automated GitHub actions that 
are run for every commit.

You can invoke the pre-configured task via `poethepoet` to perform a dry-run:

```bash
poetry run poe test
```

This performs a number of tests, involving different scenarios in which a user 
may use snakedwi.
