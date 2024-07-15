# dbCSEP

<img src="https://i.postimg.cc/W4HsKdn3/dbcsep-logo.png" width="320"> 


**A database application to store, manage and evaluate seismicity forecasts**

# Table of Contents

* [Installation](#installing-dbcsep)
* [Useful Links](#important-links)
* [Roadmap/Issues](#roadmap-and-known-issues)
* [Contributing](#contributing)
* [License](#license)


# Installing dbCSEP

* Clone the repository and access the top directory:

```shell
git clone https://github.com/cseptesting/dbcsep
cd dbcsep
```

* Using `conda` environment: https://conda.io - check out Anaconda, Miniconda or Miniforge (recommended):

```shell
conda env create -f environment.yml
conda activate dbcsep
pip install .
```

If using Miniforge, replace `conda` with `mamba` for faster installations.

* Using `pip` environments:
```shell
python -m venv venv
source venv/bin/activate
pip install .
```

# Important Links

* [CSEP Website](https://cseptesting.org)
* `pyCSEP` [Github](https://github.com/sceccode/pycsep)
* `pyCSEP` [Documentation](https://docs.cseptesting.org/)

# Roadmap and Known Issues

* Define a roadmap

# Contributing

We encourage all types of contributions, from reporting bugs, suggesting enhancements, adding new features and more. Please refer to the [Contribution Guidelines](https://github.com/cseptesting/dbcsep/blob/main/CONTRIBUTING.md) and the [Code of Conduct](https://github.com/cseptesting/dbcsep/blob/main/CODE_OF_CONDUCT.md) for more information

# License

The `dbcsep` software is distributed under the GNU Affero GPL 3-Clause open-source license. Please see the [license file](https://github.com/cseptesting/dbcsep/blob/main/LICENSE) for more information.