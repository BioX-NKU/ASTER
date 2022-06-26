[![PyPI](https://img.shields.io/pypi/v/epiaster.svg)](https://pypi.org/project/epiaster)
[![Documentation Status](https://readthedocs.org/projects/aster/badge/?version=latest)](https://aster.readthedocs.io/en/latest/?badge=stable)
[![Downloads](https://pepy.tech/badge/epiaster)](https://pepy.tech/project/epiaster)

# ASTER: accurately estimating the number of cell types in single-cell chromatin accessibility data

<div align=center>
<img src = "docs/source/logo.png" width = 65% height = 65%>
</div>  

## Installation
ASTER is available on PyPI [here](https://pypi.org/project/epiaster/) and can be installed via

```
pip install epiaster
```

You can also install ASTER from GitHub via
```
git clone git://github.com/BioX-NKU/ASTER.git
cd ASTER
python setup.py install
```
The dependencies will be automatically installed along with ASTER.


## Quick Start

### Input

* **adata**:       AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
* **search_list**: List of optional numbers of cell types for the estimation.

### Output

* **estimated_k**: Estimated number of cell types.

ASTER can also be seamlessly integrated with [EpiScanpy](https://episcanpy.readthedocs.io/en/stable/), a widely-used Python library for epigenomics single cell analysis:
```python
    import episcanpy.api as epi
    import epiaster as aster
    # Load the single-cell chromatin accessibility data as an AnnData object (adata)
    # Run ASTER
    estimated_k = aster.estimate_k(adata, search_list)
```

### Find more details on [the Documentation of ASTER](https://aster.readthedocs.io/en/latest/index.html).