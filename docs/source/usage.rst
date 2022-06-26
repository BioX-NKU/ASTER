Usage Principles
----------------

Input
^^^^^^^^
- **adata**: AnnData object of shape ``n_obs`` × ``n_vars``. Rows correspond to cells and columns to genes.
- **search_list**: List of optional numbers of cell types for the estimation.

Output
^^^^^^^^
- **estimated_k**: Estimated number of cell types.

ASTER can also be seamlessly integrated with `EpiScanpy <https://episcanpy.readthedocs.io/en/stable/>`_, a widely-used Python library for epigenomics single cell analysis::

    import episcanpy.api as epi
    import epiaster as aster
    # Load the single-cell chromatin accessibility data as an AnnData object (adata)
    # Run ASTER
    estimated_k = aster.estimate_k(adata, search_list)


AnnData
^^^^^^^^^
ASTER supports :mod:`episcanpy` and :mod:`anndata`, which provides the :class:`~anndata.AnnData` class.

.. image:: http://falexwolf.de/img/scanpy/anndata.svg
   :width: 300px

At the most basic level, an :class:`~anndata.AnnData` object `adata` stores
a data matrix `adata.X`, annotation of observations
`adata.obs` and variables `adata.var` as `pd.DataFrame` and unstructured
annotation `adata.uns` as `dict`. Names of observations and
variables can be accessed via `adata.obs_names` and `adata.var_names`,
respectively. :class:`~anndata.AnnData` objects can be sliced like
dataframes, for example, `adata_subset = adata[:, list_of_gene_names]`.
For more, see this `blog post`_.

.. _blog post: http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/

To read a data file to an :class:`~anndata.AnnData` object, call::

    import episcanpy.api as epi
    adata = epi.read(filename)

to initialize an :class:`~anndata.AnnData` object. Possibly add further annotation using, e.g., `pd.read_csv`::

    import pandas as pd
    anno = pd.read_csv(filename_sample_annotation)
    adata.obs['cell_groups'] = anno['cell_groups']  # categorical annotation of type pandas.Categorical
    adata.obs['time'] = anno['time']                # numerical annotation of type float
    # alternatively, you could also set the whole dataframe
    # adata.obs = anno

To write, use::

    adata.write(filename)
    adata.write_csvs(filename)
    adata.write_loom(filename)


.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/