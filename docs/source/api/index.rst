.. module:: epiaster
.. automodule:: epiaster
   :noindex:
API
====


Import ASTER::

   import epiaster as aster
   
ensemble
------------
.. module::epiaster.ensemble
.. currentmodule::epiaster

.. autosummary::
    :toctree: .

    ensemble.estimate_k

estimators
-----------
.. module::epiaster.estimators
.. currentmodule::epiaster


.. autosummary::
    :toctree: .

    estimators.ssd_knee_est
    estimators.davies_bouldin_est
    estimators.silhouette_est
    
utils
----------
.. module::epiaster.utils
.. currentmodule::epiaster

.. autosummary::
    :toctree: .
    
    utils.setup_seed
    utils.getNClusters
