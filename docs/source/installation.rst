Installation
------------
ASTER is available on PyPI here_ and can be installed via::

    pip install epiaster


You can also install ASTER from GitHub via::

    git clone git://github.com/BioX-NKU/ASTER.git
    cd ASTER
    python setup.py install

The dependencies will be automatically installed along with ASTER.

Anaconda
~~~~~~~~
If you do not have a working installation of Python 3.8.13 (or later), consider installing Miniconda_ (see `Installing Miniconda`_). Once Anaconda has been installed, you can create and activate a Python 3.8.13 environment via::

    conda create -n py38 python=3.8.13
    conda activate py38

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~
After downloading Miniconda_, in a unix shell (Linux, Mac), run::

    cd DOWNLOAD_DIR
    chmod +x Miniconda3-latest-VERSION.sh
    ./Miniconda3-latest-VERSION.sh

and accept all suggestions.
Either reopen a new terminal or `source ~/.bashrc` on Linux/ `source ~/.bash_profile` on Mac.
The whole process takes just a couple of minutes.

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _here: https://pypi.org/project/epiaster/

