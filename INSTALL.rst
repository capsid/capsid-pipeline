CaPSID Pipeline
=================
CaPSID's pipline is available as a package in the `Python Package Index <http://pypi.python.org/pypi/capsid/>`_.

Dependencies
------------
In order to run CaPSID's pipeline, you will need to install the following dependencies:

- Python_ - Python 2.6+ supported. No support for Python 3 at the moment.
- BioPython_ - Tools for biological computation.
- PyMongo_ - Python module needed for working with MongoDB.
- pysam_ -  Python module for reading and manipulating Samfiles.
- bx-python_ - Python library and associated set of scripts to allow for rapid implementation of genome scale analyses.

.. _Python: http://www.python.org
.. _BioPython: http://biopython.org/wiki/Main_Page
.. _PyMongo: http://api.mongodb.org/python/current/
.. _pysam: http://code.google.com/p/pysam/
.. _bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/Home

For more detailed instructions on how to install the dependencies, pick your operating system from this list:

- [[Installing capsid dependencies on Ubuntu|Installing on Ubuntu]]

Installing with Pip
-------------------
Pip_ is the preferred installation method on platforms other than Windows::

    $ sudo pip install capsid

To get a specific version of capsid::

    $ pip install capsid==1.0

To upgrade using pip::

    $ pip install --upgrade capsid

.. _Pip: http://www.pip-installer.org/en/latest/index.html

Installing with easy_install
----------------------------
If you must install capsid using `setuptools <http://pypi.python.org/pypi/setuptools>`_ do::

    $ easy_install capsid

To upgrade do::

    $ easy_install -U capsid

Building from sources
---------------------
If you'd rather install directly from the source (i.e. to stay on the bleeding edge), check out the latest source from github::

    $ git clone git://github.com/capsid/capsid-pipeline.git capsid
    $ cd capsid
    $ python setup.py install --user   # installs to your home directory -- requires Python >= 2.6
    or
    $ python setup.py build
    $ sudo python setup.py install --prefix=/usr/local   # installs to /usr/local

