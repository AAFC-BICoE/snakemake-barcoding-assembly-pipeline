==============
Bold Retriever
==============

|Pypi index| |Build Status| |Cover alls| |Dependencies status| |supported-versions|

This script accepts FASTA files containing COI sequences. It queries the BOLD
database http://boldsystems.org/ in order to get the taxa identification
based on the sequences.
 
Run this way
------------
* clone repository::

    git clone https://github.com/carlosp420/bold_retriever.git

* install dependencies (python2.7)::

    cd bold_retriever
    pip install -r requirements.txt

* run software

You have to choose one of the databases available from BOLD
http://www.boldsystems.org/index.php/resources/api?type=idengine
and enter it as argument:

* COX1_SPECIES
* COX1
* COX1_SPECIES_PUBLIC
* COX1_L640bp

For example::

    python bold_retriever.py -f ZA2013-0565.fasta -db COX1_SPECIES

* output::

    OtuID   bold_id       similarity  division  class       order       family        species                collection_country
    OTU_99  FBNE064-11    1           animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius pini        Germany
    OTU_99  NEUFI079-11   1           animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius pini        Finland
    OTU_99  FBNE172-13    0.9937      animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius atrifrons   Germany
    OTU_99  FBNE162-13    0.9936      animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius contumax    Austria
    OTU_99  TTSOW138-09   0.9811      animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius ovalis      Canada
    OTU_99  CNPAH380-13   0.9811      animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius             Canada
    OTU_99  CNKOF1602-14  0.9811      animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius pinidumus   Canada
    OTU_99  NRAS173-11    0.9748      animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius conjunctus  Canada
    OTU_99  SSBAE2911-13  0.9748      animal    Collembola  None        None          Collembola             Canada
    OTU_99  CNPAQ117-13   0.9686      animal    Insecta     Neuroptera  Hemerobiidae  Hemerobius humulinus   Canada

Speed
-----
**bold_retriever** uses the library Twisted for performing asynchronous calls.
This speeds up the total processing time.

Citation
--------
The citation should be our MolEco paper:

Vesterinen, E. J., Ruokolainen, L., Wahlberg, N., Peña, C., Roslin, T., Laine, V. N., Vasko, V., Sääksjärvi, I. E., Norrdahl, K., and Lilley, T. M. (2016) What you need is what you eat? Prey selection by the bat Myotis daubentonii. Molecular Ecology, 25(7), 1581–1594. doi:10.1111/mec.13564


Full documentation
------------------
See the full documentation at http://bold-retriever.readthedocs.org

.. |Pypi index| image:: https://badge.fury.io/py/bold_retriever.svg
   :target: http://badge.fury.io/py/bold_retriever
.. |Build Status| image:: https://travis-ci.org/carlosp420/bold_retriever.png?branch=master
   :target: https://travis-ci.org/carlosp420/bold_retriever
.. |Cover alls| image:: https://img.shields.io/coveralls/carlosp420/bold_retriever.svg
   :target: https://coveralls.io/r/carlosp420/bold_retriever?branch=master
.. |Dependencies status| image:: https://gemnasium.com/carlosp420/bold_retriever.svg
   :target: https://gemnasium.com/carlosp420/bold_retriever
.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/bold_retriever.svg?style=flat
   :alt: Supported versions
   :target: https://pypi.python.org/pypi/bold_retriever

