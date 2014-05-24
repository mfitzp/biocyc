
BioCyc Web API for Python
=========================

This package is a Python interface to the BioCyc Web API. While
incomplete the API offers access to most basic attributes for
metabolites, proteins, reactions, pathways and organisms in the
database. The Python interface comes with an disk-based caching
mechanism under ``~/.biocyc`` that greatly reduces the delay (and load)
for BioCyc servers. The interface supports multiple + configurable
caches, so it's possible to share the cache across multiple machines.

Basic initialisation
--------------------

Import the ``biocyc`` object from the ``biocyc`` module. This object
provides the base access to the database for the initial get. You can
set the organism using ``set_organism`` and one of the standard BioCyc
database identifiers. Note that this only affects the organism-database
used for direct requests on the biocyc object. Sub-requests on existing
objects will use the same database as that object (otherwise things
would be very confusing indeed).

.. code:: python

    import os
    from biocyc import biocyc
    os.environ['http_proxy'] = '' # Set your proxy if neccessary
    biocyc.set_organism('meta')

Making a request
----------------

To get an database object (of any type) simply using the unique BioCyc
identifiers for it. Here we request ``L-Lactate``. Note that if you do
this from within an IP[y] Notebook you get a nice table output of all
associated attributes for an object. This includes direct links to the
BioCyc database and other database annotations.

.. code:: python

    o=biocyc.get('L-LACTATE')
    o



==================  ===============================================================================================================================================================================================================================================================================================================================================================================================================
Name                (:sup:`S`)-lactate
BioCyc ID           `L-LACTATE <http://www.biocyc.org/META/NEW-IMAGE?object=L-LACTATE>`__
Org ID              META
Synonyms            L-lactate, L(+)-lactate
INCHI               InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1/t2-/m0/s1
Molecular weight    89.071
Gibbs 0             -72.55646
Parents             L-2-hydroxyacids, Lactate
Reactions           TRANS-RXN-104, RXN-12165, RXN-12096, LACTALDDEHYDROG-RXN, RXN0-5269, D-LACTATE-2-SULFATASE-RXN, TRANS-RXN-104, L-LACTDEHYDROGFMN-RXN, LACTATE-MALATE-TRANSHYDROGENASE-RXN, LACTATE-2-MONOOXYGENASE-RXN, L-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN, L-LACTATE-DEHYDROGENASE-RXN, RXN-9067, RXN-8076, PROPIONLACT-RXN, LACTATE-RACEMASE-RXN, LACTATE-ALDOLASE-RXN
Database links      CAS: `79-33-4 <http://www.commonchemistry.org/ChemicalDetail.aspx?ref=79-33-4>`__, PUBCHEM: `5460161 <http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=5460161>`__, LIGAND-CPD: `C00186 <http://www.genome.ad.jp/dbget-bin/www_bget?C00186>`__, CHEMSPIDER: `4573803 <http://www.chemspider.com/4573803>`__, CHEBI: `16651 <http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16651>`__, BIGG: 34179
==================  ===============================================================================================================================================================================================================================================================================================================================================================================================================



Exploring further
-----------------

Now we have an object we can perform sub-queries by accessing fields. If
you access the ``o.reactions`` field you will trigger a dynamic request
for all entities in that list. Connections to the BioCyc server are
throttled at 1/second, so this may take a little while on long lists.
However, retrieved data is cached under ``~/.biocyc`` so subsequent
requests will be much quicker. By default the cache is set to expire
objects after ~6 months, and the cache folder can be shared between
multiple machines.

*Note: If you just want access to the identifiers, you can use the
``o._reactions`` field to access these without triggering a request*

.. code:: python

    r = o.reactions
    r[0]


==================  ==============================================================================
BioCyc ID           `TRANS-RXN-104 <http://www.biocyc.org/META/NEW-IMAGE?object=TRANS-RXN-104>`__
Org ID              META
Parents             Small-Molecule-Reactions, TR-12
==================  ==============================================================================


.. code:: python

    r[1]



==================  ======================================================================
Name                NADP :sup:`+` L-lactaldehyde dehydrogenase
BioCyc ID           `RXN-12165 <http://www.biocyc.org/META/NEW-IMAGE?object=RXN-12165>`__
Org ID              META
Parents             Chemical-Reactions, Small-Molecule-Reactions
Pathways            PWY-6713
==================  ======================================================================


You can access sub-entities and manipulate objects using standard Python
list processing.

.. code:: python

    ps = [r.pathways for r in o.reactions]
    p = [p for sl in ps for p in sl]
    p



.. parsed-literal::

    [L-rhamnose degradation II,
     L-rhamnose degradation III,
     L-rhamnose degradation II,
     methylglyoxal degradation V,
     lactate biosynthesis (archaea),
     L-lactaldehyde degradation (aerobic),
     L-lactaldehyde degradation (aerobic),
     methylglyoxal degradation V,
     pyruvate fermentation to lactate,
     glucose and xylose degradation,
     Bifidobacterium shunt,
     heterolactic fermentation,
     factor 420 biosynthesis]



.. code:: python

    p[0]


==================  ====================================================================
Name                L-rhamnose degradation II
BioCyc ID           `PWY-6713 <http://www.biocyc.org/META/NEW-IMAGE?object=PWY-6713>`__
Org ID              META
Synonyms            aldolase pathway
Parents             L-rhamnose-Degradation
Species             TAX-5580, ORG-6176, TAX-95486, TAX-284592, TAX-322104
Taxonomic range     TAX-2, TAX-4751
==================  ====================================================================



Finally
-------

That's all for now! Hopefully this shows how Python (and IPython
notebook) access to the BioCyc Web API may be useful. Support for
additional attributes, API calls etc. is planned for the future. If you
have specific requests, get in touch!
