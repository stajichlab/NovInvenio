NovInvenio
==========

NovInvenio identifies lineage-specific genes: proteins present in a defined
ingroup (>= N% of members) but absent from all outgroup proteomes. It uses
pairwise protein searches (phmmer/diamond/blast), self-vs-self paralog
calibration, mmseqs2 clustering, family-profile HMM search, TBLASTN
validation against outgroup genomes, functional annotation (Pfam +
SwissProt), and model-organism gene-name lookup to produce per-species
novelty candidate tables.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   overview
   method_description
   hexA_filtering
   make_samples_intermediate

.. toctree::
   :maxdepth: 1
   :caption: Architecture decisions

   adr/0001-cluster-candidates-before-tblastn
   adr/0002-family-profile-search-pathway

.. toctree::
   :maxdepth: 1
   :caption: Agent workflows

   agents/domain
   agents/issue-tracker
   agents/triage-labels

Indices
=======

* :ref:`genindex`
* :ref:`search`
