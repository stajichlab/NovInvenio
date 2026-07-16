# NovInvenio

A pipeline for identifying lineage-specific genes (novelty candidates) in fungal genomes by searching for proteins present in an ingroup clade but absent from outgroup species.

## Language

### Core concepts

**Candidate**:
A protein that passes the pairwise-search and presence-matrix filters: present in ≥ N% of ingroup proteomes and absent from all outgroup proteomes. Has not yet been validated against outgroup genomes.
_Avoid_: Novelty (that's the next stage), hit, match

**Novelty candidate**:
A candidate that also passes genomic validation — no significant TBLASTN hit in any outgroup genome. The terminal output of the pipeline (`novelties.<SHORT>.tsv`).
_Avoid_: Novel gene, lineage-specific gene (those are interpretations, not pipeline outputs)

**Paralog**:
The best non-self match of a protein within its own proteome, as detected by self-vs-self search. Used to set a per-protein e-value cutoff: a cross-species hit must be better than the within-proteome paralog to count as presence.
_Avoid_: Homolog, duplicate (those are broader; paralog specifically means within-proteome)

**Ortholog**:
A cross-species protein match that passes the paralog cutoff and competition filters — i.e., one that establishes presence. General biology defines orthologs as genes related by speciation; in this pipeline, orthology is operationalised through the filtering steps.
_Avoid_: Homolog (too broad), hit (too raw)

**Presence**:
A protein from one proteome has a significant, paralog-filtered pairwise hit against some protein in the target proteome. The protein's own source proteome is present by definition. A raw hit that fails the paralog tests does not establish presence. Presence is directional: protein A may be present in proteome Y while Y's best match to A is absent from A's proteome. Asymmetric presence is meaningful and should be detected and noted.
_Avoid_: Hit, match (those are raw pairwise results; presence is post-filter)

### Filtering

**Paralog competition**:
A cross-species hit is disqualified when the query's own paralog hits the same target proteome with a better e-value. The hit is then better explained by the conserved domain shared with the paralog.
_Avoid_: Competition (too generic)

**Paralog cutoff**:
The e-value of the rank-2 self-hit (best within-proteome non-self match). A cross-species hit must have a better (lower) e-value than this cutoff to establish presence. Proteins with no detectable paralog fall back to a global default e-value.
_Avoid_: Threshold, e-value cutoff (those are generic; paralog cutoff is specific to self-vs-self)

### Validation

**Cluster**:
A group of similar candidate proteins, as determined by mmseqs2 clustering. The cluster is the unit of genomic validation: a TBLASTN hit on the representative disqualifies all members. Clustering simplifies the overall pattern of genes and builds consensus.
_Avoid_: Protein family (synonymous but less precise), Group, family (too generic)

**Genomic validation (TBLASTN)**:
Cluster representatives are searched against outgroup genomes via TBLASTN. Any hit disqualifies the entire cluster. This catches genes present in outgroup genomes but missed by proteome annotation.
_Avoid_: TBLASTN search (that's the method; genomic validation is the purpose)

**Cluster representative**:
The mmseqs2-designated representative sequence of a cluster of candidate proteins. TBLASTN validation is performed on representatives; hits are expanded to all cluster members.
_Avoid_: Rep (too informal)

### Organization

**Source proteome**:
The ingroup proteome from which a candidate protein was queried. Every ingroup species is a source proteome. The novelties output is organized by source proteome: `novelties.<SHORT>.tsv` contains only proteins with `source_proteome == <SHORT>`.
_Avoid_: Query proteome (that's the search role; source proteome is the organizational role)

### Annotation

**Conserved domain**:
A protein domain (e.g., Pfam, CAZy, MEROPs) detected at or above a significance threshold. Conserved domains are shared across species and paralogs; the paralog competition filter exists to strip hits that merely detect a conserved domain rather than a true ortholog.
_Avoid_: Domain (too generic), motif (different concept)

**Model organism**:
A species with a curated gene name lookup table, used to annotate candidate proteins with gene names and product descriptions. Does not need to be an ingroup species, though being ingroup allows direct protein ID matching.
_Avoid_: Reference organism (too generic)

**Annotation**:
The process of attaching functional information (gene name, product description, protein domains) to candidate proteins. All annotation sources are recorded; the product description is filled from the first available source in priority order: model organism → Pfam → SwissProt.
_Avoid_: Functional assignment (too generic)

### Taxonomic categories

**Ingroup**:
The monophyletic clade of interest — the species whose lineage-specific genes we're trying to find. All ingroup proteomes are searched against each other and against outgroups.
_Avoid_: Reference set, query set

**Outgroup**:
Species phylogenetically outside the ingroup clade, used as a contrast to determine absence. Outgroup proteomes are searched against (to detect absence) and outgroup genomes are searched via TBLASTN (to validate absence at the DNA level).
_Avoid_: Background, control set

**Taxonomic group**:
A biological classification label applied to a sample (e.g., Pezizomycotina, Taphrinomycotina). The ingroup is typically a single taxonomic group; the outgroup may span several.
_Avoid_: Clade (overloaded), lineage (too vague)
