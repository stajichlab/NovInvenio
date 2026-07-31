# =============================================================================
# NovInvenio — all bioinformatics tools + Python runtime
# Registry: ghcr.io/stajichlab/novinvenio
# Image tag : 0.5.0  (mirrors pixi.toml version)
# Build    : docker build -t ghcr.io/stajichlab/novinvenio:0.5.0 .
# Push     : docker push ghcr.io/stajichlab/novinvenio:0.5.0
# =============================================================================
# This image contains only the tool layer.  The pipeline source (bin/, lib/,
# modules/, workflows/) is supplied by the cloned repository; Nextflow stages
# it into each task working directory automatically.  User data files (proteome
# FASTAs, genomic DNA, annotation databases, config CSVs) are bind-mounted at
# runtime — see README.md "Running with Docker / Singularity".
# =============================================================================

ARG MAMBAFORGE_VERSION=24.11.2
FROM mambaforge/mambaforge:${MAMBAFORGE_VERSION}

# Install all pipeline tools in one pass.
# Channels: conda-forge (primary) + bioconda (bioinformatics tools).
# Strict channel priority ensures deterministic solving.
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict && \
    mamba install --yes --name base \
        python=3.12 \
        hmmer=3.4 \
        famsa \
        diamond=2.2 \
        blast=2.17 \
        mmseqs2=18.8 \
        openmpi=4.1 \
        pandas \
        biopython \
        pyyaml \
        matplotlib

# OpenMPI needs a writable /tmp.
ENV TMPDIR=/tmp
ENV OMPI_MCA_tmpdir_base=/tmp

# Entry point is intentionally unset — Nextflow manages per-process execution.
# All tools are already on PATH via /opt/conda/bin.
