import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))
from build_master_pool import render_master_pool  # noqa: E402

BFD_SAMPLES_CSV = (
    "ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,"
    "ORDER,FAMILY,GENUS,SPECIES,TRANSL_TABLE,LOCUSTAG\n"
    "GCA_1,Mucor circinelloides,1006PhL,PRJ1,36698,fungi_odb12,Mucoromycota,Mucoromycotina,"
    "Mucoromycetes,,Mucorales,Mucoraceae,Mucor,Mucor circinelloides,1,Mucci\n"
)

REPR_TSV = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Mucor_circinelloides_1006PhL\tMucor circinelloides\tTrue\tMucor_circinelloides_1006PhL\t100.0\tTrue\n"
)


def _make_annotation_dir(tmp_path):
    d = tmp_path / 'genome_annotation' / 'Mucor_circinelloides_1006PhL' / 'predict_results'
    d.mkdir(parents=True)
    (d / 'Mucor_circinelloides_1006PhL.proteins.fa').write_text('>x\nMKV\n')
    (d / 'Mucor_circinelloides_1006PhL.scaffolds.fa').write_text('>x\nACGT\n')
    return tmp_path / 'genome_annotation'


def test_render_master_pool_joins_representative_pick_and_keeps_absolute_paths(tmp_path):
    bfd = tmp_path / 'bfd_samples.csv'
    bfd.write_text(BFD_SAMPLES_CSV)
    repr_tsv = tmp_path / 'repr_assignments.tsv'
    repr_tsv.write_text(REPR_TSV)
    annotation_dir = _make_annotation_dir(tmp_path)

    rows = render_master_pool(bfd, annotation_dir, repr_tsv)

    assert len(rows) == 1
    row = rows[0]
    assert row['Species'] == 'Mucor circinelloides'
    assert row['Strain'] == '1006PhL'
    assert row['NCBI_TaxID'] == '36698'
    assert row['Lineage'] == 'Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Mucoraceae;Mucor'
    assert row['ProteinPath'] == str(annotation_dir / 'Mucor_circinelloides_1006PhL' / 'predict_results' / 'Mucor_circinelloides_1006PhL.proteins.fa')
    assert Path(row['ProteinPath']).is_absolute()


def test_render_master_pool_errors_on_missing_annotation_dir(tmp_path):
    bfd = tmp_path / 'bfd_samples.csv'
    bfd.write_text(BFD_SAMPLES_CSV)
    repr_tsv = tmp_path / 'repr_assignments.tsv'
    repr_tsv.write_text(REPR_TSV)
    empty_annotation_dir = tmp_path / 'no_such_dir'

    with pytest.raises(SystemExit, match='Mucor circinelloides'):
        render_master_pool(bfd, empty_annotation_dir, repr_tsv)
