"""
Run IgBLAST and output the AIRR-formatted result table
"""

# The code in here was copied from IgDiscover 0.15, https://github.com/NBISweden/IgDiscover-legacy
# and streamlined a bit to work stand-alone. Relevant files:
# - src/igdiscover/cli/igblastwrap.py
# - src/igdiscover/igblast.py
# - src/igdiscover/species.py
# - src/igdiscover/utils.py
# - src/igdiscover/dna.py

# Inline dependencies. Makes this script runnable with, e.g., 'pipx run'
# /// script
# dependencies = ["dnaio"]
# ///

import csv
import errno
import logging
import multiprocessing
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from argparse import ArgumentParser
from contextlib import ExitStack
from dataclasses import dataclass
from io import StringIO
from itertools import islice
from pathlib import Path
from typing import Dict, Optional

import dnaio

GENETIC_CODE = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    # 'TAA': stop
    "TAC": "Y",
    # 'TAG': stop,
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    # 'TGA': stop
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
}

logger = logging.getLogger(__name__)


def nt_to_aa(s, _get=GENETIC_CODE.get):
    """Translate a nucleotide sequence to an amino acid sequence"""
    return "".join([_get(s[i : i + 3], "*") for i in range(0, len(s), 3)])


class SerialPool:
    """
    An alternative to multiprocessing.Pool that runs things in serial for
    easier debugging
    """

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def imap(self, func, iterable, chunksize):
        for i in iterable:
            yield func(i)


def available_cpu_count():
    """
    Number of available virtual or physical CPUs on this system
    """
    if sys.platform != "linux":
        try:
            import multiprocessing

            return multiprocessing.cpu_count()
        except (ImportError, NotImplementedError):
            return 1
    return len(os.sched_getaffinity(0))


def escape_shell_command(command):
    return " ".join(shlex.quote(arg) for arg in command)


def run_igblast(sequences, blastdb_dir, species, sequence_type, penalty=None) -> str:
    """
    Run the igblastn command-line program.

    sequences -- list of Sequence objects
    blastdb_dir -- directory that contains BLAST databases. Files in that
    directory must be databases created by the makeblastdb program and have
    names V, D, and J.

    Return IgBLAST’s raw output as a string.
    """
    if sequence_type not in ("Ig", "TCR"):
        raise ValueError('sequence_type must be "Ig" or "TCR"')
    variable_arguments = []
    for gene in "V", "D", "J":
        variable_arguments += [
            "-germline_db_{gene}".format(gene=gene),
            os.path.join(blastdb_dir, "{gene}".format(gene=gene)),
        ]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        env = None
        # Work around internal_data/ not allowed to contain spaces on Windows
        if sys.platform == "win32" and "IGDATA" not in os.environ:
            if copy_internal_data(tmp_path):
                env = {**os.environ, "IGDATA": tmpdir}

        # An empty .aux suppresses a warning from IgBLAST. /dev/null does not work.
        empty_aux_path = tmp_path / "empty.aux"
        empty_aux_path.write_text("\n")
        arguments = []
        if penalty is not None:
            arguments += ["-penalty", str(penalty)]
        if species is not None:
            arguments += ["-organism", species]

        arguments += [
            "-auxiliary_data",
            str(empty_aux_path),
            "-ig_seqtype",
            sequence_type,
            "-num_threads",
            "1",
            "-domain_system",
            "imgt",
            "-num_alignments_V",
            "1",
            "-num_alignments_D",
            "1",
            "-num_alignments_J",
            "1",
            "-outfmt",
            "19",  # AIRR format
            "-query",
            "-",
        ]
        fasta_str = "".join(">{}\n{}\n".format(r.name, r.sequence) for r in sequences)

        # For some reason, it has become unreliable to let IgBLAST 1.10 write its result
        # to standard output using "-out -". The data becomes corrupt. This does not occur
        # when calling igblastn in the shell and using the same syntax, only with
        # subprocess.check_output. As a workaround, we write the output to a temporary file.
        path = os.path.join(tmpdir, "igblast.txt")
        command = ["igblastn"] + variable_arguments + arguments + ["-out", path]
        logger.debug("Running %s", escape_shell_command(command))
        output = subprocess.check_output(
            command, input=fasta_str, universal_newlines=True, env=env
        )
        assert output == ""
        with open(path) as f:
            return f.read()


def copy_internal_data(tmp_path) -> bool:
    """
    Guess where igblastn’s internal_data/ directory is located and copy it
    recursively to the target directory.

    igblastn on Windows does not recognize internal_data/ directories whose
    absolute path contains spaces (even though the standard installation
    location is a subdirectory of C:\\Program Files\\).

    This works only if the target (tmp_path) does not contain spaces.

    Return True iff the internal_data/ directory was found and copied.
    """
    igblastn_path = Path(shutil.which("igblastn"))
    installation_path = igblastn_path.parent.parent
    internal_data_path = installation_path / "internal_data"
    if not internal_data_path.exists():
        return False

    shutil.copytree(internal_data_path, tmp_path / "internal_data")
    logger.debug("Copied %s to %s", internal_data_path, tmp_path / "internal_data")
    return True


def chunked(iterable, chunksize: int):
    """
    Group the iterable into lists of length chunksize
    >>> list(chunked('ABCDEFG', 3))
    [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
    """
    chunk = []
    for it in iterable:
        if len(chunk) == chunksize:
            yield chunk
            chunk = []
        chunk.append(it)
    if chunk:
        yield chunk


class RawRunner:
    """
    This is the target of a multiprocessing pool. The target needs to
    be pickleable, and because nested functions cannot be pickled,
    we need this separate class.

    It runs IgBLAST and returns raw AIRR-formatted output
    """

    def __init__(self, blastdb_dir, species, sequence_type, penalty, database, cache):
        self.blastdb_dir = blastdb_dir
        self.species = species
        self.sequence_type = sequence_type
        self.penalty = penalty
        self.database = database
        self.cache = cache

    def __call__(self, sequences):
        """
        Return raw IgBLAST output
        """
        return run_igblast(
            sequences, self.blastdb_dir, self.species, self.sequence_type, self.penalty
        )


class MakeBlastDbError(subprocess.CalledProcessError):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self):
        return (
            f"Running '{escape_shell_command(self.cmd)}' failed with "
            f"exit code {self.returncode}. "
            f"Standard output:\n{self.output.decode()}\n"
            f"Standard error:\n{self.stderr.decode()}"
        )


def makeblastdb(fasta, database_name, prefix=""):
    """
    prefix -- prefix to add to sequence ids
    """
    n = 0
    with dnaio.open(fasta) as fr, open(database_name + ".fasta", "w") as db:
        for record in fr:
            name = prefix + record.name.split(maxsplit=1)[0]
            db.write(">{}\n{}\n".format(name, record.sequence))
            n += 1
    if n == 0:
        raise ValueError("FASTA file {} is empty".format(fasta))

    command = [
        "makeblastdb",
        "-parse_seqids",
        "-dbtype",
        "nucl",
        "-in",
        database_name + ".fasta",
        "-out",
        database_name,
    ]
    logger.debug("Running %s", escape_shell_command(command))
    try:
        process_output = subprocess.check_output(command, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise MakeBlastDbError(
            e.returncode, command, output=e.output, stderr=e.stderr
        ) from None

    if b"Error: " in process_output:
        raise MakeBlastDbError(0, command, stderr=process_output) from None


class Database:
    def __init__(self, path, sequence_type):
        """path -- path to database directory with V.fasta, D.fasta, J.fasta"""
        self.path = path
        self.sequence_type = sequence_type
        self._v_records = self._read_fasta(os.path.join(path, "V.fasta"))
        self.v = self._records_to_dict(self._v_records)
        self._d_records = self._read_fasta(os.path.join(path, "D.fasta"))
        self.d = self._records_to_dict(self._d_records)
        self._j_records = self._read_fasta(os.path.join(path, "J.fasta"))
        self.j = self._records_to_dict(self._j_records)
        self._cdr3_starts = dict()
        self._cdr3_ends = dict()
        for locus in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"):
            self._cdr3_starts[locus] = {
                name: cdr3_start(s, locus) for name, s in self.v.items()
            }
            self._cdr3_ends[locus] = {
                name: cdr3_end(s, locus) for name, s in self.j.items()
            }
        self.v_regions_nt, self.v_regions_aa = self._find_v_regions()

    @staticmethod
    def _read_fasta(path):
        records = []
        with dnaio.open(path) as sr:
            for record in sr:
                record.name = record.name.split(maxsplit=1)[0]
                records.append(record)
        return records

    @staticmethod
    def _records_to_dict(records):
        return {record.name: record.sequence.upper() for record in records}

    def v_cdr3_start(self, gene, locus):
        return self._cdr3_starts[locus][gene]

    def j_cdr3_end(self, gene, locus):
        return self._cdr3_ends[locus][gene]

    def _find_v_regions(self):
        """
        Run IgBLAST on the V sequences to determine the nucleotide and amino-acid sequences of the
        FR1, CDR1, FR2, CDR2 and FR3 regions
        """
        v_regions_nt = dict()
        v_regions_aa = dict()
        for record in igblast_records(self.path, self._v_records, self.sequence_type):
            nt_regions = dict()
            aa_regions = dict()
            for region in ("FR1", "CDR1", "FR2", "CDR2", "FR3"):
                nt_seq = record.region_sequence(region)
                if nt_seq is None:
                    break
                if len(nt_seq) % 3 != 0:
                    logger.warning(
                        "Length %s of %s region in %r is not divisible by three; region "
                        "info for this gene will not be available",
                        len(nt_seq),
                        region,
                        record.query_name,
                    )
                    # not codon-aligned, skip entire record
                    break
                nt_regions[region] = nt_seq
                try:
                    aa_seq = nt_to_aa(nt_seq)
                except ValueError as e:
                    logger.warning(
                        "The %s region could not be converted to amino acids: %s",
                        region,
                        str(e),
                    )
                    break
                if "*" in aa_seq:
                    logger.warning(
                        "The %s region in %r contains a stop codon (%r); region info "
                        "for this gene will not be available",
                        region,
                        record.query_name,
                        aa_seq,
                    )
                    break
                aa_regions[region] = aa_seq
            else:
                v_regions_nt[record.query_name] = nt_regions
                v_regions_aa[record.query_name] = aa_regions

        return v_regions_nt, v_regions_aa


def igblast_records(
    database, sequences, sequence_type, species=None, penalty=None, use_cache=False
):
    """
    Run IgBLAST, parse results and yield RegionRecords objects.

    database -- Path to database directory with V./D./J.fasta files
    sequences -- an iterable of Sequence objects
    sequence_type -- 'Ig' or 'TCR'
    """
    # Create the BLAST databases in a temporary directory
    with tempfile.TemporaryDirectory() as blastdb_dir:
        make_vdj_blastdb(blastdb_dir, database)
        igblast_result = run_igblast(
            sequences, blastdb_dir, species, sequence_type, penalty
        )
        sio = StringIO(igblast_result)
        for record in parse_region_records(sio):
            yield record


def parse_region_records(file):
    csv.register_dialect(
        "airr",
        delimiter="\t",
        lineterminator="\n",
        strict=True,
    )
    reader = csv.DictReader(file, dialect="airr")
    for record in reader:
        yield RegionsRecord(query_name=record["sequence_id"], fields=record)


@dataclass
class RegionsRecord:
    query_name: str
    fields: Dict[str, str]

    COLUMNS_MAP = {
        "CDR1": "cdr1",
        "CDR2": "cdr2",
        "FR1": "fwr1",
        "FR2": "fwr2",
        "FR3": "fwr3",
    }

    def region_sequence(self, region: str) -> str:
        """
        Return the nucleotide sequence of a named region. Allowed names are:
        CDR1, CDR2, CDR3, FR1, FR2, FR3. Sequences are extracted from the full read
        using begin and end coordinates from IgBLAST’s "alignment summary" table.
        """
        if region not in self.COLUMNS_MAP:
            raise KeyError(f"Region '{region}' not allowed")
        return self.fields[self.COLUMNS_MAP[region]]


def igblast_parallel_chunked(
    database,
    sequences,
    sequence_type,
    species=None,
    threads=1,
    penalty=None,
    cache=None,
):
    """
    Run IgBLAST on the input sequences and yield AIRR-formatted results.

    The input is split up into chunks of 1000 sequences and distributed to
    *threads* number of IgBLAST instances that run in parallel.

    database -- Path to database directory with V./D./J.fasta files
    sequences -- an iterable of Sequence objects
    sequence_type -- 'Ig' or 'TCR'
    threads -- number of threads.
    """
    with ExitStack() as stack:
        # Create the three BLAST databases in a temporary directory
        blastdb_dir = stack.enter_context(tempfile.TemporaryDirectory())
        make_vdj_blastdb(blastdb_dir, database)

        chunks = chunked(sequences, chunksize=1000)
        runner = RawRunner(
            blastdb_dir,
            species=species,
            sequence_type=sequence_type,
            penalty=penalty,
            database=database,
            cache=cache,
        )
        pool = stack.enter_context(
            multiprocessing.Pool(threads) if threads > 1 else SerialPool()
        )
        for igblast_output in pool.imap(runner, chunks, chunksize=1):
            yield igblast_output


def make_vdj_blastdb(blastdb_dir, database_dir):
    """Run makeblastdb for all {V,D,J}.fasta in the database_dir"""

    for gene in ["V", "D", "J"]:
        # Without adding the "%" prefix, IgBLAST reports record names that look like GenBank
        # ids as "gb|original_name|".
        makeblastdb(
            os.path.join(database_dir, gene + ".fasta"),
            os.path.join(blastdb_dir, gene),
            prefix="%",
        )


# Regular expressions for CDR3 detection
#
# The idea comes from D’Angelo et al.: The antibody mining toolbox.
# http://dx.doi.org/10.4161/mabs.27105
# The heavy-chain regex was taken directly from there, but the difference
# is that we express everything in terms of amino acids, not nucleotides.
# This simplifies the expressions and makes them more readable.
#
_CDR3_REGEX = {
    # Heavy chain
    "IGH": re.compile(
        """
        [FY] [FHVWY] C
        (?P<cdr3>
            [ADEGIKMNRSTV] .{3,31}
        )
        W[GAV]
        """,
        re.VERBOSE,
    ),
    # Light chain, kappa
    "IGK": re.compile(
        """
        [FSVY] [CFHNVY] [CDFGLSW]
        (?P<cdr3>
            .{4,15}
        )
        [FLV][GRV]
        """,
        re.VERBOSE,
    ),
    # Light chain, lambda
    "IGL": re.compile(
        """
        # the negative lookahead assertion ensures that the rightmost start is found
        [CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]
        (?P<cdr3>
            .{4,15}
        )
        [FS]G
        """,
        re.VERBOSE,
    ),
}

_CDR3_VH_ALTERNATIVE_REGEX = re.compile(
    """
        C
        (?P<cdr3> . [RK] .{3,30})
        [WF]G.G
""",
    re.VERBOSE,
)


def find_cdr3(sequence, locus):
    """
    Find the CDR3 in the given sequence, assuming it comes from the given locus
    (chain type). If the locus is not one of 'IGH', 'IGK', 'IGL', return None.

    Return a tuple (start, stop) if found, None otherwise.
    """
    try:
        regex = _CDR3_REGEX[locus]
    except KeyError:
        return None
    matches = []
    for offset in 0, 1, 2:
        aa = nt_to_aa(sequence[offset:])
        match = regex.search(aa)
        if not match and locus == "IGH":
            match = _CDR3_VH_ALTERNATIVE_REGEX.search(aa)
        if match:
            start, stop = match.span("cdr3")
            matches.append((start * 3 + offset, stop * 3 + offset))
    return min(matches, default=None)


# The following code is used for detecting CDR3 start sites within V
# reference sequences and CDR3 end sites within J reference sequences.


# Matches the start of the CDR3 within the end of a VH sequence
_CDR3_START_VH_REGEX = re.compile(
    """
    [FY] [FHVWY] C
    (?P<cdr3_start>
        [ADEGIKMNRSTV*] | $
    )
    """,
    re.VERBOSE,
)


_CDR3_START_VH_ALTERNATIVE_REGEX = re.compile(
    """
    C
    (?P<cdr3_start> . [RK])
    """,
    re.VERBOSE,
)


_CDR3_START_REGEXES = {
    "IGK": re.compile("[FSVY][CFHNVY][CDFGLSW]"),
    "IGL": re.compile("[CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]"),
    "TRG": re.compile("[YFH]C"),  # TODO test whether this also works for alpha and beta
    "TRD": re.compile("[YFH]C"),
}


def _cdr3_start_heavy(aa):
    head, tail = aa[:-15], aa[-15:]
    match = _CDR3_START_VH_REGEX.search(tail)
    if not match:
        match = _CDR3_START_VH_ALTERNATIVE_REGEX.search(tail)
    if not match:
        return None
    return len(head) + match.start("cdr3_start")


def cdr3_start(nt, locus):
    """
    Find CDR3 start location within a V gene (Ig or TCR)

    nt -- nucleotide sequence of the gene
    locus -- one of "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"
    """
    aa = nt_to_aa(nt)
    if locus == "IGH":
        start = _cdr3_start_heavy(aa)
        if start is None:
            return None
        return 3 * start
    if locus in ("IGK", "IGL", "TRG", "TRD"):
        head, tail = aa[:-15], aa[-15:]
        match = _CDR3_START_REGEXES[locus].search(tail)
        if match:
            return 3 * (len(head) + match.end())
        else:
            return None
    elif locus in ("TRA", "TRB"):
        head, tail = aa[:-8], aa[-8:]
        pos = tail.find("C")
        if pos == -1:
            return None
        else:
            return 3 * (len(head) + pos + 1)


# Matches after the end of the CDR3 within a J sequence
_CDR3_END_REGEXES = {
    "IGH": re.compile("W[GAV]"),
    "IGK": re.compile("FG"),
    "IGL": re.compile("FG"),
    "TRA": re.compile("FG"),
    "TRB": re.compile("FG"),
    "TRG": re.compile("FG"),
    "TRD": re.compile("FG"),
}


def cdr3_end(nt, locus):
    """
    Find the position of the CDR3 end within a J sequence

    nt -- nucleotide sequence of the J gene
    locus -- one of "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"
    """
    regex = _CDR3_END_REGEXES[locus]
    for frame in 0, 1, 2:
        aa = nt_to_aa(nt[frame:])
        match = regex.search(aa)
        if match:
            return match.start() * 3 + frame
    return None


# When searching for the CDR3, start this many bases to the left of the end of
# the V match.
CDR3_SEARCH_START = 30


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = ArgumentParser(description=__doc__)
    arg = parser.add_argument
    arg(
        "--threads",
        "-t",
        "-j",
        type=int,
        default=1,
        help="Number of threads. Default: 1. Use 0 for no. of available CPUs.",
    )
    arg(
        "--penalty",
        type=int,
        choices=(-1, -2, -3, -4),
        default=None,
        help="BLAST mismatch penalty (default: -1)",
    )
    arg(
        "--species",
        default=None,
        help="Tell IgBLAST which species to use. Note that this setting does "
        "not seem to have any effect since we provide our own database to "
        "IgBLAST. Default: Use IgBLAST’s default",
    )
    arg(
        "--sequence-type",
        default="Ig",
        choices=("Ig", "TCR"),
        help="Sequence type. Default: %(default)s",
    )
    arg("--limit", type=int, metavar="N", help="Limit processing to first N records")

    arg("database", help="Database directory with V.fasta, D.fasta, J.fasta.")
    arg("fasta", help="File with original reads")

    args = parser.parse_args()
    run_igblastwrap(**vars(args))


def run_igblastwrap(
    threads: int,
    penalty: Optional[int],
    species: str,
    sequence_type: str,
    limit: Optional[int],
    database: str,
    fasta: str,
):
    if threads == 0:
        threads = available_cpu_count()
    start_time = time.time()
    last_status_update = 0
    with ExitStack() as stack:
        sequences = stack.enter_context(dnaio.open(fasta))
        sequences = islice(sequences, 0, limit)

        n = 0  # number of records processed so far
        for record in igblast_parallel_chunked(
            database,
            sequences,
            sequence_type=sequence_type,
            species=species,
            threads=threads,
            penalty=penalty,
        ):
            lines = record.splitlines()
            try:
                if n == 0:
                    print(*lines, sep="\n")
                else:
                    print(*lines[1:], sep="\n")
            except IOError as e:
                if e.errno == errno.EPIPE:
                    sys.exit(1)
                raise
            n += len(lines) - 1
            if n % 1000 == 0:
                elapsed = time.time() - start_time
                if elapsed >= last_status_update + 60:
                    logger.info(
                        "Processed {:10,d} sequences at {:.3f} ms/sequence".format(
                            n, elapsed / n * 1e3
                        )
                    )
                    last_status_update = elapsed
    elapsed = time.time() - start_time
    logger.info(
        "Processed {:10,d} sequences at {:.1f} ms/sequence".format(n, elapsed / n * 1e3)
    )


if __name__ == "__main__":
    main()