from . import pgdirect as pg
from .vcf import parse_chroms
import logging
import pandas as pd
from collections import defaultdict
import lzma

default_filter = {
    "deam_only": False,
    "pos_in_read_cutoff": 2,
    "min_length": 35,
    "max_length": 1000,
    "minq": 25,
}


class AdmixfrogInput(pg.ExtCoverage):
    
    class Obs:
        def __init__(self):
            # 
            self.n_ref = 0
            self.n_alt = 0
            self.n_deam = 0
            self.n_other = 0
            self.ref_err = []
            self.alt_err = []
            # self.base_pos = 0
            # self.ref_base = ''
            # self.alt_base = ''
    def __init__(
        self,
        outfile,
        deam_cutoff,
        length_bin_size=None,
        random_read_sample=False,
        report_alleles=False,
        max_reads=100,
        error_dict = None,
        flat_error = "0.002",
        **kwargs,
    ):
        self.outfile = outfile
        self.deam_cutoff = deam_cutoff
        self.length_bin_size = length_bin_size
        self.random_read_sample = random_read_sample
        self.report_alleles = report_alleles
        self.max_reads = max_reads
        self.error_dict = error_dict
        self.flat_error = flat_error
        if random_read_sample:
            raise NotImplementedError
        try:
            self.min_length = kwargs["min_length"]
        except KeyError:
            self.min_length = 35
        self.kwargs = kwargs

    def preprocess(self, sampleset):
        self.f = lzma.open(self.outfile, "wt")
        if self.error_dict is not None:    # reformat the admixfrog-input file header
            print(
                "chrom",
                "pos",
                "lib",
                "tref",
                "talt",
                "tdeam",
                "tother",
                "ref_err",    # to make error 
                "alt_err",    # to make error pattern
                sep=",",
                file=self.f,
            )
        else:
            print(
                "chrom",
                "pos",
                "lib",
                "tref",
                "talt",
                "tdeam",
                "tother",
                sep=",",
                file=self.f,
            )
    def process_snp(self, block, snp):    # error dic to add

        flat_error = self.flat_error
        flat_error = float(flat_error)
        reads = snp.reads(**self.kwargs)
        D = defaultdict(lambda: self.Obs())
        # n_ref, n_alt, n_deam, n_other = 0, 0, 0, 0
        i = 0
        for r in reads:
            i += 1
            if i > self.max_reads:
                print(
                    f"Warning at {snp.chrom}:{snp.pos+1} more than {i-1} reads found ({len(reads)}). Skipping the rest..."
                )
                break

            DEAM = (
                "deam"
                if (r.deam[0] < self.deam_cutoff or r.deam[1] < self.deam_cutoff)
                and r.deam[0] >= 0
                else "nodeam"
            )
            if self.length_bin_size is None:
                LEN = 0
            else:
                LEN = (r.len - self.min_length) // self.length_bin_size
            pos_terminal = r.pos_in_read
            if r.pos_in_read >= 15: 
                pos_terminal = 15 + (r.len - r.pos_in_read)
                if pos_terminal >= 30:
                    pos_terminal = 30
            if self.error_dict is not None:
                D[r.RG, DEAM, LEN].base_pos = pos_terminal
                D[r.RG, DEAM, LEN].ref_base = snp.ref
                D[r.RG, DEAM, LEN].alt_base = snp.alt
                if r.base == snp.ref:
                    D[r.RG, DEAM, LEN].n_ref += 1
                    if snp.ref == "C" and snp.alt == "T":  # r.base == C
                        D[r.RG, DEAM, LEN].ref_err.append(flat_error)
                    elif snp.ref == "T" and snp.alt == "C":  # r.base = T
                        if not r.is_reverse:
                            D[r.RG, DEAM, LEN].ref_err.append(self.error_dict[r.RG]['CT'][pos_terminal]+flat_error)
                        else:
                            D[r.RG, DEAM, LEN].ref_err.append(flat_error)
                    elif snp.ref == "G" and snp.alt == "A":  # r.base == G
                        D[r.RG, DEAM, LEN].ref_err.append(flat_error)
                    elif snp.ref == "A" and snp.alt == "G":  # r.base == A
                        if r.is_reverse:
                            D[r.RG, DEAM, LEN].ref_err.append(self.error_dict[r.RG]['GA'][pos_terminal]+flat_error)
                        else:
                            D[r.RG, DEAM, LEN].ref_err.append(flat_error)
                    else:
                        D[r.RG, DEAM, LEN].ref_err.append(flat_error)
                elif r.base == snp.alt:
                    D[r.RG, DEAM, LEN].n_alt += 1
                    if snp.ref == "C" and snp.alt == "T":    # r.base == T
                        if not r.is_reverse:
                            D[r.RG, DEAM, LEN].alt_err.append(self.error_dict[r.RG]['CT'][pos_terminal]+flat_error)
                        else:
                            D[r.RG, DEAM, LEN].alt_err.append(flat_error)
                    elif snp.ref == "T" and snp.alt == "C":  # r.base = C
                        D[r.RG, DEAM, LEN].alt_err.append(flat_error)
                    elif snp.ref == "G" and snp.alt == "A":  # r.base == A
                        if r.is_reverse:
                            D[r.RG, DEAM, LEN].alt_err.append(self.error_dict[r.RG]['GA'][pos_terminal]+flat_error)
                        else:
                            D[r.RG, DEAM, LEN].alt_err.append(flat_error)
                    else:
                        D[r.RG, DEAM, LEN].alt_err.append(flat_error)
                elif r.base == "T" and not r.is_reverse and "C" in (snp.ref, snp.alt):
                    D[r.RG, DEAM, LEN].n_deam += 1
                elif r.base == "A" and r.is_reverse and "G" in (snp.ref, snp.alt):
                    D[r.RG, DEAM, LEN].n_deam += 1
                elif r.base != "N":
                    D[r.RG, DEAM, LEN].n_other += 1
            else:
                D[r.RG, DEAM, LEN].base_pos = pos_terminal
                D[r.RG, DEAM, LEN].ref_base = snp.ref
                D[r.RG, DEAM, LEN].alt_base = snp.alt
                if r.base == snp.ref:
                    D[r.RG, DEAM, LEN].n_ref += 1
                elif r.base == snp.alt:
                    D[r.RG, DEAM, LEN].n_alt += 1
                elif r.base == "T" and not r.is_reverse and "C" in (snp.ref, snp.alt):
                    D[r.RG, DEAM, LEN].n_deam += 1
                elif r.base == "A" and r.is_reverse and "G" in (snp.ref, snp.alt):
                    D[r.RG, DEAM, LEN].n_deam += 1
                elif r.base != "N":
                    D[r.RG, DEAM, LEN].n_other += 1
        if self.error_dict is not None:               
            for (rg, deam, len_), r in D.items():
                if self.length_bin_size is None:
                    lib = f"{rg}_{deam}"
                else:
                    lib = f"{rg}_{len_}_{deam}"
                if self.report_alleles:
                    alleles = "".join(sorted(snp.ref + snp.alt))
                    lib = f"{lib}_{alleles}"
                print(
                    snp.chrom,
                    snp.pos + 1,
                    lib,
                    r.n_ref,
                    r.n_alt,
                    r.n_deam,
                    r.n_other,
                    " ".join([str(x) for x in r.ref_err]),
                    " ".join([str(x) for x in r.alt_err]),
                    file=self.f,
                    sep=",",
                )
        else:
            for (rg, deam, len_), r in D.items():
                if self.length_bin_size is None:
                    lib = f"{rg}_{deam}"
                else:
                    lib = f"{rg}_{len_}_{deam}"
                if self.report_alleles:
                    alleles = "".join(sorted(snp.ref + snp.alt))
                    lib = f"{lib}_{alleles}"
                print(
                    snp.chrom,
                    snp.pos + 1,
                    lib,
                    r.n_ref,
                    r.n_alt,
                    r.n_deam,
                    r.n_other,
                    file=self.f,
                    sep=",",
                )


class AdmixfrogInput2(pg.ExtCoverage):
    def __init__(
        self,
        outfile,
        deam_cutoff,
        length_bin_size=None,
        random_read_sample=False,
        report_alleles=False,
        **kwargs,
    ):
        self.outfile = outfile
        self.deam_cutoff = deam_cutoff
        self.length_bin_size = 1
        self.random_read_sample = random_read_sample
        self.report_alleles = report_alleles
        if random_read_sample:
            raise NotImplementedError
        try:
            self.min_length = kwargs["min_length"]
        except KeyError:
            self.min_length = 35
        self.kwargs = kwargs

    def preprocess(self, sampleset):
        self.f = lzma.open(self.outfile, "wt")
        print(
            "chrom",
            "pos",
            "tref",
            "talt",
            "lib",
            "len",
            "deam",
            "dmgsite",
            sep=",",
            file=self.f,
        )

    def process_snp(self, block, snp):
        reads = snp.reads(**self.kwargs)
        D = defaultdict(lambda: self.Obs())
        # n_ref, n_alt, n_deam, n_other = 0, 0, 0, 0
        for r in reads:
            DEAM = (
                "deam"
                if (r.deam[0] < self.deam_cutoff or r.deam[1] < self.deam_cutoff)
                and r.deam[0] >= 0
                else "nodeam"
            )
            if self.length_bin_size is None:
                LEN = 0
            else:
                LEN = r.len // self.length_bin_size

            DEAM = min((r.deam[0], r.deam[1]))

            if r.base == snp.ref:
                if r.base == "T" and snp.alt == "C" and not r.is_reverse:
                    dmgsite = 1
                elif r.base == "A" and snp.alt == "G" and r.is_reverse:
                    dmgsite = 1
                else:
                    dmgsite = 0
                D[r.RG, DEAM, LEN, dmgsite].n_ref += 1
            elif r.base == snp.alt:
                if r.base == "T" and snp.ref == "C" and not r.is_reverse:
                    dmgsite = 1
                elif r.base == "A" and snp.ref == "G" and r.is_reverse:
                    dmgsite = 1
                else:
                    dmgsite = 0
                D[r.RG, DEAM, LEN, dmgsite].n_alt += 1

        for (rg, deam, len_, dmgsite_), r in D.items():
            print(
                snp.chrom,
                snp.pos + 1,
                r.n_ref,
                r.n_alt,
                rg,
                len_,
                deam,
                dmgsite_,
                file=self.f,
                sep=",",
            )


class RefIter:
    def __init__(self, ref):
        self.ref = pd.read_csv(ref, dtype={"chrom": "str"})
        self.bed = self.ref

    def __iter__(self):
        for ix, row in self.ref.iterrows():
            yield 0, (row.chrom, row.pos - 1, row.ref, row.alt)


def process_bam(
    outfile,
    bamfile,
    ref,
    deam_cutoff,
    length_bin_size,
    error_file,
    random_read_sample=False,
    max_reads=100,
    chroms=None,
    **kwargs,
):
    """generate input file from bam-file"""
    blocks = RefIter(ref)
    chroms = parse_chroms(chroms)
    sampleset = pg.CallBackSampleSet.from_file_names(
        [bamfile], blocks=blocks, chroms=chroms
    )
    
    """if error_file is not None: read error rates into a dict
    """

    default_filter.update(kwargs)
    logging.info("Filter is %s", default_filter)
    if error_file is not None:
        error_dict = defaultdict(lambda: defaultdict(dict))
        for error_file_ in error_file:
            with open(error_file_, "r") as ef:
                for i_, line in enumerate(ef):
                    if i_ == 0:
                        continue
                    lib_name, pos, CT_error, GA_error = line.strip().split()
                    error_dict[lib_name]['CT'][int(pos)] = float(CT_error)
                    error_dict[lib_name]['GA'][int(pos)] = float(GA_error)
        for i in error_dict:
            print(f"Lib {i} has position-based error rates loaded.")
            print("pos\tCT_error\tGA_error")
            for pos in range(0,31):
                print(f"{pos}\t{error_dict[i]['CT'].get(pos,'NA')}\t{error_dict[i]['GA'].get(pos,'NA')}")
        cov = AdmixfrogInput(
            **default_filter,
            random_read_sample=random_read_sample,
            length_bin_size=length_bin_size,
            deam_cutoff=deam_cutoff,
            outfile=outfile,
            max_reads=max_reads,
            error_dict=error_dict,
        )
    else:
        cov = AdmixfrogInput(
            **default_filter,
            random_read_sample=random_read_sample,
            length_bin_size=length_bin_size,
            deam_cutoff=deam_cutoff,
            outfile=outfile,
            max_reads=max_reads,
        )
    sampleset.add_callback(cov)
    sampleset.run_callbacks()


def process_bam2(
    outfile,
    bamfile,
    ref,
    deam_cutoff,
    length_bin_size,
    random_read_sample=False,
    chroms=None,
    **kwargs,
):
    """generate 2nd generation input file from bam-file
    file will have output cols
    chrom, pos, tref, talt, lib, len, deam, score
    """
    blocks = RefIter(ref)
    chroms = parse_chroms(chroms)
    sampleset = pg.CallBackSampleSet.from_file_names(
        [bamfile], blocks=blocks, chroms=chroms
    )

    default_filter.update(kwargs)
    logging.info("Filter is %s", default_filter)
    cov = AdmixfrogInput2(
        **default_filter,
        random_read_sample=random_read_sample,
        length_bin_size=length_bin_size,
        deam_cutoff=deam_cutoff,
        outfile=outfile,
    )
    sampleset.add_callback(cov)
    sampleset.run_callbacks()
