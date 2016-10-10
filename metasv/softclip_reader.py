import logging
import sys

import vcf

from sv_interval import SVInterval

logger = logging.getLogger(__name__)

tool_name = "SoftClip"
source = set([tool_name])


class SoftClipRecord:
    def __init__(self, vcf_record):
        self.name = tool_name
        info_fields = vcf_record.INFO
        self.chromosome = vcf_record.CHROM
        self.start = vcf_record.POS
        self.end = info_fields["END"]
        self.score = info_fields["FRAC_DISCORDANT_SUPPORT"]
        self.filter = "LowQual" if "LowQual" in vcf_record.FILTER else "PASS"
        self.sc_support = info_fields["NUM_SC_SUPPORT"]
        self.pe_support = info_fields["NUM_PE_SUPPORT"]
        self.read_support = info_fields["NUM_READ_SUPPORT"]
        self.left_support = info_fields["NUM_LEFT_SUPPORT"]
        self.right_support = info_fields["NUM_RIGHT_SUPPORT"]
        self.discordant_support = info_fields["NUM_DISCORDANT_SUPPORT"]
        self.coverage = info_fields["COVERAGE"]
        self.frac_discordant_support= info_fields["FRAC_DISCORDANT_SUPPORT"]
        self.imprecise= True if "IMPRECISE" in info_fields else False
        self.ciend= info_fields["CIEND"]
        self.cipos= info_fields["CIPOS"]
        self.sv_type= info_fields["SVTYPE"]
        self.sv_len = abs(info_fields["SVLEN"][0])
        self.info = {
            "SC_START": self.start,
            "SC_END": self.end,
            "SC_SCORE": self.score,
            "SC_FILTER": self.filter,
            "SC_NUM_SC_SUPPORT": self.sc_support,
            "SC_NUM_PE_SUPPORT": self.pe_support,
            "SC_NUM_READ_SUPPORT": self.read_support,
            "SC_NUM_LEFT_SUPPORT": self.left_support,
            "SC_NUM_RIGHT_SUPPORT": self.right_support,
            "SC_NUM_DISCORDANT_SUPPORT": self.discordant_support,
            "SC_COVERAGE": self.coverage,
            "SC_FRAC_DISCORDANT_SUPPORT": self.frac_discordant_support,
        }
        if "IMPRECISE" in info_fields: 
            self.info["SC_IMPRECISE"]="."   

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return "<" + self.__class__.__name__ + " " + str(self.__dict__) + ">"

    def to_sv_interval(self):
        if self.sv_type not in SoftClipReader.svs_supported:
            return None

        return SVInterval(self.chromosome,
                          self.start,
                          self.end,
                          name=self.name,
                          sv_type=self.sv_type,
                          length=self.sv_len,
                          sources=source,
                          cipos=self.cipos,
                          ciend=self.ciend,
                          info=self.info,
                          native_sv=self)

    def to_vcf_record(self, sample):
        alt = ["<%s>" % self.sv_type]
        sv_len = -self.sv_len if self.sv_type == "DEL" else self.sv_len
        info = {"SVLEN": sv_len, "SVTYPE": self.sv_type, "END": self.end}

        info.update(self.info)

        return vcf.model._Record(self.chromosome,
                                 self.start,
                                 ".",
                                 "N",
                                 alt,
                                 ".",
                                 ".",
                                 info,
                                 "GT",
                                 [0],
                                 [vcf.model._Call(None, sample, vcf.model.make_calldata_tuple("GT")(GT="1/1"))])


class SoftClipReader:
    svs_supported = set(["DEL", "INS", "DUP", "INV"])

    def __init__(self, file_name, reference_handle=None, svs_to_report=None):
        logger.info("File is " + str(file_name))
        self.file_fd = open(file_name) if file_name is not None else sys.stdin
        self.reference_handle = reference_handle
        self.svs_supported = SoftClipReader.svs_supported
        if svs_to_report is not None:
            self.svs_supported &= set(svs_to_report)

        self.vcf_reader = vcf.Reader(self.file_fd)


    def __iter__(self):
        return self

    def next(self):
        while True:
            vcf_record = self.vcf_reader.next()
            if vcf_record:
                record = SoftClipRecord(vcf_record)
                if record.filter == "LowQual":
                    continue
                if record.sv_type in self.svs_supported:
                    return record
            else:
                continue
