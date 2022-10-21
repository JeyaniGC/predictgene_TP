   
"""Tests for gene prediction"""
import pytest
import os
import re
from .context import gpred
from gpred import *


def test_read_fasta():
    sequence = read_fasta(os.path.abspath(os.path.join(os.path.dirname(__file__), "genome.fasta")))
    assert(sequence == "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGAGGAGGTAACTCAAACCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA")


def test_find_start():
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    seq_with_start = "AACGGCGTGAAACC"
    seq_without_start = "AAAAAAAAAAACCCCCCCCCCC"
    res_start = find_start(start_regex, seq_with_start, 0, len(seq_with_start))
    res_no_start = find_start(start_regex, seq_without_start, 0, len(seq_without_start))
    assert(res_start == 6)
    assert(res_no_start == None)


def test_find_stop():
    stop_regex = re.compile('TA[GA]|TGA')
    seq_with_stop = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    seq_without_stop = seq_with_stop[:-3]
    res_stop = find_stop(stop_regex, seq_with_stop, 3)
    res_no_stop = find_stop(stop_regex, seq_without_stop, 3)
    assert(res_stop == 63)
    assert(res_no_stop == None)


def test_has_shine_dalgarno():
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    seq_without_sd = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    seq_with_sd = "AGGAGGTAACTCAAACC" + seq_without_sd
    seq_with_sd_too_close = "AGGAGGTAACTC" + seq_without_sd
    seq_with_sd_too_far = "AGGAGGTAACTCAAACCGG" + seq_without_sd
    res_sd = has_shine_dalgarno(shine_regex, seq_with_sd, 17, 16)
    res_no_sd = has_shine_dalgarno(shine_regex, seq_without_sd, 15, 16)
    res_sd_too_close = has_shine_dalgarno(shine_regex, seq_with_sd_too_close, 12, 16)
    res_sd_too_far = has_shine_dalgarno(shine_regex, seq_with_sd_too_far, 19, 16)
    assert(res_sd == True)
    assert(res_no_sd == False)
    assert(res_sd_too_close == False)
    assert(res_sd_too_far == False)


def test_predict_genes():
    sequence = read_fasta(os.path.abspath(os.path.join(os.path.dirname(__file__), "genome.fasta")))
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    res = predict_genes(sequence, start_regex, stop_regex, shine_regex, 50, 16, 40)
    assert(len(res) == 1)
    assert(res[0][0] == 337)
    assert(res[0][1] == 483)
    res_neg = predict_genes(sequence[0:400], start_regex, stop_regex, shine_regex, 50, 16, 40)
    assert(len(res_neg) == 0)
