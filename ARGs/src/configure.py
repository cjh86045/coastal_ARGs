import pathlib
import numpy as np
from Bio import SeqIO
import pandas as pd
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed  # 用于多线程处理
from tqdm import tqdm
from pathlib import Path

# 软件线程配置
THREAD_CONFIG = {
    #preprocessing 流程线程配置
    "fastp": 16, # fastp软件固定最多只能使用16个线程
    "fastp_pool": 10,

    "megahit": 25,   # megahit 线程数
    "megahit_pool": 5 ,

    #"prodigal": 1,   #
    "prodigal_pool": 30,# prodigal 线程数（每个任务）

    "diamond": 9,
    "diamond_pool": 19,

    "bowtie2": 10,
    "samtools": 10,
    "get_depth_pool":15,

    "maxbin":2,
    "maxbin_pool":40,

    "metabat":2,
    "metabat_pool":40,

    "concoct_pool": 40,
    "concoct": 3,

    "plasme":30,
    "plasme_pool":1,

    "dastool": 30,
    "dastool_pool": 1,

    "checkm":6,
    "checkm_pool":5,

    "gtdb":20,
    "gtdb_pool":3,

    "sylph": 10,
    "sylph_pool": 18,

    "mafft": 15,
    "mafft_pool": 10,

    "fasttree_pool": 30,



}

SARG_DB = pathlib.Path("")
SARG_STRUCTURE = pathlib.Path("")

KO30_DB = pathlib.Path("")
KO30_STRUCTURE = pathlib.Path("")

CARD_DB = pathlib.Path("")
CARD_STRUCTURE = pathlib.Path("")


NCBI_NR_DB = pathlib.Path("")
NCBI_NR_DB_fa = pathlib.Path("")
NCBI_NR_DB_index = pathlib.Path("")

VFDB_PRO = pathlib.Path("")


SYLPH_DB = ""

valid_suffixes = [".fa", ".fastq", ".fa.gz", ".fasta.gz", ".fq", ".fastq", ".fq.gz", ".fastq.gz"]

MGEs_keywords = ["mobilization", "transposon", "recombination", "recombinase",
                 "integrase", "conjugative", "transposase", "integron", "plasmid",
                 "conjugal", "invertase", "relaxase", "resolvase"]

diamond_outfmt = []