#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/3/31 6:17 PM
# @File    : aa_index.py
# @annotation    : define properties of amino acids

ALPHABET = "QWERTYIPASDFGHKLCVNM"

hydrophobic_index = {
    'R': -0.9,
    'K': -0.889,
    'D': -0.767,
    'E': -0.696,
    'N': -0.674,
    'Q': -0.464,
    'S': -0.364,
    'G': -0.342,
    'H': -0.271,
    'T': -0.199,
    'A': -0.171,
    'P': 0.055,
    'Y': 0.188,
    'V': 0.331,
    'M': 0.337,
    'C': 0.508,
    'L': 0.596,
    'F': 0.646,
    'I': 0.652,
    'W': 0.9
}

volume_index = {
    'G': -0.9,
    'A': -0.677,
    'S': -0.544,
    'C': -0.359,
    'T': -0.321,
    'P': -0.294,
    'D': -0.281,
    'N': -0.243,
    'V': -0.232,
    'E': -0.058,
    'Q': -0.02,
    'L': -0.009,
    'I': -0.009,
    'M': 0.087,
    'H': 0.138,
    'K': 0.163,
    'F': 0.412,
    'R': 0.466,
    'Y': 0.541,
    'W': 0.9
}

helix_tendency = {
    'P': -0.9,
    'G': -0.9,
    'C': -0.652,
    'S': -0.466,
    'T': -0.403,
    'N': -0.403,
    'Y': -0.155,
    'D': -0.155,
    'V': -0.031,
    'H': -0.031,
    'I': 0.155,
    'F': 0.155,
    'W': 0.279,
    'K': 0.279,
    'Q': 0.528,
    'R': 0.528,
    'M': 0.652,
    'L': 0.714,
    'E': 0.9,
    'A': 0.9
}

sheet_tendency = {
    'G': -0.9,
    'D': -0.635,
    'E': -0.582,
    'N': -0.529,
    'A': -0.476,
    'R': -0.371,
    'Q': -0.371,
    'K': -0.265,
    'S': -0.212,
    'H': -0.106,
    'L': -0.053,
    'M': -0.001,
    'P': 0.106,
    'T': 0.212,
    'F': 0.318,
    'C': 0.476,
    'Y': 0.476,
    'W': 0.529,
    'I': 0.688,
    'V': 0.9
}

class_type_dict = {
    '_small': 'GAVSTC',
    '_large': 'FYWKRHQE',
    '_neg': 'DE',
    '_pos': 'RK',
    '_polar': 'YTSHKREDQN',
    '_non_charged_polar': 'YTSNQH',
    '_hydrophobic': 'FILVAGMW',
    '_cys': "C",
    '_pro': 'P',
    '_scan': 'ARNDCQEGHILKMFPSTWYV'
}
