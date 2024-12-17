import re
import pandas as pd
import numpy as np
from itertools import product
from collections import defaultdict, Counter
from sigmap.const import dict_pstnp, list_best_features
from sigmap.setup import make_dimers_dict, make_trimers_dict

class FeatureExtraction:

    def __init__(self, df_seq:pd.DataFrame):
        
        self.input = df_seq

        self.list_seq   = self.input['Sequence'].tolist()
        self.list_seqID = self.input['Sequence_ID'].tolist()
        
        self.cdk = self.motif_features(self.input)
        self.cdk = self.property_features(self.cdk)
        
    # def End
        

    def allp(self, k):
        nucleotides = ['A', 'T', 'C', 'G']
        return [''.join(p) for p in product(nucleotides, repeat=k)]

    def kmer(self, k, seq):
        return [seq[i:i+k] for i in range(len(seq) - k + 1)]

    def property_features(self, cdk:pd.DataFrame):
        
        # DNA_Di_Prop.txt file contains all the properties along with its value.
        dinuc_data = pd.DataFrame({
            'Physicochemical properties': [
                f'p{i}' for i in range(1, 39)
            ],
            'GA': [-0.654, -0.14400000000000002, 1.112, 0.0, 1.023, -0.27399999999999997, 0.27, -0.5, -0.21899999999999997, -0.11800000000000001, 0.8270000000000001, 0.15, -0.19, -0.5, -1.105, -1.251, 0.0, -0.413, -0.495, 1.023, -0.402, 0.0, 0.17, -0.026000000000000002, -0.036000000000000004, -0.434, 0.516, 0.6509999999999999, -0.23600000000000002, 0.08800000000000001, 0.191, 0.0, -0.081, 0.502, 0.266, 0.126, -0.39399999999999996, 0.711],
            'GC': [-2.455, -0.301, 0.7859999999999999, 1.369, 0.322, 0.47200000000000003, -1.232, 1.4040000000000001, 2.353, 0.6659999999999999, -0.22399999999999998, 0.588, 1.5190000000000001, 1.4040000000000001, 1.35, 1.162, 1.369, -1.7069999999999999, 1.4080000000000001, 0.322, -1.7069999999999999, 1.369, 1.9969999999999999, 2.354, 2.1, 0.076, 2.517, 2.45, 1.3530000000000002, 1.04, 0.8440000000000001, 1.369, -0.081, 0.215, 1.331, -0.348, 0.6459999999999999, 1.585],
            'GG': [-0.07, 0.355, -0.055999999999999994, 1.369, -1.1909999999999998, 1.3969999999999998, -1.232, 1.4040000000000001, 0.67, 2.076, -0.496, -0.579, -0.498, 1.4040000000000001, 1.306, 1.1440000000000001, 1.369, -0.485, 1.4080000000000001, -1.1909999999999998, -0.488, 1.369, 0.86, -0.726, -1.276, -0.789, 0.49700000000000005, 0.068, 0.956, 1.264, 1.2930000000000001, 1.369, 0.063, 1.077, 0.08900000000000001, 0.56, -0.8220000000000001, 0.242],
            'GT': [-0.9179999999999999, -0.831, -0.653, 0.0, -1.359, -0.156, 0.27, -0.8809999999999999, 1.107, -0.11800000000000001, -1.041, 1.025, 0.259, -0.8809999999999999, -0.703, -0.556, 0.0, -1.276, -0.887, -1.359, -1.278, 0.0, 0.10300000000000001, 0.604, 0.6759999999999999, 0.852, 0.971, 0.915, -0.23600000000000002, 0.424, 0.6409999999999999, 0.0, 1.5019999999999998, 0.502, 0.799, 0.126, 1.2890000000000001, 1.044],
            'AA': [1.0190000000000001, -0.644, -0.002, -1.369, 0.995, -1.8869999999999998, 0.833, -0.11900000000000001, -0.8420000000000001, -0.9009999999999999, 0.36, 0.515, -0.9329999999999999, -0.11900000000000001, 0.45799999999999996, 0.6679999999999999, -1.369, 0.593, -0.132, 0.995, 0.593, -1.369, -0.81, 0.46399999999999997, 0.8340000000000001, -0.7, -0.77, -1.02, -0.831, -0.361, -0.16899999999999998, -1.369, 0.063, 0.502, 0.266, 1.587, 0.111, -0.109],
            'AC': [-0.9179999999999999, -0.831, -0.653, 0.0, -1.359, -0.156, 0.27, -0.8809999999999999, 1.107, -0.11800000000000001, -1.041, 1.025, 0.259, -0.8809999999999999, -0.703, -0.556, 0.0, -1.276, -0.887, -1.359, -1.278, 0.0, 0.10300000000000001, 0.604, 0.6759999999999999, 0.852, 0.971, 0.915, -0.23600000000000002, 0.424, 0.6409999999999999, 0.0, 1.5019999999999998, 0.502, 0.799, 0.126, 1.2890000000000001, 1.044],
            'AG': [0.488, -0.894, -1.3319999999999999, 0.0, -0.799, -0.436, 0.27, -0.5, 0.016, -0.11800000000000001, -0.885, 0.15, -0.99, -0.5, -0.12300000000000001, 0.083, 0.0, 0.23399999999999999, -0.495, -0.799, 0.233, 0.0, -0.498, -1.147, -1.1440000000000001, -0.5670000000000001, -0.612, -0.489, -0.23600000000000002, -1.145, -1.406, 0.0, 0.7829999999999999, 0.359, 0.08900000000000001, 0.679, -0.24100000000000002, -0.623],
            'AT': [0.5670000000000001, -1.05, 2.089, -1.369, -0.098, -0.7509999999999999, 1.396, -1.3880000000000001, -0.5760000000000001, -1.371, -1.896, 1.973, 1.03, -0.627, 0.23399999999999999, 0.65, -1.369, -0.485, -0.615, -0.098, -0.488, -1.369, -1.456, -0.866, -0.43200000000000005, 3.159, -0.669, -0.568, -1.4269999999999998, -1.705, -1.676, -1.369, 1.071, 0.215, 0.621, -1.0190000000000001, 2.513, 1.171],
            'CA': [0.5670000000000001, 1.51, 0.596, 0.0, 1.1909999999999998, 0.98, -0.106, -0.11900000000000001, -0.915, -0.11800000000000001, 1.216, -1.38, 0.45399999999999996, -0.11900000000000001, -1.015, -1.361, 0.0, 1.0959999999999999, -0.132, 1.1909999999999998, 1.091, 0.0, -0.008, -0.23600000000000002, -0.3, 0.032, -0.043, -0.568, 0.161, -0.249, -0.371, 0.0, -1.376, -1.364, -0.266, -0.861, -0.623, -1.254],
            'CC': [-0.07, 0.355, -0.055999999999999994, 1.369, -1.1909999999999998, 1.3969999999999998, -1.232, 1.4040000000000001, 0.67, 2.076, -0.496, -0.579, -0.498, 1.4040000000000001, 1.306, 1.1440000000000001, 1.369, -0.485, 1.4080000000000001, -1.1909999999999998, -0.488, 1.369, 0.86, -0.726, -1.276, -0.789, -0.762, 0.068, 0.956, 1.264, 1.2930000000000001, 1.369, 0.063, 1.077, 0.08900000000000001, 0.56, -0.8220000000000001, 0.242],
            'CG': [-0.579, 2.229, -1.1420000000000001, 1.369, -0.266, 0.799, -2.17, 2.039, 0.187, 0.6659999999999999, 0.7490000000000001, -1.818, 2.36, 2.039, 1.7069999999999999, 1.3630000000000002, 1.369, 0.665, 2.042, -0.266, 0.662, 1.369, 1.5730000000000002, 1.6540000000000001, 1.335, -0.41200000000000003, -0.762, 0.606, 2.346, 1.768, 1.4280000000000002, 1.369, -1.6640000000000001, -1.22, -0.444, -0.8220000000000001, -0.287, -1.389],
            'CT': [0.488, -0.894, -1.3319999999999999, 0.0, -0.799, -0.436, 0.27, -0.5, 0.016, -0.11800000000000001, -0.885, 0.15, -0.99, -0.5, -0.12300000000000001, 0.083, 0.0, 0.23399999999999999, -0.495, -0.799, 0.233, 0.0, -0.498, -1.147, -1.1440000000000001, -0.5670000000000001, 0.49700000000000005, -0.489, -0.23600000000000002, -1.145, -1.406, 0.0, 0.7829999999999999, 0.359, 0.08900000000000001, 0.679, -0.24100000000000002, -0.623],
            'TA': [1.6030000000000002, 0.418, -1.061, -1.369, 0.322, 0.233, 1.396, -0.627, -1.598, -1.371, 1.41, -0.506, -1.114, -1.3880000000000001, -0.9259999999999999, -0.629, -1.369, 2.031, -1.37, 0.322, 2.036, -1.369, -1.746, -1.006, -0.511, 0.387, -1.486, -1.6030000000000002, -1.4269999999999998, -1.145, -0.956, -1.369, -1.2329999999999999, -2.3680000000000003, -0.444, -2.2430000000000003, -1.511, -1.389],
            'TC': [-0.654, -0.14400000000000002, 1.112, 0.0, 1.023, -0.27399999999999997, 0.27, -0.5, -0.21899999999999997, -0.11800000000000001, 0.8270000000000001, 0.15, -0.19, -0.5, -1.105, -1.251, 0.0, -0.413, -0.495, 1.023, -0.402, 0.0, 0.17, -0.026000000000000002, -0.036000000000000004, -0.434, 0.516, 0.6509999999999999, -0.23600000000000002, 0.08800000000000001, 0.191, 0.0, -0.081, 0.502, 0.266, 0.126, -0.39399999999999996, 0.711],
            'TG': [0.5670000000000001, 1.51, 0.596, 0.0, 1.1909999999999998, 0.98, -0.106, -0.11900000000000001, -0.915, -0.11800000000000001, 1.216, -1.38, 0.45399999999999996, -0.11900000000000001, -1.015, -1.361, 0.0, 1.0959999999999999, -0.132, 1.1909999999999998, 1.091, 0.0, -0.008, -0.23600000000000002, -0.3, 0.032, -0.612, -0.568, 0.161, -0.249, -0.371, 0.0, -1.376, -1.364, -0.266, -0.861, -0.623, -1.254],
            'TT': [1.0190000000000001, -0.644, -0.002, -1.369, 0.995, -1.8869999999999998, 0.833, -0.11900000000000001, -0.8420000000000001, -0.9009999999999999, 0.36, 0.515, -0.9329999999999999, -0.11900000000000001, 0.45799999999999996, 0.6679999999999999, -1.369, 0.593, -0.132, 0.995, 0.593, -1.369, -0.81, 0.46399999999999997, 0.8340000000000001, -0.7, -0.77, -1.02, -0.831, -0.361, -0.16899999999999998, -1.369, 0.063, 0.502, -3.284, 1.587, 0.111, -0.109],
        })

        # Not used feature, will be removed
        self.trinuc_data = pd.DataFrame({
            'Physicochemical properties': ['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12'],
            'AAA': [-2.0869999999999997, -2.745, -1.732, -2.349, -2.7439999999999998, -2.7439999999999998, 2.274, 2.1180000000000003, -1.0, -1.0, -2.342, 2.386],
            'AAC': [-1.5090000000000001, -1.354, -0.5770000000000001, -0.561, -1.3630000000000002, -1.3630000000000002, 1.105, 1.516, -1.0, -1.0, -0.555, 0.5479999999999999],
            'AAG': [-0.506, -0.257, -0.5770000000000001, 0.155, -0.26, -0.26, 0.193, 0.493, -1.0, -1.0, 0.16899999999999998, -0.179],
            'AAT': [-2.126, -2.585, -1.732, -1.9909999999999999, -2.591, -2.591, 2.141, 2.158, -1.0, -1.0, -2.004, 2.032],
            'ACA': [0.111, 0.171, -0.5770000000000001, 0.155, 0.16399999999999998, 0.16399999999999998, -0.153, -0.12300000000000001, 1.0, 1.0, 0.16899999999999998, -0.179],
            'ACC': [-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.071, 0.071, -0.078, 0.107, 1.0, 1.0, 0.266, -0.275],
            'ACG': [-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.065, 0.065, -0.07400000000000001, 0.107, 1.0, 1.0, 0.266, -0.275],
            'ACT': [-1.354, -0.685, -0.5770000000000001, 0.45299999999999996, -0.6759999999999999, -0.6759999999999999, 0.536, 1.357, 1.0, 1.0, 0.45899999999999996, -0.466],
            'AGA': [0.381, -0.15, -0.5770000000000001, -0.74, -0.158, -0.158, 0.109, -0.389, 1.0, 1.0, -0.748, 0.743],
            'AGC': [0.304, 0.92, 0.5770000000000001, 1.287, 0.9109999999999999, 0.9109999999999999, -0.753, -0.313, 1.0, 1.0, 1.28, -1.272],
            'AGG': [-0.313, -0.07, 0.5770000000000001, 0.27399999999999997, -0.07, -0.07, 0.039, 0.3, 1.0, 1.0, 0.266, -0.275],
            'AGT': [-1.354, -0.685, -0.5770000000000001, 0.45299999999999996, -0.6759999999999999, -0.6759999999999999, 0.536, 1.357, 1.0, 1.0, 0.45899999999999996, -0.466],
            'ATA': [1.615, 0.5720000000000001, -1.732, -0.978, 0.584, 0.584, -0.491, -1.585, -1.0, -1.0, -0.99, 0.988],
            'ATC': [-0.737, -0.391, -0.5770000000000001, 0.214, -0.397, -0.397, 0.307, 0.727, -1.0, -1.0, 0.217, -0.22699999999999998],
            'ATG': [1.229, 1.348, -0.5770000000000001, 0.87, 1.358, 1.358, -1.112, -1.215, -1.0, -1.0, 0.893, -0.894],
            'ATT': [-2.126, -2.585, -1.732, -1.9909999999999999, -2.591, -2.591, 2.141, 2.158, -1.0, -1.0, -2.004, 2.032],
            'CAA': [0.265, -0.231, -0.5770000000000001, -0.74, -0.226, -0.226, 0.166, -0.275, -1.0, -1.0, -0.748, 0.743],
            'CAC': [0.496, 0.7859999999999999, 0.5770000000000001, 0.81, 0.773, 0.773, -0.6459999999999999, -0.503, -1.0, -1.0, 0.797, -0.8],
            'CAG': [1.5759999999999998, 0.92, 0.5770000000000001, -0.322, 0.92, 0.92, -0.762, -1.5490000000000002, -1.0, -1.0, -0.314, 0.304],
            'CAT': [1.229, 1.348, -0.5770000000000001, 0.87, 1.358, 1.358, -1.112, -1.215, -1.0, -1.0, 0.893, -0.894],
            'CCA': [-1.8559999999999999, -1.14, 0.5770000000000001, 0.27399999999999997, -1.139, -1.139, 0.917, 1.876, 1.0, 1.0, 0.266, -0.275],
            'CCC': [0.07200000000000001, 0.358, 1.732, 0.5720000000000001, 0.345, 0.345, -0.3, -0.084, 1.0, 1.0, 0.555, -0.562],
            'CCG': [-0.9690000000000001, -0.7120000000000001, 1.732, -0.084, -0.705, -0.705, 0.5579999999999999, 0.9620000000000001, 1.0, 1.0, -0.07200000000000001, 0.062],
            'CCT': [-0.313, -0.07, 0.5770000000000001, 0.27399999999999997, -0.07, -0.07, 0.039, 0.3, 1.0, 1.0, 0.266, -0.275],
            'CGA': [0.111, 1.0, 0.5770000000000001, 1.645, 1.012, 1.012, -0.8340000000000001, -0.12300000000000001, 1.0, 1.0, 1.666, -1.646],
            'CGC': [-0.46799999999999997, 0.385, 1.732, 1.287, 0.379, 0.379, -0.326, 0.455, 1.0, 1.0, 1.28, -1.272],
            'CGG': [-0.9690000000000001, -0.7120000000000001, 1.732, -0.084, -0.705, -0.705, 0.5579999999999999, 0.9620000000000001, 1.0, 1.0, -0.07200000000000001, 0.062],
            'CGT': [-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.065, 0.065, -0.07400000000000001, 0.107, 1.0, 1.0, 0.266, -0.275],
            'CTA': [0.882, -0.09699999999999999, -0.5770000000000001, -1.276, -0.09699999999999999, -0.09699999999999999, 0.062, -0.88, -1.0, -1.0, -1.28, 1.285],
            'CTC': [0.419, 0.43799999999999994, 0.5770000000000001, 0.27399999999999997, 0.42700000000000005, 0.42700000000000005, -0.365, -0.42700000000000005, -1.0, -1.0, 0.266, -0.275],
            'CTG': [1.5759999999999998, 0.92, 0.5770000000000001, -0.322, 0.92, 0.92, -0.762, -1.5490000000000002, -1.0, -1.0, -0.314, 0.304],
            'CTT': [-0.506, -0.257, -0.5770000000000001, 0.155, -0.26, -0.26, 0.193, 0.493, -1.0, -1.0, 0.16899999999999998, -0.179],
            'GAA': [-0.159, -0.605, -0.5770000000000001, -0.9179999999999999, -0.6, -0.6, 0.474, 0.146, -1.0, -1.0, -0.893, 0.89],
            'GAC': [0.034, 0.171, 0.5770000000000001, 0.27399999999999997, 0.17800000000000002, 0.17800000000000002, -0.165, -0.046, -1.0, -1.0, 0.266, -0.275],
            'GAG': [0.419, 0.43799999999999994, 0.5770000000000001, 0.27399999999999997, 0.42700000000000005, 0.42700000000000005, -0.365, -0.42700000000000005, -1.0, -1.0, 0.266, -0.275],
            'GAT': [-0.737, -0.391, -0.5770000000000001, 0.214, -0.397, -0.397, 0.307, 0.727, -1.0, -1.0, 0.217, -0.22699999999999998],
            'GCA': [0.7659999999999999, 0.8390000000000001, 0.5770000000000001, 0.5720000000000001, 0.8420000000000001, 0.8420000000000001, -0.7020000000000001, -0.767, 1.0, 1.0, 0.555, -0.562],
            'GCC': [1.036, 2.097, 1.732, 2.479, 2.089, 2.089, -1.6869999999999998, -1.0290000000000001, 1.0, 1.0, 2.487, -2.4330000000000003],
            'GCG': [-0.46799999999999997, 0.385, 1.732, 1.287, 0.379, 0.379, -0.326, 0.455, 1.0, 1.0, 1.28, -1.272],
            'GCT': [0.304, 0.92, 0.5770000000000001, 1.287, 0.9109999999999999, 0.9109999999999999, -0.753, -0.313, 1.0, 1.0, 1.28, -1.272],
            'GGA': [0.265, -0.09699999999999999, 0.5770000000000001, -0.501, -0.10300000000000001, -0.10300000000000001, 0.066, -0.275, 1.0, 1.0, -0.507, 0.499],
            'GGC': [1.036, 2.097, 1.732, 2.479, 2.089, 2.089, -1.6869999999999998, -1.0290000000000001, 1.0, 1.0, 2.487, -2.4330000000000003],
            'GGG': [0.07200000000000001, 0.358, 1.732, 0.5720000000000001, 0.345, 0.345, -0.3, -0.084, 1.0, 1.0, 0.555, -0.562],
            'GGT': [-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.071, 0.071, -0.078, 0.107, 1.0, 1.0, 0.266, -0.275],
            'GTA': [0.342, -0.07, -0.5770000000000001, -0.561, -0.062, -0.062, 0.031, -0.35100000000000003, -1.0, -1.0, -0.555, 0.5479999999999999],
            'GTC': [0.034, 0.171, 0.5770000000000001, 0.27399999999999997, 0.17800000000000002, 0.17800000000000002, -0.165, -0.046, -1.0, -1.0, 0.266, -0.275],
            'GTG': [0.496, 0.7859999999999999, 0.5770000000000001, 0.81, 0.773, 0.773, -0.6459999999999999, -0.503, -1.0, -1.0, 0.797, -0.8],
            'GTT': [-1.5090000000000001, -1.354, -0.5770000000000001, -0.561, -1.3630000000000002, -1.3630000000000002, 1.105, 1.516, -1.0, -1.0, -0.555, 0.5479999999999999],
            'TAA': [0.6890000000000001, -0.284, -1.732, -1.395, -0.275, -0.275, 0.20600000000000002, -0.6920000000000001, -1.0, -1.0, -1.376, 1.3840000000000001],
            'TAC': [0.342, -0.07, -0.5770000000000001, -0.561, -0.062, -0.062, 0.031, -0.35100000000000003, -1.0, -1.0, -0.555, 0.5479999999999999],
            'TAG': [0.882, -0.09699999999999999, -0.5770000000000001, -1.276, -0.09699999999999999, -0.09699999999999999, 0.062, -0.88, -1.0, -1.0, -1.28, 1.285],
            'TAT': [1.615, 0.5720000000000001, -1.732, -0.978, 0.584, 0.584, -0.491, -1.585, -1.0, -1.0, -0.99, 0.988],
            'TCA': [1.73, 1.348, -0.5770000000000001, 0.27399999999999997, 1.348, 1.348, -1.103, -1.696, 1.0, 1.0, 0.266, -0.275],
            'TCC': [0.265, -0.09699999999999999, 0.5770000000000001, -0.501, -0.10300000000000001, -0.10300000000000001, 0.066, -0.275, 1.0, 1.0, -0.507, 0.499],
            'TCG': [0.111, 1.0, 0.5770000000000001, 1.645, 1.012, 1.012, -0.8340000000000001, -0.12300000000000001, 1.0, 1.0, 1.666, -1.646],
            'TCT': [0.381, -0.15, -0.5770000000000001, -0.74, -0.158, -0.158, 0.109, -0.389, 1.0, 1.0, -0.748, 0.743],
            'TGA': [1.73, 1.348, -0.5770000000000001, 0.27399999999999997, 1.348, 1.348, 4.522, -1.696, 1.0, 1.0, 0.266, -0.275],
            'TGC': [0.7659999999999999, 0.8390000000000001, 0.5770000000000001, 0.5720000000000001, 0.8420000000000001, 0.8420000000000001, -0.7020000000000001, -0.767, 1.0, 1.0, 0.555, -0.562],
            'TGG': [-1.8559999999999999, -1.14, 0.5770000000000001, 0.27399999999999997, -1.139, -1.139, 0.917, 1.876, 1.0, 1.0, 0.266, -0.275],
            'TGT': [0.111, 0.171, -0.5770000000000001, 0.155, 0.16399999999999998, 0.16399999999999998, -0.153, -0.12300000000000001, 1.0, 1.0, 0.16899999999999998, -0.179],
            'TTA': [0.6890000000000001, -0.284, -1.732, -1.395, -0.275, -0.275, 0.20600000000000002, -0.6920000000000001, -1.0, -1.0, -1.376, 1.3840000000000001],
            'TTC': [-0.159, -0.605, -0.5770000000000001, -0.9179999999999999, -0.6, -0.6, 0.474, 0.146, -1.0, -1.0, -0.893, 0.89],
            'TTG': [0.265, -0.231, -0.5770000000000001, -0.74, -0.226, -0.226, 0.166, -0.275, -1.0, -1.0, -0.748, 0.743],
            'TTT': [-2.0869999999999997, -2.745, -1.732, -2.349, -2.7439999999999998, -2.7439999999999998, -2.615, 2.1180000000000003, -1.0, -1.0, -2.342, 2.386],
        })

        # Generate rs and initialize sq, phy
        self.rs  = self.allp(2)
        self.phy = list(dinuc_data['Physicochemical properties'])

        # Create the temp dictionary
        self.temp = defaultdict(list)
        for p in self.phy:
            for k in self.rs:
                for idx, prop_value in enumerate(dinuc_data['Physicochemical properties']):
                    if prop_value == p:
                        self.temp[k].append(dinuc_data[k][idx])
        
        # extract features from sequence
        cdk = self.DACC(cdk)
        cdk = self.DCC(cdk)
        cdk = self.DDON(cdk)
        cdk = self.MAC(cdk)
        cdk = self.NMBAC(cdk)
        cdk = self.pcpsetnc(cdk)
        
        return cdk

    # def End
    
    def DACC(self, cdk:pd.DataFrame, lagvalue:int=2):
        
        # Precompute lengths and caches for performance
        seq_lengths = [len(seq) for seq in self.list_seq]
        phy_len = len(self.phy)
        lag_range = range(1, lagvalue + 1)
    
        # DAC and DCC calculation (combined optimization)
        dac = defaultdict(list)
        dcc = defaultdict(list)
        selected_feature = [1460, 1556, 1938, 2013, 2798, 2795]

        # Precompute the PU (mean of temp) for each sequence
        pu_dict = {}
        for seq, seq_len in zip(self.list_seq, seq_lengths):
            pu = [sum(self.temp[k][i] for k in self.temp.keys()) / seq_len for i in range(phy_len)]
            pu_dict[seq] = pu

        # Calculate DAC and DCC values
        for seq, seq_len in zip(self.list_seq, seq_lengths):
            pu = pu_dict[seq]
            mer = self.kmer(2, seq)

            # DAC Calculation
            count = 0
            for l in lag_range:
                for p in range(phy_len):
                    count += 1
                    if count in selected_feature:
                        av = sum(self.temp[mer[i]][p] for i in range(seq_len - l - 1)) / (seq_len - l - 1)  # Use seq_len - l - 1
                        sum_dacc = sum((self.temp[mer[i]][p] - av) * self.temp[mer[i + l]][p] for i in range(seq_len - l - 1)) / (seq_len - l - 1)
                        dac[f"DACC_{count}"].append(sum_dacc)
                    else: 
                        continue

            # DCC Calculation
            count1 = count
            for l in lag_range:
                for p1 in range(phy_len):
                    for p2 in range(phy_len):
                        if p1 != p2:
                            count1 += 1
                            if count1 in selected_feature:
                                av1 = sum(self.temp[mer[i]][p1] for i in range(seq_len - l - 1)) / (seq_len - l - 1)
                                av2 = sum(self.temp[mer[i]][p2] for i in range(seq_len - l - 1)) / (seq_len - l - 1)
                                sum_dcc = sum((self.temp[mer[i]][p1] - av1) * (self.temp[mer[i + l]][p2] - av2) for i in range(seq_len - l - 1)) / (seq_len - l - 1)
                                dcc[f"DACC_{count1}"].append(sum_dcc)
                            else:
                                continue
                            
        # Append DAC and DCC values to DataFrame
        dac_df = pd.DataFrame.from_dict(dac)
        dcc_df = pd.DataFrame.from_dict(dcc)
        dacc_df = pd.concat([dac_df, dcc_df], axis=1)[['DACC_1460', 'DACC_1556','DACC_1938','DACC_2013','DACC_2798','DACC_2795']]
        cdk = pd.concat([cdk, dacc_df], axis=1)

        return cdk

    def DCC(self, cdk:pd.DataFrame, lagvalue:int=2):
        
        # calculating dcc value
        dcc= defaultdict(list)
        selected_feature = [1384, 1580, 1862, 1937, 2135, 2719, 2722]
        
        for seq in self.list_seq:
            pu= []
            for i in range(len(self.phy)):
                mean=0.0
                for j in self.temp.keys():
                    mean += self.temp[j][i]
                mean = mean/len(seq)
                pu.append(mean)  
            mer=self.kmer(2,seq)   
            count=0
            for l in range(1,lagvalue+1):
                for p1 in range(len(self.phy)):
                    for p2 in range(len(self.phy)):
                        if(p1!=p2):
                            count = count+1
                            if count not in selected_feature:
                                continue
                            
                            else:
                                av1=0.0
                                av2=0.0
                                for i in range(len(seq)-l-1):
                                    av1 += self.temp[mer[i]][p1]
                                    av2 += self.temp[mer[i]][p2]
                                av1= av1/len(seq)
                                av2=av2/len(seq)
                                sum=0.0
                                st = "DCC_"+ str(count)
                                for i in range(len(seq)-l-1):
                                    pu1 = self.temp[mer[i]][p1]-av1
                                    pu2 = self.temp[mer[i+l]][p2]-av2
                                    sum += (pu1*pu2)
                                sum = sum/(len(seq)-l-1)
                                dcc[st].append(sum)
                            
        # Appending final output to the output file
        dcc_df = pd.DataFrame.from_dict(dcc)
        cdk = pd.concat([cdk, dcc_df], axis=1)
        
        return cdk



    def DDON(self, cdk:pd.DataFrame):
        
        ddon = defaultdict(list)
        # Original code extract ATGC, but used only feature from G
        # alphabet=['A','C','G','T']
        alphabet=['G',]
            
        for seq in self.list_seq:
            for i in alphabet:
                c=0
                s1=0
                s2=0
                j=0
                o=0
                while(j<(len(seq))):
                    if(seq[j]!=i):
                        c=c+1
                    elif(seq[j]==i):
                        o=o+1
                        s1 = (c**2)+s1
                        c=0
                    j=j+1
                if(c!=0):
                    s1 = s1+(c**2)
                    c=0
                s2 = (len(seq)-o+1)
                t = s1/s2
                if(o==0):
                    ddon["DDON_"+i].append(0.0)
                else:
                    ddon["DDON_"+i].append(t)  

        ddon_df = pd.DataFrame.from_dict(ddon)
        cdk = pd.concat([cdk, ddon_df], axis=1)
            
        return cdk


    def MAC(self, cdk:pd.DataFrame, lagvalue:int=2, kvalue:int=2):
  
        mac= defaultdict(list)
        k = kvalue
        selected_feature = [54]
        
        for seq in self.list_seq:
            pu= []
            for i in range(len(self.phy)):
                mean=0.0
                for j in self.temp.keys():
                    mean += self.temp[j][i]
                mean = mean/len(seq)
                pu.append(mean)  
            mer=self.kmer(k,seq)   
            count=0
            for l in range(1,lagvalue+1):
                motif_len = len(seq)-l-k+1

                for p in range(len(self.phy)):
                    count = count+1
                    if count not in selected_feature:
                        continue

                    av=0.0
                    su=0.0
                    su2=0.0
                    
                    for i in range(motif_len):
                        av += self.temp[mer[i]][p]
                    av= av/len(self.list_seq[0]) # ????? maybe mistake?
                    st = "MAC_"+ str(count)
                    
                    for i in range(motif_len):
                        pu1 = self.temp[mer[i]][p]-av
                        pu2 = self.temp[mer[i+l]][p]-av
                        su += (pu1*pu2)
                    for i in range(len(seq)-k+1):
                        pu3 = self.temp[mer[i]][p]-av
                        su2 +=(pu3*pu3) 
                    su = ((1/(motif_len))*su)/((1/(len(seq)-k+1))*(su2))
                    mac[st].append(su)
        
        mac_df = pd.DataFrame.from_dict(mac)
        cdk    = pd.concat([cdk, mac_df], axis=1)
        
        return cdk



    def NMBAC(self, cdk:pd.DataFrame, lagvalue:int=2, kvalue:int=2):
  
        nmbac= defaultdict(list)
        k = kvalue
        selected_feature = [47, 64]
        
        for seq in self.list_seq:
            pu= []
            for i in range(len(self.phy)):
                mean=0.0
                for j in self.temp.keys():
                    mean += self.temp[j][i]
                mean = mean/len(seq)
                pu.append(mean)  
            mer   = self.kmer(k,seq)   
            count = 0
            for l in range(1,lagvalue+1):

                motif_len = len(seq)-l-k+1

                for p in range(len(self.phy)):
                    count = count+1
                    if count not in selected_feature:
                        continue
                    
                    av=0.0
                    su=0.0
                    for i in range(motif_len):
                        av += self.temp[mer[i]][p]
                    av = av/len(seq)
                    st = "NMBAC_"+ str(count)
                    for i in range(motif_len):
                        pu1 = self.temp[mer[i]][p]
                        pu2 = self.temp[mer[i+l]][p]
                        su += (pu1*pu2)**2
                    
                    su = (su)/(motif_len)
                    nmbac[st].append(su)
                    
        nmbac_df = pd.DataFrame.from_dict(nmbac)
        cdk      = pd.concat([cdk, nmbac_df], axis=1)
            
        return cdk
    
    def pcpsetnc(self, cdk:pd.DataFrame, kvalue:int=3, wvalue:int=0.05, lmvalue:int=1):

        k  = kvalue
        rs = self.allp(3) 
        rs.sort()
        ph_v = defaultdict(list)
        res  = defaultdict(list)
        data = self.trinuc_data
        selected_feature = [19, 42, 46]
        
        for km in rs:
            ph_v[km] = list(data[km])
                
        for s in self.list_seq:
            if len(s) < k or lmvalue + k > len(s):
                continue
            mer = self.kmer(3,s)
            fre = [mer.count(i) for i in rs]
            fre_sum = sum(fre)
            fre = [(f/fre_sum) for f in fre]
            
            theta =[]
            pyn = len(ph_v[rs[0]])
            for i in range(1,lmvalue+1):
                temp_sum =0.0
                for j in range(len(s)-k-i+1):
                    n1 = mer[j]
                    n2 = mer[j+i]
                    temp = 0.0
                    for y in range(pyn):
                        temp += ((ph_v[n1][y]-ph_v[n2][y])**2)
                    temp = temp/pyn
                    temp_sum += temp
                temp_sum = (temp_sum/(len(s)-i-k+1))
                theta.append(temp_sum)

            t_sum = sum(theta)
            dm = 1 + wvalue*t_sum
            temp_vec =[f/dm for f in fre]

            for i in theta:
                temp_vec.append(wvalue*i/dm)

            for i in range(1,len(temp_vec)+1):
                if i not in selected_feature: continue
                st = "PC_PTNC_"+str(i)
                res[st].append(temp_vec[i-1])
                
        res_df = pd.DataFrame.from_dict(res)
        cdk    = pd.concat([cdk, res_df], axis=1)
        
        return cdk


    def motif_features(self, cdk:pd.DataFrame):
            
        #creating the first dataframe for initial columns (features)
        list_of_col_names=["f"+""+str(i) for i in range(0, 1598)]
        df_first=pd.DataFrame(columns=list_of_col_names)

        for i in range(len(self.list_seq)):
            list_to_add=[]
            nuc=str(self.list_seq[i])
            
            # Count all nucleotide occurrences
            nucleotide_counts = Counter(nuc)
            list_to_add.extend([nucleotide_counts.get(base, 0) for base in 'ATGC'])
            
            # Count dimer and trimer motif
            list_dim_wg  = self.calculate_dimer_counts(nuc)
            list_trim_wg = self.calculate_trimer_counts(nuc)
            list_dimer_gaps = [self.calculate_dimer_counts(nuc, gap=g) for g in range(1, 6)]
            list_to_add += list_dim_wg + list_trim_wg + sum(list_dimer_gaps, [])
            
            # Calculate trimer counts for rgap and lgap in one loop each
            gap_values = [1, 2, 3, 7, 8, 9, 10, 15, 16, 17]
            list_trimer_rgap = [self.calculate_trimer_counts(nuc, gap=gap, frame=1) for gap in gap_values]
            list_trimer_lgap = [self.calculate_trimer_counts(nuc, gap=gap, frame=2) for gap in gap_values]
            list_to_add += sum(list_trimer_rgap + list_trimer_lgap, [])

            # Define patterns to count
            patterns = [
                "TTGAC", "TATAAT", "TTATAA", "TTGACA", "AACGAT", "ACAGTT", "AGGAGG", "TAAAAT", "TTGATT"
            ]

            # Calculate counts for patterns and skew
            pattern_counts = [self.pattern(nuc, pat) for pat in patterns]
            skew_counts = [self.GC_skew(nuc), self.AT_skew(nuc)]
            list_to_add.extend(pattern_counts + skew_counts)

            list_new  = self.pstnp(nuc)
            list_new2 = self.eiip(nuc)
            list_to_add+=list_new+list_new2

            final_list=np.array(list_to_add)
            df_first.loc[i]=final_list

        features_list2=["f"+str(i) for i in list_best_features if i< 1598]

        df_first = df_first[features_list2]
        
        return df_first
    
    def calculate_dimer_counts(self, nuc, gap=0):
        """
        Calculate dimer counts for a given nucleotide sequence.
        Supports both contiguous dimers and gapped dimers.

        Parameters:
        - nuc (str): Nucleotide sequence.
        - gap (int): Gap between the two nucleotides in a dimer. Default is 0 (contiguous).

        Returns:
        - list: A list of counts for each dimer.
        """
        dim = [0] * 16
        dimers_dict = make_dimers_dict()

        for i in range(len(nuc) - (gap + 1)):
            dimer = nuc[i] + nuc[i + gap + 1]
            pos = dimers_dict[dimer]
            dim[pos] += 1

        return dim

    def calculate_trimer_counts(self, nuc, gap=0, frame=1):
        """
        Calculate trimer counts for a given nucleotide sequence.
        Supports contiguous and gapped trimers.

        Parameters:
        - nuc (str): Nucleotide sequence.
        - gap (int): Gap between the nucleotides in a trimer. Default is 0 (contiguous).
        - frame (int): Frame of position of trimer for making gap. Default is 1. 1 or 2

        Returns:
        - list: A list of counts for each trimer.
        """
        tri = [0] * 64
        trimers_dict = make_trimers_dict()
        
        if frame == 1: 
            for i in range(len(nuc) - (gap + 2)):
                trimer = nuc[i] + nuc[i + 1] + nuc[i + gap + 2]
                pos = trimers_dict[trimer]
                tri[pos] += 1
                
        elif frame == 2: 
            for i in range(len(nuc) - (gap + 2)):
                trimer = nuc[i] + nuc[i + gap + 1] + nuc[i + gap + 2]
                pos = trimers_dict[trimer]
                tri[pos] += 1
        else:
            raise ValueError('Not available frame.')

        return tri

    def pattern(self, nuc,pattern):
        return len(re.findall(pattern, nuc))

    def GC_skew(self, nuc):
        g=nuc.count('G')
        c=nuc.count('C')
        return (c-g)/(c+g)

    def AT_skew(self, nuc):
        a=nuc.count('A')
        t=nuc.count('T')
        return (a-t)/(a+t)

    def pstnp(self, seq):
        return [dict_pstnp[seq[i:i+3]][i] for i in range(79)]

    def eiip(self, seq):
        
        dict_val={'A': 0.1260,'C': 0.1335, 'G': 0.0806,'T': 0.1340}
        
        dict_tmp={}
        dict_out={}
        dict_pos_tri = make_trimers_dict()
        
        for pat in dict_pos_tri.keys():
            dict_tmp[pat]=0
            dict_out[pat]=0
        for i in range(79):
            var = seq[i:i+3]
            dict_tmp[var]+=1
        for i in range(79):
            var = seq[i:i+3]
            etemp=dict_val[seq[i]]+dict_val[seq[i+1]]+dict_val[seq[i+2]]

            res=etemp*(dict_tmp[var])/79
            dict_out[var]=res
            
        return [dict_out[k] for k in dict_out]
