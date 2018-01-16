# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:34:56 2016

@author: daniellefever
"""

#!/usr/bin/python

from __future__ import division

import numpy 
import os

from prody import parsePDB, GNM
from numpy import sqrt, diff, sign, where
from collections import OrderedDict
from itertools import izip
from Bio.SeqUtils import GC



def gene_id_finder(feature_seqrec):
    for db_xref_item in feature_seqrec.qualifiers['db_xref']:
        if 'GeneID' in db_xref_item:
            gene_id_num = db_xref_item.split(':')[1]
    return int(gene_id_num)


def full_fp_list(input_dir):
    file_list = []
    for file in os.listdir(input_dir):
        full_fp = os.path.join(input_dir, file)
        file_list.append(full_fp)
    return file_list



def seq_streams(parsed_gb, feature_seqrec, consec_bases = 750):
    start_chrome_ind = feature_seqrec.location.start.position
    end_chrome_ind = feature_seqrec.location.end.position

    start_upseq = parsed_gb[(start_chrome_ind - consec_bases):start_chrome_ind]
    start_downseq = parsed_gb[(start_chrome_ind + 3):(start_chrome_ind + consec_bases + 3)]
    end_upseq = parsed_gb[(end_chrome_ind - consec_bases):end_chrome_ind]
    end_downseq = parsed_gb[(end_chrome_ind + 3):(end_chrome_ind + consec_bases + 3)]

    return(start_upseq, start_downseq, end_upseq, end_downseq)



def feature_getter(parsed_gb, feature_string='CDS', spec_gene=None):
    cds_seqs = []
    for index, cur_feature in enumerate(parsed_gb.features):
        if cur_feature.type == feature_string:
            if spec_gene is None:
                extracted_feature = parsed_gb.features[index]
                feature_seqrec = extracted_feature.extract(parsed_gb.seq)
                #cds_seqs.append((extracted_feature,feature_seqrec))
                cds_seqs.append(extracted_feature)
            if spec_gene is not None:
                gene_id = cur_feature.qualifiers['gene']
                if spec_gene in gene_id[0]:
                    extracted_feature = parsed_gb.features[index]
                    feature_seqrec = extracted_feature.extract(parsed_gb.seq)
                    cds_seqs.append(extracted_feature)
                    #cds_seqs.append((extracted_feature,feature_seqrec))
            
    return(cds_seqs) 
    
    


def seq_parser(seqfile, file_type = 'genbank'):
    with open(seqfile, 'rU') as f:
        parsed_seqfile = SeqIO.read(seqfile, file_type)
    return parsed_seqfile 


def file_to_dict(file_name):
    '''
    Function to turn a file with a key value
    pair per line
    '''
    cur_dict = {}
    with open(file_name, 'rU') as f:
        for line in f:
            k, v = line.split('\t')
            cur_dict[k.strip(',')] = v.strip('\n')
        return cur_dict   

def file_len(fname):
    '''
    Taken from 
    http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    '''
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def all_hinge_setter(hingedict):
    hinge_set = set()
    for cur_set in hingedict.itervalues():
        for val in cur_set:
            hinge_set.add(val)
    return hinge_set

def mode_str_maker(res_index, hingedict):
        mode_id = ''
        for mode_num, hinge_set in hingedict.iteritems():
            if res_index in hinge_set:
                mode_id = mode_id + ' ' + str(mode_num)
        return mode_id
                   
def res_set_maker(parsed_pdb):
    res_index_set = set()
    for i in parsed_pdb.getResindices():
        res_index_set.add(i)
    return res_index_set

def hingedict_maker(gnm):
        hinge_dict = OrderedDict()
        for mode in gnm[:3]:
            mode_str = str(mode).split()
            hinge_dict[mode_str[0] + '_' + mode_str[1]] = set(hinger(mode.getArray()))
        return hinge_dict

def hinger(gnm_array):
    cur_hinges = []
    hinge_loc_list = list(where(diff(sign(gnm_array)))[0])
    for index, hinge_loc in enumerate(hinge_loc_list):
        try:
            consec_hinge_loc = hinge_loc_list[(index + 1)]
            if (consec_hinge_loc - hinge_loc) > 4:
                cur_hinges.append(hinge_loc)
                cur_hinges.append((hinge_loc + 1))
        except IndexError:
            cur_hinges.append(hinge_loc)
            cur_hinges.append((hinge_loc + 1))

    return cur_hinges

def lig_dist(pdb_prot, pdb_lig):
    '''
    Calculates the distance between each atom 
    of the ligand to each atom of each residue
    '''
    dist_dict = ordered_distdict(pdb_prot)
    for prot_coords, prot_index in izip(pdb_prot.getCoords(), pdb_prot.getResindices()):
        dist_list = []
        for lig_coords in pdb_lig.getCoords():
            coord_diff = lig_coords - prot_coords
            coord_dist = round(sqrt(sum(coord_diff**2)), 3)
            dist_list.append(coord_dist)
        
        min_dist = min(dist_list)
        dist_dict[prot_index].append(min_dist)

    closest_distdict = closest_atom_dist(dist_dict, pdb_prot)
    return(closest_distdict)
 

def closest_atom_dist(dist_dict, pdb_prot):
    closest_distdict = ordered_distdict(pdb_prot)
    for prot_index, dist_list in dist_dict.iteritems():
        closest_distdict[prot_index] = min(dist_list)
    return closest_distdict

def ordered_distdict(pdb_prot):
    dist_dict = OrderedDict()
    for i in pdb_prot.getResindices():
        dist_dict[i] = []
    return dist_dict

def gnm_maker(protein_name, parsed_pdb):
    gnm = GNM(protein_name)
    gnm.buildKirchhoff(parsed_pdb)
    gnm.calcModes()
    return gnm

def ligbound_resfinder(closest_distdict):
    bound_residues = []
    for res_index, dist_val in closest_distdict.iteritems():
        if dist_val <= 4.5:
            bound_residues.append(res_index)
    return bound_residues
