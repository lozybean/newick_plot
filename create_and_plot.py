#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2016-01-08'
This create newick format from qiime tax assignment file and otu profile
"""

from __future__ import division
import re
import sys
import os
import math
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace, TreeNode, ClusterTree, PieChartFace
import pandas as pd
import numpy as np
import scipy as sp
import argparse
from collections import OrderedDict
from tree import TreeNode as Node
from tree import Tree as TaxTree

COLS_BREWER = ['#00447E', '#F34800', '#64A10E', '#930026', '#464E04', '#049a0b', '#4E0C66', '#D00000', '#FF6C00',
               '#FF00FF', '#c7475b', '#00F5FF', '#BDA500', '#A5CFED', '#f0301c', '#2B8BC3', '#FDA100', '#54adf5',
               '#CDD7E2', '#9295C1']

LEVELS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

COLOR_DICT = OrderedDict()
for ind, level in enumerate(LEVELS):
    COLOR_DICT[level] = COLS_BREWER[ind % 20]

TAX_DICT = {
    'k': 'kingdom',
    'p': 'phylum',
    'c': 'class',
    'o': 'order',
    'f': 'family',
    'g': 'genus',
    's': 'species',
}

SEARCH_REG = {
    'k': re.compile('(k__[^;]+)'),
    'p': re.compile('(p__[^;]+)'),
    'c': re.compile('(c__[^;]+)'),
    'o': re.compile('(o__[^;]+)'),
    'f': re.compile('(f__[^;]+)'),
    'g': re.compile('(g__[^;]+)'),
    's': re.compile('(s__[^;]+)'),
}


def read_params(args):
    parser = argparse.ArgumentParser(description='plot tax tree| v1.0 at 2015/01/08 by liangzb')
    parser.add_argument('--profile', dest='profile', metavar='FILE', type=str, required=True,
                        help="set the profile file")
    parser.add_argument('--tax_ass', dest='tax_ass', metavar='FILE', type=str, required=True,
                        help="set the tax_assignment.txt")
    parser.add_argument('-g', '--group_file', dest='group', metavar='FILE', type=str, default=None,
                        help="set the group file if need")
    parser.add_argument('--top', dest='top', metavar='INT', type=int, default=20,
                        help="set the top num, [default is 20]")
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='DIR', type=str, required=True,
                        help="set the outdir")
    parser.add_argument('--with_leaf_pie', dest='with_leaf_pie', action='store_true',
                        help="plot with leaf pie chart, [default is false]")
    parser.set_defaults(with_leaf_pie=False)
    parser.add_argument('--with_branch_circle', dest='with_branch_circle', action='store_true',
                        help="plot with branch circle, [default is false]")
    parser.set_defaults(with_branch_circle=False)

    args = parser.parse_args()
    params = vars(args)
    params['group'] = parse_group_file(params['group'])
    return params


def table_uniform(table):
    for col_name, col in table.iteritems():
        table[col_name] = col.map(lambda s: s / col.sum())
    return table


def parse_group_file(file):
    if file is None:
        return None
    group = OrderedDict()
    with open(file) as g:
        for line in g:
            tabs = line.strip().split('\t')
            if len(tabs) >= 2:
                group[tabs[0]] = tabs[1]
            else:
                group[tabs[0]] = tabs[0]
    return group



def mkdir(dirname):
    if os.path.isdir(dirname):
        return None
    if dirname == '':
        return None
    os.mkdir(dirname)
    return None


def float_trans(string):
    try:
        num = float(string)
    except ValueError:
        return string
    if num >= 0.01:
        value = format(num, '.2f')
    else:
        value = format(num, '.2e')
    if num == 0:
        value = 0
    return value


def get_tax(level=None):
    def f(string):
        if level is None:
            return np.nan
        try:
            result = SEARCH_REG[level[0].lower()].search(string).group(1)
        except AttributeError:
            return np.nan
        return result

    return f


def get_level_profile(otu_profile, level):
    level_profile = otu_profile.loc[level.index]
    level_profile.index = level
    result = {}
    for p in level_profile.index.unique():
        if level_profile.loc[p].__class__ == pd.DataFrame:
            s = level_profile.loc[p].sum()
            s.name = p
        else:
            s = level_profile.loc[p]
        result[p] = s
    result = pd.DataFrame(result).T
    return result


def pick_high_abundance(level, otu_profile, top):
    def filt_fun(n):
        if n in picked_level:
            return n
        else:
            return np.nan

    genus_profile = get_level_profile(otu_profile, level)
    genus_profile = genus_profile.T.sum()
    picked_level = genus_profile.sort(inplace=False)[:-top - 1:-1].index
    level = level.map(filt_fun).dropna()
    return otu_profile.loc[level.index]


def read_profile(tax_ass, profile, top):
    tax_table = pd.DataFrame.from_csv(tax_ass, sep='\t', header=None)
    tax_table.columns = pd.Index(['tax_line', 'confidence'])
    tax_line = tax_table['tax_line']
    level_dict = {}
    for level in LEVELS:
        level_dict[level] = tax_line.map(get_tax(level=level)).dropna()
    otu_profile = pd.DataFrame.from_csv(profile, sep='\t')
    otu_profile = pick_high_abundance(level_dict['genus'], otu_profile, top)
    level_profile = {}
    total_profile = pd.DataFrame()
    for level in LEVELS:
        level_profile[level] = get_level_profile(otu_profile, level=level_dict[level]).dropna()
        level_profile[level] = table_uniform(level_profile[level])
        total_profile = pd.concat([total_profile, level_profile[level]])
    return level_profile, total_profile


def modify_tax_ass(tax_ass, outfile):
    with open(tax_ass) as fp, open(outfile, 'w') as out:
        for line in fp:
            tabs = line.rstrip().split('\t')
            tax_line = tabs[1]
            tax_list = tax_line.split(';')
            for ind, level in enumerate(LEVELS):
                try:
                    tax_name = tax_list[ind]
                except IndexError:
                    break
                if level[0] != tax_name[0]:
                    tax_list.insert(ind, '%s%s_unidentified' % (level[0], tax_name[1:13]))
            out.write('%s\t%s\t%s\n' % (tabs[0], ';'.join(tax_list), tabs[2]))


def read_tax(tax_ass, profile_dict):
    tree = TaxTree()
    with open(tax_ass) as fp:
        for line in fp:
            tabs = line.rstrip().split('\t')
            tax_line = tabs[1]
            tax_list = tax_line.split(';')
            for ind, tax_name in enumerate(tax_list):
                if tax_name not in profile_dict:
                    continue
                tax_level = TAX_DICT[tax_name[0]]
                if tax_name not in tree.nodes:
                    tree.nodes[tax_name] = Node(level=tax_level, name=tax_name)
                node_present = tree.nodes[tax_name]
                if ind > 0:
                    node_previous = tree.nodes[tax_list[ind - 1]]
                else:
                    node_previous = tree.root
                node_previous.add_child(child=node_present)
                node_present.profile = profile_dict[node_present.name]
                node_present.tree = tree
    return tree


def get_tree_style():
    ts = TreeStyle()
    # ts.mode = 'c'
    ts.margin_top = 10
    ts.margin_bottom = 10
    ts.margin_left = 10
    ts.margin_right = 10
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_scale = False
    title = TextFace("     Tax Assignment Tree", fsize=10)
    title.hz_align = 2
    title.vt_align = 2
    ts.title.add_face(TextFace(" "), column=0)
    ts.title.add_face(TextFace(" "), column=0)
    ts.title.add_face(title, column=0)
    return ts


def set_node_default(tree, node_dict):
    root_node = tree & 'root'
    for node in tree.traverse():
        if node is root_node:
            continue
        node = node_dict[node.name]
        nstyle = NodeStyle()
        node.style = nstyle
        nstyle["shape"] = "circle"
        nstyle["hz_line_type"] = 3
        nstyle["hz_line_color"] = "black"
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "black"
        nstyle["vt_line_width"] = 1


def add_node_circle(tree, node_dict):
    root_node = tree & 'root'
    for node in tree.traverse():
        if node is root_node:
            continue
        nstyle = node_dict[node.name].style
        nstyle["size"] = node_dict[node.name].size
        try:
            nstyle["fgcolor"] = COLOR_DICT[node_dict[node.name].level]
        except KeyError:
            nstyle["fgcolor"] = "grey"


def add_branch_text(tree, tree_style, node_dict):
    root_node = tree & 'root'
    for node in tree.traverse():
        if node is root_node:
            continue
        T1 = TextFace(node.name, ftype="Monaco", fsize=10)
        T1.hz_align = 0
        node.add_face(T1, 0, 'branch-top')
        T2 = TextFace('%s%%' % float_trans(node_dict[node.name].profile * 100), ftype="Monaco", fsize=10)
        T2.hz_align = 0
        node.add_face(T2, 0, 'branch-bottom')
        # print node_dict[node.name].size
        node.dist = node_dict[node.name].branch_length
        tree_style.scale = 1


def set_node_style(tree, node_dict):
    root_node = tree & 'root'
    for node in tree.traverse():
        if node is root_node:
            node.set_style(NodeStyle(size=0))
            continue
        node.set_style(node_dict[node.name].style)


def remove_node_circle(tree, node_dict):
    root_node = tree & 'root'
    for node in tree.traverse():
        if node is root_node:
            continue
        else:
            node_dict[node.name].style['size'] = 0
            node_dict[node.name].min_size = 0
            node_dict[node.name].size = 0


def add_pie_face(tree, ts, total_profile, group):
    for leaf_name in tree.iter_leaf_names():
        node = tree & leaf_name
        profile_list = total_profile[node.name]
        if group is not None:
            profile_list_grouped = OrderedDict()
            for sample_name, profile in profile_list.iteritems():
                profile_list_grouped[group[sample_name]] = profile
            profile_list = pd.Series(profile_list_grouped)
        col_num = len(profile_list)
        times = int(math.ceil(col_num / len(COLS_BREWER)))
        cols = (COLS_BREWER * times)[:col_num]
        summary = sum(profile_list)
        percents = map(lambda s: s / summary * 100, profile_list)
        P = PieChartFace(percents=percents, width=50, height=50, colors=cols)
        node.add_face(P, 0, 'aligned')
    for ind, g in enumerate(profile_list.index):
        ts.legend.add_face(TextFace(" "), 0)
        ts.legend.add_face(TextFace(" "), 1)
        T = TextFace('  %s  ' % g)
        ts.legend.add_face(T, 0)
        C = CircleFace(radius=10, color=cols[ind], style="circle")
        ts.legend.add_face(C, 1)
        ts.legend_position = 1


if __name__ == '__main__':
    params = read_params(sys.argv)
    mkdir(params['outdir'])
    modified_tax = '%s/tax_ass_modified.txt' % params['outdir']
    newick_file = '%s/tax_tree.nwk' % params['outdir']
    modify_tax_ass(params['tax_ass'], modified_tax)
    level_profile, total_profile = read_profile(modified_tax, params['profile'], params['top'])
    tree = read_tax(modified_tax, total_profile.T.mean())
    tree.adjust_profile()
    with open(newick_file, 'w') as out:
        out.write(str(tree))
    t = Tree(newick_file, format=1)
    ts = get_tree_style()

    set_node_default(t, node_dict=tree.nodes)
    suf = ''
    if params['with_branch_circle']:
        add_node_circle(t, node_dict=tree.nodes)
        suf += '_circle'
    if params['with_leaf_pie'] and not params['with_branch_circle']:
        remove_node_circle(t, node_dict=tree.nodes)
    if params['with_leaf_pie']:
        add_pie_face(t, ts, total_profile.T, group=params['group'])
        suf += '_pie'
    add_branch_text(t, tree_style=ts, node_dict=tree.nodes)
    set_node_style(t, node_dict=tree.nodes)

    pdf_file = '%s/tax_tree%s.pdf' % (params['outdir'], suf)
    t.render(pdf_file, tree_style=ts, dpi=100)
