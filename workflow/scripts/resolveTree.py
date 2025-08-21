from ete3 import Tree
from argparse import ArgumentParser
import numpy as np
import time

parser = ArgumentParser(description='Convert a tree with polytomies into a fully binary tree.')
parser.add_argument('in_tree', help='Input Tree')
parser.add_argument('out_tree', help='Output Tree')
args = parser.parse_args()

t = Tree(args.in_tree, format=1)
t.resolve_polytomy(recursive=True)
t.write(format=5, outfile=args.out_tree)