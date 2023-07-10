import sys
from ete3 import Tree

phylogeneticTree = sys.argv[1] 
outGroup =sys.argv[2]   #'Cya_NS01_5_2B_1_CK_Cya_NS01_01838_1666461_1666877_1_CK_00001561_null'

t= Tree(phylogeneticTree, format=1)

if outGroup != '':
    #find full name node
    setancestor = ''
    for node in t.traverse():
        if outGroup in node.name:
            setancestor = node.name
    
    if setancestor == '':
        raise ValueError("Outgroup not found in tree")
    
    t.set_outgroup(setancestor)

else:
    ancestor = t.get_midpoint_outgroup()
    t.set_outgroup(ancestor)


file_extension = len(phylogeneticTree.split('.')[-1]) 
out_name = phylogeneticTree[: -file_extension]+ 'rooted.treefile'

t.write(format=1, outfile=out_name)
