#!/usr/bin/python3

import os
import sys
from ete3 import EvolTree #, Tree, TreeStyle, TextFace
from ete3.treeview.layouts import evol_clean_layout

os.mkdir("pamlwd")
os.mkdir("plots")

name = sys.argv[1]
tree_file_path = sys.argv[2]
alignment_file_path = sys.argv[3]
models=sys.argv[4]

#evol_models = ('M0',"M1","M2") #,'fb' runs forever??

evol_models = models.split(',')

def mylayout(node):
        if node.is_leaf():
                node.img_style["size"] = 8
                node.img_style["shape"] = "circle"


def tree_plot(tree, model):
        #ts = TreeStyle()

        # ts.scale = 100
        # ts.title.add_face(TextFace(name, fsize=20), column=0)
        # ts.layout_fn = mylayout

        modname = model.replace(".", "_")
        hist=model.split(".",1)

        image_name = modname+"_dnds.pdf"
        plot_filename = os.path.join("plots", image_name)

        #tree.render(plot_filename, w=18000, tree_style=ts, layout=evol_clean_layout, histfaces=[model])
        if hist[0]=="fb":
            tree.render(plot_filename, layout=evol_clean_layout)
        elif hist[0]=="M0":
            tree.render(plot_filename, layout=evol_clean_layout)
        else:
            # tree.render(plot_filename,tree_style=ts, layout=evol_clean_layout)
            tree.render(plot_filename,layout=evol_clean_layout, histfaces=[model])
        print("Plot EXECUTED WITH SUCCESS: "+ plot_filename)


def dnds_ete3(tree_file_path, alignment_file_path, mod):
        tree = EvolTree(tree_file_path, format=1)
        tree.link_to_alignment(alignment_file_path)
        print('alignment linked')

        tree.workdir =  'pamlwd'

        #run model
        tree.run_model(mod)
        print("model " +mod+ " done")

        #wdir=output_directory+'/'+mod+'/out'
        #tree.link_to_evol_model(wdir,evmodel)

        tree_plot(tree, mod)


for evol_model in evol_models:
    model_name = evol_model+'.'+name
    dnds_ete3(tree_file_path, alignment_file_path, model_name)