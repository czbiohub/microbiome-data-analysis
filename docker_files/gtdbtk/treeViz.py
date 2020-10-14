#!/usr/bin/env python3

import sys

# import ete3
from ete3 import PhyloTree, faces, AttrFace, TreeStyle, NodeStyle
import pandas as pd
import numpy as np
import colorsys
import logging


def _get_colors(num_colors):
    """Generate N distinct colors.
    Source: https://stackoverflow.com/a/9701141/7609803

    Args:
        num_colors (int): number of colors

    Returns:
        list: list of 'num_color' color hex values
    """
    colors = []
    for i in np.arange(0, 360, 360 / num_colors):
        hue = i / 360
        lightness = (50 + np.random.rand() * 10) / 100
        saturation = (90 + np.random.rand() * 10) / 100
        # rgb_colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
        colors.append(
            "".join(
                [
                    "%0.2X" % int(c * 255)
                    for c in colorsys.hls_to_rgb(hue, lightness, saturation)
                ]
            )
        )
    return [f"#{hex}" for hex in colors]


def parse_summary_file(summary_file, out_prefix, taxa_rank=None):
    """Parse the gtdb summary file to extract list of genomes to prune the tree,
    genomes per taxa rank and generate colors per taxa rank for the tree.

    Args:
        summary_file (str): path to summary file from GTDBtk classify_wf
        out_prefix (str): prefix for the parsed summary dataframe.
        taxa_rank ([type], optional): [description]. Defaults to phylum.

    Returns:
        tuple: containing 3 elements
            list:   list of all genomes in this subset
            dict:   taxa rank --> background color hex value
            dict:   list of genomes per taxa rank. taxa rank --> [genome1, genome2, ... , genomeN]
    """

    taxa_ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    if taxa_rank is None:
        taxa_rank = "phylum"

    taxa_rank = taxa_rank.lower()
    assert (
        taxa_rank in taxa_ranks
    ), f"Taxonomy rank must be one of [{','.join(taxa_ranks)}]"

    # Parse the summary file
    summary_df = pd.read_table(
        summary_file, header=0, usecols=["user_genome", "classification"]
    )
    summary_df[taxa_ranks] = summary_df["classification"].apply(
        lambda x: pd.Series(str(x).split(";"))
    )
    summary_df.to_csv(f"{out_prefix}.taxa_levels.csv", index=False)

    # Get the column number for the taxa rank for use in the itertuple for loop later
    taxa_rank_col_num = (
        taxa_ranks.index(taxa_rank) + 3
    )  # add 2 columns that already exist in the df

    genomes = sorted(summary_df["user_genome"].unique())
    logging.info(f"Final tree will contain {len(genomes)} genomes.")

    t_rank = sorted(summary_df[taxa_rank].unique())
    num_taxa_ranks = len(t_rank)
    if num_taxa_ranks > 1:
        colors = _get_colors(num_taxa_ranks)
    elif num_taxa_ranks == 1:
        colors = ["#FFFFFF"]  # White
    else:
        logging.error(
            f"Could not find any taxonomic levels at {taxa_rank}. Did the summary file {summary_file} come from GTDBtk classify_wf output? Exiting ..."
        )
        sys.exit(1)

    logging.info(f"Found {num_taxa_ranks} levels at {taxa_rank} rank.")

    # rank --> color
    color_dict = dict(zip(t_rank, colors))

    # rank --> [genome1, genome2, ... , genomeN]
    common_ancestor = dict()
    for row in summary_df.itertuples():
        taxa = row[taxa_rank_col_num]
        if taxa in common_ancestor:
            common_ancestor[taxa].append(row.user_genome)
        else:
            common_ancestor[taxa] = [row.user_genome]
    return genomes, color_dict, common_ancestor


def generate_tree(tree_file, genomes, color_dict, common_ancestor, out_prefix):
    """Read the tree from GTDBtk
      - prune it to the list of genomes provided.
      - color the nodes based on the taxa rank chosen
      - write a tree file in newick format and a draw a linear tree in pdf.

    Args:
        tree_file (str): large tree file from GTDBtk classify_wf.
        genomes (list): list of genomes to prune the large tree file.
        color_dict (dict): taxa rank --> background color hex value
        common_ancestor (dict): list of genomes per taxa rank. taxa rank --> [genome1, genome2, ... , genomeN]
        out_prefix (str): prexif for pruned tree file and pdf
    """
    output_tree = f"{out_prefix}.pruned.tree"
    output_tree_fig = f"{out_prefix}.pruned.pdf"
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#node-backgrounds
    def layout(node):
        if node.is_leaf():
            N = AttrFace("name", fsize=20)
            faces.add_face_to_node(N, node, 0, position="aligned")

    tree = PhyloTree(tree_file, format=1, quoted_node_names=True)

    # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#node-backgrounds
    for phyla, genome_list in common_ancestor.items():
        node_style = NodeStyle()
        node_style["bgcolor"] = color_dict[phyla]
        tree.get_common_ancestor(genome_list).set_style(node_style)

    logging.info(f"Pruning tree to {len(genomes)} organisms ...")
    tree.prune(genomes, preserve_branch_length=True)

    logging.info(f"Saving pruned tree in newick format ...")
    tree.write(format=1, outfile=output_tree)

    logging.info(f"Drawing pruned tree ...")
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.show_scale = True

    # For circular tree
    # ts.mode = "c"
    # ts.root_opening_factor = 1

    # tree.show(tree_style=ts) # for interactive tree manipulation
    logging.info(f"Saving pruned tree image ...")
    tree.render(output_tree_fig, w=183, units="mm", tree_style=ts)
    return


def main(tree_file, summary_file, out_prefix):
    genomes, color_dict, common_ancestor = parse_summary_file(summary_file, out_prefix)
    generate_tree(tree_file, genomes, color_dict, common_ancestor, out_prefix)
    logging.info(f"All done! Huzzah!")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    # tree_file = (
    #     "/Users/sunit.jain/Research/Min/figures/db_Min_v1_1/gtdbtk.bac120.classify.tree"
    # )
    # summary_file = (
    #     "/Users/sunit.jain/Research/Min/figures/db_Min_v1_1/gtdbtk.bac120.summary.tsv"
    # )

    tree_file = sys.argv[1]
    summary_file = sys.argv[2]
    prefix = sys.argv[3]

    main(tree_file, summary_file, prefix)