from __future__ import print_function
import sys
import os
from sys import exit
from argparse import ArgumentParser, ArgumentError
from collections import OrderedDict
from random import choice
import ete3
from ete3 import PhyloTree, faces, TreeStyle, NCBITaxa

#from ete3.treeview.faces import StaticItemFace, QGraphicsEllipseItem

__DESCRIPTION__ = """
Plotting gene trees
"""

#class TriangleFace(StaticItemFace):
    
# lay not used for plotting
def lay(node):
    node.img_style['size'] = 0
    if node.is_leaf():
        if node.name not in seqid2name:
            name = faces.TextFace("%s [%s]" % (node.name, node.sci_name))
        elif node.taxid:
            name = faces.TextFace("%s (%s) [%s]" % (seqid2name[node.name], node.name, node.sci_name))
        else:
            print("skipped")
            name = faces.TextFace("%s" % node.name)
        faces.add_face_to_node(name, node, 0, aligned=True)
        if node.name in seqid2descr:
            descr = faces.TextFace("%s" % seqid2descr[node.name])
            faces.add_face_to_node(descr, node, 1, aligned=True)
        if 'color' in node.features:
            node.img_style['bgcolor'] = node.color

        if node.taxid == 9606:
            node.img_style['bgcolor'] = 'red'
        elif node.taxid == 7227:
            node.img_style['bgcolor'] = 'blue'
        elif node.taxid == 6239:
            node.img_style['bgcolor'] = 'purple'
        if 'sequence' in node.features:
            seqFace = faces.SeqMotifFace(node.sequence, seq_format="seq")
            faces.add_face_to_node(seqFace, node, 2, aligned=True)
    if 'evoltype' in node.features and node.evoltype == 'D':
        node.img_style['size'] = 8
        node.img_style['fgcolor'] = 'red'

def layout(node):
    node.img_style['size'] = 0
    node.img_style['vt_line_width'] = 3
    node.img_style['hz_line_width'] = 3
    #if "color" in node.features:
    #    node.img_style['bgcolor'] = node.color
    if args.mark_clades:
        if "markcolor" in node.features:
            node.img_style['vt_line_width'] = 15
            node.img_style['hz_line_width'] = 15
            node.img_style['hz_line_color'] = node.markcolor
            node.img_style['vt_line_color'] = node.markcolor
        if "cladename" in node.features:
            name = faces.TextFace("%s\n%s" % (node.cladename, node.sci_name), fsize=12)
            faces.add_face_to_node(name, node, 1, position="branch-top")
            
    if node.is_leaf():
        if not args.no_names:
            n = 1
            #name = faces.TextFace(".".join(node.name.split(".")[1:]), fsize=12, fstyle='normal')
            #name.margin_left = 6
            #faces.add_face_to_node(name, node, 0, aligned=False)
            if (args.seqid2name or args.seqid2descr) and node.name in seqid2name:
                seqname = seqid2name[node.name]
                #seqname = node.name
            else:
                #seqname = node.name.split(".")[1]
                seqname = node.name
            if "aep_" in seqname:
                seqname = seqname.split("aep_")[0]
            if node.taxid:
                name = faces.TextFace("%s {%s}" % (seqname, node.sci_name), fsize=14, fstyle='italic')
                #name = faces.TextFace("%s" % seqname, fsize=14, fstyle='italic')
                
            else:
                name = faces.TextFace("%s" % seqname, fsize=14, fstyle='italic')
            name.margin_left = 2
            if "color" in node.features:
                name.background.color = node.color
            faces.add_face_to_node(name, node, n, aligned=False)
            # test triangle
            #print("triangle")
            #diamond = faces.QGraphicsTriangleItem(10,7,orientation=1)
            #faces.add_face_to_node(triangle, node, 1, aligned=False)
            n += 1

        if args.tag and "tag" in node.features:
            tag = faces.TextFace("(%s)" % node.tag, fsize=14, fstyle='bold')
            faces.add_face_to_node(tag, node, n, aligned=True)

            n += 1

            
        if args.aln:
            if "sequence" in node.features:
                # seqFace = faces.SeqMotifFace(node.sequence[578:579], seq_format="seq")
                if args.aln_cols:
                    indices = [x-1 for x in map(int, args.aln_cols.split(","))]
                    #print("plotting aln positions %s" % args.aln_cols)
                    seq = "".join([node.sequence[i] for i in indices])
                    node.add_feature("motif", seq)
                else:
                    seq = node.sequence
                seqFace = faces.SeqMotifFace(seq, seq_format="seq", scale_factor=4, height=22)
                faces.add_face_to_node(seqFace, node, n, aligned=True)
                n += 1

        if args.seqlength:
            bar = faces.SeqMotifFace(seq=None, motifs=[[0,  int(node.seqlength/2.), "[]", None, 10, "gray", "gray", None]], gap_format="blank")
            bar.margin_left = 10
            faces.add_face_to_node(bar, node, n, aligned=True)
            n += 1

        if args.seqid2descr and node.name in seqid2descr:
            descr = faces.TextFace("(%s)" % seqid2descr[node.name], fsize=14)
            if node.taxid in [9606,10090,7227,6239,6359]:
                descr.background.color = node.color
            faces.add_face_to_node(descr, node, n, aligned=True)            
    else:
        if args.dups and "evoltype" in node.features and node.evoltype == "D":
            sphere = faces.CircleFace(radius=12, color='red', style='sphere')
            faces.add_face_to_node(sphere, node, 0, position="branch-right")

        if not args.no_support:
            support = faces.TextFace(int(round(node.support)), fgcolor="red", fsize=16, fstyle='normal')
            faces.add_face_to_node(support, node, 0, position="branch-bottom")

if __name__ == "__main__":
    parser = ArgumentParser(description=__DESCRIPTION__)
        # name or flags - Either a name or a list of option strings, e.g. foo or -f, --foo.
        # action - The basic type of action to be taken when this argument is encountered at the command line. (store, store_const, store_true, store_false, append, append_const, version)
        # default - The value produced if the argument is absent from the command line.
        # type - The type to which the command-line argument should be converted.
        # choices - A container of the allowable values for the argument.
        # required - Whether or not the command-line option may be omitted (optionals only).
        # help - A brief description of what the argument does.
        # dest - The name of the attribute to be added to the object returned by parse_args().


    parser.add_argument("-i", "--infile", dest="infile",
                        type=str,
                        required=True,
                        help="""full path of taxids infile""")


    parser.add_argument("--mode", dest="mode",
                        type=str,
                        help="""
                        "r" rectangular or "c" circular mode
                        """)

    parser.add_argument("--taxoncolors", dest="taxoncolors",
                        type=str,
                        help="""
                        path of color dictionary
                        """)

    parser.add_argument("--save", dest="save",
                        type=str,
                        help="""
                        If provided, tree rendered to file
                        """)

    parser.add_argument("--no_internal_names", dest="no_internal_names",
                        action="store_true",
                        help="""
                        If true, internal names are not plotted
                        """)

    parser.add_argument("--no_intermediate_nodes", dest="no_intermediate_nodes",
                        action="store_true",
                        help="""
                        If true, internal nodes are removed
                        """)    

    parser.add_argument("--no_names", dest="no_names",
                        action="store_true",
                        help="""
                        If true, leaf names are not plotted
                        """)

    parser.add_argument("--no_support", dest="no_support",
                        action="store_true",
                        help="""
                        If true, support not shown
                        """)        

    parser.add_argument("--ultrametric", dest="ultrametric",
                        type=int,
                        help="""
                        If true, tree is transformed to ultrametic                        
                        """)

    parser.add_argument("--scale", dest="scale",
                        type=int,
                        help="""
                        If true, tree is scaled to int
                        """)
    
    parser.add_argument("--swap", dest="swap",
                        type=str,
                        #default="all",
                        help="""
                        If string=='', all branches are swaped. Otherwise node taxids shall be defined - for visual inspection
                        """)

    parser.add_argument("--aln", dest="aln",
                        type=str,
                        help="""
                        path of alignment
                        """)

    parser.add_argument("--aln_cols", dest="aln_cols",
                        type=str,
                        help="""
                        comma separated column indeces to keep
                        """)

    parser.add_argument("--hmmsearch_domains", dest="hmmsearch_domains",
                        type=str,
                        help="""
                        add hmmsearch outfile domain annotation
                        """)  
    
    parser.add_argument("--seqid2name", dest="seqid2name",
                        type=str,
                        help="""
                        seqID to name dictionary
                        """)

    parser.add_argument("--seqid2descr", dest="seqid2descr",
                        type=str,
                        help="""
                        seqID to description dictionary
                        """)    

    parser.add_argument("--root", dest="root",
                        type=str,
                        help="""
                        Comma separated leaf names, their ancestor used as root
                        """)

    parser.add_argument("--keep", dest="keep",
                        type=str,
                        help="""
                        Comma separated leaf names, their ancestor used as clade to keep
                        """)    

    parser.add_argument("--root_in_group", dest="root_in_group",
                        type=int,
                        help="""
                        Taxonomic id used as root
                        """)

    parser.add_argument("--dups", dest="dups",
                        action="store_true",
                        help="""
                        Show duplications, species overlap
                        """)    

    parser.add_argument("--itol", dest="itol",
                        type=str,
                        help="""
                        path of itol annotation files (path+prefix)
                        """)

    parser.add_argument("--write", dest="write",
                        type=str,
                        help="""
                        path of newly written newick file (prefix)
                        """)

    parser.add_argument("--write_leaf_names", dest="write_leaf_names",
                        type=str,
                        help="""
                        path to write leaf names
                        """)

    parser.add_argument("--mark_clades", dest="mark_clades",
                        type=str,
                        help="""
                        mark the clade/outortholog/paralog selection. Path,family
                        """)
    
    parser.add_argument("--tag", dest="tag",
                        type=str,
                        help="""
                        tag leaves with any string
                        """)

    parser.add_argument("--species_selection", dest="species_selection",
                        type=str,
                        help="""
                        comma-delimited species taxids 2 keep
                        """)

    parser.add_argument("--extract_names", dest="extract_names",
                        action="store_true",
                        help="""
                        print leaf names
                        """)

    parser.add_argument("--seqlength", dest="seqlength",
                        type=str,
                        help="""
                        length profile based on fasta
                        """)


    
    
    args = parser.parse_args()
    infile = args.infile
    mode = args.mode

    t = PhyloTree(args.infile, sp_naming_function=lambda name: name.split(".")[0])
    #t = PhyloTree(args.infile, sp_naming_function=lambda name: name.split(".")[0], format=1)
    
    # IQTREE mods
    #for node in t.get_descendants():
    #    if not node.is_leaf():
    #        try:
    #            node.support = float(node.name.split("/")[1])
    #            node.name = ""
    #        except IndexError:
    #            node.support = 0

    if args.root:
        root = args.root.split(",")
        if len(root) == 1:
            t.set_outgroup(t.search_nodes(name=root[0])[0])
        else:
            try:
                print("rooting first in midpoint")
                outgroup = t.get_midpoint_outgroup()
                t.set_outgroup(outgroup)
                print("rooting in common ancestor")
                root_ancestor = t.get_common_ancestor(root)
                t.set_outgroup(root_ancestor)
            except:
                print("root exception")
                print("outgroup not monophyletic, first mid")
                t.set_outgroup(t.get_midpoint_outgroup())
                t.set_outgroup(t.get_common_ancestor(root))
                
    elif args.root_in_group:
        t.annotate_ncbi_taxa()
        root_seqids = []
        for node in t.get_leaves():
            for l in node.lineage[::-1]:
                if l == args.root_in_group:
                    root_seqids.append(node.name)
        if len(root_seqids) == 1:
            outgroup = t.search_nodes(name=list(root_seqids)[0])[0]
        else:
            outgroup = t.get_common_ancestor(root_seqids)
        if not outgroup.is_root():
            t.set_outgroup(outgroup)
    else:
        print("midpoint")
        outgroup = t.get_midpoint_outgroup()
        t.set_outgroup(outgroup)


    if args.keep:
        keep = args.keep.split(",")
        print("finding common ancestor")
        keep_ancestor = t.get_common_ancestor(keep).detach()
        t = keep_ancestor
        
    t.annotate_ncbi_taxa()
    for node in t.traverse():
        node.sci_name = str(node.sci_name)

    #t.ladderize()
    #print("ladderize")
        
    if args.species_selection:
        sps = args.species_selection.split(",")
        for leaf in t.get_leaves():
            if leaf.species not in sps or leaf.species in sps and leaf.name.startswith("6087.XP_"):
                leaf.delete(preserve_branch_length=True)
        for node in t.get_descendants():
            if len(node.children) == 1:
                node.delete(preserve_branch_length=True)
        if len(t.children) == 1:
            t = t.children[0]
            t.up.delete()
                

            
    if args.tag:
        seqid2tag = {}
        for line in open(args.tag):
            f = line.strip().split()
            seqid2tag[f[0]] = f[1]

        for leaf in t.get_leaves():
            print(leaf.name)
            if leaf.name.startswith("6087."):
                leaf.name = leaf.name.split("_")[0]
            if leaf.name in seqid2tag:
                print("yes")
                leaf.add_feature("tag", seqid2tag[leaf.name])
    
    if args.seqlength:
        from ete3.parser.fasta import read_fasta
        F = read_fasta(args.seqlength)
        for leaf in t.get_leaves():
            leaf.add_feature("seqlength", len(F.get_seq(leaf.name)))
    
    fast = False
    if max([node.support for node in t.traverse()]) <= 1:
        fast = True
    if fast:
        for node in t.traverse():
            node.support = node.support * 100

    if args.aln:
        t.link_to_alignment(args.aln)        

    if args.seqid2name:
        pass
        #seqid2name = dict([line.split()[0], line.split()[1].split("|")[2]] for line in open(args.seqid2name))
        
    if args.seqid2descr:
        #pass
        seqid2descr = dict([line.split()[0], " ".join(line.split("OS=")[0].split()[2:])] for line in open(args.seqid2descr))
        seqid2name = dict([line.split()[0], line.split()[1].split("|")[2]] for line in open(args.seqid2descr))

    if args.extract_names:
        names = [seqid2name[seqid] if seqid in seqid2name else seqid for seqid in t.get_leaf_names() if not seqid.startswith("6087.XP")]
        #print(",".join(names))
        print("\n".join(names))        
        sys.exit("names printed and exit")
        
    if args.taxoncolors:
        #taxon2color = dict([int(line.split()[0]), line.split()[1]] for line in open(args.taxoncolors))
        taxon2color = dict([int(line.split()[0]), line.split()[1].strip()] for line in open(args.taxoncolors) if not line.startswith("#"))

        for node in t.get_descendants():
            for l in node.lineage[::-1]:
                if l in taxon2color:
                    node.add_feature("color", taxon2color[l])
                    break

    if args.mark_clades:
        profile, selection = args.mark_clades.split(",")
        fam2seqids = OrderedDict()
        for line in open(profile):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.strip().split()
            family = f[0]
            subset = f[1]
            fam = "%s.%s" % (family, subset)
            if family == selection:
                cat = f[2]
                if 'clade' in cat:
                    color = 'indianred'
                elif 'out' in cat:
                    color = 'gray'
                elif 'paralog' in cat:
                    color = 'gainsboro'
                name = f[3]
                tree = PhyloTree(f[4], sp_naming_function=lambda name: name.split(".")[0])
                seqids = tree.get_leaf_names()
                for leaf in t.get_leaves():
                    if leaf.name in seqids:
                        leaf.add_feature("mark", name)

                for monoclade in t.get_monophyletic([name],'mark'):
                    monoclade.add_feature("cladename", name)
                    for n in monoclade.traverse():
                        n.add_feature("markcolor", color)
            
    if args.dups:
        try:
            events = t.get_descendant_evol_events()
        except TypeError:
            t.resolve_polytomy()
            #for node in t.get_descendants():
            #    if len(node.children) == 1:
            #        node.delete()
            if len(t.children) == 1:
                t = t.children[0]
            events = t.get_descendant_evol_events()
    
    if args.swap:
        if args.swap == "all":
            t.swap_children()
            for node in t.get_descendants():
                if not node.is_leaf():
                    node.swap_children()
        elif args.swap == "base":
            t.swap_children()        
        else:
            swaps = args.swap.split("-")
            # for taxid in toswap:
            #    t.search_nodes(name=taxid)[0].swap_children()
            for swap in swaps:
                toswap = swap.split(",")
                if len(toswap) == 1:
                    for node in t.search_nodes(taxid=int(toswap[0])):
                        node.swap_children()
                else:
                    t.get_common_ancestor(toswap).swap_children()

    if args.ultrametric:
        print("ultrametric tree")
        t.convert_to_ultrametric(tree_length=args.ultrametric)

    if args.hmmsearch_domains:
        domains = []
        seqid2dom = {}

        for line in open(args.hmmsearch_domains):
            if line.startswith("#"): continue
            f = line.strip().split()
            start = int(f[19])
            end = int(f[20])

            if not f[4].startswith("PF"):
                seqid = f[3]
                dom_acc = f[1]
                dom_name = f[0]
            else:
                seqid = f[0]
                dom_acc = f[4]
                dom_name = f[3]

            seqid2dom.setdefault(seqid, [])
            seqid2dom[seqid].append([start, end, dom_name])

            domains.append(dom_name)

        domains = set(domains)
        dom2color = {}
        # Random color per domain
        for dom in domains:
            color = choice(list(ete3.SVG_COLORS))
            dom2color[dom] = color
            #dom2color[dom] = "darkgreen"
            #dom2color[dom] = "darkred"

        for leaf in t.get_leaves():
            # nameFace = faces.TextFace("%s{%s}" % (leaf.name, leaf.sci_name))
            # leaf.add_face(nameFace, 0, "branch-right")
            seqid = leaf.name
            if not seqid in seqid2dom:
                continue
            motifs = []
            for domain in seqid2dom[seqid]:
                start, end, dom_name = domain
                color = dom2color[dom_name]
                #motif = [start,  end, "()", None, 10, "white", "rgradient:%s" % color, "arial|2|black|%s" % dom_name]
                motif = [start,  end, "()", None, 16, "white", color, "arial|13|black|%s" % dom_name]
                motifs.append(motif)
            domFace = faces.SeqMotifFace(seq=None, motifs=motifs, gap_format="line")
            leaf.add_face(domFace, 1, "aligned")
            #leaf.add_face(domFace, 2, "branch-right")
            
    S = TreeStyle()
    #S.allow_face_overlap = True
    S.show_leaf_name = False
    S.guiding_lines_color = 'white'
    if args.scale:
        S.scale = args.scale
    #S.scale = 100
    #S.draw_aligned_faces_as_table = True
    #S.aligned_table_style = 0
    #S.min_leaf_separation = 1
    if args.mode == 'r':
        S.mode = 'r'
    elif args.mode == 'c':
        S.extra_branch_line_color = 'white'
        S.mode = 'c'

    #S.arc_span = 300
    #S.allow_face_overlap = True
    #S.optimal_scale_level = 'full'
    #S.rotation = 225
    
    #S.scale = None
    #S.show_branch_support = True

    if args.save:
        t.render(file_name=args.save, layout=layout, tree_style=S)#, dpi=120)
    elif args.itol:
        out = open("%s.newick.itol.txt" % args.itol, "w")
        print(t.write(), file=out)
        out.close()
        
        out = open("%s.colors.itol.txt" % args.itol, "w")
        
        print("TREE_COLORS\nSEPARATOR COMMA\nDATA", file=out)
        colors = set([leaf.color for leaf in t.get_leaves()])
        for color in colors:
            for clade in t.get_monophyletic(values=[color], target_attr="color"):
                if len(clade) == 1:
                    print("%s,branch,%s,normal,1" % (clade.get_leaves()[0].name, clade.get_leaves()[0].color), file=out)
                else:
                    print("%s|%s,clade,%s,normal,1" % (clade.children[0].get_leaves()[0].name, clade.children[1].get_leaves()[0].name, color), file=out)

        out.close()
        out = open("%s.labels.itol.txt" % args.itol, "w")
                    
        print("LABELS\nSEPARATOR COMMA\nDATA", file=out)
        for leaf in t.get_leaves():
            if (args.seqid2name and leaf.name in seqid2name) or (args.seqid2descr and leaf.name in seqid2descr):
                print("%s,%s (%s) {%s} [%s]" % (leaf.name, leaf.name, seqid2name[leaf.name], leaf.sci_name, seqid2descr[leaf.name]), file=out)
                #print("%s,%s" % (leaf.name, seqid2name[leaf.name]), file=out) # For mouse only dataset
            else:                
                print("%s,%s {%s}" % (leaf.name, leaf.name, leaf.sci_name), file=out)

        out.close()
    elif args.write_leaf_names:
        out = open(args.write_leaf_names, "w")
        print("\n".join(t.get_leaf_names()), file=out)
        out.close()
    else:
        print("showing")
        t.show(layout=layout, tree_style=S)


        
