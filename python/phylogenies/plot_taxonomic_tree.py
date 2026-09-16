import sys
import os
from sys import exit
from argparse import ArgumentParser, ArgumentError
from collections import OrderedDict, Counter
#from string import strip
from ete3 import PhyloTree, faces, NCBITaxa, TreeStyle
#import scipy.cluster.hierarchy as sch
from numpy import mean

__DESCRIPTION__ = """
Plotting of NCBI taxonomy and species trees
and a phylogenetic profile associated to them.
"""

def remove_single_inter_nodes(tNCBI):
    for d in tNCBI.get_descendants():
        if len(d.children) == 1:
            d.delete(prevent_nondicotomic=False)
    return tNCBI

def get_modified_NCBI(tNCBI):
    tax2node = dict([node.taxid, node] for node in tNCBI.get_cached_content())
    euk = tax2node[2759]
    # get 4 major eukaryotic groups as structure
    available_taxids = [3, 4, 5, 8]
    template = euk.copy()
    for child in template.get_descendants(): child.detach()
    # Unikonts
    opisthokonta = euk.search_nodes(sci_name='Opisthokonta')[0].detach()

    holozoa = template.copy()
    holomycota = template.copy()
    
    metazoa = opisthokonta.search_nodes(name='33208')[0].detach()
    metazoa.search_nodes(name='6072')[0].add_child(metazoa.search_nodes(name='6073')[0].detach())

    trichoplax = metazoa.search_nodes(name='10226')[0].detach()
    metazoa.add_child(trichoplax) # Trichoplax
    
    ctenophora = metazoa.search_nodes(name='10197')[0].detach()
    metazoa.add_child(ctenophora)

    porifera = metazoa.search_nodes(name='6040')[0].detach()
    metazoa.add_child(porifera) # Porifera
    

    choano = opisthokonta.search_nodes(name='28009')[0].detach()
    filasterea = opisthokonta.search_nodes(name='2687318')[0].detach()
    ichthyosporea = opisthokonta.search_nodes(name='127916')[0].detach()

    fungi = opisthokonta.search_nodes(name='4751')[0].detach()
    rotosphaerida = opisthokonta.search_nodes(name='2686024')[0].detach()

    a = template.copy()
    a.add_child(metazoa)
    a.add_child(choano)
    b = template.copy()
    b.add_child(a)
    b.add_child(filasterea)
    holozoa.add_child(b)
    holozoa.add_child(ichthyosporea)

    for ch in [fungi, rotosphaerida]:
        holomycota.add_child(ch)

    opisthokonta.add_child(holozoa)
    opisthokonta.add_child(holomycota)

    try:
        obazoa = template.copy()
        obazoa.add_child(opisthokonta)
        obazoa.add_child(euk.search_nodes(name='554296')[0].detach()) # Apusozoa
    except:
        pass

    try:
        amorphea = template.copy()
        amorphea.add_child(obazoa)
        amorphea.add_child(euk.search_nodes(name='554915')[0].detach()) # Amoebozoa
        euk.add_child(amorphea)
    except:
        pass
        
    try:
        # Archaeplastida
        viridiplantae = euk.search_nodes(sci_name='Viridiplantae')[0].detach()
        rhodophyta = euk.search_nodes(sci_name='Rhodophyta')[0].detach()

        archaeplastida = template.copy()
        archaeplastida.add_child(viridiplantae)
        archaeplastida.add_child(rhodophyta)
        euk.add_child(archaeplastida)
    except:
        pass    

    try:
        sar = euk.search_nodes(name='2698737')[0].detach()
        euk.add_child(sar)
    except:
        pass

    try:
        haptista = euk.search_nodes(name='2608109')[0].detach()
        euk.add_child(haptista)
    except:
        pass

    try:
        excavates = template.copy()
        discoba  = euk.search_nodes(name='2611352')[0].detach()
        metamonads = euk.search_nodes(name='2611341')[0].detach()
        excavates.add_child(discoba)
        excavates.add_child(metamonads)
        euk.add_child(excavates)
    except:
        pass

    try:
        cryptos = euk.search_nodes(name='3027')[0].detach()
        euk.add_child(cryptos)
    except:
        pass

    
    tNCBI = remove_single_inter_nodes(tNCBI)
    return tNCBI, tax2node

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

    parser.add_argument("--phylome_trees", dest="phylome_trees", default=".",
                        help="directory of gene trees for the duplication-counting mode")
    parser.add_argument("--tsv", dest="tsv",
                        action="store_true",
                        help="""parse tsv file""")    

    parser.add_argument("--newick", dest="newick",
                        type=str,
                        help="""
                        The newick of a species tree + code2taxid in infile
                        """)    

    parser.add_argument("--mode", dest="mode",
                        type=str,
                        help="""
                        "r" rectangular or "c" circular mode
                        """)

    parser.add_argument("--inter", dest="inter",
                        action="store_true",
                        help="""
                        If true, intermediate nodes are kept
                        """)

    parser.add_argument("--taxoncolors", dest="taxoncolors",
                        type=str,
                        help="""
                        path of color dictionary
                        """)

    parser.add_argument("--color_clades", dest="color_clades",
                        action="store_true",
                        help="""
                        If true, clades are coloured
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

    parser.add_argument("--no_strains", dest="no_strains",
                        action="store_true",
                        help="""
                        If true, terminal nodes are removed if internal present
                        """)    

    parser.add_argument("--no_names", dest="no_names",
                        action="store_true",
                        help="""
                        If true, leaf names are not plotted
                        """)    

    parser.add_argument("--ultrametric", dest="ultrametric",
                        type=int,
                        help="""
                        If true, tree is transformed to ultrametic                        
                        """)

    parser.add_argument("--scale", dest="scale",
                        type=int,
                        help="""
                        If true, tree is scaled
                        """)    

    parser.add_argument("--swap", dest="swap",
                        type=str,
                        #default="all",
                        help="""
                        If string=='', all branches are swaped. Otherwise node taxids shall be defined - for visual inspection
                        """)

    parser.add_argument("--swap_up", dest="swap_up",
                        type=str,
                        #default="all",
                        help="""
                        Ordered taxids to bring up
                        """)            

    parser.add_argument("--profile", dest="profile",
                        type=str,
                        help="""
                        Phylogenetic profile, tab delimited, taxid per line, column per COG, if T/F presence/absence, else number
                        """)

    parser.add_argument("--profile_shape", dest="profile_shape",
                        type=str,
                        help="""
                        'box' (b) or 'circle' (c) 
                        """)    

    parser.add_argument("--bubbles", dest="bubbles",
                        type=str,
                        help="""
                        Plots bubbles in nodes based on taxid 2 value
                        tab file.
                        """)    

    parser.add_argument("--mcl_clusters", dest="mcl_clusters",
                        type=str,
                        help="""
                        MCL output file, cluster per line
                        """)

    parser.add_argument("--neurotrans_clusters", dest="neurotrans_clusters",
                        type=str,
                        help="""
                        Neurotransmitters clades file (profile)
                        """)

    parser.add_argument("--barplots", dest="barplots",
                        action="store_true",
                        help="""
                        If true, profile plotted as barplots
                        """)    
    
    parser.add_argument("--per_neurotrans", dest="per_neurotrans",
                        type=str,
                        help="""
                        neurotrans clusters profile
                        """)    

    parser.add_argument("--map_dups", dest="map_dups",
                        type=str,
                        help="""
                        Map dups of gene tree(s). Path of tree per line
                        """)

    parser.add_argument("--map_dups_tree_collection", dest="map_dups_tree_collection",
                        action="store_true",
                        help="""
                        tree collection dup analysis
                        """)

    parser.add_argument("--dups_type", dest="dups_type",
                        type=str,
                        help="""
                        internal, external or both
                        """)    

    parser.add_argument("--mcl_cluster_colors", dest="mcl_cluster_colors",
                        type=str,
                        help="""
                        color per line
                        """)    

    parser.add_argument("--hmmsearch_profile", dest="hmmsearch_profile",
                        type=str,
                        help="""
                        profile from hmmsearch outfile
                        """)

    parser.add_argument("--counts", dest="counts",
                        type=str,
                        help="""
                        barplot of counts per species
                        """)        

    parser.add_argument("--seqid2taxid", dest="seqid2taxid",
                        type=str,
                        help="""
                        seqID to taxid dictionary
                        """)

    parser.add_argument("--colorbar", dest="colorbar",
                        action="store_true",
                        help="""
                        Colorbar for the heatmap with matplotlib
                        """)

    parser.add_argument("--colorbar_save", dest="colorbar_save",
                        type=str,
                        help="""
                        save path of Colorbar for the heatmap with matplotlib
                        """)    

    parser.add_argument("--coverage", dest="coverage",
                        type=str,
                        help="""
                        path of cogs file to find coverage of species
                        """)    

    parser.add_argument("--root", dest="root",
                        type=str,
                        help="""
                        Comma separated leaf names, their ancestor used as root
                        """)

    args = parser.parse_args()
    infile = args.infile
    mode = args.mode
    newick = args.newick
    mcl = args.mcl_clusters
    neurotrans = args.neurotrans_clusters
    domtbl = args.hmmsearch_profile
    
    if newick:
        t = PhyloTree(newick, sp_naming_function=lambda name: name.split(".")[0])
    elif infile:
        if args.tsv:
            taxid2data = dict([line.split("\t")[2], line.split("\t")[4]] for line in open(infile).readlines()[1:] if line.split("\t")[4])
            taxids = taxid2data.keys()
        else:
            taxids = map(str.strip, open(infile))
        taxids = set(taxids)
        ncbi = NCBITaxa()
        t = ncbi.get_topology(taxids, intermediate_nodes=True)
        
        
                
        
    #t = t.collapse_lineage_specific_expansions()

    if len(t.children) == 1:
        t = t.children[0]
    
    if args.root:
        t.set_outgroup(t.get_common_ancestor(args.root.split(",")))
    
    if args.taxoncolors:
        #taxon2color = dict([int(line.split()[0]), line.split()[1]] for line in open(args.taxoncolors))
        taxon2color = dict([int(line.split()[0]), line.split()[1].strip()] for line in open(args.taxoncolors) if not line.startswith("#"))

        for node in t.traverse():
            for l in node.lineage[::-1]:
                if "taxid" in node.features and node.taxid in taxon2color:
                    node.add_feature("bgcolor", taxon2color[node.taxid])
                    for n in node.get_descendants():
                        n.add_feature("bgcolor", taxon2color[node.taxid])

                        
    t, tax2node = get_modified_NCBI(t)
    
    t = remove_single_inter_nodes(t)
        
    commons = []
    for node in t.traverse():
        node.sci_name = str(node.sci_name)
    
    if args.no_strains:
        print("Tree size before prunning: %s" % len(t))
        for node in t.get_descendants():
            if not node.is_leaf():
                if node.name in taxids:
                    for d in node.get_children():
                        d.detach()
                        #if d.is_leaf():
                        #    node.up.add_child(d)
                # elif node.rank == 'genus':
                #     # get leaf with most non-zero profile postions
                #     most_complete = sorted(node.get_leaves(), key= lambda l: len([p for p in l.profile if p != 0]), reverse=True)[0]
                #     most_complete.detach()
                #     for d in node.get_children():
                #         d.detach()
                #     node.add_child(most_complete)
                elif len(node.children) == 1:
                    node.delete()
        print("Tree size after prunning: %s" % len(t))

    if domtbl:
        clustering = False
        container2contained = False
        #container2contained = {'PF00001.21' : ['PF13853.6']}
        #clcolors = ['black']
        clcolors = ['darkgray']
        cl2seqids = {}
        cl2sp = {}
        cl2max = {}
        for line in open(domtbl):
            if not line.startswith("#"):
                f = line.strip().split()
                seqid = f[0]
                taxid = int(seqid.split(".")[0])
                domid = f[3]
                cl2seqids.setdefault(domid, set())
                cl2seqids[domid].add(seqid)

        if container2contained:
            for container in container2contained:
                for contained in container2contained[container]:
                    cl2seqids[container] = cl2seqids[container] - cl2seqids[contained]

        for cl in cl2seqids:                
            cl2sp.setdefault(cl, Counter())
            cl2sp[cl].update([seqid.split(".")[0] for seqid in cl2seqids[cl]])
            
        leaves = t.get_leaves()
        doms = sorted(cl2sp.keys())
        
        for cl in doms:
            clmax = float(max(cl2sp[cl].values()))
            cl2max[cl] = clmax
            
        # doms = sorted(cl2sp, key=lambda l: len(cl2sp[l]), reverse=True)
        if clustering:
            # Clustering
            X = []

            for cl in doms:
                profile = []
                for leaf in leaves:
                    profile.append(cl2sp[cl][leaf.taxid])
                X.append(profile)

            print('clustering')

            # Y = sch.linkage(X, method='ward')
            Y = sch.linkage(X, method='single', metric='hamming')
            Z = sch.dendrogram(Y)
            index = Z['leaves']
            print('done')
        
            for leaf in leaves:
                leaf.add_feature("profile", [cl2sp[doms[i]][leaf.taxid] for i in index])

        else:
            # assign profiles to species (leaves)
            for leaf in leaves:
                leaf.add_feature("profile", [(cl2sp[cl][leaf.name], cl2sp[cl][leaf.name]/cl2max[cl]) for cl in doms])

        if not args.barplots:
            # headers
            clnames = []
            headers = []
            for cl in doms:
                clnames.append(cl)
                header = faces.TextFace("%s [size: %s]" % (cl, sum(cl2sp[cl].values())), fsize=12, bold=True, fgcolor="#aa0000")
                header.rotation = -45
                header.vt_align = 0
                headers.append(header)
        else:
            # headers
            clnames = []
            headers = []
            for cl in doms:
                clnames.append(cl)
                header = faces.TextFace("%s [size: %s, max: %s]" % (cl, sum(cl2sp[cl].values()), int(cl2max[cl])), fsize=30, bold=True, fgcolor="#aa0000")
                header.rotation = -45 
                header.hz_align = 1
                header.margin_left = 50
                headers.append(header)

    if args.counts:
        sp2count = {}
        for line in open(args.counts).readlines()[1:]:
            sp = line.split("\t")[7]
            count = line.split("\t")[11]
            #sp = line.split()[0]
            #count = line.split()[1].strip()
            if count:
                count = int(count)
                sp2count.setdefault(sp, [])
                sp2count[sp].append(count)

        for sp in sp2count.keys():
            sp2count[sp] = max(sp2count[sp])
        
        sp2data = dict([line.split("\t")[7], line.split("\t")[9]] for line in open(args.counts).readlines()[1:])
        sp2new = dict([line.split("\t")[7], line.split("\t")[2]] for line in open(args.counts).readlines()[1:])
        sp2sc = dict([line.split("\t")[7], line.split("\t")[3].strip()] for line in open(args.counts).readlines()[1:])
        
        #sp2data = dict([line.split()[0], line.split()[1].strip()] for line in open(args.counts).readlines()[1:])
        vmax = float(max([v for v in sp2count.values() if v]))
        
        leaves = t.get_leaves()
            
        # assign profiles to species (leaves)
        for leaf in leaves:
            if leaf.name in sp2count:
                leaf.add_feature("profile", [(sp2count[leaf.name], sp2count[leaf.name]/vmax)])
            if leaf.name in sp2new:
                leaf.add_feature("status", sp2new[leaf.name].split("-")[0].strip())
                #print(leaf.sci_name, leaf.status)
            if sp2sc[leaf.name] == "YES":
                leaf.add_feature("sc", True)

        for node in t.get_monophyletic(values=["NEW"], target_attr="status"):
            node.add_feature("new", True)
            #print(node.get_ascii(attributes=['sci_name']))
            #print(node.sci_name)
            for d in node.get_descendants():
                d.add_feature("new", True)
                
        headers = []
        header = faces.TextFace("# predicted peptides", fsize=30, bold=True, fgcolor="#aa0000")
        # header = faces.TextFace("# species-specific clusters", fsize=30, bold=True, fgcolor="#aa0000")
        #header = faces.TextFace("times among orthologs", fsize=30, bold=True, fgcolor="#aa0000")
        #header.rotation = -45 
        header.hz_align = 1
        #header.margin_left = 10
        headers.append(header)                
                
        
    if neurotrans:
        def cumulative(l):
            # Cumulative list
            newl = []
            for i, n in enumerate(l):
                newl.append(n+sum(l[:i]))
            return newl
        
        clcolors = ['darkred','darkgreen','darkgray','indianred', 'blue', 'orange']
        #cl2sp = dict([n+1, Counter([int(seqid.split(".")[0]) for seqid in line.split()])] for n, line in enumerate(open(mcl)))


        
        f2fam = OrderedDict()
        f2fam["biosynthesis_a"] = ['Carn_acetyltransf.0','Biopterin_H.0','Pyridoxal_deC.1','Pyridoxal_deC.2','DOMON.0','Glutaminase.0','PNMT.0']
        f2fam["vesicular_transport"] = ['MFS_SLC18.0','MFS_SLC17.0','VGAT_aa.0']
        f2fam["metabolism"] = ['AChE_tetra.0','HNMT.0']
        f2fam["reuptake"] = ['SSF.0','SNF.1','SNF.2','SDF.0']                               
        f2fam["metabotropic_1"] = ['7tm_1.0','7tm_1.1','7tm_1.2','7tm_1.3','7tm_1.4','7tm_1.5','7tm_1.6','7tm_1.7']
        f2fam["metabotropic_3"] = ['7tm_3.1','7tm_3.2']
        f2fam["ionotropic_1"] = ['Neur_chan.1','Neur_chan.2','Neur_chan.3','Neur_chan.4']
        f2fam["ionotropic_2"] = ['Lig_chan.1','Lig_chan.2']

        categories = ['clade','clade.1','clade.2','clade.3','clade.4','clade.5','clade.6','paralog','paralog.1','out','out.1']

        fam2set2name = {}
        fam2set = OrderedDict()
        for line in open(neurotrans):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.strip().split()
            family = f[0]
            subset = f[1]
            cat = f[2]
            name = f[3]
            tree = PhyloTree(f[4], sp_naming_function=lambda name: name.split(".")[0])
            tree.annotate_ncbi_taxa()
            species = set([int(leaf.name.split(".")[0]) for leaf in tree.get_leaves() if 2759 in leaf.lineage])

            fam = "%s.%s" % (family, subset)
            fam2set.setdefault(fam, {})
            fam2set[fam][cat] = species

            fam2set2name.setdefault(fam, {})
            fam2set2name[fam][cat] = name
            
        gaps = []
        headers = []
        for func in f2fam:
            n = 0
            for fam in f2fam[func]:
                for cat in sorted(fam2set[fam], key= lambda l: categories.index(l)):
                    if cat.startswith("clade"):
                        name = fam2set2name[fam][cat]
                        header = faces.TextFace("%s" % name, fsize=12, bold=True, fgcolor="#aa0000")
                        header.rotation = -90
                        header.vt_align = 1
                        #header.hz_align = 1
                        headers.append(header)
                        n += 1

                        for leaf in t.get_leaves():
                            if 'profile' not in leaf.features:
                                leaf.add_feature("profile", [])
                            
                            if cat == 'clade':
                                if leaf.taxid in fam2set[fam]['clade']:
                                    leaf.profile.append(1)
                                elif 'out' in fam2set[fam] and leaf.taxid in fam2set[fam]['out']:
                                    leaf.profile.append(2)                                
                                elif ( 'paralog' in fam2set[fam] and leaf.taxid in fam2set[fam]['paralog'] ) or\
                                     ( 'paralog.1' in fam2set[fam] and leaf.taxid in fam2set[fam]['paralog.1'] ):
                                    leaf.profile.append(3)
                                else:
                                    leaf.profile.append(0)
                            elif cat == 'clade.1':
                                if leaf.taxid in fam2set[fam]['clade.1']:
                                    leaf.profile.append(1)
                                elif ( 'out' in fam2set[fam] and leaf.taxid in fam2set[fam]['out'] ) or\
                                     ( 'out.1' in fam2set[fam] and leaf.taxid in fam2set[fam]['out.1'] ):
                                    leaf.profile.append(2)                                    
                                elif ( 'paralog' in fam2set[fam] and leaf.taxid in fam2set[fam]['paralog'] ) or\
                                     ( 'paralog.1' in fam2set[fam] and leaf.taxid in fam2set[fam]['paralog.1'] ):
                                    leaf.profile.append(3)
                                else:
                                    leaf.profile.append(0)
                            else:
                                if leaf.taxid in fam2set[fam][cat]:
                                    leaf.profile.append(1)
                                elif 'out' in fam2set[fam] and leaf.taxid in fam2set[fam]['out']:
                                    leaf.profile.append(2)
                                elif ( 'paralog' in fam2set[fam] and leaf.taxid in fam2set[fam]['paralog'] ) or\
                                     ( 'paralog.1' in fam2set[fam] and leaf.taxid in fam2set[fam]['paralog.1'] ):
                                    leaf.profile.append(3)
                                else:
                                    leaf.profile.append(0)

            
            gaps.append(n)
        print(gaps)
        gaps = cumulative(gaps)

    if args.mcl_clusters:
        clustering = False
        leaves = t.get_leaves()
        cl2sp_ordered = OrderedDict()
        cl2sp = {}
        for n, line in enumerate(open(args.mcl_clusters).readlines()):
            seqids = line.strip().split("\t")
            cl2sp[n+1] = Counter([seqid.split(".")[0] for seqid in seqids])

        if clustering:
            # Clustering
            print('clustering')
            X = []

            cl2max = {}
            for cl in cl2sp.keys():                
                clmax = float(max(cl2sp[cl].values()))
                cl2max[cl] = clmax
                profile = []
                for leaf in leaves:
                    profile.append(cl2sp[cl][leaf.name]/clmax)
                X.append(profile)

            #Y = sch.linkage(X, method='ward', optimal_ordering=True)
            Y = sch.linkage(X, method='single', optimal_ordering=True)
            #Y = sch.linkage(X, method='single', metric='hamming', optimal_ordering=True)
            Z = sch.dendrogram(Y)
            index = Z['leaves']
            for i in index:
                cl = cl2sp.keys()[i]
                cl2sp_ordered[cl] = cl2sp[cl]
            print('done')

        else:
            cl2max = {}
            for cl in sorted(cl2sp, key= lambda l: sum(cl2sp[l].values()), reverse=True):
                cl2max[cl] = float(max(cl2sp[cl].values()))
                cl2sp_ordered[cl] = cl2sp[cl]

        #(cl, cl2sp[cl]) for cat in catorder for cl in cat2cl[cat])
        cl2sp = cl2sp_ordered
        cls = cl2sp.keys()
        #cls = [cl for cl in cls if set(['27923','34499']) & set(cl2sp[cl].keys())] # Filter for ctenos
        cls = cls#[:10]
        

        if not args.barplots:
            # headers
            clnames = []
            headers = []
            for cl in cls:
                clnames.append(cl)
                header = faces.TextFace("cl%02d [size: %s]" % (int(cl), sum(cl2sp[cl].values())), fsize=12, bold=True, fgcolor="#aa0000")
                header.rotation = -90
                header.vt_align = 0
                headers.append(header)
        else:
            # headers
            clnames = []
            headers = []
            for cl in cls:
                clnames.append(cl)
                header = faces.TextFace("cl%02d [size: %s, max: %s]" % (int(cl), sum(cl2sp[cl].values()), int(cl2max[cl])), fsize=30, bold=True, fgcolor="#aa0000")
                #header.rotation = -90
                header.hz_align = 1
                header.margin_left = 50
                headers.append(header)


        # assign profiles to species (leaves)
        for leaf in leaves:
            leaf.add_feature("profile", [(cl2sp[cl][leaf.name], cl2sp[cl][leaf.name]/cl2max[cl]) for cl in cls])


    if args.map_dups_tree_collection:
        from goatools.obo_parser import GODag
        from glob import glob
        for node in t.traverse():
            node.add_feature('duprate', 0)
        godag = GODag("Documents/CGenomics/Projects/Phylome/01.Data/go-basic.obo", load_obsolete=True)
        slimdag = GODag("Documents/CGenomics/Projects/Phylome/01.Data/goslim_generic.obo")
        godag_keys = set(godag.keys())
        seqid2gos = {}
        for line in open("Documents/CGenomics/Projects/Phylome/02.Analysis/Phylome_Arun/singlecell_data_jake/data2jake/singlecell_emapper_v2.1.0.emapper.annotations"):
            if not line.startswith("#"):
                f = line.strip().split("\t")
                seqid = f[0]
                gos = f[12].split(",")
                if not gos == '-':
                    seqid2gos[seqid] = gos
                else:
                    seqid2gos[seqid] = []
                    
        # seed_go = 'GO:0045202' # synapse (C)
        # seed_go = 'GO:0048699' # generation of neurons
        seed_go = 'GO:0050877' # nervous system process
        # seed_go = 'GO:0007268' # chemical synaptic transmission
        # seed_go = 'GO:0035264' # multicellular organism growth
        # seed_go = 'GO:0097722' # sperm motility
        # seed_go = 'GO:0007613' # memory
        # seed_go = 'GO:0007049' # cell cycle
        # seed_go = 'GO:0006936' # muscle contraction

        seqids_with_go = set([seqid for seqid in seqid2gos if seed_go in seqid2gos[seqid]])
        seed_seqid2gos = {}
        slim_mode = False
        for seqid in seqids_with_go:
            gos = seqid2gos[seqid]
            if not gos:
                seed_seqid2gos[seqid] = seqid2gos[seqid]
                continue
            if slim_mode:
                slims = []
                for go in gos:
                    for slim in mapslim.mapslim(go, godag, slimdag)[0]:
                        slims.append(slim)
                slims = list(set(slims))
                seed_seqid2gos[seqid] = slims
            else:
                seed_seqid2gos[seqid] = gos

        tree_count = 0
        seed_go = False
        taxa2dupcounts = Counter()

        for n, treefile in enumerate(glob(os.path.join(args.phylome_trees, "*"))):
            if (n+1) % 100:
                #print("%s trees seen" % (n+1), end="\r")
                print("%s trees seen" % (n+1))
            treename = treefile.split("/")[-1]
            tree = PhyloTree(treefile, sp_naming_function= lambda seqid: seqid.split(".")[0])
            # Continue only with seed go containing trees
            if seed_go:
                common = seqids_with_go & set(tree.get_leaf_names())
            if seed_go and common:
                tree.set_outgroup(tree.get_midpoint_outgroup())
                root = tree.get_common_ancestor(common)
                if root != tree:
                    tree.set_outgroup(root)
                    root.get_sisters()[0].detach()
                    if len(tree.children) == 1:
                        tree.children[0].delete()
                        
                root = root.collapse_lineage_specific_expansions()
                root.resolve_polytomy() # for fasttree polytomies
                if len(tree) < 3:
                    continue
                tree_count += 1
                try:
                    root.get_descendant_evol_events()
                except:
                    continue
                for node in root.traverse():
                    if 'evoltype' in node.features and node.evoltype == 'D':
                        t.get_common_ancestor(node.get_species()).duprate += 1

            else:                
                tree.set_outgroup(tree.get_midpoint_outgroup())
                tree = tree.collapse_lineage_specific_expansions()
                tree.resolve_polytomy() # for fasttree polytomies
                if len(tree) < 3:
                    continue
                tree_count += 1
                tree.get_descendant_evol_events()
                for node in tree.traverse():
                    if 'evoltype' in node.features and node.evoltype == 'D':
                        t.get_common_ancestor(node.get_species()).duprate += 1
                        #taxa2dupcounts.update([",".join(sorted(node.get_species()))])
                
        print("tree storage and processing done")

        if seed_go:
            print("INFO: %s trees contain %s sequences annotated with '%s' (%s) and %s duplications" % (tree_count, len(seqids_with_go), godag[seed_go].name, seed_go, [node.duprate for node in t.traverse()]))
        else:
            print("INFO: %s trees contain %s duplications" % (tree_count, sum([node.duprate for node in t.traverse()]) ) )

        #for taxa in taxa2dupcounts:
        #    dupnode = t.get_common_ancestor(taxa.split(","))
        #    if 'duprate' not in dupnode.features:
        #        dupnode.add_feature('duprate', 0)
        #    dupnode.duprate += 1

        max_rate = max([node.duprate for node in t.traverse() if 'duprate' in node.features])
        
    if args.map_dups:
        tree = PhyloTree(args.map_dups, sp_naming_function=lambda name: name.split(".")[0])
        #outgroup = tree.get_midpoint_outgroup()
        #tree.set_outgroup(outgroup)
        try:
            events = tree.get_descendant_evol_events()
        except TypeError:
            print("polytomy")
            tree.resolve_polytomy()
            events = tree.get_descendant_evol_events()

        tree.annotate_ncbi_taxa()
        species_count = Counter([name.split(".")[0] for name in tree.get_leaf_names()])
        dups = [node for node in tree.traverse() if not node.is_leaf() and node.evoltype == 'D']
        #node2dups = Counter([node.taxid for node in dups])
        
        node2dups = Counter()
        occurences = Counter()
        for dup in dups:
            species = dup.get_species()
            if len(species) == 1:
                if args.dups_type in ['external','both']:
                    occurences.update(species)
                    common = t.search_nodes(name=list(species)[0])[0]
                else:
                    continue
            else:
                if args.dups_type in ['internal','both']:
                    occurences.update(species)
                    try:
                        common = t.get_common_ancestor(species & taxids)
                    except:
                        common = t
                else:
                    continue    
                
            node2dups.update([common])
            
        
            
    if args.swap:
        if args.swap == "all":
            print('swap all')
            t.swap_children()
            for node in t.get_descendants():
                if not node.is_leaf():
                    node.swap_children()
        elif args.swap == "ladder":
            print('ladderize')
            t.ladderize(direction=0)
        else:
            taxids2swap = args.swap.split()
            for taxid in taxids2swap:
                t.search_nodes(name=taxid)[0].swap_children()

    if args.swap_up:
        taxids2swapup = map(int, args.swap_up.split(","))
        for taxid in taxids2swapup:
            n = t.search_nodes(taxid=taxid)[0]
            up = n.up
            void = n.detach()
            up.add_child(n)

                
    if args.coverage:
        species2cov = {}
        n = 0
        for line in open(args.coverage):
            n += 1
            seqids = line.strip().split()
            for seqid in seqids:
                taxid = seqid.split(".")[0]
                species2cov.setdefault(taxid, 0)
                species2cov[taxid] += 1

    if args.ultrametric:
        print("ultrametric tree")
        t.convert_to_ultrametric(tree_length=args.ultrametric)

    if args.bubbles:
        taxid2value = dict([line.split("\t")[0], float(line.split("\t")[1])] for line in open(args.bubbles))
        max_bubble = max(taxid2value.values())
        
    def layout(node):
        node.img_style['size'] = 0
        node.img_style['vt_line_width'] = 2
        node.img_style['hz_line_width'] = 2
        node.img_style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
        node.img_style["hz_line_type"] = 0

        if "new" in node.features:
            #print(node)
            node.img_style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
            node.img_style["hz_line_type"] = 0
            node.img_style['vt_line_width'] = 6
            node.img_style['hz_line_width'] = 6
        
        if args.color_clades and "bgcolor" in node.features:
            descendants = node.get_descendants()
            if descendants and len(set([d.bgcolor for d in descendants])) == 1:
                for d in descendants:
                    d.features.remove('bgcolor')
            
            #node.img_style['size'] = 10
            node.img_style['bgcolor'] = node.bgcolor        

        if args.map_dups_tree_collection:
            if 'duprate' in node.features:
                #print(node.name, node.sci_name, node.duprate, "\n")
                sphere = faces.CircleFace( (node.duprate / max_rate) * 100, 'indianred', 'sphere')
                sphere.opacity = 0.5
                faces.add_face_to_node(sphere, node, 0, position="branch-right")
            
        if args.map_dups:
            if node in node2dups:
                print(node.name, node.sci_name, node2dups[node])
                sphere = faces.CircleFace(node2dups[node] / 3., 'indianred', 'sphere')
                sphere.opacity = 0.5
                faces.add_face_to_node(sphere, node, 0, position="float")
        if node.is_leaf():
            if args.coverage:
                present = (species2cov[node.name] / float(n)) * 100
                absent = 100 - present
                pie = faces.PieChartFace(percents=[present, absent], width=10, height=10, colors=["black", "white"])
                faces.add_face_to_node(pie, node, 0, position="branch-right")
            if not args.no_names:
                #name = faces.TextFace(".".join(node.name.split(".")[1:]), fsize=12, fstyle='normal')
                #name.margin_left = 6
                #faces.add_face_to_node(name, node, 0, aligned=False)
                penwidth = 0
                fsize = 12
                bold = False
                if "status" in node.features:
                    if node.status in ["NEW","UPDATE"]:
                        #print(node)
                        penwidth = 10
                        fsize = 14
                        bold = True
                                
                species = faces.TextFace("%s" % node.sci_name, fsize=fsize, fstyle='italic', penwidth=penwidth, bold=bold)
                
                
                if not args.color_clades and "bgcolor" in node.features:
                    species.background.color = node.bgcolor

                species.margin_left = 10
                faces.add_face_to_node(species, node, 0, aligned=False)
            if "sc" in node.features:
                c = faces.CircleFace(6, 'black', 'sphere')
                c.margin_left = 20
                faces.add_face_to_node(c, node, 1, position="branch-right")
            if args.map_dups:
                if occurences[node.name]:                  
                    
                    # pass if count is 0
                
                    occ = faces.TextFace(str(occurences[node.name]), fsize=12)
                    occ.inner_background.color = 'indianred'
                    occ.opacity = 0.5 + (occurences[node.name] / float(max(occurences.values()))) / 2.
                    faces.add_face_to_node(occ, node, 1, aligned=True)
                    
                if species_count[node.name]:
                    counts = faces.TextFace(str(species_count[node.name]), fsize=12)
                    counts.inner_background.color = 'steelblue'
                    counts.opacity = 0.5 + (species_count[node.name] / float(max(species_count.values()))) / 2.
                    counts.margin_left = 10
                    faces.add_face_to_node(counts, node, 2, aligned=True)
            else:
                pass
                #species = faces.TextFace(" ".center(50), fsize=12)
                #species.margin_left = 3
                #if "bgcolor" in node.features:
                #    species.background.color = node.bgcolor
                #faces.add_face_to_node(species, node, 0, aligned=True)

            """
            if domtbl:
                for i, p in enumerate(node.profile):
                    #bar = faces.BarChartFace(values=[p], width=10, height=1000, max_value=clmax[i], colors=['gray'], label_fsize=0, scale_fsize=0)
                    #bar.inner_background.color = 'white'
                    #bar.rotation = 90
                    #if p!=0:
                    #    p = p / 2.
                    bar = faces.SeqMotifFace(seq=None, motifs=[[0,  p, "[]", None, 10, "black", "black", None]], gap_format="blank")
                    bar.margin_left = 50
                    faces.add_face_to_node(bar, node, i, aligned=True)
            """
            
            if neurotrans:
                if args.mode == 'c':
                    size = 16
                else:
                    size = 10
                for i, p in enumerate(node.profile):
                    if p != 0:
                        text = ' '.center(8)
                        if p == 1:
                            color = 'indianred'
                        elif p == 2:
                            color = 'gray'
                        elif p == 3:
                            color = 'gainsboro'
                    else:
                        text = ' '.center(8)
                        color = 'white'
                                                
                    if args.profile_shape and args.profile_shape == 'c':
                        cell = faces.TextFace(size, color)
                    else:
                        cell = faces.TextFace(text, fsize=size)
                        cell.inner_background.color = color
                        
                    if i in gaps:
                        if args.mode == 'c':
                            cell.margin_left = 500
                        else:
                            cell.margin_left = 50
                        
                    else:
                        if args.mode == 'c':
                            cell.margin_left = 30
                        else:
                            cell.margin_left = 1
                            
                    #cell.border.color = 'lightgray'
                    #cell.border.width = .5

                    faces.add_face_to_node(cell, node, i, aligned=True)

            if mcl or domtbl or args.counts:
                abslt = False
                if args.mode == 'c':
                    size = 16
                else:
                    size = 10
                if "profile" in node.features:
                    for i, p in enumerate(node.profile):
                        absolute = p[0]
                        relative = p[1]

                        if args.barplots:
                            if abslt:                        
                                val = absolute
                            else:
                                if args.counts:
                                    val = int(round(relative * 1000))
                                else:
                                    val = int(round(relative * 1000))

                            if args.counts and node.name in sp2data:
                                if sp2data[node.name].startswith("genome"):
                                    color = "#2980b9"
                                elif sp2data[node.name].startswith("transcriptome"):
                                    color = "indianred"
                                elif sp2data[node.name].startswith("REMOVE"):
                                    color = "#cacfd2"                                    
                                elif node.name in sp2data:
                                    color = "darkgray"
                                else:
                                    color = "black"
                            else:
                                #color = "black"
                                color = "#383838"
                                #color = "darkgray"
                            if val != 0:
                                bar = faces.SeqMotifFace(seq=None, motifs=[[0,  val, "[]", None, 10, color, color, None]], gap_format="blank")
                                bar.margin_left = 50
                                faces.add_face_to_node(bar, node, i+1, aligned=True)

                        else:
                            if absolute != 0:
                                text = str(absolute).center(5)
                                color = 'indianred'
                                cell = faces.TextFace(text, fsize=size)
                                cell.inner_background.color = color
                                #cell.background.color = color
                                cell.margin_left = 3
                                cell.opacity = 0.5 + relative / 2.
                            else:
                                text = " ".center(5)
                                color = 'white'
                                cell = faces.TextFace(text, fsize=size)
                                cell.inner_background.color = color
                                #cell.background.color = color
                                cell.margin_left = 3

                            faces.add_face_to_node(cell, node, i, aligned=True)    

                
        else:
            pass
            #support = faces.TextFace(node.support, fgcolor="red", fsize=8, fstyle='normal')
            #faces.add_face_to_node(support, node, 0, position="branch-bottom")

            #    name = faces.TextFace(node.sci_name, fsize=10, fstyle='italic')
            #    faces.add_face_to_node(name, node, 0, position='float')

    # Manual mods
    #if t.search_nodes(name='140493'):
    #    t.search_nodes(name='33208')[0].add_child(t.search_nodes(name='140493')[0].detach())
    # t = t.search_nodes(name='10197')[0] # Ctenophore
    
    S = TreeStyle()
    S.show_scale = False
    #S.allow_face_overlap = True
    S.show_leaf_name = False
    #S.scale = 100
    #S.draw_aligned_faces_as_table = True
    #S.aligned_table_style = 0
    #S.min_leaf_separation = 1
    if args.mode == 'r':
        S.mode = 'r'
    elif args.mode == 'c':
        S.mode = 'c'

    if args.scale:
        S.scale = args.scale
    #S.show_branch_support = True
    if neurotrans:
        for i, header in enumerate(headers):
            if i in gaps:
                if args.mode == 'c':
                    header.margin_left = 500
                else:
                    header.margin_left = 50
                        
            else:
                if args.mode == 'c':
                    header.margin_left = 30
                else:
                    header.margin_left = 1
            S.aligned_header.add_face(header, column=i)
        S.draw_aligned_faces_as_table = True

    if mcl or domtbl and args.barplots:
        for i, header in enumerate(headers):
            #pass
            S.aligned_header.add_face(header, i)
        S.draw_aligned_faces_as_table = True

    if args.counts:
        for i, header in enumerate(headers):
            S.aligned_header.add_face(header, i)
        S.draw_aligned_faces_as_table = True        
    
    if args.save:
        t.render(file_name=args.save, layout=layout, tree_style=S)
    else:
        print("showing")
        t.show(layout=layout, tree_style=S)


        
