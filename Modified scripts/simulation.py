import random
import os
import dendropy
import csv
import pandas as pd

"""
Script for phylogenetic simulations of pangenomes with selection
"""

def write_seqs(seqs, fname):
    """
    Writes binary gene presence/absence sequences to a file.

    Parameters:
        seqs (dict or list of dict): A dictionary or a list of dictionaries where each key is a gene ID
            and each value is a list of integers (0 or 1) indicating gene presence/absence.
        fname (str): The output filename where the formatted sequences will be written.
    """

    outf = open(fname, 'w')
    if isinstance(seqs, dict):
        seqs = [seqs]
    genomes = [k for k in seqs[0]]
    for s in seqs[0]:
        lens = []
        for aseqs in seqs:
            l = "".join([str(i) for i in aseqs[s]])
            outf.write(l)
            lens.append(aseqs[s].count(1))
        outf.write(" {}".format(s))
        for l in lens:
            outf.write(" {}".format(l))
        outf.write("\n")


def write_pa(seqs, fname, remove_constant=False):
    """
    Writes gene presence/absence (PA) data to a CSV file.

    Parameters:
        seqs (dict or list of dict): A dictionary (for one partition) or a list of dictionaries (multiple partitions),
            where each key is a genome ID and each value is a list of binary values indicating gene presence (1)
            or absence (0) across that genome.
        fname (str): Path to the output CSV file.
        remove_constant (bool): If True, genes with constant presence/absence across all genomes are excluded
            (i.e., genes that are present in all or none of the genomes).

    Output Format:
        The output CSV has the header: Gene,<genome_1>,<genome_2>,...,<genome_n>
        Each row: gene_<index>,<0 or 1>,<0 or 1>,...,<0 or 1>
    """

    outf = open(fname, 'w')
    if isinstance(seqs, dict):
        # only one partition
        seqs = [seqs]
    genomes = [k for k in seqs[0] if k.startswith("Node_")]
    outf.write("Gene,{}\n".format(",".join(genomes)))
    off = 0
    for aseqs in seqs:
        genes = len(aseqs[genomes[0]])
        for i in range(genes):
            pa = [aseqs[genomes[g]][i] for g in range(len(genomes))]
            if remove_constant and (sum(pa) == 0 or sum(pa) == len(pa)):
                continue
            outf.write("gene_{},".format(off + i + 1))
            outf.write("{}\n".format(",".join([str(p) for p in pa])))
        off += genes


def write_params(outf, psize, gsize, mut_rate, scale_events, num_pairs):
    """
    Writes simulation parameter values to an output file in a tab-separated format.

    Parameters:
        outf (file-like object): An open writable file or file handle.
        psize (int): Number of gene in pangenome.
        gsize (int): Number of genes per genome.
        mut_rate (float): Mutation rate used in the simulation.
        scale_events (float): Scaling factor applied to evolutionary events (e.g., gains/losses).
        num_pairs (int): Number of genome pairs used in pairwise comparisons.

    Output Format:
        Each parameter is written on its own line in the format: <parameter_name>\t<value>.
    """

    outf.write("psize\t{}\ngsize\t{}\nmut_rate\t{}\nscale_events\t{}\nnum_pairs\t{}\n".format(psize, gsize, mut_rate,
                                                                                              scale_events, num_pairs))


def write_events(events, fname):
    """
    Writes lists of event data to a file, one value per line.

    Parameters:
        events (list of int or list of list of int): A list of event values, or a list of event lists.
        fname (str): Path to the output file where event data will be written.

    Notes:
        - Each event is written on a new line, they will be ordered by the partitions
    """

    outf = open(fname, 'w')
    if isinstance(events[0], int):
        events = [events]
    for event in events:
        outf.write("{}\n".format("\n".join([str(e) for e in event])))

def simulate(adir, num_taxa, psize, gsize, lamb, mut_rate=1, scale_events=0, num_pairs=0, selection=0, mu=0):
    """Main simulation function

    :param adir: directory to create and write
    :param num_taxa: number of taxa
    :param psize: total pangenome size, integer for one rate or [int,int] for 2 rates
    :param gsize: number of genes per genome, integer for one rate or [int,int] for 2 rates
    :param lamb: lambda rate for birth death tree, set to -1 for coalescent tree
    :param mut_rate: turnover rate, expected rate is mut_rate*gsize for a branch length of 1, float for one rate or [float,float] for 2 rates
    :param scale_events: determine the mutation rate automatically to expect this number of events per gene, integer for one rate or [int,int] for 2 rates
    :param num_pairs: number of coevolving pairs, set to 0 for neutral simulation, integer for one rate or [int,int] for 2 rates (pairs only coevolve within rate class)
    :param selection: rate for genes in a pair if their mate is different, should be > 1 (positive association) or < -1 (negative association), set to 1 for co-transfer, set to -1 for co-avoidance
    :param mu: mu rate for birth death tree
    """

    def generate_ancestor(psize, gsize):
        """
         Generates a binary presence/absence vector representing an ancestral genome.

         Parameters:
             psize (int): Total number of possible gene families (pangenome size).
             gsize (int): Number of genes present in the ancestral genome.

         Returns:
             list of int: A binary list of length `psize`, with `gsize` ones and the rest zeros.
         """

        anc = [0] * (psize - gsize) + [1] * gsize
        return anc

    def generate_descendent(anc, d, psize, gsize, mut_rate, rate, num_pairs, selection, pairs, real_events):
        """
        Simulates the evolution of a descendent genome from an ancestral genome under a mutation model with optional selection.

        Parameters:
            anc (list of int): Binary vector representing the ancestral genome (1 = gene present, 0 = absent).
            d (float): Branch length (i.e., time available for evolution).
            psize (int): Total number of possible genes (length of the presence/absence vector).
            gsize (int): Number of genes per genome (used in mutation rate calculations).
            mut_rate (float): mutation rate.
            rate (list of float): Baseline mutation rates for each gene state.
            num_pairs (int): Number of gene pairs under selection.
            selection (float): Selection coefficient.
            pairs (list of tuple): List of gene index pairs under selection.
            real_events (list of int): A counter list to record the number of mutations per gene.

        Returns:
            list of int: A binary vector representing the descendent genome.

        This function simulates binary gene presence/absence changes (0<>1) over a branch of length `d`,
        starting from the given ancestral state `anc`. Mutations occur at an exponential rate proportional
        to the number of genes per genome and the mutation rate. Selection can modify mutation rates for
        certain gene pairs. Selection modifies mutation probabilities by scaling relevant gene rates. The
        simulation stops when the next mutation time exceeds the remaining branch length.
        """

        abs_sel = abs(selection)
        while True:
            t = random.expovariate(gsize * mut_rate)
            if t > d:
                # print ("{} changes on a bl of {}, expected mut: {}".format(mut, orgd,orgd*gsize*mut_rate))
                return anc
            rates = [rate[i] for i in anc]
            if num_pairs and selection > 1:
                for i in range(num_pairs):
                    for j in [0, 1]:
                        # only increase the rate for the 0-genes where the mate is 1
                        if anc[pairs[i][j]] == 0 and anc[pairs[i][1 - j]] == 1:
                            rates[pairs[i][j]] *= abs_sel
            if num_pairs and selection < 1:
                for i in range(num_pairs):
                    # only increase the rate of 11 pairs, then only increase every one by half, so the total rate to get out of the pair is abs_sel
                    if anc[pairs[i][0]] == 1 and anc[pairs[i][1]] == 1:
                        rates[pairs[i][0]] *= abs_sel / 2
                        rates[pairs[i][1]] *= abs_sel / 2
            i = random.choices(range(psize), rates)[0]
            anc[i] = 1 - anc[i]
            real_events[i] += 1
            d = d - t

    def simulate_one_psize(tree, psize, gsize, mut_rate, scale_events, num_pairs, selection):

        """
        Simulates the evolution of gene presence/absence patterns over a phylogenetic tree for a single partition.

        Parameters:
            tree (DendroPy Tree): The phylogenetic tree along which to simulate evolution.
            psize (int): The total number of genes in the simulated pangenome.
            gsize (int): The number of genes per genome.
            mut_rate (float): Mutation rate.
            scale_events (float or None): If provided, rescales the mutation rate.
            num_pairs (int): Number of gene pairs to apply selection to.
            selection (float): Selection coefficient for gene pairs.

        Returns:
            tuple:
                - mut_rate (float): The possibly rescaled mutation rate.
                - pairs (list): List of gene index pairs under selection.
                - seqs (dict): Dictionary mapping node labels to gene presence/absence vectors.
                - real_events (list): Total number of mutation events that occurred per gene across the tree.
                - obs_events (list): Total number of observed state changes per gene (descendent and ancestral node differs).
        """
        print("simulate one", psize, gsize, mut_rate, scale_events, num_pairs, selection)

        psize_org = psize  # the original pangenome size including the linked genes
        rate = [gsize / psize, 1 - gsize / psize]  # rates for leaving 0,1
        abs_sel = abs(selection)

        pairs = []
        if num_pairs:
            if selection==0:
                raise ValueError("Pairs are set but no selection coefficient.")
            if abs_sel != 1:
                in_pairs = random.sample(range(psize), num_pairs * 2)
                pairs = [[in_pairs[i], in_pairs[i + 1]] for i in range(0, num_pairs * 2, 2)]
            else:  # in this case, just do a neutral simulation and add the linked genes later
                psize = psize - num_pairs
                pairs = random.sample(range(psize), num_pairs)
            print(pairs)

        if scale_events:
            l = tree.length()
            mut_rate = scale_events * psize_org / l / gsize
            print(l, mut_rate)

        i = 1
        seqs = {}
        real_events = [0] * psize  # count the number of events that actually happened
        obs_events = [0] * psize  # number of observed events (node and ancestor differ)
        for node in tree.preorder_node_iter():
            if not node.parent_node:
                label = "root"
                node.label = label
                seq = generate_ancestor(psize, gsize)
                seq = generate_descendent(seq, 100, psize, gsize, mut_rate, rate, num_pairs, selection, pairs,
                                          real_events)  # simulate equilibrium at root
                real_events = [0] * psize  # reset them after root simulation
            else:
                if not node.taxon:
                    label = "inner_{}".format(i)
                    node.label = label
                    i += 1
                else:
                    label = node.taxon.label
                dseq = seqs[node.parent_node.label]
                seq = generate_descendent(dseq[:], node.edge_length, psize, gsize, mut_rate, rate, num_pairs, selection,
                                          pairs, real_events)
                for g in range(psize):
                    if seq[g] != dseq[g]: obs_events[g] += 1
            seqs[label] = seq

        # Simulation finished
        # add the linked pairs
        if num_pairs and abs_sel == 1:
            for a in pairs:
                obs_events.append(obs_events[a])
                real_events.append(obs_events[a])
                for i in seqs:
                    if selection == 1:
                        seqs[i].append(seqs[i][a])
                    else:
                        seqs[i].append(1 - seqs[i][a])
            pairs = [[pairs[i], psize + i] for i in range(0, num_pairs)]

        return mut_rate, pairs, seqs, real_events, obs_events

    print(adir, num_taxa, psize, gsize, lamb, mut_rate, num_pairs, selection, mu)

    taxa = dendropy.TaxonNamespace(["Node_" + str(i + 1) for i in range(num_taxa)])
    if lamb == -1:
        tree = dendropy.simulate.treesim.pure_kingman_tree(taxon_namespace=taxa)
    else:
        tree = dendropy.simulate.treesim.birth_death_tree(birth_rate=lamb, death_rate=mu, num_extant_tips=num_taxa,
                                                          taxon_namespace=taxa)
    # remove 0 branch lengths
    for node in tree.leaf_node_iter():
        if node.edge.length == 0:
            node.edge.length = 1e-6

    os.makedirs(adir, exist_ok=True)
    os.mkdir("{}/simulation".format(adir))
    outf = open("{}/tree.nwk".format(adir), 'w')
    ts = tree.as_string(unquoted_underscores=True, suppress_internal_node_labels=True, schema="newick")
    outf.write(ts.replace("[&R] ", ""))
    outf = open("{}/simulation/tree_inner.nwk".format(adir), 'w')
    outf.write(tree.as_string(unquoted_underscores=True, schema="newick").replace("[&R] ", ""))

    outf = open("{}/simulation/params.txt".format(adir), 'w')
    outf.write(
        "lambda\t{}\nmu\t{}\nnum_taxa\t{}\n".format(lamb, mu, num_taxa))
    if num_pairs:
        outf.write("selection\t{}\n".format(selection))

    if isinstance(psize, int):
        # one partition
        mut_rate, pairs, seqs, real_events, obs_events = simulate_one_psize(tree, psize, gsize, mut_rate, scale_events,
                                                                            num_pairs, selection)
        write_params(outf, psize, gsize, mut_rate, scale_events, num_pairs)

    else:
        # multiple partitions
        seqs = []
        real_events = []
        obs_events = []
        pairs = []
        for p in range(len(psize)):
            outf.write("Partition {}\n".format(p + 1))
            tmut_rate = mut_rate
            if isinstance(mut_rate, list): tmut_rate = mut_rate[p]
            amut_rate, apairs, aseqs, areal_events, aobs_events = simulate_one_psize(tree, psize[p], gsize[p],
                                                                                     tmut_rate, scale_events[p],
                                                                                     num_pairs[p], selection)
            write_params(outf, psize[p], gsize[p], amut_rate, scale_events[p], num_pairs[p])

            seqs.append(aseqs)
            real_events.append(areal_events)
            obs_events.append(aobs_events)
            pairs.append(apairs)

    write_pa(seqs, "{}/gene_presence_absence_all.csv".format(adir))
    write_pa(seqs, "{}/gene_presence_absence.csv".format(adir), True)
    gpa_to_coinfinder(adir)
    write_seqs(seqs, "{}/simulation/pg.txt".format(adir))
    write_events(real_events, "{}/simulation/real_events.txt".format(adir))
    write_events(obs_events, "{}/simulation/obs_events.txt".format(adir))

    if num_pairs:
        outf = open("{}/simulation/pairs.txt".format(adir), 'w')
        if isinstance(psize, int):
            apsize = [psize]
            apairs = [pairs]
            aseqs = [seqs]
            areal_events = [real_events]
            aobs_events = [obs_events]
        else:
            apsize = psize
            apairs = pairs
            aseqs = seqs
            areal_events = real_events
            aobs_events = obs_events

        off = 0
        for i in range(0, len(apsize)):
            for p in apairs[i]:
                outf.write("gene_{}\tgene_{}\n".format(p[0] + 1 + off, p[1] + 1 + off))
            off += apsize[i]

        outf = open("{}/simulation/pg_sel.txt".format(adir), 'w')
        print(apairs)
        for s in aseqs[0]:
            l = []
            for i in range(len(aseqs)):
                for p in apairs[i]:
                    l.append(str(aseqs[i][s][p[0]]) + str(aseqs[i][s][p[1]]))
            l = " ".join(l)
            outf.write("{} {}\n".format(l, s))

        outf = open("{}/simulation/real_events_sel.txt".format(adir), 'w')
        for i in range(len(apairs)):
            for p in apairs[i]:
                outf.write("{}\t{}\n".format(areal_events[i][p[0]], areal_events[i][p[1]]))
        outf = open("{}/simulation/obs_events_sel.txt".format(adir), 'w')
        for i in range(len(apairs)):
            for p in apairs[i]:
                outf.write("{}\t{}\n".format(aobs_events[i][p[0]], aobs_events[i][p[1]]))


def gpa_to_coinfinder(adir):
    """
    Script to generate the alpha beta format for coinfinder

    :param adir: simulation directory
    :return: None
    """

    print("convert format {}".format(adir))
    outf = open("{}/alpha_beta.tab".format(adir), 'w')
    header = []
    for line in open("{}/gene_presence_absence.csv".format(adir)):
        spl = line.rstrip().split(",")
        if not header:
            header = spl
            continue
        for i in range(1, len(spl)):
            if int(spl[i]):
                outf.write("{}\t{}\n".format(spl[0], header[i]))
