__author__ = 'Polina'
import os.path
import copy
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import operator
import scipy.stats as stats
from collections import defaultdict

def get_GT_cost(TS, sources, sinks, unreported):
    infected = set()
    cost = 0.0
    interactions = []
    causality = []
    for i in TS:
        t, n1, n2, inf1, inf2, rep1, rep2, ext1, ext2 = i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8]
        if inf1 == 1:
            infected.add(n1)
            if n2 not in infected and inf2 == 1:
                penalty_n1 = 0.5*np.abs((sources[n1]-t).total_seconds()) if n1 in sinks else 0.5*np.abs((unreported[n1]-t).total_seconds())
                penalty_n2 = 0.5*np.abs((sources[n2]-t).total_seconds()) if n2 in sinks else 0.5*np.abs((unreported[n2]-t).total_seconds())
                cost += penalty_n1 + penalty_n2
                interactions.append((t, n1, n2))
                infected.add(n2)
                causality.append((n1, n2))

    return cost, interactions, causality

def get_out_cost(output_paths, sources, sinks, unreported):
    cost = 0.0
    interactions = []
    causality = []
    for p in output_paths:
        for step in p:
            t, n1, n2 = step[0], step[1], step[2]
            if n1 != n2 and step not in interactions:
                interactions.append(step)
                causality.append((n1, n2))
                penalty_n1 = 0.5*np.abs((sources[n1]-t).total_seconds()) if n1 in sinks else 0.5*np.abs((unreported[n1]-t).total_seconds())
                penalty_n2 = 0.5*np.abs((sources[n2]-t).total_seconds()) if n2 in sinks else 0.5*np.abs((unreported[n2]-t).total_seconds())
                cost += penalty_n1 + penalty_n2
    return cost, interactions, causality


def validatePaths(paths, moment_of_infection, reported_infected):
    sorted_gt = [i[0] for i in sorted(moment_of_infection.items(), key=operator.itemgetter(1))]
    infected_gt = set(sorted_gt)
    res, lengths = [], []
    infected = set()
    for p in paths:
        order = []
        for step in p:
            infected.add(step[1])
            infected.add(step[2])
            if not order or order[-1] != step[1]:
                order.append(step[1])
            if not order or order[-1] != step[2]:
                order.append(step[2])
        order_gt = [i for i in sorted_gt if i in order]
        order = [i for i in order if i in order_gt]
        tau = stats.kendalltau(order, order_gt)
        res.append(tau[0])
        lengths.append(len(order))

    print 'GT and output infection intersection', len(infected.intersection(infected_gt))
    print 'output as infected, but not in GT', len(infected.difference(infected_gt))
    print 'output as infected, but not in reported', len(infected.difference(reported_infected))
    print 'output as infected, in GT, but not in reported', len((infected.intersection(infected_gt)).difference(reported_infected))
    print len(infected_gt.difference(infected))
    print len(reported_infected.difference(infected)), len(infected.difference(reported_infected))
    print len(infected), len(infected_gt), len(reported_infected)
    return res, lengths

def get_upstream(TS, sinks):
    out = set()

    from_to = defaultdict(set)
    to_from = defaultdict(set)
    from_to_idx = defaultdict(set)
    to_from_idx = defaultdict(set)

    for ind in xrange(len(TS)):
        j = TS[ind]
        t, n1, n2, inf1, inf2, rep1, rep2, ext1, ext2 = j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8]

        from_to[n1].add(n2)
        from_to_idx[(n1, n2)].add(ind)

        for i in to_from[n1]:
            from_to[i].add(n2)
            to_from[n2].add(i)

            from_to_idx[(i,n2)].add(ind)
            to_from_idx[(n2,i)].update(to_from_idx[(n1,i)])

        to_from[n2].add(n1)
        to_from_idx[(n2, n1)].add(ind)

        if n1 in sinks and sinks[n1] >= t:
            for elem in to_from[n1]:
                out.update(to_from_idx[(n1, elem)])
        if n2 in sinks and sinks[n2] >= t:
            for elem in to_from[n2]:
                out.update(to_from_idx[(n2, elem)])

    return out

def get_snapshots(TS, m, num_nodes, threshold = 0.75, maxTicks = False, upstream_ind = set()):
    #coloring, local_coloring = {}, {}
    moment_of_infection = {}
    infected = set()
    recovered = set()
    upstream = set()
    seeds = set()
    snapshots = {}
    output_snapshots = {}
    first_report = 0

    if maxTicks == True:
        ticks = [0] + range(first_report, len(TS))
    else:
        if m < len(TS):
            ticks = [0] + range(first_report, len(TS), (len(TS)-first_report)/m) + [len(TS)-1]
        else:
            ticks = range(len(TS))

    for i in xrange(0, len(ticks)-1):
        for j_ind in xrange(ticks[i], ticks[i+1]+1):
            j = TS[j_ind]
        #for j in TS[ticks[i]:ticks[i+1]+1]:
            t, n1, n2, inf1, inf2, rep1, rep2, ext1, ext2 = j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8]
            if j_ind in upstream_ind:
                upstream.add(n1)
                upstream.add(n2)


            if inf1 == 1:
                if n1 not in infected:
                    moment_of_infection[inf1] = t
                infected.add(n1)
                #local_coloring[j[1]] = local_coloring.get(j[1], 0) + 1

            if inf2 == 1:
                if n2 not in infected:
                    moment_of_infection[n2] = t
                    infected.add(n2)

            if ext1 == 1:
                seeds.add(n1)
        output_snapshots['infected'] = copy.deepcopy(infected)
        output_snapshots['upstream'] = copy.deepcopy(upstream)
        output_snapshots['recovered'] = copy.deepcopy(recovered)
        output_snapshots['seeds'] = copy.deepcopy(seeds)
        snapshots[TS[ticks[i+1]][0]] = copy.deepcopy(output_snapshots)
        #coloring[TS[ticks[i+1]][0]] = copy.deepcopy(local_coloring)

    return snapshots, moment_of_infection, seeds, ticks[1:]