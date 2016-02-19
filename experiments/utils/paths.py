__author__ = 'Polina'

import sys
#import weighting as w
#import GMC
import numpy as np
from datetime import datetime, timedelta
import accuracy as acc
import copy

def get_preinfected_penalty(t, infected_t, parameter = 1.0):
    if infected_t > t:
        return parameter*(infected_t - t).total_seconds()
    else:
        return 0.0

def get_prerecovered_penalty(t, recovered_t, beta = 0.0):
    if recovered_t > t:
        return beta/(recovered_t - t).total_seconds()
    else:
        return 0.0

def shortestPath_noweight(TS, sources, sinks, immuned, unreported, parameter = 1.0): #src: node -> time, sink: node-> time, immuned: node -> time

    shortest_paths = {i: {} for i in sources.keys()}
    out_paths = {i: {} for i in sources.keys()}
    c = 0
    for interaction in TS:
        c += 1
        if c % 100 == 0:
            print c
        t, n1, n2 = interaction[0], interaction[1], interaction[2]
        for src_node, t_start in sources.iteritems():

            if (n1 not in immuned or t < immuned[n1]):
                prerecovered_penalty = 0.0

                if n1 == src_node:
                    preinfected_penalty = get_preinfected_penalty(t, sources[n1], parameter) if n1 in sources else 0.0
                    dn1 = 0.0
                    if src_node not in shortest_paths[src_node] or shortest_paths[src_node][n1][0] >= dn1:
                        shortest_paths[src_node][src_node] = (dn1, t, [(t, n1, n1)])

                if n1 in shortest_paths[src_node].keys():

                    prime = shortest_paths[src_node][n1]
                    dn2 = prime[0] + 1.0
                    pathn2 = prime[2] + [(t, n1, n2)]
                    an2 = t
                    if n2 not in immuned or t < immuned[n2]: # if node is not immune yet
                        if n2 not in shortest_paths[src_node] or shortest_paths[src_node][n2][0] >= dn2:
                            shortest_paths[src_node][n2] = (dn2, an2, pathn2)
                    if n2 in sinks and sinks[n2] >= t and (n2 not in out_paths[src_node] or out_paths[src_node][n2][0] >= shortest_paths[src_node][n2][0]):
                        out_paths[src_node][n2] = copy.deepcopy(shortest_paths[src_node].get(n2, (np.Inf, -1, [])))
                        #c = c
                    if n1 in sinks and sinks[n1] >= t and (n1 not in out_paths[src_node] or out_paths[src_node][n1][0] >= shortest_paths[src_node][n1][0]):
                        out_paths[src_node][n1] = copy.deepcopy(shortest_paths[src_node].get(n1, (np.Inf, -1, [])))
                        #c = c
    SP = {}
    for i in sources.keys():
        SP[i] = {}
        for j in sinks.keys():
            p = out_paths[i].get(j, (np.Inf, -1, []))
            if len(p[-1]) > 1:
                p = (p[0], p[1], p[2][1:])
            SP[i][j] = p    
    return SP

def shortestPath1(TS, sources, sinks, immuned, unreported, parameter = 1.0): #src: node -> time, sink: node-> time, immuned: node -> time

    shortest_paths = {i: {} for i in sources.keys()}
    out_paths = {i: {} for i in sources.keys()}
    c = 0
    for interaction in TS:
        c += 1
        if c % 100 == 0:
            print c
        t, n1, n2 = interaction[0], interaction[1], interaction[2]
        for src_node, t_start in sources.iteritems():

            #if (n1 not in immuned or t < immuned[n1]):
            if n1 == src_node:
                dn1 = 0.0
                if src_node not in shortest_paths[src_node] or shortest_paths[src_node][n1][0] >= dn1:
                    shortest_paths[src_node][src_node] = (dn1, t, [(t, n1, n1)])

            if n1 in shortest_paths[src_node].keys():

                prime = shortest_paths[src_node][n1]

                if n1 in sinks:                        
                    penalty_n1 = parameter*0.5*np.abs((sources[n1]-t).total_seconds())
                else:
                    penalty_n1 = parameter*0.5*np.abs((unreported[n1]-t).total_seconds())
                if n2 in sinks:
                    penalty_n2 = parameter*0.5*np.abs((sources[n2]-t).total_seconds())
                else:
                    penalty_n2 = parameter*0.5*np.abs((unreported[n2]-t).total_seconds())
                
                dn2 = prime[0] + penalty_n1 + penalty_n2
                pathn2 = prime[2] + [(t, n1, n2)]
                an2 = t
                #if n2 not in immuned or t < immuned[n2]: # if node is not immune yet
                if n2 not in shortest_paths[src_node] or shortest_paths[src_node][n2][0] >= dn2:
                    shortest_paths[src_node][n2] = (dn2, an2, pathn2)
                        
                if n2 in sinks and sinks[n2] >= t and (n2 not in out_paths[src_node] or out_paths[src_node][n2][0] >= shortest_paths[src_node][n2][0]):
                    out_paths[src_node][n2] = copy.deepcopy(shortest_paths[src_node].get(n2, (np.Inf, -1, [])))                        
                if n1 in sinks and sinks[n1] >= t and (n1 not in out_paths[src_node] or out_paths[src_node][n1][0] >= shortest_paths[src_node][n1][0]):
                    out_paths[src_node][n1] = copy.deepcopy(shortest_paths[src_node].get(n1, (np.Inf, -1, [])))
                    
    SP = {}
    for i in sources.keys():
        SP[i] = {}
        for j in sinks.keys():
            p = out_paths[i].get(j, (np.Inf, -1, []))
            if len(p[-1]) > 1:
                p = (p[0], p[1], p[2][1:])
            SP[i][j] = p
    #return {i: {j: out_paths[i].get(j, (np.Inf, -1, [])) for j in sinks.keys()} for i in sources.keys()}
    return SP

def add_shortestPath1(new_TS, SP, sources, sinks, immuned, unreported, parameter = 1.0): #src: node -> time, sink: node-> time, immuned: node -> time

    #shortest_paths = {i: {} for i in sources.keys()}
    shortest_paths = SP
    #out_paths = {i: {} for i in sources.keys()}
    out_paths = SP
    c = 0
    for interaction in new_TS:
        c += 1
        if c % 100 == 0:
            print c
        t, n1, n2 = interaction[0], interaction[1], interaction[2]
        for src_node, t_start in sources.iteritems():

            #if (n1 not in immuned or t < immuned[n1]):
            if n1 == src_node:
                dn1 = 0.0
                if src_node not in shortest_paths[src_node] or shortest_paths[src_node][n1][0] >= dn1:
                    shortest_paths[src_node][src_node] = (dn1, t, [(t, n1, n1)])

            if n1 in shortest_paths[src_node].keys():

                prime = shortest_paths[src_node][n1]

                if n1 in sinks:
                    penalty_n1 = parameter*0.5*np.abs((sources[n1]-t).total_seconds())
                else:
                    penalty_n1 = parameter*0.5*np.abs((unreported[n1]-t).total_seconds())
                if n2 in sinks:
                    penalty_n2 = parameter*0.5*np.abs((sources[n2]-t).total_seconds())
                else:
                    penalty_n2 = parameter*0.5*np.abs((unreported[n2]-t).total_seconds())

                dn2 = prime[0] + penalty_n1 + penalty_n2
                pathn2 = prime[2] + [(t, n1, n2)]
                an2 = t
                #if n2 not in immuned or t < immuned[n2]: # if node is not immune yet
                if n2 not in shortest_paths[src_node] or shortest_paths[src_node][n2][0] >= dn2:
                    shortest_paths[src_node][n2] = (dn2, an2, pathn2)

                if n2 in sinks and sinks[n2] >= t and (n2 not in out_paths[src_node] or out_paths[src_node][n2][0] >= shortest_paths[src_node][n2][0]):
                    out_paths[src_node][n2] = copy.deepcopy(shortest_paths[src_node].get(n2, (np.Inf, -1, [])))
                if n1 in sinks and sinks[n1] >= t and (n1 not in out_paths[src_node] or out_paths[src_node][n1][0] >= shortest_paths[src_node][n1][0]):
                    out_paths[src_node][n1] = copy.deepcopy(shortest_paths[src_node].get(n1, (np.Inf, -1, [])))

    SP = {}
    for i in sources.keys():
        SP[i] = {}
        for j in sinks.keys():
            p = out_paths[i].get(j, (np.Inf, -1, []))
            if len(p[-1]) > 1:
                p = (p[0], p[1], p[2][1:])
            SP[i][j] = p
    #return {i: {j: out_paths[i].get(j, (np.Inf, -1, [])) for j in sinks.keys()} for i in sources.keys()}
    return SP

def shortestPath2(TS, sources, sinks, immuned, unreported, parameter = 1.0): #src: node -> time, sink: node-> time, immuned: node -> time

    shortest_paths = {i: {} for i in sources.keys()}
    out_paths = {i: {} for i in sources.keys()}
    c = 0
    for interaction in TS:
        c += 1
        if c % 100 == 0:
            print c
        t, n1, n2 = interaction[0], interaction[1], interaction[2]
        for src_node, t_start in sources.iteritems():

            if (n1 not in immuned or t < immuned[n1]):
                prerecovered_penalty = 0.0

                if n1 == src_node:
                    preinfected_penalty = get_preinfected_penalty(t, sources[n1], parameter) if n1 in sources else 0.0
                    #if n1 not in sinks:
                    #    preinfected_penalty = len(TS)
                    dn1 = 0.0 + preinfected_penalty + prerecovered_penalty
                    #dn1 = 0.0
                    if src_node not in shortest_paths[src_node] or shortest_paths[src_node][n1][0] >= dn1:
                        shortest_paths[src_node][src_node] = (dn1, t, [(t, n1, n1)])

                if n1 in shortest_paths[src_node].keys():

                    prime = shortest_paths[src_node][n1]

                    if n1 in sinks:
                        penalty_n1 = get_preinfected_penalty(sources[n1], t, parameter) if prime[1] < sources[n1] else get_preinfected_penalty(prime[1], t, parameter)
                        #penalty_n1 = get_preinfected_penalty(sources[n1], t, parameter) if prime[1] < sources[n1] else get_preinfected_penalty(prime[1], t, parameter)
                        #penalty_n1 = parameter*0.5*np.abs((sources[n1]-t).total_seconds())
                    else:
                        penalty_n1 = get_preinfected_penalty(unreported[n1], t, parameter) if prime[1] < unreported[n1] else get_preinfected_penalty(prime[1], t, parameter)
                        #penalty_n1 = parameter*0.5*np.abs((unreported[n1]-t).total_seconds())
                    if n2 in sinks:
                        penalty_n2 = get_preinfected_penalty(t, sources[n2], parameter) if t < sources[n2] else get_preinfected_penalty(sources[n2], t, parameter)
                        #penalty_n2 = parameter*0.5*np.abs((sources[n2]-t).total_seconds())
                    else:
                        penalty_n2 = get_preinfected_penalty(t, unreported[n2], parameter) if t < unreported[n2] else get_preinfected_penalty(unreported[n2], t, parameter)
                        #penalty_n2 = parameter*0.5*np.abs((unreported[n2]-t).total_seconds())
                    #dn2 = prime[0] + 1.0 + penalty_n1 + penalty_n2
                    dn2 = prime[0] + penalty_n1 + penalty_n2
                    pathn2 = prime[2] + [(t, n1, n2)]
                    an2 = t
                    if n2 not in immuned or t < immuned[n2]: # if node is not immune yet
                        if n2 not in shortest_paths[src_node] or shortest_paths[src_node][n2][0] >= dn2:
                            shortest_paths[src_node][n2] = (dn2, an2, pathn2)
                    if n2 in sinks and sinks[n2] >= t and (n2 not in out_paths[src_node] or out_paths[src_node][n2][0] >= shortest_paths[src_node][n2][0]):
                        out_paths[src_node][n2] = copy.deepcopy(shortest_paths[src_node].get(n2, (np.Inf, -1, [])))
                        #c = c
                    if n1 in sinks and sinks[n1] >= t and (n1 not in out_paths[src_node] or out_paths[src_node][n1][0] >= shortest_paths[src_node][n1][0]):
                        out_paths[src_node][n1] = copy.deepcopy(shortest_paths[src_node].get(n1, (np.Inf, -1, [])))
                        #c = c
    SP = {}
    for i in sources.keys():
        SP[i] = {}
        for j in sinks.keys():
            p = out_paths[i].get(j, (np.Inf, -1, []))
            if len(p[-1]) > 1:
                p = (p[0], p[1], p[2][1:])
            SP[i][j] = p
    #return {i: {j: out_paths[i].get(j, (np.Inf, -1, [])) for j in sinks.keys()} for i in sources.keys()}
    return SP

def shortestPath3(TS, sources, sinks, immuned, unreported, parameter = 1.0): #src: node -> time, sink: node-> time, immuned: node -> time
    L_N = get_L_N(TS[-1][0]-TS[0][0])
    shortest_paths = {i: {} for i in sources.keys()}
    out_paths = {i: {} for i in sources.keys()}
    infLog = np.log(0.5)
    c = 0
    for interaction in TS:
        c += 1
        if c % 100 == 0:
            print c
        t, n1, n2 = interaction[0], interaction[1], interaction[2]
        for src_node, t_start in sources.iteritems():

            if (n1 not in immuned or t < immuned[n1]):
                #prerecovered_penalty = 0.0

                if n1 == src_node:
                    preinfected_penalty = get_preinfected_penalty(t, sources[n1], parameter) if n1 in sources else 0.0
                    #dn1 = np.log(preinfected_penalty + 1.0) + prerecovered_penalty
                    dn1 = np.log(preinfected_penalty + 1.0)
                    if src_node not in shortest_paths[src_node] or shortest_paths[src_node][n1][0] > dn1:
                        shortest_paths[src_node][src_node] = (dn1, t, [(t, n1, n1)])

                if n1 in shortest_paths[src_node].keys():

                    prime = shortest_paths[src_node][n1]

                    if n1 in sinks:
                        #penalty_n1 = get_preinfected_penalty(sources[n1], t, parameter) if prime[1] < sources[n1] else get_preinfected_penalty(prime[1], t, parameter)
                        penalty_n1 = get_preinfected_penalty(prime[1], t, parameter)
                    else:
                        #penalty_n1 = get_preinfected_penalty(unreported[n1], t, parameter) if prime[1] < unreported[n1] else get_preinfected_penalty(prime[1], t, parameter)
                        penalty_n1 = get_preinfected_penalty(prime[1], t, parameter)
                    if n2 in sinks:
                        penalty_n2 = get_preinfected_penalty(t, sources[n2], parameter) if t < sources[n2] else get_preinfected_penalty(sources[n2], t, parameter)
                        #penalty_n2 = get_preinfected_penalty(t, sources[n2], parameter) if t < sources[n2] else get_preinfected_penalty(sources[n2], t, parameter)
                    else:
                        penalty_n2 = get_preinfected_penalty(t, unreported[n2], parameter) if t < unreported[n2] else get_preinfected_penalty(unreported[n2], t, parameter)

                    dn2 = prime[0] + infLog + infLog*np.log(penalty_n1 + 1.0) + L_N[penalty_n2] + 1.0
                    pathn2 = prime[2] + [(t, n1, n2)]
                    an2 = t
                    if n2 not in immuned or t < immuned[n2]: # if node is not immune yet
                        if n2 not in shortest_paths[src_node] or shortest_paths[src_node][n2][0] > dn2:
                            shortest_paths[src_node][n2] = (dn2, an2, pathn2)
                    if n2 in sinks and sinks[n2] >= t and (n2 not in out_paths[src_node] or out_paths[src_node][n2][0] >= shortest_paths[src_node][n2][0]):
                        out_paths[src_node][n2] = copy.deepcopy(shortest_paths[src_node].get(n2, (np.Inf, -1, [])))
                        #c = c
                    if n1 in sinks and sinks[n1] >= t and (n1 not in out_paths[src_node] or out_paths[src_node][n1][0] >= shortest_paths[src_node][n1][0]):
                        out_paths[src_node][n1] = copy.deepcopy(shortest_paths[src_node].get(n1, (np.Inf, -1, [])))
                        #c = c
    SP = {}
    for i in sources.keys():
        SP[i] = {}
        for j in sinks.keys():
            p = out_paths[i].get(j, (np.Inf, -1, []))
            if len(p[-1]) > 1:
                p = (p[0], p[1], p[2][1:])
            SP[i][j] = p
    #return {i: {j: out_paths[i].get(j, (np.Inf, -1, [])) for j in sinks.keys()} for i in sources.keys()}
    return SP



if __name__ == "__main__":

    p = 0.7
    type = 'ER'
    M = 100
    TS, snapshots, G = g.generateTS(n = 10,
                                  p = 0.5,
                                  seed = 1.0,
                                  st = datetime(2000, 01, 01, 00, 00, 00),
                                  m = M,
                                  srcN = 1,
                                  reportingP = 1.0,
                                  infectionP = 0.9,
                                  recoveringP = 0.1,
                                  type = type)
    #TS, snapshots, G,_ = g.readFile('generated.txt')
    sources, immuned = get_sinks_and_sources(TS)
    print sources
    print immuned
    #exit()
    SP = shortestPath(TS, sources, sources, immuned)
    #print shortestPath(TS, (1, datetime(2000, 01, 01, 00, 00, 00)), [(6, datetime(2000, 01, 01, 00, 00, 07))])
    exit()
    #print len(terminals)
    weights = w.get_flow_weights(TS, p)
    for i in weights:
        print i, weights[i]

    costs = w.get_GMC_costs(weights)

    m = 3
    unreachable = True
    while unreachable == True:
        covering, _ = GMC.cover(weights, costs, m)
        unreachable = False
        print covering
        for k,v in covering.iteritems():
            if v[-1] == 0.0:
                unreachable = True
        m += 1
    print m-1
    print unreachable
    print covering
    #exit()
    selectedBins = {}
    for item, bin in covering.iteritems():
        #print item, '||' , bin
        selectedBins[bin[0]] = selectedBins.get(bin[0],[]) + [item]

    print 'picked bins'
    print len(selectedBins)
    for i in selectedBins:
        print i, selectedBins[i]

    paths = {}
    #print 'debug'
    print shortestPath(TS, (1, datetime(2000, 01, 01, 00, 00, 00)), [(6, datetime(2000, 01, 01, 00, 00, 07))])


    for bin, items in selectedBins.iteritems():
        paths[bin] = shortestPath(TS, bin, items)

    #print acc.snapshot_accuracy(TS, paths, n)
    print 'paths'
    for i in paths:
        print i, paths[i]
    #print acc.getPrecision(TS, paths)
    #print acc.accuracy_interactions(TS, paths, infected_interactions)

    # print 'paths'
    # output = set(weights.keys())
    # for src, covered in paths.iteritems():
    #     #print p, paths[p]
    #     for target, path in covered.iteritems():
    #         for p in path[1]:
    #             output.add((p[1],p[0]))
    #             output.add((p[2],p[0]))
    #             #print p
    #         #exit()
    #         #output.add()
    # print output



