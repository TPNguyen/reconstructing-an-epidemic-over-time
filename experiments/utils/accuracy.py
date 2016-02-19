__author__ = 'Polina'
import os.path
import copy
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import operator
import scipy.stats as stats
from itertools import groupby


def get_infection_paths_noise(TS, paths, seeds):
    out = {}
    seed_time = {}
    for p in paths:
        for e in p:
            if e[0] not in out:
                out[e[0]] = set()
            out[e[0]].add((e[1], e[2]))

            if e[1] in seeds and (e[1] not in seed_time or seed_time[e[1]] > e[0]):
                seed_time[e[1]] = e[0]
    out_sorted_keys = sorted(out.keys())
    found_edge_list = set()

    H = nx.DiGraph()
    last_instance = {}
    reported = {}
    found = {}
    inf_moment = {}
    seed_nodelist = []
    counter = 0
    for idx in xrange(len(TS)):
        i = TS[idx]
        t, n1, n2, inf1, inf2, rep1, rep2, ext1, ext2 = i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8]
        t = idx
        H.add_node((n1, t), inf = inf1, rep = rep1, ext = ext1)
        H.add_node((n2, t), inf = inf2, rep = rep2, ext = ext2)
        H.add_edge((n1, t), (n2, t))
        if inf1 == 1 and n1 not in inf_moment:
            inf_moment[n1] = t
        if inf2 == 1 and n2 not in inf_moment:
            inf_moment[n2] = t
            H[(n1, t)][(n2, t)]['inf'] = 1

        if n1 in seed_time and i[0] == seed_time[n1]:
            seed_nodelist.append((n1,t))


        if rep1:
            reported[n1] = 1
        if rep2:
            reported[n2] = 1
        if n1 in last_instance:
            H.add_edge((n1, last_instance[n1]), (n1, t))
            H.node[(n1, t)]['rep'] = reported.get(n1,0)
        if n2 in last_instance:
            H.add_edge((n2, last_instance[n2]), (n2, t))
            H.node[(n2, t)]['rep'] = reported.get(n2,0)
        last_instance[n1] = t
        last_instance[n2] = t

        if counter < len(out_sorted_keys) and out_sorted_keys[counter] == TS[idx][0]:
            for i in out[out_sorted_keys[counter]]:
                found_edge_list.add(((i[0], t), (i[1], t)))
                found[i[0]] = 1
                found[i[1]] = 1
            counter += 1

        H.node[(n1, t)]['found'] = found.get(n1, 0)
        H.node[(n2, t)]['found'] = found.get(n2, 0)
    print 'H', H.number_of_nodes(), H.number_of_edges()
    pos = {i: (i[0], len(TS)-i[1]) for i in H.nodes_iter()}

    node_color = {}
    node_alpha = {}
    linewidths = []
    node_shape = []
    true_seeds = []
    for i in H.nodes_iter(data = True):
        attr = i[1]
        if attr['inf'] == 1:
            c = 'red'
            if attr['rep'] == 1:
                c = 'cyan'
            if attr['ext'] == 1:
                c = 'magenta'
                true_seeds.append(i[0])
                if attr['rep'] == 1:
                    c = 'pink'
        else:
            c = 'grey'

        if attr['found'] == 1:
            linewidths.append(2)
            node_shape.append('8')
            alpha = 1.0
        else:
            linewidths.append(1)
            node_shape.append('o')
            alpha = 0.5

        #node_color.append(c)
        node_color[i[0]] = c
        node_alpha[i[0]] = alpha

    edge_color = {}
    for i in H.edges_iter(data = True):
        if i[-1] and i[-1]['inf'] == 1:
            edge_color[(i[0],i[1])] = 'black'

            #edge_color.append('black')
        else:
            edge_color[(i[0],i[1])] = 'grey'
            #edge_color.append('grey')

    print H.number_of_nodes()
    print len([i[0] for i in H.nodes_iter(data=True) if i[1]['found'] == 1])


    nx.draw_networkx_edges(H, pos = pos, edge_color = 'grey', width = 1.5);
    nx.draw_networkx_edges(H, pos = pos, edgelist = list(found_edge_list),
                           edge_color = 'black', width = 2.0)
    nx.draw_networkx_edges(H, pos = pos, edgelist = [(i[0],i[1]) for i in H.edges_iter(data = True) if i[-1] and i[-1]['inf'] == 1],
                           edge_color = 'red', width = 1.5)
    correct_edges = found_edge_list.intersection(set([(i[0],i[1]) for i in H.edges_iter(data = True) if i[-1] and i[-1]['inf'] == 1]))
    nx.draw_networkx_edges(H, pos = pos, edgelist = list(correct_edges),
                           edge_color = 'green', width = 1.5)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = [i[0] for i in H.nodes_iter(data=True) if i[1]['found'] != 1],
                           node_color = [node_color[i[0]] for i in H.nodes_iter(data=True) if i[1]['found'] != 1],
                           linewidths = 1.0, alpha = 0.5)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = [i[0] for i in H.nodes_iter(data=True) if i[1]['found'] == 1],
                           node_color = [node_color[i[0]] for i in H.nodes_iter(data=True) if i[1]['found'] == 1],
                           linewidths = 2.0, alpha = 1.0)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = seed_nodelist,
                           node_color = 'orange',
                           linewidths = 2.0, alpha = 1.0)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = set(seed_nodelist).intersection(set(true_seeds)),
                           node_color = 'green',
                           linewidths = 2.0, alpha = 1.0)

    plt.show()
    #exit()
    #nx.draw(H, pos = pos, node_color = node_color, edge_color = edge_color, node_shape = node_shape)

    #plt.show()
    return H

def get_infection_paths(TS, paths, seeds):
    out = {}
    seed_time = {}
    for p in paths:
        for e in p:
            if e[0] not in out:
                out[e[0]] = set()
            out[e[0]].add((e[1], e[2]))

            if e[1] in seeds and (e[1] not in seed_time or seed_time[e[1]] > e[0]):
                seed_time[e[1]] = e[0]
    out_sorted_keys = sorted(out.keys())
    found_edge_list = set()

    H = nx.DiGraph()
    last_instance = {}
    reported = {}
    found = {}
    inf_moment = {}
    seed_nodelist = []
    counter = 0
    for idx in xrange(len(TS)):
        i = TS[idx]
        t, n1, n2, inf1, inf2, rep, ext = i[0], i[1], i[2], i[3], i[4], i[5], i[6]
        t = idx
        H.add_node((n1, t), inf = inf1, rep = rep, ext = ext)
        H.add_node((n2, t), inf = inf2, rep = 0, ext = 0)
        H.add_edge((n1, t), (n2, t))
        if inf1 == 1 and n1 not in inf_moment:
            inf_moment[n1] = t
        if inf2 == 1 and n2 not in inf_moment:
            inf_moment[n2] = t
            H[(n1, t)][(n2, t)]['inf'] = 1

        if n1 in seed_time and i[0] == seed_time[n1]:
            seed_nodelist.append((n1,t))


        if rep:
            reported[n1] = 1
        if n1 in last_instance:
            H.add_edge((n1, last_instance[n1]), (n1, t))
            H.node[(n1, t)]['rep'] = reported.get(n1,0)
        if n2 in last_instance:
            H.add_edge((n2, last_instance[n2]), (n2, t))
            H.node[(n2, t)]['rep'] = reported.get(n2,0)
        last_instance[n1] = t
        last_instance[n2] = t

        if counter < len(out_sorted_keys) and out_sorted_keys[counter] == TS[idx][0]:
            for i in out[out_sorted_keys[counter]]:
                found_edge_list.add(((i[0], t), (i[1], t)))
                found[i[0]] = 1
                found[i[1]] = 1
            counter += 1

        H.node[(n1, t)]['found'] = found.get(n1, 0)
        H.node[(n2, t)]['found'] = found.get(n2, 0)
    print H.number_of_nodes(), H.number_of_edges()
    pos = {i: (i[0], len(TS)-i[1]) for i in H.nodes_iter()}
    node_color = {}
    node_alpha = {}
    linewidths = []
    node_shape = []
    true_seeds = []
    for i in H.nodes_iter(data = True):
        attr = i[1]
        if attr['inf'] == 1:
            c = 'red'
            if attr['rep'] == 1:
                c = 'cyan'
            if attr['ext'] == 1:
                c = 'magenta'
                true_seeds.append(i[0])
                if attr['rep'] == 1:
                    c = 'pink'
        else:
            c = 'grey'

        if attr['found'] == 1:
            linewidths.append(2)
            node_shape.append('8')
            alpha = 1.0
        else:
            linewidths.append(1)
            node_shape.append('o')
            alpha = 0.5

        #node_color.append(c)
        node_color[i[0]] = c
        node_alpha[i[0]] = alpha

    edge_color = {}
    for i in H.edges_iter(data = True):
        if i[-1] and i[-1]['inf'] == 1:
            edge_color[(i[0],i[1])] = 'black'

            #edge_color.append('black')
        else:
            edge_color[(i[0],i[1])] = 'grey'
            #edge_color.append('grey')

    print H.number_of_nodes()
    print len([i[0] for i in H.nodes_iter(data=True) if i[1]['found'] == 1])


    nx.draw_networkx_edges(H, pos = pos, edge_color = 'grey', width = 1.5);
    nx.draw_networkx_edges(H, pos = pos, edgelist = list(found_edge_list),
                           edge_color = 'black', width = 2.0)
    nx.draw_networkx_edges(H, pos = pos, edgelist = [(i[0],i[1]) for i in H.edges_iter(data = True) if i[-1] and i[-1]['inf'] == 1],
                           edge_color = 'red', width = 1.5)
    correct_edges = found_edge_list.intersection(set([(i[0],i[1]) for i in H.edges_iter(data = True) if i[-1] and i[-1]['inf'] == 1]))
    nx.draw_networkx_edges(H, pos = pos, edgelist = list(correct_edges),
                           edge_color = 'green', width = 1.5)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = [i[0] for i in H.nodes_iter(data=True) if i[1]['found'] != 1],
                           node_color = [node_color[i[0]] for i in H.nodes_iter(data=True) if i[1]['found'] != 1],
                           linewidths = 1.0, alpha = 0.5)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = [i[0] for i in H.nodes_iter(data=True) if i[1]['found'] == 1],
                           node_color = [node_color[i[0]] for i in H.nodes_iter(data=True) if i[1]['found'] == 1],
                           linewidths = 2.0, alpha = 1.0)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = seed_nodelist,
                           node_color = 'orange',
                           linewidths = 2.0, alpha = 1.0)

    nx.draw_networkx_nodes(H, pos = pos, nodelist = set(seed_nodelist).intersection(set(true_seeds)),
                           node_color = 'green',
                           linewidths = 2.0, alpha = 1.0)

    plt.show()
    #exit()
    #nx.draw(H, pos = pos, node_color = node_color, edge_color = edge_color, node_shape = node_shape)

    #plt.show()
    return H

def get_infection_tree(TS):
    H = nx.DiGraph()
    last_instance = {}
    reported = {}
    inf_moment = {}
    for idx in xrange(len(TS)):
        i = TS[idx]
        t, n1, n2, inf1, inf2, rep, ext = i[0], i[1], i[2], i[3], i[4], i[5], i[6]
        t = idx
        H.add_node((n1, t), inf = inf1, rep = rep, ext = ext)
        H.add_node((n2, t), inf = inf2, rep = 0, ext = 0)
        H.add_edge((n1, t), (n2, t))
        if inf1 == 1 and n1 not in inf_moment:
            inf_moment[n1] = t
        if inf2 == 1 and n2 not in inf_moment:
            inf_moment[n2] = t
            H[(n1, t)][(n2, t)]['inf'] = 1


        if rep:
            reported[n1] = 1
        if n1 in last_instance:
            H.add_edge((n1, last_instance[n1]), (n1, t))
            H.node[(n1, t)]['rep'] = reported.get(n1,0)
        if n2 in last_instance:
            H.add_edge((n2, last_instance[n2]), (n2, t))
            H.node[(n2, t)]['rep'] = reported.get(n2,0)
        last_instance[n1] = t
        last_instance[n2] = t
    print H.number_of_nodes(), H.number_of_edges()
    pos = {i: (i[0], len(TS)-i[1]) for i in H.nodes_iter()}
    node_color = []
    for i in H.nodes_iter():
        attr = H.node[i]
        if attr['inf'] == 1:
            c = 'red'
            if attr['rep'] == 1:
                c = 'cyan'
            if attr['ext'] == 1:
                c = 'magenta'
                if attr['rep'] == 1:
                    c = 'pink'
        else:
            c = 'grey'

        node_color.append(c)

    edge_color = []
    for i in H.edges_iter(data = True):

        if i[-1] and i[-1]['inf'] == 1:
            edge_color.append('red')
        else:
            edge_color.append('black')
    nx.draw(H, pos=pos, node_color = node_color, edge_color = edge_color)

    plt.show()
    return H

def plotSnapshots(snapshots, snapshots_GT, G, mode, folder, pos, lb_snapshots = [], coloring = []):
    print mode
    maxColoring = 0.0
    if coloring:
        for i in coloring.values():
            for j in i.values():
                if j > maxColoring:
                    maxColoring = float(j)
    nodelist = G.nodes()
    index = {nodelist[i]: i for i in xrange(len(nodelist))}
    node_color = ['grey' for i in nodelist]
    node_size = [1 for i in nodelist]
    sorted_snapshots = sorted(snapshots.keys())
    #folder = time.strftime(mode+"_%Y%m%d-%H%M%S")
    #os.mkdir(folder)
    counts = {}
    for t in sorted_snapshots:
        for n in snapshots[t]['infected']:
            if coloring:
                node_color[index[n]] = str(1.0 - coloring[t][n]/maxColoring)
            else:
                node_color[index[n]] = 'grey'
            if lb_snapshots and n not in lb_snapshots[t]['infected']:
                node_color[index[n]] = 'royalblue'
            if n not in snapshots_GT[t]['infected']:
                node_color[index[n]] = 'orange'
            node_size[index[n]] = 50

        for n in snapshots_GT[t]['infected']:
            if n not in snapshots[t]['infected']:
                node_color[index[n]] = 'cyan'
            node_size[index[n]] = 50

        for n in snapshots[t]['recovered']:
            node_color[index[n]] = 'white'
            node_size[index[n]] = 50

        for n in snapshots[t]['seeds']:
            if n in snapshots_GT[t]['seeds']:
                node_color[index[n]] = 'r'
                if lb_snapshots and n not in lb_snapshots[t]['infected']:
                    node_color[index[n]] = 'magenta'
            else:
                node_color[index[n]] = 'black'
            node_size[index[n]] = 50

        #pos = nx.spectral_layout(G)
        #pos = nx.spring_layout(G)
        pos = dict(zip(G,G))
        plt.figure()
        nx.draw_networkx_nodes(G, pos, nodelist = nodelist, node_size = node_size, labels = nodelist, node_color = node_color, title = 'mode')
        nx.draw_networkx_edges(G, pos, edge_color = 'lightgray')

        filename = os.path.join(folder, mode+ '_' + t.strftime("%Y-%m-%d-%H-%M-%S")+'.png')
        plt.axis('off')
        plt.savefig(filename, bbox_inches="tight")
        plt.close()
        #plt.show()
    return

def getOutput(paths):
    output = {}
    for src, covered in paths.iteritems():
        #print p, paths[p]
        for target, path in covered.iteritems():
            for p in path[1]:
                output[p[0]] = output.get(p[0], set())
                output[p[0]].add(p[1])
                output[p[0]].add(p[2])
                #print p
            #exit()
            #output.add()
    return output

def getPrecision(TS, paths):
    output = getOutput(paths)
    TP = 0.0
    for i in TS:
        t = i[0]
        if t in output.keys():
            if i[1] in output[t]:
                TP += 1.0
            if i[2] in output[t]:
                TP += 1.0
    P = 0.0
    for k, v in output.iteritems():
        P += len(v)

    print P
    precision = TP/P
    print precision
    return precision


def how_far_seeds(G, seeds, found_seeds):
    res = []
    while seeds:
        minSP_global, closest_from, closest_to = np.Inf, -1, -1
        for u in seeds:
            minSP, closest = np.Inf, -1
            for v in found_seeds:
                if nx.has_path(G, u, v):
                    s = nx.shortest_path_length(G, u, v)
                    if s < minSP:
                        minSP, closest = s, v
            if minSP < minSP_global:
                minSP_global, closest_from, closest_to = minSP, u, closest
        if closest_from != -1:
            seeds.remove(closest_from)
            found_seeds.remove(closest_to)
            res.append(minSP_global)
        else:
            break
    print len(res)
    return res

def how_far_snapshots(snapshots_gt, snapshots_out):
    sorted_idx = sorted(snapshots_gt.keys())
    gt_moment_of_infection = {}
    out_moment_of_infection = {}
    n1, n2 = 0.0, 0.0
    for i in sorted_idx:
        #n1_l, n2_l = 0.0, 0.0
        for j in snapshots_gt[i]['infected']:
            if j not in gt_moment_of_infection.keys():
                gt_moment_of_infection[j] = i
                #n1_l += 1.0
        for j in snapshots_out[i]['infected']:
            if j not in out_moment_of_infection.keys():
                out_moment_of_infection[j] = i
                #n2_l += 1.0
        # if n1_l > 1:
        #     n1 += n1_l*(n1_l-1)/2.0
        # if n2_l > 1:
        #     n2 += n2_l*(n2_l-1)/2.0

    s1_out = [i for i in sorted(out_moment_of_infection.items(), key=operator.itemgetter(1)) if i[0] in gt_moment_of_infection]
    sorted1_out = [i[0] for i in s1_out]
    s1_ties = [len(list(group)) for key, group in groupby([i[1] for i in s1_out])]
    for i in s1_ties:
        if i > 1:
            n1 += i*(i-1)/2.0

    s2_out = [i for i in sorted(gt_moment_of_infection.items(), key=operator.itemgetter(1)) if i[0] in out_moment_of_infection]
    sorted2_out = [i[0] for i in s2_out]
    s2_ties = [len(list(group)) for key, group in groupby([i[1] for i in s2_out])]
    for i in s2_ties:
        if i > 1:
            n2 += i*(i-1)/2.0

    #sorted_out = [i[0] for i in sorted(out_moment_of_infection.items(), key=operator.itemgetter(1)) if i[0] in gt_moment_of_infection]

    #sorted_gt = [i[0] for i in sorted(gt_moment_of_infection.items(), key=operator.itemgetter(1)) if i[0] in out_moment_of_infection]
    tau, _ = stats.kendalltau(sorted1_out, sorted2_out)
    n = len(sorted1_out)
    #tau = (tau*n*(n-1)/2.0)/(np.sqrt((n-n1)*(n-n2)))

    return tau

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
        #print order_gt, order, p
        tau = stats.kendalltau(order, order_gt)
        res.append(tau[0])
        lengths.append(len(order))
    # plt.figure('path lengths')
    # plt.title('path length')
    # plt.plot(sorted(lengths))
    # plt.show()
    print 'sus+-'
    print 'GT and output infection intersection', len(infected.intersection(infected_gt))
    print 'output as infected, but not in GT', len(infected.difference(infected_gt))
    print 'output as infected, but not in reported', len(infected.difference(reported_infected))
    print 'output as infected, in GT, but not in reported', len((infected.intersection(infected_gt)).difference(reported_infected))
    print len(infected_gt.difference(infected))
    print len(reported_infected.difference(infected)), len(infected.difference(reported_infected))
    print len(infected), len(infected_gt), len(reported_infected)
    return res, lengths

def how_far_intime(paths, moment_of_infection, mode = 'abs'):
    res = []
    out_moment_of_infection = {}
    for p in paths:
        for step in p:
            if step[1] in out_moment_of_infection:
                out_moment_of_infection[step[1]] = min(out_moment_of_infection[step[1]], step[0])
            else:
                out_moment_of_infection[step[1]] = step[0]

            if step[2] in out_moment_of_infection:
                out_moment_of_infection[step[2]] = min(out_moment_of_infection[step[2]], step[0])
            else:
                out_moment_of_infection[step[2]] = step[0]
    sorted_out = [i[0] for i in sorted(out_moment_of_infection.items(), key=operator.itemgetter(1)) if i[0] in moment_of_infection]
    sorted_gt = [i[0] for i in sorted(moment_of_infection.items(), key=operator.itemgetter(1)) if i[0] in out_moment_of_infection]

    for k, v in out_moment_of_infection.iteritems():
        if k in moment_of_infection:
            if v > moment_of_infection[k]:
                t = (v - moment_of_infection[k]).total_seconds()
            else:
                t = -(moment_of_infection[k] - v).total_seconds()
                if mode == 'abs':
                    t = np.abs(t)
            res.append(t)
    #return res, stats.kendalltau(sorted_gt, sorted_out), stats.pearsonr(sorted_gt, sorted_out)
    try:
        tau = stats.kendalltau(sorted_gt, sorted_out)
    except:
        tau = 0.0
    return res, tau


def get_snapshots(TS, m, num_nodes, threshold = 0.75, maxTicks = False):
    coloring, local_coloring = {}, {}
    moment_of_infection = {}
    infected = set()
    recovered = set()
    seeds = set()
    snapshots = {}
    output_snapshots = {}
    for i in xrange(len(TS)):
        if TS[i][5] == 1:
            first_report = i
            break
    if maxTicks == True:
        ticks = [0] + range(first_report, len(TS))
    else:
        ticks = [0] + range(first_report, len(TS), (len(TS)-first_report)/m) + [len(TS)-1]

    for i in xrange(0, len(ticks)-1):
        for j in TS[ticks[i]:ticks[i+1]+1]:
            if j[3] == 1:
                if j[1] not in infected:
                    moment_of_infection[j[1]] = j[0]
                infected.add(j[1])
                local_coloring[j[1]] = local_coloring.get(j[1], 0) + 1
            elif j[3] == -1:
                recovered.add(j[1])
                if j[1] in infected:
                    infected.remove(j[1])
                if j[1] in seeds:
                    seeds.remove(j[1])
            # take into account only acting nodes
            if j[4] == 1:
                if j[2] not in infected:
                    moment_of_infection[j[2]] = j[0]
                    infected.add(j[2])
            # elif j[4] == -1:
            #     recovered.add(j[2])
            #     if j[2] in infected:
            #         infected.remove(j[2])
            #     if j[2] in seeds:
            #         seeds.remove(j[2])
            if j[6] == 1:
                seeds.add(j[1])
        output_snapshots['infected'] = copy.deepcopy(infected)
        output_snapshots['recovered'] = copy.deepcopy(recovered)
        output_snapshots['seeds'] = copy.deepcopy(seeds)
        snapshots[TS[ticks[i+1]][0]] = copy.deepcopy(output_snapshots)
        coloring[TS[ticks[i+1]][0]] = copy.deepcopy(local_coloring)
        #if float(len(infected))/num_nodes > threshold:
        #    print 'everything is infected by ', i
        #    break
    return snapshots, coloring, moment_of_infection, seeds, ticks[1:]


def get_added_snapshots(added, snapshots, immuned):
    output_infected = set()
    output_recovered = set()
    output_seeds = set()
    found_seeds = {}
    output_snapshots = {}
    output = {}
    output_interactions = set()
    node_activity = {}

    sorted_output = sorted(added.items(), key=operator.itemgetter(1))

    snapshots_time = sorted(snapshots.keys())
    iter_output = 0
    to_recover = set()
    for i in xrange(len(snapshots_time)):
        t1 = snapshots_time[i]
        while iter_output < len(sorted_output) and sorted_output[iter_output][1] <= t1:

            n1 = sorted_output[iter_output][0]

            output_infected.add(n1)

            if n1 in immuned:
                if sorted_output[iter_output][1] >= immuned[n1]:
                    output_recovered.add(n1)
                    if n1 in output_infected:
                        output_infected.remove(n1)

            iter_output += 1

        output_snapshots['infected'] = copy.deepcopy(output_infected)
        output_snapshots['recovered'] = copy.deepcopy(output_recovered)
        output[t1] = copy.deepcopy(output_snapshots)

    return output


def get_output_snapshots_no_recov_pred(pths, snapshots, immuned):# pths - list of lists
    output_infected = set()
    output_recovered = set()
    output_seeds = set()
    found_seeds = {}
    output_snapshots = {}
    output = {}
    output_interactions = set()
    node_activity = {}

    for p in pths:
        for step in p:
            output_interactions.add(step)
            # if step[1] in immuned:
            #     node_activity[step[1]] = node_activity.get(step[1], []) + [step[0]]
            # if step[2] in immuned:
            #     node_activity[step[2]] = node_activity.get(step[2], []) + [step[0]]
        if p:
            found_seeds[p[0][1]] = min(p[0][0], found_seeds.get(p[0][1], p[0][0]))


    # for node in node_activity:
    #     node_activity[node].sort()

    sorted_output = sorted(list(output_interactions))

    snapshots_time = sorted(snapshots.keys())
    iter_output = 0
    to_recover = set()
    for i in xrange(len(snapshots_time)):
        t1 = snapshots_time[i]
        #output_recovered.update(to_recover)
        #output_infected.difference_update(output_recovered)
        #output_seeds.difference_update(output_recovered)
        while iter_output < len(sorted_output) and sorted_output[iter_output][0] <= t1:
            #output_recovered.update(to_recover)
            #output_infected.difference_update(output_recovered)
            #output_seeds.difference_update(output_recovered)

            n1 = sorted_output[iter_output][1]
            n2 = sorted_output[iter_output][2]
            output_infected.add(n1)
            output_infected.add(n2)

            if n1 in immuned:
                if sorted_output[iter_output][0] >= immuned[n1]:
                    output_recovered.add(n1)
                    if n1 in output_infected:
                        output_infected.remove(n1)
                    if n1 in output_seeds:
                        output_seeds.remove(n1)

            if n2 in immuned:
                if sorted_output[iter_output][0] >= immuned[n2]:
                    output_recovered.add(n2)
                    if n2 in output_infected:
                        output_infected.remove(n2)
                    if n2 in output_seeds:
                        output_seeds.remove(n2)

            if n1 in found_seeds:
                if sorted_output[iter_output][0] >= found_seeds[n1]:
                    output_seeds.add(n1)

            if n2 in found_seeds:
                if sorted_output[iter_output][0] >= found_seeds[n2]:
                    output_seeds.add(n2)
            iter_output += 1

            # if n1 in found_seeds.keys() and found_seeds[n1] <= t1:
            #     output_seeds.add(n1)

        output_snapshots['seeds'] = copy.deepcopy(output_seeds)
        #if iter_output >= len(sorted_output):
        #    output_recovered.update(to_recover)
        #    output_infected.difference_update(output_recovered)
        output_snapshots['infected'] = copy.deepcopy(output_infected)
        output_snapshots['recovered'] = copy.deepcopy(output_recovered)
        output[t1] = copy.deepcopy(output_snapshots)
    return output, found_seeds


def get_output_snapshots(pths, snapshots, immuned):# pths - list of lists
    output_infected = set()
    output_recovered = set()
    output_seeds = set()
    found_seeds = set()
    output_snapshots = {}
    output = {}
    output_interactions = set()
    node_activity = {}


    for p in pths:
        for step in p:
            output_interactions.add(step)
            if step[1] in immuned:
                node_activity[step[1]] = node_activity.get(step[1], []) + [step[0]]
            if step[2] in immuned:
                node_activity[step[2]] = node_activity.get(step[2], []) + [step[0]]
        #found_seeds[p[0][1]] = found_seeds[p[0][0]]
        found_seeds[p[0][1]] = min(p[0][0], found_seeds.get(p[0][1], p[0][0]))


    for node in node_activity:
        node_activity[node].sort()

    sorted_output = sorted(list(output_interactions))

    snapshots_time = sorted(snapshots.keys())
    iter_output = 0
    to_recover = set()
    for i in xrange(len(snapshots_time)):
        t1 = snapshots_time[i]
        output_recovered.update(to_recover)
        output_infected.difference_update(output_recovered)
        while iter_output < len(sorted_output) and sorted_output[iter_output][0] <= t1:
            output_recovered.update(to_recover)
            output_infected.difference_update(output_recovered)

            n1 = sorted_output[iter_output][1]
            n2 = sorted_output[iter_output][2]
            output_infected.add(n1)
            output_infected.add(n2)
            if n1 in immuned:
                if sorted_output[iter_output][0] == node_activity[n1][-1]:
                    to_recover.add(n1)
                if sorted_output[iter_output][0] == immuned[n1]:
                    output_recovered.add(n1)
                    output_infected.remove(n1)
                    if n1 in output_seeds:
                        output_seeds.remove(n1)

            if n2 in immuned:
                if sorted_output[iter_output][0] == node_activity[n2][-1]:
                    to_recover.add(n2)
                if sorted_output[iter_output][0] == immuned[n2]:
                    output_recovered.add(n2)
                    output_infected.remove(n2)
                    if n2 in output_seeds:
                        output_seeds.remove(n2)
            iter_output += 1


        if n1 in found_seeds.keys() and found_seeds[n1] <= t1:
            output_seeds.add(n1)

        output_snapshots['seeds'] = copy.deepcopy(output_seeds)
        #if iter_output >= len(sorted_output):
        #    output_recovered.update(to_recover)
        #    output_infected.difference_update(output_recovered)
        output_snapshots['infected'] = copy.deepcopy(output_infected)
        output_snapshots['recovered'] = copy.deepcopy(output_recovered)
        output[t1] = copy.deepcopy(output_snapshots)
    return output


def get_lb_snapshots(sources, immune, reported, snapshots):
    #sorted_input = sorted([(i[1], i[0]) for i in sources.keys()])
    sorted_infected = sorted([(t, src) for src, t in sources.iteritems() if src in reported['infected']])
    sorted_immune = sorted([(t, node) for node, t in immune.iteritems() if src in reported['recovered']])
    snapshots_time = sorted(snapshots.keys())
    input_infected = set()
    input_immune = set()
    input_snapshots = {}
    output_snapshots = {}
    input_seeds = set()
    iter_infected, iter_immune = 0, 0
    for i in xrange(len(snapshots_time)):
        t1 = snapshots_time[i]
        while iter_infected < len(sorted_infected) and sorted_infected[iter_infected][0] <= t1:
            input_infected.add(sorted_infected[iter_infected][1])
            iter_infected += 1
        while iter_immune < len(sorted_immune) and sorted_immune[iter_immune][0] <= t1:
            input_immune.add(sorted_immune[iter_immune][1])
            input_infected.remove(sorted_immune[iter_immune][1])
            iter_immune += 1

        output_snapshots['infected'] = copy.deepcopy(input_infected)
        output_snapshots['recovered'] = copy.deepcopy(input_immune)
        output_snapshots['seeds'] = copy.deepcopy(input_seeds)
        input_snapshots[t1] = copy.deepcopy(output_snapshots)
    return input_snapshots


def get_ub_snapshots(TS, snapshots, sinks):
    interactions, causality = [], []
    snapshots_time = sorted(snapshots.keys())
    infected = set()
    recovered = set()
    input_seeds = set()
    #reported = set()
    snapshots = {}
    output_snapshots = {}
    iter = 0
    for i in xrange(len(snapshots_time)):
        t1 = snapshots_time[i]
        while iter < len(TS) and TS[iter][0] <= t1:
            record = TS[iter]
            t, n1, n2, inf1, inf2, rep1, rep2 = record[0], record[1], record[2], record[3], record[4], record[5], record[6]
            if n1 in sinks and sinks[n1] <= t:
                if n2 not in infected:
                    interactions.append((t, n1, n2))
                    causality.append((n1, n2))
                    infected.add(n2)
                infected.add(n1)
                #reported.add(n1)
            if n2 in sinks and sinks[n2] <= t:
                infected.add(n2)
                #reported.add(n2)
            iter += 1
        output_snapshots['infected'] = copy.deepcopy(infected)
        output_snapshots['recovered'] = copy.deepcopy(recovered)
        output_snapshots['seeds'] = copy.deepcopy(input_seeds)
        snapshots[t1] = copy.deepcopy(output_snapshots)
        #snapshots[t1] = copy.deepcopy(infected)
    print 'len of ubs', len(interactions), len(causality)
    return snapshots, interactions, causality

def get_ub_snapshots_cascade(TS, snapshots, sinks):
    interactions, causality = [], []
    snapshots_time = sorted(snapshots.keys())
    infected = set()
    recovered = set()
    input_seeds = set()
    snapshots = {}
    output_snapshots = {}
    iter = 0
    for i in xrange(len(snapshots_time)):
        t1 = snapshots_time[i]
        while iter < len(TS) and TS[iter][0] <= t1:
            if TS[iter][1] in sinks and sinks[TS[iter][1]] <= TS[iter][0]:
                if TS[iter][3] == 1:
                    infected.add(TS[iter][1])
                elif TS[iter][3] == -1:
                    recovered.add(TS[iter][1])
                    if TS[iter][1] in infected:
                        infected.remove(TS[iter][1])
            if TS[iter][2] in sinks and sinks[TS[iter][2]] <= TS[iter][0]:
                if TS[iter][4] == 1:
                    if TS[iter][2] not in infected:
                        interactions.append((TS[iter][0], TS[iter][1], TS[iter][2]))
                        causality.append((TS[iter][1], TS[iter][2]))
                    infected.add(TS[iter][2])

            if TS[iter][1] in infected and TS[iter][1] not in recovered:
                infected.add(TS[iter][2])
            iter += 1
        output_snapshots['infected'] = copy.deepcopy(infected)
        output_snapshots['recovered'] = copy.deepcopy(recovered)
        output_snapshots['seeds'] = copy.deepcopy(input_seeds)
        snapshots[t1] = copy.deepcopy(output_snapshots)
        #snapshots[t1] = copy.deepcopy(infected)
    return snapshots, interactions, causality


def snapshot_accuracy_uppstr(GT_snapshots, output_snapshots, pths, immuned, sources, reported, TS, G, num_nodes, mode = 'main', pred_recover = False, draw = False, folder = ''):

    precision_infected = []
    recall_infected = []
    abs_values_TP = []
    gt_values = []
    abs_values_T = []
    set_T = []
    set_gt_T = []
    MCC = []
    MCC_upstr = []
    F1 = []
    for k in sorted(GT_snapshots.keys()):
        #print k, output_snapshots
        TP_infected = float(len(GT_snapshots[k]['infected'] & output_snapshots[k]['infected']))
        TN_ = float(num_nodes - len(GT_snapshots[k]['infected'] | output_snapshots[k]['infected']))
        FP_ = float(len(output_snapshots[k]['infected'] - GT_snapshots[k]['infected']))
        FN_ = float(len(GT_snapshots[k]['infected'] - output_snapshots[k]['infected']))

        TP_infected_upstr = float(len(GT_snapshots[k]['infected'] & output_snapshots[k]['infected'] & GT_snapshots[k]['upstream']))
        TN_upstr = float(len(GT_snapshots[k]['upstream']) - len((GT_snapshots[k]['infected'] | output_snapshots[k]['infected']) & GT_snapshots[k]['upstream']))
        FP_upstr = float(len((output_snapshots[k]['infected'] - GT_snapshots[k]['infected']) & GT_snapshots[k]['upstream']))
        FN_upstr = float(len((GT_snapshots[k]['infected'] - output_snapshots[k]['infected']) & GT_snapshots[k]['upstream']))

        print num_nodes, len(GT_snapshots[k]['infected']), len(output_snapshots[k]['infected'])
        print len(GT_snapshots[k]['infected'] | output_snapshots[k]['infected'])

        #float(num_nodes - len(GT_snapshots[k]['infected'] | output_snapshots[k]['infected']))

        precision_infected.append(np.divide(TP_infected, len(output_snapshots[k]['infected'])))
        recall_infected.append(np.divide(TP_infected, len(GT_snapshots[k]['infected'])))

        abs_values_TP.append(TP_infected)
        gt_values.append(len(GT_snapshots[k]['infected']))
        abs_values_T.append(len(output_snapshots[k]['infected']))
        set_T.append(output_snapshots[k]['infected'])
        set_gt_T.append(GT_snapshots[k]['infected'])
        print 'accuracies',k, mode, TP_infected, FP_, TP_infected, FN_, TN_, FP_, TN_, FN_
        MCC.append(np.divide((TP_infected*TN_ - FP_*FN_), np.sqrt((TP_infected+FP_)*(TP_infected+FN_)*(TN_+FP_)*(TN_+FN_))))
        MCC_upstr.append(np.divide((TP_infected_upstr*TN_upstr - FP_upstr*FN_upstr), np.sqrt((TP_infected_upstr+FP_upstr)*(TP_infected_upstr+FN_upstr)*(TN_upstr+FP_upstr)*(TN_upstr+FN_upstr))))

        F1.append(np.divide(2.0*precision_infected[-1]*recall_infected[-1], (precision_infected[-1]+recall_infected[-1])))
        precision_recovered = []
    recall_recovered = []

    return precision_infected, recall_infected, precision_recovered, recall_recovered, abs_values_TP, gt_values, abs_values_T, set_T, set_gt_T, MCC, F1, MCC_upstr

def snapshot_accuracy(GT_snapshots, output_snapshots, pths, immuned, sources, reported, TS, G, num_nodes, mode = 'main', pred_recover = False, draw = False, folder = ''):

    precision_infected = []
    recall_infected = []
    abs_values_TP = []
    gt_values = []
    abs_values_T = []
    set_T = []
    set_gt_T = []
    MCC = []
    F1 = []
    for k in sorted(GT_snapshots.keys()):
        #print k, output_snapshots
        TP_infected = float(len(GT_snapshots[k]['infected'] & output_snapshots[k]['infected']))
        TN_ = float(num_nodes - len(GT_snapshots[k]['infected'] | output_snapshots[k]['infected']))
        FP_ = float(len(output_snapshots[k]['infected'] - GT_snapshots[k]['infected']))
        FN_ = float(len(GT_snapshots[k]['infected'] - output_snapshots[k]['infected']))

        print num_nodes, len(GT_snapshots[k]['infected']), len(output_snapshots[k]['infected'])
        print len(GT_snapshots[k]['infected'] | output_snapshots[k]['infected'])

        #float(num_nodes - len(GT_snapshots[k]['infected'] | output_snapshots[k]['infected']))

        precision_infected.append(np.divide(TP_infected, len(output_snapshots[k]['infected'])))
        recall_infected.append(np.divide(TP_infected, len(GT_snapshots[k]['infected'])))

        abs_values_TP.append(TP_infected)
        gt_values.append(len(GT_snapshots[k]['infected']))
        abs_values_T.append(len(output_snapshots[k]['infected']))
        set_T.append(output_snapshots[k]['infected'])
        set_gt_T.append(GT_snapshots[k]['infected'])
        print 'accuracies',k, mode, TP_infected, FP_, TP_infected, FN_, TN_, FP_, TN_, FN_
        MCC.append(np.divide((TP_infected*TN_ - FP_*FN_), np.sqrt((TP_infected+FP_)*(TP_infected+FN_)*(TN_+FP_)*(TN_+FN_))))

        F1.append(np.divide(2.0*precision_infected[-1]*recall_infected[-1], (precision_infected[-1]+recall_infected[-1])))
        precision_recovered = []
    recall_recovered = []

    return precision_infected, recall_infected, precision_recovered, recall_recovered, abs_values_TP, gt_values, abs_values_T, set_T, set_gt_T, MCC, F1


def snapshot_accuracy_lastonly(GT_snapshots, output_snapshots, pths, immuned, sources, reported, TS, G, num_nodes, mode = 'main', pred_recover = False, draw = False, folder = ''):

    precision_infected = []
    recall_infected = []
    abs_values_TP = []
    gt_values = []
    abs_values_T = []
    set_T = []
    set_gt_T = []
    MCC = []
    F1 = []
    last = sorted(GT_snapshots.keys())[-1]
    for k in sorted(GT_snapshots.keys()):
        TP_infected = float(len(GT_snapshots[last]['infected'] & output_snapshots[k]['infected']))

        TN_ = float(num_nodes - len(GT_snapshots[last]['infected'] | output_snapshots[k]['infected']))
        FP_ = float(len(output_snapshots[k]['infected'] - GT_snapshots[last]['infected']))
        FN_ = float(len(GT_snapshots[last]['infected'] - output_snapshots[k]['infected']))
        if len(output_snapshots[k]['infected']) == 0:
            precision_infected.append(0.0)
        else:
            precision_infected.append(TP_infected/len(output_snapshots[k]['infected']))
        if len(GT_snapshots[last]['infected']) == 0:
            recall_infected.append(0.0)
        else:
            recall_infected.append(TP_infected/len(GT_snapshots[last]['infected']))

        abs_values_TP.append(TP_infected)
        gt_values.append(len(GT_snapshots[last]['infected']))
        abs_values_T.append(len(output_snapshots[k]['infected']))
        set_T.append(output_snapshots[k]['infected'])
        set_gt_T.append(GT_snapshots[last]['infected'])
        MCC.append((TP_infected*TN_ - FP_*FN_)/np.sqrt((TP_infected+FP_)*(TP_infected+FN_)*(TN_+FP_)*(TN_+FN_)))
        F1.append(2.0*precision_infected[-1]*recall_infected[-1]/(precision_infected[-1]+recall_infected[-1]))

    precision_recovered = []
    recall_recovered = []
    #MCC = []
    for k in sorted(GT_snapshots.keys()):
        TP_recovered = float(len(GT_snapshots[last]['recovered'] & output_snapshots[k]['recovered']))
        #TN = float(num_nodes - len(GT_snapshots[k] | output_snapshots[k]))
        #FP = float(len(output_snapshots[k] - GT_snapshots[k]))
        #FN = float(len(GT_snapshots[k] - output_snapshots[k]))
        if len(output_snapshots[k]['recovered']) == 0:
            precision_recovered.append(0.0)
        else:
            precision_recovered.append(TP_recovered/len(output_snapshots[k]['recovered']))

        if len(GT_snapshots[last]['recovered']) == 0:
            recall_recovered.append(0.0)
        else:
            recall_recovered.append(TP_recovered/len(GT_snapshots[last]['recovered']))
        #MCC.append((TP*TN - FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
    return precision_infected, recall_infected, precision_recovered, recall_recovered, abs_values_TP, gt_values, abs_values_T, set_T, set_gt_T, MCC, F1


def snapshot_accuracy_shifted(GT_snapshots, output_snapshots, pths, immuned, sources, reported, TS, G, num_nodes, mode = 'main', pred_recover = False, draw = False, folder = ''):

    precision_infected = []
    recall_infected = []
    abs_values_TP = []
    gt_values = []
    abs_values_T = []
    set_T = []
    set_gt_T = []
    MCC = []
    F1 = []
    #last = sorted(GT_snapshots.keys())[-1]
    sorted_keys = sorted(GT_snapshots.keys())
    for i in xrange(len(sorted_keys)-1):
        k = sorted_keys[i]
        next = sorted_keys[i+1]
        TP_infected = float(len(GT_snapshots[next]['infected'] & output_snapshots[k]['infected']))

        TN_ = float(num_nodes - len(GT_snapshots[next]['infected'] | output_snapshots[k]['infected']))
        FP_ = float(len(output_snapshots[k]['infected'] - GT_snapshots[next]['infected']))
        FN_ = float(len(GT_snapshots[next]['infected'] - output_snapshots[k]['infected']))
        if len(output_snapshots[k]['infected']) == 0:
            precision_infected.append(0.0)
        else:
            precision_infected.append(TP_infected/len(output_snapshots[k]['infected']))
        if len(GT_snapshots[next]['infected']) == 0:
            recall_infected.append(0.0)
        else:
            recall_infected.append(TP_infected/len(GT_snapshots[next]['infected']))

        abs_values_TP.append(TP_infected)
        gt_values.append(len(GT_snapshots[next]['infected']))
        abs_values_T.append(len(output_snapshots[k]['infected']))
        set_T.append(output_snapshots[k]['infected'])
        set_gt_T.append(GT_snapshots[next]['infected'])
        MCC.append((TP_infected*TN_ - FP_*FN_)/np.sqrt((TP_infected+FP_)*(TP_infected+FN_)*(TN_+FP_)*(TN_+FN_)))
        F1.append(2.0*precision_infected[-1]*recall_infected[-1]/(precision_infected[-1]+recall_infected[-1]))

    precision_recovered = []
    recall_recovered = []
    #MCC = []
    for k in sorted(GT_snapshots.keys()):
        TP_recovered = float(len(GT_snapshots[next]['recovered'] & output_snapshots[k]['recovered']))
        #TN = float(num_nodes - len(GT_snapshots[k] | output_snapshots[k]))
        #FP = float(len(output_snapshots[k] - GT_snapshots[k]))
        #FN = float(len(GT_snapshots[k] - output_snapshots[k]))
        if len(output_snapshots[k]['recovered']) == 0:
            precision_recovered.append(0.0)
        else:
            precision_recovered.append(TP_recovered/len(output_snapshots[k]['recovered']))

        if len(GT_snapshots[next]['recovered']) == 0:
            recall_recovered.append(0.0)
        else:
            recall_recovered.append(TP_recovered/len(GT_snapshots[next]['recovered']))
        #MCC.append((TP*TN - FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
    return precision_infected, recall_infected, precision_recovered, recall_recovered, abs_values_TP, gt_values, abs_values_T, set_T, set_gt_T, MCC, F1



def snapshot_acc_getsets(GT_snapshots, output_snapshots, G):
    k = sorted(GT_snapshots.keys())[-1]

    out_P = output_snapshots[k]['infected']
    out_TP = GT_snapshots[k]['infected'] & output_snapshots[k]['infected']
    out_TN = set(G.nodes()) - (GT_snapshots[k]['infected'] | output_snapshots[k]['infected'])
    out_FP = output_snapshots[k]['infected'] - GT_snapshots[k]['infected']
    out_FN = GT_snapshots[k]['infected'] - output_snapshots[k]['infected']

    return out_TP, out_TN, out_FP, out_FN


def snapshot_acc_barplot(GT_snapshots, output_snapshots, G):

    out_TP, out_TN, out_FP, out_FN = snapshot_acc_getsets(GT_snapshots, output_snapshots, G)

    degree_TP = {i: G.degree(i) for i in out_TP}
    degree_TN = {i: G.degree(i) for i in out_TN}
    degree_FP = {i: G.degree(i) for i in out_FP}
    degree_FN = {i: G.degree(i) for i in out_FN}

    return degree_TP, degree_TN, degree_FP, degree_FN

def snapshot_acc_barplot_time(GT_snapshots, output_snapshots, G, TS):
    sorted_snapshots = sorted(GT_snapshots.keys())
    snapshot_num = 1

    out_TP, out_TN, out_FP, out_FN = snapshot_acc_getsets(GT_snapshots, output_snapshots, G)

    FN_infection_moments = {}
    for itr in xrange(len(TS)):
        i = TS[itr]
        if i[0] > sorted_snapshots[snapshot_num]:
            snapshot_num += 1
        if i[1] in out_FN and i[3] == 1 and i[1] not in FN_infection_moments.keys():
            FN_infection_moments[i[1]] = snapshot_num
        if i[2] in out_FN and i[4] == 1 and i[2] not in FN_infection_moments.keys():
            FN_infection_moments[i[2]] = snapshot_num


    FN_degree = {i: G.degree(i) for i in out_FN}
    return FN_infection_moments, FN_degree


def accuracy_interactions(TS, paths, infected_interactions):
    count = 0.0
    output_interactions = set()
    print infected_interactions
    for src, covered in paths.iteritems():
        #print p, paths[p]
        for target, path in covered.iteritems():
            for p in path[1]:
                print p
                output_interactions.add(p)
                if p in infected_interactions:
                    count += 1.0
    print 'recall:', count/len(infected_interactions)
    print 'precision:',count/len(output_interactions)
    return count