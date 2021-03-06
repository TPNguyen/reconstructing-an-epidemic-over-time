__author__ = 'Polina'
import copy
import numpy as np
import operator
import scipy.stats as stats

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
