import os
import sys
import uuid

from utils.generator_noise import *
from utils.visualization import *
from utils.paths import *
from utils.postprocessing import *
from utils.sinks_sources import *
from utils.greedyKcover import *
from utils.get_path_stats import *
import pickle

srcN = 9
#reportingP = float(sys.argv[1])
reportingP = 0.5
K = srcN
real = True

#filepath = os.path.join('..', 'Data', 'lipstick_empty.txt')
dataset = 'flixter_10000_7K'
#dataset = 'flixter_1000_9K'
filepath = os.path.join('..', 'Data', dataset+'.txt')

TS, G, infection_order, infection_track, seeds = readFile(filepath, mode = 'general')

print len(TS), len(seeds), len(infection_order), len(infection_track)

dt_sec = 100
sources, immuned, sinks, reported, unreported, sources_TnI, sinks_TnI, unreported_TnI = get_sinks_and_sources_shifted(
    TS, G=nx.Graph(), mode='all', dt_sec=dt_sec, rep_prob=reportingP)

upstream = get_upstream(TS, sinks)

SP = shortestPath1(TS, sources, sinks, immuned, unreported)

cover, output_paths, cover_cost, legal_alpha = greedyBS(SP, len(sinks), K)
gt_cost, gt_interactions, gt_causality = get_GT_cost(TS, sources, sinks, unreported)
print 'cost of GT', gt_cost
print cover
print output_paths

out_cost, out_interactions, out_causality = get_out_cost(output_paths, sources, sinks, unreported)
print 'cost of our solution', out_cost
#exit()

print 'correct interactions', len(set(gt_interactions).intersection(set(out_interactions)))
print 'intesection precision', np.divide(1.0*len(set(gt_interactions).intersection(set(out_interactions))), len(set(gt_interactions)))
print 'intesection recall', np.divide(1.0*len(set(gt_interactions).intersection(set(out_interactions))), len(set(out_interactions)))

print 'correct of causality', len(set(gt_causality).intersection(set(out_causality)))
print 'causality precision', np.divide(1.0*len(set(gt_causality).intersection(set(out_causality))), len(set(gt_causality)))
print 'causality recall', np.divide(1.0*len(set(gt_causality).intersection(set(out_causality))), len(set(out_causality)))
#exit()
#get_infection_paths_noise(TS, output_paths, cover.keys())
N = G.number_of_nodes()
M = len(TS)

ticksN = 100
GT_snapshots, moment_of_infection, _, _ = get_snapshots(TS, ticksN, N, upstream_ind = upstream)

print 'accuracy'
draw = False
ticksN = 100

lb_snapshots = get_lb_snapshots(sources, immuned, reported, GT_snapshots)
#ub_snapshots_cascade = get_ub_snapshots_cascade(TS, GT_snapshots, sinks)
ub_snapshots_cascade, ub_interactions, ub_causality= get_ub_snapshots(TS, GT_snapshots, sinks)

print 'intesection precision ub', np.divide(1.0*len(set(gt_interactions).intersection(set(ub_interactions))),len(set(gt_interactions)))
print 'intesection recall ub', np.divide(1.0*len(set(gt_interactions).intersection(set(ub_interactions))),len(set(ub_interactions)))

print 'causality precision ub', np.divide(1.0*len(set(gt_causality).intersection(set(ub_causality))),len(set(gt_causality)))
print 'causality recall ub', np.divide(1.0*len(set(gt_causality).intersection(set(ub_causality))),len(set(ub_causality)))


output_snapshots, found_seeds = get_output_snapshots_no_recov_pred(output_paths, GT_snapshots, immuned)
print set(found_seeds.keys()).intersection(set(seeds))

print 'how far in time'
slack, tau = how_far_intime(output_paths, moment_of_infection)
print 'total tau', tau
print 'slack', np.nanmean(slack), stats.nanmedian(slack)

folder = ''
pred_recover = False
prec_infected, recall_infected, prec_recovered, recall_recovered, abs_values_tp, gt_positive, abs_values_p, set_nodes, set_nodes_gt, MCC, F1, MCC_upstr \
    = snapshot_accuracy_uppstr(GT_snapshots, output_snapshots, output_paths, immuned, sources, reported, TS, G, N, 'main', pred_recover = pred_recover,draw = draw, folder = folder)
prec_lb_infected, recall_lb_infected, prec_lb_recovered, recall_lb_recovered, abs_values_lb_tp, _,_,_,_,MCC_lb,F1_lb, MCC_lb_upstr \
    = snapshot_accuracy_uppstr(GT_snapshots, lb_snapshots, output_paths, immuned, sources, reported, TS, G, N,'lb', pred_recover = pred_recover, draw = draw, folder = folder)
prec_ub_infected, recall_ub_infected, prec_ub_recovered, recall_ub_recovered, abs_values_ub_tp,_,abs_values_ub_p,_,_,MCC_ub,F1_ub, MCC_ub_upstr \
    = snapshot_accuracy_uppstr(GT_snapshots, ub_snapshots_cascade, output_paths, immuned, sources, reported, TS, G, N, 'ub', pred_recover = pred_recover, draw = draw, folder = folder)
degrees_out = [G.degree(i) for i in set_nodes[-1]]
degrees_gt = [G.degree(i) for i in set_nodes_gt[-1]]

print stats.nanmedian(prec_infected), stats.nanmedian(recall_infected), stats.nanmedian(F1)
print stats.nanmedian(prec_lb_infected), stats.nanmedian(recall_lb_infected), stats.nanmedian(F1_lb)
print stats.nanmedian(prec_ub_infected), stats.nanmedian(recall_ub_infected), stats.nanmedian(F1_ub)

unique_out = str(reportingP).replace('.', '-')+'_'+ dataset + '_' + str(uuid.uuid4())
pickle.dump([abs_values_tp, abs_values_p, abs_values_lb_tp, gt_positive, abs_values_ub_tp, abs_values_ub_p,
             prec_infected, prec_lb_infected, prec_ub_infected,
             recall_infected, recall_lb_infected, recall_ub_infected,
             MCC, MCC_lb, MCC_ub,
             F1, F1_lb, F1_ub,
             MCC_upstr, MCC_lb_upstr, MCC_ub_upstr
             ], open(unique_out+".p", "wb"))

title = 'MCC'
plt.figure('MCC')
plt.plot(MCC)
plt.plot(MCC_lb)
plt.plot(MCC_ub, 'k:')

plt.xlabel('snapshots')
plt.ylim(ymax = 1.01, ymin = -0.1)
plt.legend(['MCC', 'MCC reports', 'MCC BL'], loc = 3)
#plt.title(title)
name = dataset + '_mcc'+'.pdf'
plt.tight_layout()
name = plt.savefig(name)

title = 'MCC_upstream'
plt.figure('MCC_upstream')
plt.plot(MCC_upstr, 'k-')
plt.plot(MCC_lb_upstr, 'k--')
plt.plot(MCC_ub_upstr, 'k:')

plt.xlabel('snapshots')
plt.ylim(ymax = 1.01, ymin = -0.1)
plt.legend(['MCC', 'MCC reports', 'MCC BL'], loc = 3)
#plt.title(title)
name = dataset + '_mcc_upstr'+'.pdf'
plt.tight_layout()
name = plt.savefig(name)
