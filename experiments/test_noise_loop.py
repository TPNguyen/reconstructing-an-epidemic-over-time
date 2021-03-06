import os
import sys
import uuid

#from utils.generator import *
from utils.generator_noise import *
#from utils.generator_noise_misc import *
from utils.sinks_sources import *
from utils.greedyKcover import *
from utils.visualization import *
from utils.paths import *
from utils.postprocessing import *
from utils.get_path_stats import *
import pickle
from utils import powerlaw
from collections import Counter

M = 30
N = 100
srcN = 5
infectionP = 0.3
reportingP = 0.5
recoveringP = 0.0
K = srcN
#dataset = str(sys.argv[1])
#model = str(sys.argv[2])
#parameter = float(sys.argv[3])
#type = 'ER'
dataset = 'facebook' #1.32099106109 #3.05636370627
dataset = 'twitter' #0.967784678046 #5.72835362312
dataset = 'students' #0.778465105778 #2.01248419318
dataset = 'tumblr' #1.22235609767 #2.37233199614
dataset = 'enron' #0.843106547984 # 2.37358302613
model = 'SI'
parameter = 0.0
p = 4.0/N
real = False

cost_ratio = []
correct_interaction_recall = []
correct_interaction_precision = []
correct_causality_recall = []
correct_causality_precision = []
path_tau_median = []
path_len_median = []
total_order_tau = []
mean_slack = []

prec_out = []
recall_out = []
F_out = []
mcc_out = []

prec_out_05 = []
recall_out_05 = []
F_out_05 = []

prec_lb = []
recall_lb = []
F_lb = []
mcc_lb = []

prec_lb_05 = []
recall_lb_05 = []
F_lb_05 = []

prec_ub = []
recall_ub = []
F_ub = []
mcc_ub = []

prec_ub_05 = []
recall_ub_05 = []
F_ub_05 = []

correct_causality_F = []
correct_interaction_F = []
ub_correct_causality_F = []
ub_correct_interaction_F = []

TS_len = []
inf_len = []

correct_found_seeds = []

fitted_gamma = []
#noise_range = range(0, 501, 10)
noise_range = range(0, 501, 5)
#for noise in noise_range:
for noise in xrange(100):
    #TS, G, infection_order, infection_track, seeds, nodes = get_noisy_TS_not_only_leaves(n = N, inf = 0.5, seeds_num = 5, noise = noise)
    TS, G, infection_order, infection_track, seeds, nodes = get_noisy_TS_leaves_with_P(type = dataset, n = N, inf_fraction = 0.5, seeds_num = srcN, noise = noise, reportP = 1.0, model = model, parameter = parameter)
    fitted_pl = powerlaw.Fit(G.degree().values())
    #degs = {k: v for k, v in Counter(G.degree().values()).iteritems() if k > 0}
    #gamma, b = np.polyfit(np.log(degs.keys()), np.log(degs.values()), 1.0)
    fitted_gamma.append(fitted_pl.alpha - 1.0)
    print fitted_pl.alpha
    #fitted_gamma.append(-gamma)
    #print -gamma
    #exit()
    continue

    K = len(seeds)
    sources, immuned, sinks, reported, unreported, sources_TnI, sinks_TnI, unreported_TnI = get_sinks_and_sources_noise(TS, mode = 'all')

    SP = shortestPath1(TS, sources, sinks, immuned, unreported)
    print len(sources), len(sinks)

    cover, output_paths, cover_cost, legal_alpha = greedyBS(SP, len(sinks), K)

    gt_cost, gt_interactions, gt_causality = get_GT_cost(TS, sources, sinks, unreported)
    print 'cost of GT', gt_cost
    print cover
    print output_paths

    out_cost, out_interactions, out_causality = get_out_cost(output_paths, sources, sinks, unreported)
    print 'cost of our solution', out_cost
    cost_ratio.append(1.0*out_cost/gt_cost)

    intesection_recall = np.divide(1.0*len(set(gt_interactions).intersection(set(out_interactions))), len(set(gt_interactions)))
    interaction_precision = np.divide(1.0*len(set(gt_interactions).intersection(set(out_interactions))), len(set(out_interactions)))

    correct_interaction_recall.append(intesection_recall)
    correct_interaction_precision.append(interaction_precision)
    correct_interaction_F.append(np.divide(2.0*intesection_recall*interaction_precision, (interaction_precision + intesection_recall)))

    causality_recall = np.divide(1.0*len(set(gt_causality).intersection(set(out_causality))),len(set(gt_causality)))
    causality_precision = np.divide(1.0*len(set(gt_causality).intersection(set(out_causality))),len(set(out_causality)))

    correct_causality_recall.append(causality_recall)
    correct_causality_precision.append(causality_precision)
    correct_causality_F.append(np.divide(2.0*causality_precision*causality_recall,(causality_precision + causality_recall)))

    #    print 'correct of causality', len(set(gt_causality).intersection(set(out_causality)))

    #get_infection_paths_noise(TS, output_paths, cover.keys())
    N = G.number_of_nodes()
    M = len(TS)

    ticksN = 100
    GT_snapshots, moment_of_infection, _, _ = get_snapshots(TS, ticksN, N)

    folder = time.strftime('test_' + dataset + "_%Y%m%d-%H%M%S")
    #os.mkdir(folder)

    #added_infections, counters, gt_uninf_neighbors, out_uninf_neighbors, gt_inf_neighbors, out_inf_neighbors = postPr(TS, output_paths, 0.3)


    print 'accuracy'
    draw = False
    ticksN = 100

    #GT_snapshots,_,moment_of_infection, _, _ = get_snapshots(TS, ticksN, N)
    lb_snapshots = get_lb_snapshots(sources, immuned, reported, GT_snapshots)
    ub_snapshots, ub_interactions, ub_causality = get_ub_snapshots(TS, GT_snapshots, sinks)
    output_snapshots, found_seeds = get_output_snapshots_no_recov_pred(output_paths, GT_snapshots, immuned)

    ub_interaction_recall = np.divide(1.0*len(set(gt_interactions).intersection(set(ub_interactions))), len(set(gt_interactions)))
    ub_interaction_precision = np.divide(1.0*len(set(gt_interactions).intersection(set(ub_interactions))), len(set(ub_interactions)))
    ub_correct_interaction_F.append(np.divide(2.0*ub_interaction_recall*ub_interaction_precision, (ub_interaction_precision + ub_interaction_recall)))

    ub_causality_recall = np.divide(1.0*len(set(gt_causality).intersection(set(ub_causality))), len(set(gt_causality)))
    ub_causality_precision = np.divide(1.0*len(set(gt_causality).intersection(set(ub_causality))), len(set(ub_causality)))
    ub_correct_causality_F.append(np.divide(2.0*ub_causality_precision*ub_causality_recall, (ub_causality_precision + ub_causality_recall)))

    print len(set(gt_interactions).intersection(set(ub_interactions)))
    print set(gt_causality).intersection(set(ub_causality))


    print gt_interactions
    print ub_interactions
    print gt_causality
    print ub_causality

    print ub_correct_interaction_F,ub_correct_causality_F
    print ub_interaction_precision, ub_interaction_recall
    print ub_causality_precision, ub_causality_recall

    correct_found_seeds.append(len(set(found_seeds) & set(seeds)))

    print 'how far in time'
    slack, tau = how_far_intime(output_paths, moment_of_infection)
    print 'total tau', tau
    total_order_tau.append(tau[0])
    print 'slack', np.nanmean(slack), stats.nanmedian(slack)
    mean_slack.append(np.nanmean(slack))


    pred_recover = False
    prec_infected, recall_infected, prec_recovered, recall_recovered, abs_values_tp, gt_positive, abs_values_p, set_nodes, set_nodes_gt, MCC, F1 \
        = snapshot_accuracy(GT_snapshots, output_snapshots, output_paths, immuned, sources, reported, TS, G, N, 'main', pred_recover = pred_recover,draw = draw, folder = folder)
    prec_lb_infected, recall_lb_infected, prec_lb_recovered, recall_lb_recovered, abs_values_lb_tp, _,_,_,_,MCC_lb,F1_lb \
        = snapshot_accuracy(GT_snapshots, lb_snapshots, output_paths, immuned, sources, reported, TS, G, N,'lb', pred_recover = pred_recover, draw = draw, folder = folder)
    prec_ub_infected, recall_ub_infected, prec_ub_recovered, recall_ub_recovered, abs_values_ub_tp,_,abs_values_ub_p,_,_,MCC_ub,F1_ub \
        = snapshot_accuracy(GT_snapshots, ub_snapshots, output_paths, immuned, sources, reported, TS, G, N, 'ub', pred_recover = pred_recover, draw = draw, folder = folder)
    degrees_out = [G.degree(i) for i in set_nodes[-1]]
    degrees_gt = [G.degree(i) for i in set_nodes_gt[-1]]

    print stats.nanmedian(prec_infected), stats.nanmedian(recall_infected), stats.nanmedian(F1)
    print stats.nanmedian(prec_lb_infected), stats.nanmedian(recall_lb_infected), stats.nanmedian(F1_lb)
    print stats.nanmedian(prec_ub_infected), stats.nanmedian(recall_ub_infected), stats.nanmedian(F1_ub)
    prec_out.append(prec_infected[-1])
    recall_out.append(recall_infected[-1])
    F_out.append(F1[-1])
    mcc_out.append(MCC[-1])

    len_05 = len(prec_infected)/2

    prec_out_05.append(prec_infected[len_05])
    recall_out_05.append(recall_infected[len_05])
    F_out_05.append(F1[len_05])

    prec_lb.append(prec_lb_infected[-1])
    recall_lb.append(recall_lb_infected[-1])
    F_lb.append(F1_lb[-1])
    mcc_lb.append(MCC_lb[-1])

    prec_lb_05.append(prec_lb_infected[len_05])
    recall_lb_05.append(recall_lb_infected[len_05])
    F_lb_05.append(F1_lb[len_05])

    prec_ub.append(prec_ub_infected[-1])
    recall_ub.append(recall_ub_infected[-1])
    F_ub.append(F1_ub[-1])
    mcc_ub.append(MCC_ub[-1])

    prec_ub_05.append(prec_ub_infected[len_05])
    recall_ub_05.append(recall_ub_infected[len_05])
    F_ub_05.append(F1_ub[len_05])
    
    TS_len.append(len(TS))
    inf_len.append(len(infection_order))
print 'avg', np.mean(fitted_gamma)

exit()
unique_out = dataset + '_' + model + '_' + str(parameter).replace('.','-') + str(uuid.uuid4())
pickle.dump([noise_range, correct_found_seeds, cost_ratio, correct_interaction_recall, correct_interaction_precision,
              correct_causality_recall, correct_causality_precision, total_order_tau, mean_slack,
              prec_out, prec_lb, prec_ub, recall_out, recall_lb, recall_ub,
              prec_out_05, prec_lb_05, prec_ub_05, recall_out_05, recall_lb_05, recall_ub_05,
              F_out, F_lb, F_ub,
              correct_causality_F, ub_correct_causality_F, correct_interaction_F, ub_correct_interaction_F,
              mcc_out, mcc_lb, mcc_ub, len(TS_len), len(inf_len)], open(unique_out+".p", "wb"))
exit()

plt.figure('seeds')
plt.xlabel('number of correct seeds')
plt.plot(noise_range, correct_found_seeds)
plt.xlabel('noise')
plt.title('found seeds')
#plt.show()
name = plt.savefig(dataset+'_seeds.pdf')

plt.figure('cost ratio')
plt.title('cost of found solution / cost of GT solution')
plt.xlabel('noise')
plt.plot(noise_range, cost_ratio)

#plt.show()
name = plt.savefig(dataset+'_cost_ratio.pdf')

plt.figure('correctness_order')
plt.title('correctness_order')
plt.plot(noise_range, correct_interaction_recall)
plt.plot(noise_range, correct_interaction_precision)
plt.plot(noise_range, correct_causality_recall)
plt.plot(noise_range, correct_causality_precision)
plt.plot(noise_range, total_order_tau)
plt.legend(['recall, interactions ', 'precision, interactions', 'recall, order', 'precision, order', 'total order tau'], loc = 3)
plt.xlabel('noise')
plt.tight_layout()
name = plt.savefig(dataset+'_correctness_order.pdf')

# plt.figure('path_tau_median')
# plt.plot(noise_range, path_tau_median)
# plt.tight_layout()
# name = plt.savefig('path_tau_median.pdf')
#
# plt.figure('path_len_median')
# plt.plot(noise_range, path_len_median)
# plt.tight_layout()
# name = plt.savefig('path_len_median.pdf')

plt.figure('mean_slack')
plt.title('mean slack')
plt.plot(noise_range, mean_slack)
plt.ylabel('interactions')
plt.xlabel('noise')
plt.tight_layout()
name = plt.savefig(dataset+'_mean_slack.pdf')

plt.figure('accuracy_end')
plt.title('accuracy, by the end')
plt.plot(noise_range, prec_out, 'k-')
plt.plot(noise_range, prec_lb, 'k--')
plt.plot(noise_range, prec_ub, 'k:')
plt.plot(noise_range, recall_out, 'r-')
plt.plot(noise_range, recall_lb, 'r--')
plt.plot(noise_range, recall_ub, 'r:')
plt.legend(['Prec output', 'Prec reports', 'Prec BL', 'Recall output', 'Recall reports', 'Recall BL'], loc = 3)
#plt.ylabel('F1')
plt.xlabel('noise')
plt.tight_layout()

name = plt.savefig(dataset+'_accuracy.pdf')

plt.figure('accuracy_mid')
plt.title('accuracy, middle moment')
plt.plot(noise_range, prec_out_05, 'k-')
plt.plot(noise_range, prec_lb_05, 'k--')
plt.plot(noise_range, prec_ub_05, 'k:')
plt.plot(noise_range, recall_out_05, 'r-')
plt.plot(noise_range, recall_lb_05, 'r--')
plt.plot(noise_range, recall_ub_05, 'r:')
plt.legend(['Prec output', 'Prec repotrs', 'Prec BL', 'Recall output', 'Recall reports', 'Recall BL'], loc = 3)
#plt.ylabel('F1')
plt.xlabel('noise')

plt.tight_layout()
name = plt.savefig(dataset+'_accuracy_05.pdf')

print 'here_05:'
print prec_out_05
print recall_out_05

plt.figure('F1_end')
plt.title('F1, by the end')
plt.plot(noise_range, F_out)
plt.plot(noise_range, F_lb)
plt.plot(noise_range, F_ub)
plt.legend(['F1', 'F1 reports', 'F1 BL'], loc = 3)
plt.ylabel('F1')
plt.xlabel('noise')
plt.tight_layout()
name = plt.savefig(dataset+'_F1.pdf')

plt.figure('order_F1')
plt.title('compare order, F1')
plt.plot(noise_range, correct_causality_F)
plt.plot(noise_range, ub_correct_causality_F)
plt.plot(noise_range, correct_interaction_F)
plt.plot(noise_range, ub_correct_interaction_F)
plt.legend(['order output', 'order BL', 'interactions output', 'interactions BL'], loc = 3)
plt.ylabel('F1')
plt.xlabel('noise')
plt.tight_layout()
name = plt.savefig(dataset+'_compare order.pdf')

#print MCC_out, MCC_lb, MCC_ub
plt.figure('MCC_end')
plt.title('MCC, by the end')
plt.plot(noise_range, mcc_out)
plt.plot(noise_range, mcc_lb)
plt.plot(noise_range, mcc_ub)
plt.legend(['MCC1', 'MCC1 reports', 'MCC1 BL'], loc = 3)
plt.ylabel('MCC')
plt.xlabel('noise')
plt.tight_layout()
name = plt.savefig(dataset+'_MCC1.pdf')
plt.show()
