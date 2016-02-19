__author__ = 'Polina'
#from generator import *
from paths import *
from accuracy import *
#from trees import *
import uuid

def greedy(SP, N, K = 10, alpha = 1.0):
    norm = 0.0
    for src, sinks in SP.iteritems():
        for sink, info in sinks.iteritems():
            if np.isfinite(info[0]):
                norm += info[0]
    weights = []
    for src, sinks in SP.iteritems():
        for sink, info in sinks.iteritems():
            weights.append((info[0]/norm, src, sink, info[-1]))
    weights.sort()

    #weights = [sorted([(info[0]/norm, src, sink, info[-1]) for sink, info in sinks.iteritems()]) for src, sinks in SP.iteritems()]
    covered_nodes = set()
    seeds = set()
    output_paths = []
    cover = {}

    cover_cost = 0
    for i in weights:
        if len(seeds) < K or i[1] in seeds:
            seeds.add(i[1])
            output_paths.append(SP[i[1]][i[2]][-1])
            cover_cost += i[0]
            covered_nodes.add(i[2])
            cover[i[1]] = cover.get(i[1], []) + [i[2]]
            if len(covered_nodes) == N:
                break
    cover_cost += K*alpha
    return cover, output_paths, cover_cost

def greedy_plot_KvsAlpha(SP, N, TS, ticks = 30, ua = -1):
    norm = 0.0
    maxCost = 0.0
    for src, sinks in SP.iteritems():
        for sink, info in sinks.iteritems():
            if np.isfinite(info[0]):
                norm += info[0]
                if maxCost < info[0]:
                    maxCost = info[0]

    weights = {src: sorted([(info[0]/norm, sink, info[-1]) for sink, info in sinks.iteritems()]) for src, sinks in SP.iteritems()}
    #print weights.values()

    la = 0.0
    #ua = N*maxCost
    if ua < 0:
        ua = N*maxCost/norm

    #plt.plot(alphas, Ks)
    #plt.plot(alphas, [np.ln(100)*np.ln(10000)]*len(alphas))
    #plt.xlabel('alpha')
    ##plt.ylabel('K')
    #plt.title(str(np.log(100)*np.log(21000)))
    #plt.show()


    open_cost = {s: 1.0 for s in weights}

    output_paths = []

    step = ua/ticks

    iter = 0
    alpha = 0
    alphas = []
    Ks = []
    costs = []
    cost_RHT = []
    for i in xrange(ticks+1):
        iter += 1
        covered_nodes = set()
        actually_covered = set()
        cover = {}
        total_cost = 0.0
        while len(covered_nodes) < N:
        #for k in xrange(K, 0, -1):
            best_src, best_mgain, best_cost, best_covered_nodes = -1, np.Inf, np.Inf, set()
            for src, sinks in weights.iteritems():
                #if src not in cover: #with or without repetition
                    m_profit = 0.0
                    cost = open_cost[src] * alpha
                    local_best_mgain, local_best_cost, local_best_covered_nodes = np.Inf, np.Inf, set()
                    local_covered_nodes = set()
                    for s in sinks:
                        if s[1] not in covered_nodes:
                            m_profit += 1.0
                            cost += s[0]
                            local_covered_nodes.add(s[1])
                            if local_best_mgain > np.divide(cost, 1.0 * m_profit):
                                local_best_mgain = np.divide(cost, 1.0 * m_profit)
                                local_best_covered_nodes = copy.deepcopy(local_covered_nodes)
                                local_best_cost = cost
                    if best_mgain > local_best_mgain:
                        best_src = src
                        best_mgain = local_best_mgain
                        best_cost = local_best_cost
                        best_covered_nodes = copy.deepcopy(local_best_covered_nodes)

            covered_nodes.update(best_covered_nodes)
            #N -= len(best_covered_nodes)
            cover[best_src] = cover.get(best_src, []) + list(best_covered_nodes)
            print 'best cost', best_cost
            if best_cost != np.Inf:
                actually_covered.update(best_covered_nodes)
            total_cost += best_cost

        alphas.append(alpha)
        Ks.append(len(cover))
        costs.append(total_cost)
        cost_RHT.append(total_cost - len(cover)*alpha)
        alpha += step

    # unique_name = str(uuid.uuid4())
    # plt.figure()
    # plt.plot(alphas, Ks)
    # plt.xlabel('alpha')
    # plt.ylabel('K')
    # plt.savefig("K_vs_Alpha_" + unique_name)
    # #plt.show()
    #
    # plt.figure()
    # plt.plot(Ks, costs)
    # plt.xlabel('K')
    # plt.ylabel('cost')
    # plt.savefig("cost_" + unique_name)
    # #plt.show()
    #
    # plt.figure()
    # plt.plot(Ks, [cost_RHT[i] for i in xrange(len(cost_RHT))])
    # plt.xlabel('K')
    # plt.ylabel('cost, alpha=0')
    # #plt.savefig("cost, alpha=0")
    # plt.savefig("Alpha=0_" + unique_name)
    # #plt.show()
    #
    # plt.figure()
    # plt.plot(Ks, [cost_RHT[i] + Ks[i] for i in xrange(len(cost_RHT))])
    # plt.xlabel('K')
    # plt.ylabel('cost, alpha=1')
    # plt.savefig("Alpha=1_" + unique_name)
    # #plt.show()
    #
    # plt.figure()
    # plt.plot(Ks, [cost_RHT[i] + Ks[i]*np.log(len(weights))*np.log(len(TS)) for i in xrange(len(cost_RHT))])
    # plt.xlabel('K')
    # plt.ylabel('"AIC" cost')
    # #plt.show()
    # plt.savefig("AIC_" + unique_name)

    return Ks, alphas, cost_RHT

def greedyBS(SP, N, K = 10, la = 0.0, ua = 0.0):
    norm = 0.0
    maxCost = 0.0
    for src, sinks in SP.iteritems():
        for sink, info in sinks.iteritems():
            if np.isfinite(info[0]):
                norm += info[0]
                if maxCost < info[0]:
                    maxCost = info[0]

    weights = {src: sorted([(info[0]/norm, sink, info[-1]) for sink, info in sinks.iteritems()]) for src, sinks in SP.iteritems()}
    #print weights.values()

    #ua = K*N*maxCost/norm
    ua = K*N*maxCost/norm
    #ua = N*maxCost/norm

    open_cost = {s: 1.0 for s in weights}

    output_paths = []

    alpha = (la + ua)/2.0
    legal_cover = {}
    legal_alpha = -1
    #alpha = ua
    #alpha = 1.0
    #for iter in xrange(100):
    iter = 0
    while iter < 100 or (ua - la) > 1e-8:
        iter += 1
        covered_nodes = set()
        actually_covered = set()
        cover = {}
        total_cost = 0.0
        while len(covered_nodes) < N:
        #for k in xrange(K, 0, -1):
            best_src, best_mgain, best_cost, best_covered_nodes = -1, np.Inf, np.Inf, set()
            for src, sinks in weights.iteritems():
                #if src not in cover: #with or without repetition
                    m_profit = 0.0
                    cost = open_cost[src] * alpha
                    local_best_mgain, local_best_cost, local_best_covered_nodes = np.Inf, np.Inf, set()
                    local_covered_nodes = set()
                    for s in sinks:
                        if s[1] not in covered_nodes:
                            m_profit += 1.0
                            cost += s[0]
                            local_covered_nodes.add(s[1])
                            if local_best_mgain > np.divide(cost, 1.0 * m_profit):
                                local_best_mgain = np.divide(cost, 1.0 * m_profit)
                                local_best_covered_nodes = copy.deepcopy(local_covered_nodes)
                                local_best_cost = cost
                    if best_mgain > local_best_mgain:
                        best_src = src
                        best_mgain = local_best_mgain
                        best_cost = local_best_cost
                        best_covered_nodes = copy.deepcopy(local_best_covered_nodes)

            covered_nodes.update(best_covered_nodes)
            #N -= len(best_covered_nodes)
            cover[best_src] = cover.get(best_src, []) + list(best_covered_nodes)
            print 'best cost', best_cost
            if best_cost != np.Inf:
                actually_covered.update(best_covered_nodes)
            total_cost += best_cost

        print 'iteration: ', iter, 'sinks: ',N, len(actually_covered)
        print 'K: ', K, 'cover: ', len(cover), 'alpha: ', alpha
        if len(actually_covered) < N:
            ua = alpha
        else:
            if len(cover) > K:
                la = alpha
            elif len(cover) < K:
                ua = alpha
            else:
                ua = alpha
                legal_cover = cover
                legal_alpha = alpha

                break
        alpha = (ua + la)/2.0
        print 'alpha', la, ua, alpha

    if not legal_cover:
        legal_cover = cover
        legal_alpha = alpha
    print 'num of srcs', len(legal_cover)
    for src, sinks in legal_cover.iteritems():
        for s in sinks:
            output_paths.append(SP[src][s][-1])
    print total_cost
    print len(covered_nodes), len(actually_covered)
    return legal_cover, output_paths, total_cost, legal_alpha


def get_earliest(weights, choice):
    t, best_src = -1, -1
    for src, info in weights.iteritems():
        if src in choice and info[2]:
            if t == -1 or info[2][0][0] < t:
                t = info[2][0][0]
                best_src = src
    return best_src, t


def get_next_furthest(SP, sources, taken):
    t, dist, best_src = -1, -1, -1
    for new_src in SP:
        if taken and new_src not in taken:
            closest_dist, closest_t, closest_elem = -1, -1, -1
            for selected_src in taken:
                if (closest_dist == -1 or closest_dist > SP[selected_src][new_src][0]) or (closest_dist == SP[selected_src][new_src][0] and (closest_t == -1 or closest_t > sources[new_src])):
                    closest_dist = SP[selected_src][new_src][0]
                    #closest_elem = new_src
                    closest_t = sources[new_src]
            if dist < closest_dist or (dist == closest_dist and (t == -1 or t > closest_t)):
                dist = closest_dist
                best_src = new_src
                t = closest_t
        elif not taken:
            if t == -1 or t > sources[new_src]:
                best_src = new_src
                t = sources[new_src]

    return best_src, dist, t

def greedy_furthest_first(SP, sources, sinks, K = 1):
    taken = set()
    for i in xrange(K):
        best_src, dist, t = get_next_furthest(SP, sources, taken)
        print best_src, dist, t
        taken.add(best_src)

    cover, out_paths = assign_to_sources(SP, taken, sinks)
    return cover, out_paths

def assign_to_sources(SP, taken, sinks):
    cover = {s: set() for s in taken}
    out_paths = []
    for sink in sinks:
        best_src, best_dist = -1, np.Inf
        for src in taken:
            if SP[src][sink][0] < best_dist:
                best_src, best_dist = src, SP[src][sink][0]
        if best_src != -1:
            cover[best_src].add(sink)
            out_paths.append(SP[best_src][sink][-1])
        print best_dist
    return cover, out_paths

if __name__ == "__main__":

    p = 0.7
    type = 'ER'
    M = 10000
    N = 1000
    TS, snapshots, G = generateTS(n = N,
                                  p = 0.5,
                                  seed = 1.0,
                                  st = datetime(2000, 01, 01, 00, 00, 00),
                                  m = M,
                                  srcN = 10,
                                  reportingP = 0.7,
                                  infectionP = 0.9,
                                  recoveringP = 0.1,
                                  type = type)
    #TS, snapshots, G,_ = g.readFile('generated.txt')
    sources, immuned = get_sinks_and_sources(TS)
    print sources
    print immuned
    #exit()
    SP = shortestPath(TS, sources, sources, immuned)

    K = 10
    best_cost, best_root = get_best_root(SP)
    print best_cost, best_root
    #cover, output_paths = greedy(SP, len(sources), K)
    cover, output_paths = greedy_furthest_first(SP, sources, sources, K)
    print len(cover)

    #print cover, output_paths