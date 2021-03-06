import argparse
import json
from collections import defaultdict, namedtuple
from copy import deepcopy

FactorRow = namedtuple('FactorRow', 'groundVariables, value')

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", help="path to json representation of the CPD graph")
    parser.add_argument("--evidence", nargs="+",
            help="events given as evidence; formatted event=value")
    args = parser.parse_args()
    return args.graph, args.evidence

def parse_graph(graph):
    """
    `baggage` are the elements that will make up a new factor
    when summing out a given variable.
    """
    with open(graph, 'r') as f:
        g = json.load(f)
    full_factors = g['graph']['factors']
    factors = {f['name']:
        [FactorRow(r, r.pop('value')) for r
        in f['groundVariables']] for f in full_factors}
    scopes = defaultdict(set)
    reverse_scopes = defaultdict(set)
    for f in factors:
        for v in factors[f][0].groundVariables.keys(): # peek at vars in scope of first row
            scopes[f].add(v)
            reverse_scopes[v].add(f)
    baggage = {}
    for var in reverse_scopes:
        concomitant = set()
        for f in reverse_scopes[var]:
            concomitant = concomitant.union(scopes[f])
        baggage[var] = concomitant
    return factors, scopes, reverse_scopes, baggage

def build_cluster_graph(scopes, reverse_scopes, baggage):
    cluster_nodes = {} # maps cluster names to factors
    tau_to_cluster = {} # maps resulting product (tau) to cluster names
    cluster_edges = defaultdict(set) # maps cluster names to clusters they have an edge to
    counter = 0
    while len(baggage):
        # create cluster and resulting tau
        cluster_name = "C%s" % counter
        tau_name = "T%s" % counter
        tau_to_cluster[tau_name] = cluster_name
        # pick the easiest variable to eliminate
        to_eliminate = sorted(baggage, key=lambda x : len(baggage[x]))[0]
        # get all factors that contain that variable
        pertinent_factors = reverse_scopes[to_eliminate]
        # assign these (non-tau) factors to the cluster
        cluster_nodes[cluster_name] = pertinent_factors.difference(tau_to_cluster.keys())
        tau_scope = set()
        for f in pertinent_factors:
            # get the scope of tau by "multiplying" all pertinent factors together
            tau_scope = tau_scope.union(scopes[f])
            # if this factor is a tau, draw an edge from this cluster
            # to the cluster that created said tau
            if f in tau_to_cluster.keys():
                cluster_edges[cluster_name].add(tau_to_cluster[f])
                cluster_edges[tau_to_cluster[f]].add(cluster_name)
        # sum out the variable to_eliminate
        tau_scope = tau_scope.difference(to_eliminate)
        # add tau to remaining factors
        scopes[tau_name] = tau_scope
        for var in tau_scope:
            reverse_scopes[var].add(tau_name)
            # remove the now non-existent factors
            reverse_scopes[var] = reverse_scopes[var].difference(pertinent_factors)
            # recalibrate baggage
            baggage[var] = baggage[var].union(tau_scope)
            baggage[var] = baggage[var].difference([to_eliminate])
        # remove variable from baggage
        baggage.pop(to_eliminate)
        counter += 1
    # merge clusters that are sub/supersets of each other's scopes
    subsets = {}
    for c1 in list(cluster_nodes.keys()):
        scope1 = set().union(*[scopes[f] for f in cluster_nodes[c1]])
        for c2 in cluster_edges[c1]:
            scope2 = set().union(*[scopes[f] for f in cluster_nodes[c2]])
            if c1 != c2 and scope1.issubset(scope2):
                subsets[c1] = c2
                # push factors from c1 into its superset c2
                cluster_nodes[c2] = cluster_nodes[c2].union(cluster_nodes.pop(c1))
                # redirect all edges to/from c1 to c2 (except edge from c1 to c2)
                sub_edges = cluster_edges.pop(c1)
                cluster_edges[c2] = cluster_edges[c2].union(sub_edges).difference([c2])
                for c in sub_edges:
                    cluster_edges[c].remove(c1)
                    if c != c2: cluster_edges[c].add(c2)
                break
    return cluster_nodes, cluster_edges

def make_clusters(factors, cluster_nodes):
    clusters = {} # maps cluster names to resulting factors
    for c in cluster_nodes:
        product = [factors[f] for f in cluster_nodes[c]]
        while len(product) != 1:
            product.append(make_new_factor(product))
        clusters[c] = product.pop()
    return clusters

def propogate(clusters, cluster_edges):
    messages = {} # maps (from_cluster, to_cluster) to a factor
    upstream_clusters = set() # clusters that have already passed their message this sweep
    neighbor_count = {c: len(cluster_edges[c]) for c in clusters}
    ordering = []
    while len(neighbor_count): # pass messages to appropriate neighbors
        leaves = [c for c in neighbor_count if neighbor_count[c] == 1]
        if set(leaves) == set(clusters.keys()): # all (two) clusters are a leaf
            leaves = leaves[:-1]
        for c in leaves:
            ordering.append(c)
            c_scope = set(clusters[c][0].groundVariables.keys())
            upstream_neighbors = cluster_edges[c].intersection(upstream_clusters)
            # only one downstream_neighbor, since c is a leaf
            downstream_neighbor = cluster_edges[c].difference(upstream_clusters).pop()
            product = [deepcopy(messages[(n, c)]) for n in upstream_neighbors] # product of messages received
            product.append(deepcopy(clusters[c])) # psi (deepcopy VERY IMPORTANT, else modifies `clusters`)
            while len(product) != 1:
                product.append(make_new_factor(product))
            n_scope = clusters[downstream_neighbor][0].groundVariables.keys()
            f = product.pop()
            messages[(c, downstream_neighbor)] = sum_out(f, c_scope.difference(n_scope))
            upstream_clusters.add(c)
            if neighbor_count[downstream_neighbor] == 1:
                ordering.append(downstream_neighbor)
                neighbor_count.pop(downstream_neighbor)
            else:
                neighbor_count[downstream_neighbor] -= 1
            if neighbor_count[c] == 1:
                neighbor_count.pop(c)
            else:
                neighbor_count[c] -= 1
    # pass messages in opposite order
    upstream_clusters = set()
    for c in ordering[::-1]:
        c_scope = set(clusters[c][0].groundVariables.keys())
        downstream_neighbors = cluster_edges[c].difference(upstream_clusters)
        for n in downstream_neighbors:
            # product of messages received minus messages from the receiving cluster
            product = [deepcopy(messages[(n_, c)]) for n_ in downstream_neighbors if n_ != n]
            product.append(deepcopy(clusters[c])) # psi (deepcopy VERY IMPORTANT, else modifies `clusters`)
            while len(product) != 1:
                product.append(make_new_factor(product))
            n_scope = clusters[n][0].groundVariables.keys()
            f = product.pop()
            messages[(c, n)] = sum_out(f, c_scope.difference(n_scope))
        upstream_clusters.add(c)
    return messages

def get_beliefs(clusters, cluster_edges, messages):
    beliefs = {}
    for c in clusters:
        product = [messages[(n, c)] for n in cluster_edges[c]] + [clusters[c]]
        while len(product) > 1:
            product.append(make_new_factor(product))
        beliefs[c] = product.pop()
    return beliefs

def get_marginals(beliefs):
    result = {} # maps variables to an ordered sequence of probabilities
    for c in beliefs:
        for v in beliefs[c][0].groundVariables:
            if not v in result: # we don't have these marginals yet
                marginals = sum_out(deepcopy(beliefs[c]), set(beliefs[c][0].groundVariables.keys()).difference([v]))
                marginals = renormalize(marginals)
                probabilities = []
                for r in sorted(marginals, key = lambda r : r.groundVariables.keys()):
                    probabilities.append(r.value)
                result[v] = tuple(probabilities)
    return result

def sum_out(factor, vars_to_sum_out):
    """ sum out of a factor every var_to_sum_out """
    result = []
    for r in factor:
        foundAMatch = False
        for i in range(len(result)):
            matches = True
            for k in r.groundVariables.keys():
                if not k in vars_to_sum_out and r.groundVariables[k] != result[i].groundVariables[k]:
                    matches = False
            if matches:
                foundAMatch = True
                result[i] = FactorRow(result[i].groundVariables, r.value + result[i].value)
                break
        if not foundAMatch:
            for var in vars_to_sum_out:
                r.groundVariables.pop(var)
            result.append(r)
    return result

def parse_evidence(evidence):
    if evidence:
        split_evidence = [e.split("=") for e in evidence]
        return {i[0]: int(i[1]) for i in split_evidence}
    else:
        return None

def remove_noncompatible_evidence(factors, evidence):
    for f in factors:
        remove = []
        for i in range(len(factors[f])):
            for e in evidence:
                if e in factors[f][i].groundVariables.keys() \
                    and factors[f][i].groundVariables[e] != evidence[e]:
                    remove.append(i)
                    break
        for j in remove[::-1]:
            factors[f].pop(j)
    return factors

def make_new_factor(product):
    f1 = product.pop()
    f2 = product.pop()
    shared_vars = set(f1[0].groundVariables.keys()).intersection(f2[0].groundVariables.keys())
    new_factor = []
    for r in f1:
        for r2 in f2:
            matches = True
            for var in shared_vars:
                if r.groundVariables[var] != r2.groundVariables[var]:
                    matches = False
            if matches:
                new_factor.append(multiply_rows(r, r2))
    return new_factor

def renormalize(factor):
    partition_func = sum([r.value for r in factor])
    for i in range(len(factor)):
        factor[i] = FactorRow(factor[i].groundVariables, factor[i].value / partition_func)
    return factor

def multiply_rows(r, r2):
    result = FactorRow({**r.groundVariables, **r2.groundVariables}, r.value * r2.value)
    return result

def main():
    graph, evidence = read_args()
    factors, scopes, reverse_scopes, baggage = parse_graph(graph)
    evidence = parse_evidence(evidence)
    cluster_nodes, cluster_edges = build_cluster_graph(scopes, reverse_scopes, baggage)
    clusters = make_clusters(factors, cluster_nodes)
    messages = propogate(clusters, cluster_edges)
    beliefs = get_beliefs(clusters, cluster_edges, messages)
    marginals = get_marginals(beliefs)
    print(marginals)

if __name__ == "__main__":
    main()
