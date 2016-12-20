import argparse
import json
from collections import defaultdict, namedtuple

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
    ordering = sorted(cluster_edges, key=lambda c : len(cluster_edges[c]))
    messages = {} # maps (from_cluster, to_cluster) to a factor
    upstream_clusters = set() # clusters that have already passed their message this sweep
    for c in ordering: # pass messages to appropriate neighbors
        print("calculating message of", c)
        c_scope = set(clusters[c][0].groundVariables.keys())
        print("c_scope", c_scope)
        upstream_neighbors = cluster_edges[c].intersection(upstream_clusters)
        print("upstream_neighbors", upstream_neighbors)
        downstream_neighbors = cluster_edges[c].difference(upstream_clusters)
        print("downstream_neighbors", downstream_neighbors)
        product = [messages[(n, c)] for n in upstream_neighbors]
        product.append(clusters[c])
        print("product", product)
        while len(product) != 1:
            product.append(make_new_factor(product))
        for n in downstream_neighbors:
            n_scope = clusters[n][0].groundVariables.keys()
            print("n_scope", n_scope)
            sepset = c_scope.intersection(n_scope) # sum out everything but these
            print('sepset', sepset)
            f = product.pop()
            messages[(c, n)] = sum_out(f, sepset)
            print('messages[(c, n)]', messages[(c, n)])
        upstream_clusters.add(c)

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
    print("clique tree", cluster_nodes, cluster_edges, "\n")
    #clusters = make_clusters(factors, cluster_nodes)
    #propogate(clusters, cluster_edges)

if __name__ == "__main__":
    main()
