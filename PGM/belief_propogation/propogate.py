import argparse
import json
from collections import defaultdict, namedtuple

FactorRow = namedtuple('FactorRow', 'groundVariables, value')

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", help="path to json representation of the CPD graph")
    parser.add_argument("event",
        help="which event to calcualate marginal probability for")
    parser.add_argument("--evidence", nargs="+",
            help="events given as evidence; formatted event=value")
    args = parser.parse_args()
    return args.graph, args.event, args.evidence

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

def build_cluster_graph(event, scopes, reverse_scopes, baggage):
    cluster_nodes = {} # maps cluster names to factors
    tau_to_cluster = {} # maps cluster names to resulting product (tau)
    cluster_edges = defaultdict(set) # maps cluster names to clusters they have an edge to
    baggage.pop(event)
    counter = 0
    while len(baggage):
        event_to_sum_out = sorted(baggage, key=lambda x : len(baggage[x]))[0]
        pertinent_factors = reverse_scopes[event_to_sum_out]
        for k in reverse_scopes:
            reverse_scopes[k] = reverse_scopes[k].difference(pertinent_factors)
        cluster_name = "C%s" % counter
        tau_name = "T%s" % counter
        tau_to_cluster[tau_name] = cluster_name
        cluster_nodes[cluster_name] = pertinent_factors.difference(tau_to_cluster.keys())
        tau_scope = set()
        for f in pertinent_factors:
            tau_scope = tau_scope.union(scopes[f])
            if f in tau_to_cluster:
                cluster_edges[cluster_name].add(tau_to_cluster[f])
                cluster_edges[tau_to_cluster[f]].add(cluster_name)
        for v in tau_scope.difference([event_to_sum_out] + [event]):
            scopes[tau_name].add(v)
            reverse_scopes[v].add(tau_name)
            baggage[v] = baggage[v].union(tau_scope)
            baggage[v] = baggage[v].difference([event_to_sum_out])
        baggage.pop(event_to_sum_out)
        counter += 1
    return cluster_nodes, cluster_edges

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
    graph, event, evidence = read_args()
    factors, scopes, reverse_scopes, baggage = parse_graph(graph)
    evidence = parse_evidence(evidence)
    cluster_nodes, cluster_edges = build_cluster_graph(event, scopes, reverse_scopes, baggage)

if __name__ == "__main__":
    main()
