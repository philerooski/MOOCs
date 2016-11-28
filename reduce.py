import argparse
import json
from collections import defaultdict

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
    with open(graph, 'r') as f:
        g = json.load(f)
    full_factors = g['graph']['factors']
    factors = {f['name']: f['groundVariables'] for f in full_factors}
    scopes = {k: set(factors[k][0].keys()).difference({"value"}) for k in factors}
    reverse_scopes = defaultdict(set)
    for k in factors:
        for v in factors[k][0]: # peek at vars in scope of first row
            if v != "value":
                reverse_scopes[v].add(k)
    baggage = {}
    for var in reverse_scopes:
        concomitant = set()
        for f in reverse_scopes[var]:
            concomitant = concomitant.union(scopes[f])
        baggage[var] = concomitant
    print("baggage", baggage)
    return factors, scopes, reverse_scopes, baggage

def parse_evidence(evidence):
    split_evidence = [e.split("=") for e in evidence]
    return {i[0]: int(i[1]) for i in split_evidence}

def calculate_marginal(event, factors, evidence, scopes, reverse_scopes, baggage):
    """
    We want to push out of the summation the variable with
    the least amount of baggage, since that will give us the
    smallest resulting factor after multiplying all the
    pre-existing factors in which that variable exists
    """
    #while len(factors) > 1:
    event_to_sum_out = sorted(baggage, key=lambda x : len(baggage[x]))[0]

def main():
    graph, event, evidence = read_args()
    factors, scopes, reverse_scopes, baggage = parse_graph(graph)
    evidence = parse_evidence(evidence)
    calculate_marginal(event, factors, evidence, scopes, reverse_scopes, baggage)

if __name__ == "__main__":
    main()
