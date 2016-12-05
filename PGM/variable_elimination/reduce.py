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
    return factors, scopes, reverse_scopes, baggage

def parse_evidence(evidence):
    if evidence:
        split_evidence = [e.split("=") for e in evidence]
        return {i[0]: int(i[1]) for i in split_evidence}
    else:
        return None

def calculate_marginal(event, factors, evidence, scopes, reverse_scopes, baggage):
    """
    We want to push out of the summation the variable with
    the least amount of baggage, since that will give us the
    smallest resulting factor after multiplying all the
    pre-existing factors in which that variable exists
    """
    if evidence:
        factors = remove_noncompatible_evidence(factors, evidence)
    baggage.pop(event)
    while len(baggage):
        event_to_sum_out = sorted(baggage, key=lambda x : len(baggage[x]))[0]
        pertinent_factors = reverse_scopes[event_to_sum_out]
        for k in reverse_scopes:
            reverse_scopes[k] = reverse_scopes[k].difference(pertinent_factors)
        product = [factors.pop(p) for p in pertinent_factors]
        counter = 0
        while len(product):
            if len(product) == 1: # sum out
                f = product.pop()
                tau = [] # factor table
                for r in f:
                    if not len(tau):
                        r.pop(event_to_sum_out)
                        tau.append(r)
                    else:
                        state = {k: r[k] for k in r if k != event_to_sum_out and k != 'value'}
                        foundAMatch = False
                        for r2 in tau:
                            matches = True
                            for k in state.keys():
                                if state[k] != r2[k]:
                                    matches = False
                            if matches:
                                r['value'] += r2['value']
                                foundAMatch = True
                                break
                        if not foundAMatch:
                            r.pop(event_to_sum_out)
                            tau.append(r)
                if len(baggage) == 1:
                    for f in reverse_scopes[event]:
                        if f in factors: # singelton potential for event exists
                            tau = make_new_factor([tau, factors[f]])
                    tau = renormalize(tau)
                    return tau
                new_factor_name = "T%s" % counter
                counter += 1
                factors[new_factor_name] = tau
                for var in tau[0].keys():
                    if var != 'value':
                        reverse_scopes[var].add(new_factor_name)
                        scopes[new_factor_name] = set()
                        scopes[new_factor_name].add(var)
                baggage.pop(event_to_sum_out)
            else: # factors left in product
                product.append(make_new_factor(product))

def make_new_factor(product):
    f1 = product.pop()
    f2 = product.pop()
    shared_vars = set(f1[0].keys()).intersection(f2[0].keys()).difference({'value'})
    new_factor = []
    for r in f1:
        for r2 in f2:
            matches = True
            for var in shared_vars:
                if r[var] != r2[var]:
                    matches = False
            if matches:
                new_factor.append(multiply_rows(r, r2))
    return new_factor

def remove_noncompatible_evidence(factors, evidence):
    for f in factors:
        remove = []
        for i in range(len(factors[f])):
            for e in evidence:
                if e in factors[f][i].keys() and factors[f][i][e] != evidence[e]:
                    remove.append(i)
                    break
        for j in remove[::-1]:
            factors[f].pop(j)
    return factors

def renormalize(factor):
    partition_func = sum([r['value'] for r in factor])
    for r in factor:
        r['value'] /= partition_func
    return factor

def multiply_rows(r, r2):
    value = r['value']
    value2 = r2['value']
    result = {k: r[k] for k in r.keys() if k != 'value'}
    for var in r2.keys():
        if var != 'value':
            result[var] = r2[var]
    result['value'] = value * value2
    return result

def main():
    graph, event, evidence = read_args()
    factors, scopes, reverse_scopes, baggage = parse_graph(graph)
    evidence = parse_evidence(evidence)
    marginal = calculate_marginal(event, factors, evidence, scopes, reverse_scopes, baggage)
    print(marginal)

if __name__ == "__main__":
    main()
