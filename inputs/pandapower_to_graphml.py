from enum import Enum

import pandapower as pp
from pandapower import networks as ppn
import simbench as sb
import networkx as nx
import argparse

class Sources(Enum):
    PANDA = 'Panda',
    SIMBENCH = 'Simbench',


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', help="case file name", nargs='*')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-p', '--pandacase', help='pandapower case (default)', dest='source', action='store_const', const=Sources.PANDA)
    group.add_argument('-s', '--simbenchnet', help='simbench case', dest='source', action='store_const', const=Sources.SIMBENCH)
    parser.set_defaults(source=Sources.PANDA)
    parser.add_argument('-l', '--list', help="list available case files", action='store_true')
    parser.add_argument('-a', '--all', help="use all case files in the specified set (pandapower or simbench)", action='store_true')
    parser.add_argument('-o', '--output', help="output file name (default: {case name}.graphml")
    parser.add_argument('--stats', help="print case statistics", action='store_true')
    args = parser.parse_args()

    source = args.source

    if args.list:
        if source == Sources.PANDA:
            print('\n'.join(net for net in ppn.__dict__ if net.startswith('case')))
        elif source == Sources.SIMBENCH:
            print('\n'.join(sb.collect_all_simbench_codes()))
        else:
            raise Exception('invalid source')
    else:
        if args.all:
            if source == Sources.PANDA:
                cases = [net for net in ppn.__dict__ if net.startswith('case')]
            elif source == Sources.SIMBENCH:
                cases = sb.collect_all_simbench_codes()
            else:
                raise Exception('invalid source')
        else:
            cases = args.case

        for case in cases:
            if args.output is None:
                outname = f'{case}.graphml'
            else:
                outname = args.output

            if source == Sources.PANDA:
                net = ppn.__dict__[case]()
            elif source == Sources.SIMBENCH:
                net = sb.get_simbench_net(case)
            else:
                raise Exception(f'unknown source: "{source}"')

            if args.stats:
                print(case)
                print(net)
            else:
                graph = pp.topology.create_nxgraph(net, multi=False)
                loads = {i for i in net.load['bus']}
                generators = {i for i in net.gen['bus']}
                static_generators = {i for i in net.sgen['bus']}

                load_gen = loads.union(generators).union(static_generators)

                outgraph = nx.Graph()
                outgraph.add_nodes_from(graph.nodes)
                outgraph.add_edges_from(graph.edges)
                # mark zero injection
                for n in outgraph:
                    outgraph.nodes[n]['zero_injection'] = 1 * (n not in load_gen)

                nx.write_graphml(outgraph, outname)
