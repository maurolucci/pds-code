"""Este módulo se ocupa de generar instancias con el módulo matpower.

La idea es generar algunas de las instancias de la literatura que no estan disponibles en el 
módulo pandapower. El procedimiento consiste en generar la instancia con matpower, exportarla 
a formato mpc (con extensión ".m"), levantarla con pandapower, exportarla a NetworkX y 
escribirla.

Hacer funcionar matpower no es sencillo. Si bien se puede instalar desde pip, corre sobre 
MATLAB o GNU Octave, por lo que también se necesitan tener instalado alguno de dichos softwares. 
Para hacerlo desde GNU Octave, también son necesarios algunos módulos adicionales para Python: 
oct2py y matpower[octave].

Personalmente, también tuve que cambiar la línea 46 de matpower/__init__.py a:
m.install_matpower(1,1,0,1).

Referencias:
* Manual de matpower: https://matpower.org/docs/MATPOWER-manual.pdf
"""

import os
from matpower import start_instance
import networkx as nx
import pandapower.topology as top
import pandapower.converter as pc


def generar_redes() -> list[str]:
    """Genera un listado con algunas redes de interés.

    Retorno:
    list[(str, str)]
        Una lista de pares con el nombre de la instancia y su nombre correspondiente
        dentro del paquete matpower.
    """

    return [
        ("case73rts", "case_RTS_GMLC"),
        ("case13659pegase", "case13659pegase"),
        ("case1951rte", "case1951rte"),
        ("case2868rte", "case2868rte"),
        ("case6468rte", "case6468rte"),
        ("case2383wp", "case2383wp"),
        ("case2736sp", "case2736sp"),
        ("case2737sop", "case2737sop"),
        ("case2746wop", "case2746wop"),
        ("case2746wp", "case2746wp"),
        ("case3012wp", "case3012wp"),
        ("case3375wp", "case3375wp"),
        ("case_south_carolina500", "case_ACTIVSg500"),
        # ("texas_2000", "case_ACTIVSg2000"),
        # ("western_10k", "case_ACTIVSg10k"),
        # ("northeast_25k", "case_ACTIVSg25k"),
        # ("eastern_70k", "case_ACTIVSg70k"),
        # ("usa_82k", "case_SyntheticUSA"),
    ]

def generar_gml(dir_instancias: str):
    """Genera algunas instancias provenientes de matpower.

    Argumentos:
    dir_instancias : str
        Ruta con el directorio donde escribir las instancias.

    La red de matpower de convierte a una red de pandapower con la función from_mpc.
    La red de pandapower se convirte a un grafo con la funcion create_nxgraph.
    Se ignoran aristas paralelas. Los grafos se escriben en formato GraphML.
    Se marcan los vértices pre-seleccionados, excluidos y propagativos.
    """
    for nombre, red in generar_redes():
        m = start_instance()
        mpc = m.loadcase(red)
        m.savecase(nombre + ".m", mpc)
        n = pc.from_mpc(nombre + ".m")
        graph = top.create_nxgraph(n, multi=False, respect_switches=False)

        loads = set(n.load['bus'])
        generators = set(n.gen['bus'])
        static_generators = set(n.sgen['bus'])

        load_gen = loads.union(generators).union(static_generators)

        outgraph = nx.Graph()
        outgraph.add_nodes_from(graph.nodes)
        outgraph.add_edges_from(graph.edges)

        # mark excluded, pre-selected, zero injection (propagating)
        for n in outgraph:
            outgraph.nodes[n]['propagating'] = 1 * (n not in load_gen)

        print(
            "Instancia "
            + nombre
            + " con "
            + str(graph.number_of_nodes())
            + " vértices y "
            + str(graph.number_of_edges())
            + " aristas."
        )

        nx.write_graphml(outgraph, dir_instancias + nombre + ".graphml")
        os.remove(nombre + ".m")

generar_gml("")
