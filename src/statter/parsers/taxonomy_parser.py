import logging
from networkx import DiGraph, number_of_nodes, number_of_edges, set_node_attributes
from typing import Dict

"""
Parse NCBI taxonomy into a DAG
"""
logger = logging.getLogger(__file__)


class Taxonomy:
    def __init__(self, nodes: str, names: str) -> None:
        """
        nodes: nodes.dmp, see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt
        names: names.dmp, see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt
        """
        self.nodes = nodes
        self.names = names

    def _parse_names(self) -> Dict[str, str]:
        """
        Helper function
        parse names.dmp file for all "scientific name"s
        """
        names = {}
        with open(self.names, "r") as nh:
            for n in nh:
                ndat = n.strip("\t|\n").split("\t|\t")
                if ndat[-1] != "scientific name":
                    continue
                names[ndat[0]] = ndat[1]
        if len(names) == 0:
            raise RuntimeError(
                f"Cannot parse scientific names from {names}! Check the file format"
            )
        logger.info(f"Found {len(names)} taxonomy ids and scientific names")
        return names

    def build_dag(self) -> DiGraph:
        """
        Helper function
        build taxonomy DAG from nodes.dmp file
        """
        taxonomy = DiGraph()
        ranks = {}  # taxonomy ranks
        with open(self.nodes, "r") as nh:
            for n in nh:
                ndat = n.strip("\t|\n").split("\t|\t")
                ranks[ndat[0]] = ndat[2]
                if ndat[0] == ndat[1]:  # no self loops
                    continue
                taxonomy.add_edge(ndat[0], ndat[1])
        _nedge = number_of_edges(taxonomy)
        if _nedge == 0:
            raise RuntimeError(f"Cannot create taxonomy rank tree from {self.nodes}!")
        logger.info(
            f"Found {number_of_nodes(taxonomy)} taxonomy ids and {_nedge} relationships"
        )
        # add rank
        set_node_attributes(taxonomy, ranks, "rank")
        # add names
        set_node_attributes(taxonomy, self._parse_names(), "name")
        return taxonomy
