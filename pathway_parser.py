from Bio.KEGG.KGML.KGML_parser import read


class Node:
    def __init__(self, names):
        self.gene = names
        self.weight = 0
        self.n = len(names)


class Parser:

    def __init__(self):
        self.path = {}

    def get_path(self):
        return self.path

    def parse(self, file):

        def _set_gene(way, path):
            gene_name = [e.name for e in way.entries.values() if e.type == 'gene']
            ent_id = [e.id for e in way.entries.values() if e.type == 'gene']
            name = _parse_name(gene_name)
            for i in range(len(name)):
                path[ent_id[i]] = Node(name[i])

        def _parse_name(nodes):
            names = []
            for node in nodes:
                names.append(node.split())
            for i in range(0,len(names)):
                temp = []
                for name in names[i]:
                    temp.append(name.split(':')[1])
                names[i] = temp
            return names

        def _set_weight(way, path):

            for relation in way.relations:

                if relation.entry1.type == 'gene':
                    path[relation.entry1.id].weight += 1
                if relation.entry2.type == 'gene':
                    path[relation.entry2.id].weight += 1

        def _change_keys(path):
            i = 0
            for key in list(path.keys()):
                path[i] = path.pop(key)
                i = i+1
        f = open('D:/NIR/bernmix_method/hsaFiles/'+file)
        pathway = read(f, 'r')
        _set_gene(pathway, self.path)
        _set_weight(pathway, self.path)
        _change_keys(self.path)
