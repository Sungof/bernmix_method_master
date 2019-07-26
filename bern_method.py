import pathway_parser
from scipy.special import comb
import os
def pathway_score(genes, pathway):
    score = 0
    for node in pathway:
        flag = False
        for gene in genes:
            if node.gene.count(gene) != 0:
                flag = True
                break
        if flag:
            score = score+node.weight
    return score


def num_of_genes(path):
    n = 0
    for node in path.values():
        n = node.n + n
    return n


def probability(path, n, m):
    p = []
    for node in path.values():
        prob = 1-comb(n-node.n, m)/comb(n, m)
        p.append(prob)
    return p


def num_de_inpath(pathway, deGenes):
    count=0
    for i in range(0,len(pathway)):
        node = pathway[i]
        for gene in deGenes:
            if node.gene.count(gene) != 0:
                count=count+1
    return count

def main():

    pathways = []
    files = os.listdir('D:/NIR/bernmix_method/hsaFiles')
    for file in files:
        if file.endswith('.xml'):
            p = pathway_parser.Parser()
            p.parse(file)
            pathways.append(p.get_path())

    allGenes = []
    deGenes = []
    datasets = os.listdir('D:/NIR/bernmix_method/datasets')

    for dataset in datasets:
        file = open('D:/NIR/bernmix_method/datasets/' + dataset, 'r')
        data = file.read()
        data = data.split(" ")
        if dataset.endswith('all.txt'):
            allGenes.append(data)
        if dataset.endswith('de.txt'):
            deGenes.append(data)
    probabilities=[]
    j = 0
    for pathway in pathways:
        for i in range(0, len(allGenes)):
            probabilities.append(probability(pathway,len(allGenes[i]),num_de_inpath(pathway,deGenes[i])))
        j=j+1
        print(j)
    print(probability(pathways[0], num_of_genes(pathways[0]), 20))
    '''
    path = read(open('hsa04150.xml'), 'r')
    ent = path.entries.copy()
    print([e.name for e in path.entries.values() if e.type == 'gene'])
    print([e.name for e in path.entries.values() if e.type == 'compound'])

    print("\n".join([e.name for e in path.entries.values() if e.type == 'gene']))
    arr = [e.type for e in path.entries.values()]
    '''


if __name__ == '__main__':
    main()