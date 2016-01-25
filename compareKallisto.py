import argparse, os, sys
#from HTMLChart import *

class KallistoResult:

    def __init__(self, target_id, length, eff_length, est_counts, tpm):
        self.target_id = target_id
        self.length = length
        self.eff_length = eff_length
        self.est_counts = est_counts
        self.tpm = tpm

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}'.format(self.target_id, self.length, self.eff_length, self.est_counts, self.tpm)

def parseFullGenes(name, data, tid_to_gene):
    with open('/projects/btl/egibb/CCLE/' + name + '/abundance.tsv', 'r') as f:
        header = f.readline()
        for line in f:
            k = KallistoResult(*line.strip().split('\t'))
            try:
                gene = tid_to_gene[k.target_id]
            except KeyError:
                continue
            if gene not in data:
                data[gene] = {name: k}
            if name not in data[gene]:
                data[gene][name] = k
    return None

def parseEnsemblGeneNames(ens):
    tid_to_gene = {}
    gene_to_tid = {}
    with open(ens, 'r') as f:
        line = f.readline()
        for line in f:
            tid, gene = line.strip().split()
            tid_to_gene[tid] = gene
            gene_to_tid[gene] = tid
    return tid_to_gene, gene_to_tid

def parseKallisto(tsv, data, delim=('|')):
    with open(tsv, 'r') as f:
        line = f.readline()
        for line in f:
            k = KallistoResult(*line.strip().split('\t'))
            k.target_id = os.path.basename(tsv).split('.')[0] + delim + k.target_id
            sample, kind, gene, chrom, start, end = k.target_id.split(delim)
            region = ('-').join([chrom, str(start), str(end)])
            if chrom not in data:
                data[chrom] = {gene: {region: {sample: k.tpm}}}
            if gene not in data[chrom]:
                data[chrom][gene] = {region: {sample: k.tpm}}
            if region not in data[chrom][gene]:
                data[chrom][gene][region] = {sample: k.tpm}
            else:
                data[chrom][gene][region][sample] = k.tpm
        return data

def r_bar_plot(save_path, values, labels, title):
    result = 'B <- c({})\n'.format((',').join(values))
    result += 'L <- c("{}")\n'.format(('","').join(labels))
    result += 'jpeg("{}")\n'.format(save_path)
    result += 'barplot(B, main="{}", horiz=TRUE, names.arg=L)\n'.format(title)
    result += 'dev.off()'
    return result

def r_line_graph(samples, regions, gene, max_y, title='title', xlab="xlab", ylab="ylab"):
    keys = samples.keys()
    colors = ['red','blue','green','orange','pink','purple','gold','darkslategray','darkred','darkolivegreen1', 'dodgerblue']
    plot = 'first <- c({})\n'.format((',').join(samples[keys[0]]))
    plot += 'jpeg("{}.jpg", width=1500, height=1000, units="px", res=300)\n'.format(gene)
    plot += 'plot(first, xlab="start-end coords", ylab="TPM", ylim=c(0,{}), type="o", col="{}", xaxt="n")\n'.format(max_y,colors[0])
    for i,sample in enumerate(keys[1:]):
        plot += 'samp_{} <- c({})\n'.format(i, (',').join(samples[sample]))
        plot += 'lines(samp_{}, type="o", col="{}")\n'.format(i, colors[i+1])
    plot += 'axis(1, at=1:{}, lab=c("{}"))\n'.format(len(regions),('","').join(regions))
    plot += 'title(main="{}", col.main="red")\n'.format(title)
    plot += 'dev.off()'
    return plot

def genGeneMatrix(gene_dict):
    m = []
    col_labels = gene_dict.keys()
    row_labels = gene_dict[col_labels[0]].keys()
    samples = {}
    for region in sorted(gene_dict, key = lambda x: int(x.split('-')[1])):
        for sample in gene_dict[region]:
            if sample not in samples:
                samples[sample] = [gene_dict[region][sample]]
            else:
                samples[sample].append(gene_dict[region][sample])
    for row in row_labels:
        m.append(samples[row])
    return [m, row_labels, col_labels]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare kallisto tsv files with the same regions.')
    parser.add_argument('-k', '--kallisto', nargs='+', help='Kallisto quant tsv output file')

    args = parser.parse_args()

    gene_names_file = '/home/dmacmillan/annotations/ensembl/ensemblToGeneName.original'

    tid_to_gene, gene_to_tid = parseEnsemblGeneNames(gene_names_file)
    
    data = {}

    names = []

    for k in args.kallisto:
        names.append(k.split('.')[0])
        parseKallisto(k, data)
    
    full_genes = {}

    for sample in names:
        parseFullGenes(sample, full_genes, tid_to_gene)

    for chrom in data:
        #print 'chrom: {}'.format(chrom)
        last_chrom = chrom
        last_gene = None
        for gene in data[chrom]:
            m, rows, cols = genGeneMatrix(data[chrom][gene])
            if len(m[0][0]) < 2:
                continue
            print 'gene: {}'.format(gene)
            print m
            continue
#            last_gene = gene
#            xrange = 'xrange <- {}'.format(len(data[chrom][gene]))
#            samples = {}
#            max_y = 0
#            regions = []
#            ordered_regions = sorted(data[chrom][gene], key=lambda x: int(x.split('-')[1]))
#            first = ordered_regions[0].split('-')
#            strand = None
#            if (int(first[2])-int(first[1]) == 100):
#                ordered_regions = reversed(ordered_regions)
#            for region in ordered_regions:
#                r = region.split('-')
#                #print '    region:\t{}\t{}'.format(r[1], r[2])
#                #raw_input('*')
#                regions.append(r[1] + '-' + r[2])
#                for sample in data[chrom][gene][region]:
#                    if sample not in samples:
#                        samples[sample] = [data[chrom][gene][region][sample]]
#                    else:
#                        samples[sample].append(data[chrom][gene][region][sample])
#                    temp_max = max([float(x) for x in samples[sample]])
#                    if temp_max > max_y:
#                        max_y = temp_max
#            yrange = 'yrange <- {}'.format(max_y)
#            if len(samples[samples.keys()[0]]) < 4:
#                continue
#            r = r_line_graph(samples, regions, last_gene, max_y, title=('_').join([last_chrom,last_gene]))
#            with open('./regions/{}.r'.format(last_gene), 'w') as f:
#                f.write(r)
