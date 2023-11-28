import os
import pandas as pd
import json
import click


@click.command()
@click.argument('pairtree_ssm')
@click.argument('pairtree_clusters')
@click.argument('calder_input_file')
def pairtree2calder(pairtree_ssm,
                    pairtree_clusters,
                    calder_input_file):
    pairtree_json = json.load(open(pairtree_clusters, 'r'))
    clusters = pairtree_json['clusters']
    samples = pairtree_json['samples']
    print(f"Assuming that samples are in the following longitudinal order: {samples}")
     
    ssm = pd.read_table(os.path.join(pairtree_ssm))
    
    rows = [[f'sample{i}'] for i in range(len(samples))]
    for cl in clusters:
        my_ssm = ssm[ssm.id.isin(cl)]
        my_var = my_ssm.var_reads.str.split(',', expand=True).astype(int).sum(axis = 0).to_numpy()
        my_tot = my_ssm.total_reads.str.split(',', expand=True).astype(int).sum(axis = 0).to_numpy()
        
        rows[0].extend([str(my_tot[0]), str(my_var[0])])
        rows[1].extend([str(my_tot[1]), str(my_var[1])])

    with open(os.path.join(calder_input_file), 'w') as f:
        header = [''] + [f'cluster{i//2 + 1}' for i in range(len(clusters) * 2)]
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(r) + '\n')

if __name__ == "__main__":
    pairtree2calder()