# Author: Matt Myers (created 4/4/2019)

import sys

def apply_clustering(datafile, cluster_assignments, outfile):
    data = {}
    sample_info = {}
    anatomical_site_info = {}
    mutations = {}
    with open(datafile, "r") as f:
        a = f.readline().split()[0] # number of anatomical sites
        m = int(f.readline().split()[0]) # number of samples
        n = int(f.readline().split()[0]) # number of mutations
        header = f.readline().split() # header
        for line in f:
            tkns = line.split()
            sample_idx = int(tkns[0])
            sample_label = tkns[1]
            site_idx = tkns[2]
            site_label = tkns[3]
            # Store information about each sample
            sample_info[sample_idx] = sample_label, site_idx, site_label

            mut = tkns[-3]
            mutations[mut] = True

            ref = int(tkns[-2])
            var = int(tkns[-1])
            data[sample_idx, mut] = ref, var

    print("Found %d samples and %d mutations" % (m, n))

    cluster_to_muts = {} # map from cluster index to list of mutations
    idx = 0
    with open(cluster_assignments, "r") as f:
        for line in f:
            tkns = line.strip().split(';')
            assert idx not in cluster_to_muts
            cluster_to_muts[idx] = tkns
            idx += 1

    print("Found %d clusters" % idx)
    for cluster in cluster_to_muts.values():
        for mut in cluster:
            if mut not in mutations:
                print("Error: cluster assignments include a mutation absent from the original data: %s" % mut)
                quit(1)
            if not mutations[mut]:
                print("Error: cluster assignments include a mutation multiple times: %s" % mut)
                quit(1)
            mutations[mut] = False
    if any(mutations.values()):
        print("Error: cluster assignments do not include all mutations.")
        print("Missing mutations: ", [x for x in mutations.keys() if mutations[x]])
        quit(1)

    # Create output file
    rows = []

    # Add header
    header = [""]
    for cl_idx in range(idx):
        cl_name = ";".join(cluster_to_muts[cl_idx])
        header.append(cl_name)
        header.append(cl_name)

    rows.append(header)        

    for sample_idx in range(int(m)):
        myrow = []
        myrow.append(sample_info[sample_idx][0])
        for cl_idx in range(idx):
            total_ref = sum([data[sample_idx, x][0] for x in cluster_to_muts[cl_idx]])
            total_var = sum([data[sample_idx, x][1] for x in cluster_to_muts[cl_idx]])
            myrow.append(total_ref)
            myrow.append(total_var)
        rows.append([str(x) for x in myrow])

    with open(outfile, "w") as f:
        f.write('\n'.join(['\t'.join(row) for row in rows]))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python %s [clustering input file] [cluster assignments] [outfile]" % sys.argv[0])
        quit(1)
    else:
        datafile = sys.argv[1]
        cluster_assignments = sys.argv[2]
        outfile = sys.argv[3]
        apply_clustering(datafile, cluster_assignments, outfile)