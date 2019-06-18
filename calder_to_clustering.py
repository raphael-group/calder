# Author: Matt Myers (created 6/14/2019)

import sys

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    p = read_CALDER(infile)
    write_MACHINA(p, outfile)

# Reads a CALDER input file and constructs a patient data structure
def read_CALDER(infile):
    # read in CALDER format input file
    with open(infile, "r") as f:
        tokens = f.readline().split("\t")
        header = [x.strip() for x in tokens if len(x) > 0]
        if len(header) % 2 == 1:
            print("Parsing error: found an odd number of non-empty entries in the header row.")
            sys.exit(1)
        n = int(len(header) / 2)
        mut_labels = [header[i] for i in range(0, len(header), 2)]
        sample_labels = []
        
        data = {}
        for line in f:
            tokens = line.split("\t")
            sample_labels.append(tokens[0])
            values = tokens[1:]
            if len(values) != n * 2:
                print("Parsing error: expected %d entries in line %d, found %d." % (n * 2 + 1, len(sample_labels) + 1, len(tokens)))
                sys.exit(1)
            
            for i in range(n):
                data[len(sample_labels) - 1, i] = [int(values[i * 2].strip()), int(values[i * 2 + 1].strip())]
        m = len(sample_labels)
        print("Found %d samples and %d mutations" % (int(m), int(n)))

    patient = {}
    patient['m'] = len(sample_labels)
    patient['sample_labels'] = sample_labels
    patient['mut_labels'] = mut_labels
    patient['n'] = n
    patient['data'] = data
    return patient

# Writes a "patient" data structure in MACHINA format, which is what the absence-aware clustering uses
def write_MACHINA(patient, filename):
    m = patient['m']
    n = patient['n']
    
    header = ["#sample_index","sample_label","anatomical_site_index","anatomical_site_label","character_index","character_label","ref","var"]

    rows = []
    rows.append(["1 # anatomical sites"])
    rows.append([str(m) + " # samples"])
    rows.append([str(n) + " # mutations"])
    rows.append(header)
    
    for t in range(m):
        for p in range(n):
            if (t,p) in patient['data']:
                my_row = [t, patient['sample_labels'][t]] + [0, "null"] + [p, patient['mut_labels'][p]] + patient['data'][(t, p)]
            else:
                print(t,p)
                print("Found 0-0 entry")
                my_row = [t, patient['sample_labels'][t]] + [0, "null"] + [p, patient['mut_labels'][p]] + [100, 0]
            rows.append([str(x) for x in my_row])
    
    with open(filename, 'w') as f:
        f.write('\n'.join(['\t'.join(row) for row in rows]))
        
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python %s [CALDER-formatted input file] [output filename]" % sys.argv[0])
        quit(1)
    else:
        main()