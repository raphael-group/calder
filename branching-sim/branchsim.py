import random
from math import pow
from datetime import datetime
from scipy.stats import dirichlet, binom, nbinom, poisson
from multiprocessing import Pool
import numpy as np

class Simulation:    
    seed = 0
    mut_idx = 0
        
    N_genotype = None
    driver_muts = None
    
    extant_clones = None
    edges = None
    
    WGS_depth = 80
    WGS_threshold = 0.05
    
    deepseq_depth = 3000
    deepseq_threshold = 0.01
    
    exome_depth = 200
    exome_threshold = 0.02
    
    
    def __init__(self, n_samples, seed=0, generations_between_samples = 10, timeout=1000, abundance_threshold=3000, mut_rate = 0.6, 
                 driver_rate = pow(10, -5), fitness_adv = 0.01, sampling = "exome"):
        
        assert n_samples.__class__ == int and n_samples > 0
        assert seed.__class__ == int and seed >= 0
        assert abundance_threshold > 1
        assert generations_between_samples.__class__ == int and generations_between_samples > 0
        assert timeout.__class__ == int and timeout > 0
        assert mut_rate.__class__ == float and mut_rate > 0
        assert driver_rate.__class__ == float and driver_rate > 0
        assert fitness_adv.__class__ == float and fitness_adv > 0
        assert sampling.__class__ == str and sampling in ["exome", "targeted"]
            
        self.n_samples = n_samples
            
        self.seed = seed # default 0
        self.abundance_threshold = abundance_threshold # default 3000
        self.generations_between_samples = generations_between_samples # default 10
        self.timeout = timeout # default 1000000

        self.mut_rate = mut_rate
        self.driver_rate = driver_rate
        self.fitness_adv = fitness_adv
        
        self.sampling = sampling
        
        random.seed(self.seed)
        np.random.seed(self.seed)
        self.run()


    def run(self):
        self.samples = []
        self.simulate()
        _, _, last_clones, last_muts = self.samples[-1]
        self.edges, _, _, _ = Simulation.clones_to_tree(last_clones, last_muts)
                
            
    def child_nondriver(self, genotype, N_genotype, N_phenotype, genotypes, geno_to_pheno):
        # for each mutating that isn't driver: create new genotype, increment phenotype
        new_mut = self.mut_idx
        self.mut_idx += 1

        # make my new genotype
        new_genotype = genotype + (new_mut, )
        phenotype = geno_to_pheno[genotype]
        
        N_phenotype[phenotype] += 1
        geno_to_pheno[new_genotype] = phenotype
        genotypes[new_genotype] = True
        assert new_genotype not in N_genotype
        N_genotype[new_genotype] = 1
        
    def child_driver(self, genotype, N_genotype, N_phenotype, genotypes, geno_to_pheno):
        # for each new driver genotype: create a new genotype and new phenotype at abundance 1
        new_mut = self.mut_idx
        self.driver_muts[self.mut_idx] = True
        self.mut_idx += 1

        
        # make my new genotype
        new_genotype = genotype + (new_mut, )
        new_phenotype = geno_to_pheno[genotype] + (new_mut, )

        assert new_phenotype not in N_phenotype
        N_phenotype[new_phenotype] = 1
        geno_to_pheno[new_genotype] = new_phenotype
        genotypes[new_genotype] = True
        assert new_genotype not in N_genotype
        N_genotype[new_genotype] = 1
    
    # should be called before any other functions
    def simulate(self):
        
        # mutation index for the next "novel" mutation - should be incremented when a new mutation is added to a clone
        self.mut_idx = 0

        # List of genotype tuples, each tuple representing the set of mutations in a clone (tuple -> True)
        genotypes = {}

        # (genotype tuple) -> number of cells with that genotype
        N_genotype = {}

        # genotype tuple -> phenotype tuple (convenience to just select the driver mutations)
        geno_to_pheno = {}
        
        # (phenotype tuple) -> number of cells with those driver mutations
        N_phenotype = {}

        # integer -> boolean
        self.driver_muts = {}
        
        samples = []
        generated_samples = 0
        
        t = 0
        T = self.timeout
        t_last_sample = -1 * self.generations_between_samples
        
        founder = (self.mut_idx,)

        genotypes[founder] = True
        geno_to_pheno[founder] = founder
        N_genotype[founder] = 1
        N_phenotype[founder] = 1
        
        self.driver_muts[self.mut_idx] = True
        self.mut_idx += 1
        
        while t < T:
            t += 1
            i = 0

            prev_genotypes = list(genotypes.keys())
            
            total_newdrivers = 0
            for genotype in prev_genotypes:
                ## Compute birth probability for this genotype

                phenotype = geno_to_pheno[genotype]
                
                # number of cells with my phenotype
                myN = N_phenotype[phenotype]
                i += myN
                
                # carrying capacity = 50000 * number of driver genes
                myK = 50000 * len(phenotype)
                
                
                birth_prob = .5
                for m in genotype:
                    s = self.fitness_adv if m in self.driver_muts else 0
                    assert s == 0 or m in phenotype
                    birth_prob *= (1 + s) * (1 - myN / myK)
                    
                #birth_prob = min(birth_prob, 0.99)
                
                n_cells = N_genotype[genotype]
                if n_cells == 0:
                    continue
                #print("Genotype " + str(genotype) + "has %d cells" % n_cells)

                    
                #if len(phenotype) > 1:
                    #print(phenotype, n_cells)
                
                n_replicating = sum([1 for _ in range(n_cells) if random.random() < birth_prob])
                n_mutating = sum([1 for _ in range(n_replicating) if random.random() < self.mut_rate])
                n_newdrivers = sum([1 for _ in range(n_mutating) if random.random() < self.driver_rate])
                    

                [self.child_nondriver(genotype, N_genotype, N_phenotype, genotypes, geno_to_pheno) for _ in range(n_mutating - n_newdrivers)]
                [self.child_driver(genotype, N_genotype, N_phenotype, genotypes, geno_to_pheno) for _ in range(n_newdrivers)]

                total_newdrivers += n_newdrivers

                #print("%d cells, %d replicating, %d mutating, %d new drivers" % (n_cells, n_replicating, n_mutating, n_newdrivers))
                
                # Account for population increase for cells that replicate but do not mutate
                N_genotype[genotype] += n_replicating - (n_cells - n_replicating)
                N_phenotype[phenotype] += n_replicating - (n_cells - n_replicating)

                # for mutating cells: call functions
                # for non-replicating cells: decrement genotypes
                
            #if total_newdrivers > 0:
            #    print("%d new drivers in generation %d" % (total_newdrivers, t))    
                    
            if i == 0:
                #print("All cells are dead :(")
                # all cells are dead
                break
            
            if i >= self.abundance_threshold and t > t_last_sample + self.generations_between_samples: 
                print("Taking a sample at generation %d" % t)

                # we have a detectable amount of cells
                self.N_genotype = N_genotype
                
                if self.sampling == "exome":
                    self.samples.append(self.sample())
                elif generated_samples == 0:
                    # Sequence using WGS and then immediately do targeted sequencing on "the same biological sample"
                    self.samples.append(self.sample(kind="wgs"))
                    self.samples.append(self.sample(kind="targeted"))
                else:
                    self.samples.append(self.sample(kind="targeted"))
                    
                generated_samples += 1
                if generated_samples == self.n_samples:
                    self.N_genotype = N_genotype
                    self.t = t
                    return       
                else:
                    t_last_sample = t
    
    def sample(self, kind="exome"):
        # identify the set of mutations that occur at at least detection_threashold frequency
        N = 0
        mut_prevalence = {}
            
        threshold, coverage = {"exome":(self.exome_threshold, self.exome_depth),
                              "wgs":(self.WGS_threshold, self.WGS_depth),
                              "targeted":(self.deepseq_threshold, self.deepseq_depth)}[kind]

        
        # First, aggregate the total number of cells, and the number of cells with each mutation
        for geno, ncells in self.N_genotype.items():
            # add up all the cells we have
            N += ncells
            
            for m in geno:
                # add my population to the prevalence of each of my mutations
                if m in mut_prevalence:
                    mut_prevalence[m] += ncells
                else:
                    mut_prevalence[m] = ncells
        
        if kind == "targeted":
            # Only sequence mutations that were detected in the initial WGS sample
            extant_muts = self.targeted_muts
        else:
            extant_muts = {}
            
        extant_clones = {}
        # Shrink genomes to only detectable mutations and aggregate cells with each genome (detectable clones)
        for geno, ncells in self.N_genotype.items():
            new_geno = []
            for m in geno:
                if m in extant_muts:
                    new_geno.append(m)

                elif kind != "targeted" and mut_prevalence[m] / N >= threshold:
                    new_geno.append(m)
                    extant_muts[m] = True

            if len(new_geno) > 0:
                new_geno = tuple(new_geno)
                if new_geno in extant_clones.items():
                    extant_clones[new_geno] += ncells
                else:
                    extant_clones[new_geno] = ncells
                    
        if self.sampling == "targeted" and kind == "wgs":
            self.targeted_muts = extant_muts
            
                
        geno_to_id = {}
        
        detectable_clones = {}
        n_detectable_cells = 0
        for geno, ncells in extant_clones.items():
            if ncells > 0:
                detectable_clones[geno] = ncells
                n_detectable_cells += ncells
        
        k = 0
        id_to_geno = {}
        alpha = []
        for geno, ncells in detectable_clones.items():
            id_to_geno[k] = geno
            k += 1
            alpha.append(ncells)
        mixture = dirichlet.rvs(alpha, random_state=self.seed)[0]

        # Compute true frequencies, including mutations that are not detectable (in case they are detectable in other samples)
        true_freqs = {}
        for m in mut_prevalence.keys():
            true_freqs[m] = mut_prevalence[m] / (2 * N)

        # For each mutation, compute sequencing depth and number of alternate reads
        sample = {}
        for mut in extant_muts.keys():
            my_depth = np.random.poisson(coverage)
            alt_reads = np.random.binomial(my_depth, true_freqs[mut])
            sample[mut] = int(my_depth- alt_reads), int(alt_reads)
        
        return sample, true_freqs, extant_clones, extant_muts
                        
        
    def clones_to_tree(extant_clones, extant_muts):
        """
        Using the set of extant clones, reconstructs a perfect phylogeny matrix Mbar, 
        row labels "taxa", and column labels 
        """
        
        # Create binary matrix (taxa x mutations)
        rows = []
        idx_to_mut = {}
        for clone in extant_clones.keys():
            my_row = []
            i = 0
            for mut in extant_muts.keys():
                idx_to_mut[i] = mut
                i += 1
                if isinstance(clone, int):
                    my_row.append(1 if mut == clone else 0)
                else:
                    my_row.append(1 if mut in clone else 0)
                
            rows.append(my_row)
        
        col_counts = [sum([rows[r][i] for r in range(len(rows))]) for i in range(len(rows[0]))]
        
        # Do counting sort on columns
        # list of the old column index corresponding to the new column in each position: new_ind[0] is the old index of the (new) 1st col
        new_idx = []
       
        minsum = min(col_counts)
        maxsum = max(col_counts)
        for count in range(maxsum, minsum - 1, -1):
            for i in range(len(col_counts)):
                if col_counts[i] == count:
                    new_idx.append(i)
        
               
        # Sorted matrix by columns
        Mbar = []
        for r in range(len(rows)):
            new_row = [rows[r][new_idx[i]] for i in range(len(rows[0]))]
            Mbar.append(new_row)
        
        col_counts2 = [sum([Mbar[r][i] for r in range(len(Mbar))]) for i in range(len(Mbar[0]))]
        #print(col_counts2)
        
        assert all(col_counts2[i] >= col_counts2[i+1] for i in range(len(col_counts2)-1))
         
        # vertices are labeled by the mutation on the incoming edge    
        edges = {} # source -> [dest, dest, ...]
        
        # add the root vertex to the tree
        taxon = Mbar[0]
        prev = None
        for k in range(len(taxon)):
            if taxon[k] == 1:
                if prev is not None:
                    edges[prev] = [k]
                prev = k
                 
        for idx in range(1, len(Mbar)):
            taxon = Mbar[idx]
            # all taxa should have the first mutation in common
            assert taxon[0] == 1
            prev = 0
            for k in range(1, len(taxon)):
                if taxon[k] == 1:
                    # need to traverse the edge corresponding to this mutation
                    if prev not in edges:
                        # if the vertex isn't in the tree, add it
                        edges[prev] = [k]
                    elif k not in edges[prev]:
                        # if the edge isn't in the tree, add it
                        edges[prev].append(k)
                    prev = k
        
        row_labels = []
        for encoded_geno in Mbar:
            # Make sure that each detected mutation is encoded as present or absent
            assert len(encoded_geno) == len(idx_to_mut)
            my_muts = [idx_to_mut[x] for x in range(len(encoded_geno)) if encoded_geno[x] == 1]
            row_labels.append(tuple(my_muts))

        column_labels = [idx_to_mut[x] for x in range(len(idx_to_mut))]
        
        real_edges = {}
        for source, dests in edges.items():
            real_edges[idx_to_mut[source]] = [idx_to_mut[x] for x in dests]

        return real_edges, Mbar, row_labels, column_labels

def run_sim_exome(seed):
    try:
        s = Simulation(5, seed=seed, abundance_threshold= pow(10, 7), mut_rate = 0.2, driver_rate = pow(10, -3), 
                       fitness_adv=0.02, generations_between_samples = 40)
        totaldrivers = len(s.driver_muts)
        founddrivers = len([x for x in s.driver_muts.keys() if x in s.edges])
        print("seed %d - Total muts: %d, extant muts %d, extant drivers %d" % (seed, s.mut_idx, len(s.samples[-1][-1]), founddrivers))
        return s 
    except IndexError:
        return -1

def construct_patient_exome(sim, seed, error_rate=.001):
    random.seed(seed)
    np.random.seed(seed)

    # Construct patient object from list of samples
    detected_muts = {}
    depth = []
    for sample in sim.samples:
        for mut, (ref, alt) in sample[0].items():
            if alt / (ref + alt) > 0.02:
                detected_muts[mut] = True
                
            depth.append(ref + alt)
    med_depth = np.median(depth)
    
    instance = {}
    instance['m'] = len(sim.samples)

    mut_to_idx = {}
    idx_to_mut = {}
    sample_labels = {}
    for t in range(len(sim.samples)):
        sample_labels[t] = "t%d" % t
        
    # Use all mutations present in any sample
    i = 0    
    for sample in sim.samples:
        for mut in sample[0].keys():
            if mut in detected_muts and detected_muts[mut]:
                if mut not in mut_to_idx:
                    mut_to_idx[mut] = i
                    idx_to_mut[i] = mut
                    i += 1
    instance['n'] = i
    instance['mut_labels'] = idx_to_mut
    instance['mut_to_idx'] = mut_to_idx
    instance['sample_labels'] = sample_labels
    
    data = {}
    for t in range(len(sim.samples)):
        for i in range(instance['n']):
            if idx_to_mut[i] in sim.samples[t][0] and sim.samples[t][0][idx_to_mut[i]][0] > 0:
                data[t, i] = list(sim.samples[t][0][idx_to_mut[i]])
            else:
                var = np.random.binomial(200, error_rate)
                data[t, i] = [200 - var, var]

    instance['edges'] = sim.edges            
    instance['data'] = data

    m = instance['m']
    n = instance['n']
    trueF = np.zeros((m, n))
    for i in range(n):
        mut = instance['mut_labels'][i]
        for t in range(m):
            if mut in sim.samples[t][1]:
                trueF[t][i] = sim.samples[t][1][mut]
            else:
                trueF[t][i] = 0
    instance['trueF'] = trueF

    return instance
        
def write_clustering_input(patient, filename):
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
                print("Found 0-0 entry")
                my_row = [t, patient['sample_labels'][t]] + [0, "null"] + [p, patient['mut_labels'][p]] + [100, 0]
            rows.append([str(x) for x in my_row])
    
    with open(filename, 'w') as f:
        f.write('\n'.join(['\t'.join(row) for row in rows]))

def example_loop():
    print(datetime.now())
    seeds = []
    sims = []
    for seed in range(10):
        try:
            s = Simulation(5, seed=seed, abundance_threshold= pow(10, 7), mut_rate = 0.2, driver_rate = pow(10, -3), fitness_adv=0.02, generations_between_samples = 40)
            totaldrivers = len(s.driver_muts)
            founddrivers = len([x for x in s.driver_muts.keys() if x in s.edges])
            
            if founddrivers > 2:
                print(seed)
                print("Total muts: %d, extant muts %d, extant drivers %d" % (s.mut_idx, len(s.samples[-1][-1]), founddrivers))
                sims.append(s)
                seeds.append(seed)
                if len(sims) > 5:
                    break
            else:
                print("seed %d had only %d drivers" % (seed, founddrivers))
        except IndexError:
            print("seed %d failed" % seed)
    print(datetime.now())
    return sims, seeds