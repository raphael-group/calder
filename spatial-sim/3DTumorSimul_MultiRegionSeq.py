#! /usr/bin/python

###################################################################################
## A Python script to simulate 3D tumor growth and multi-region sequencing data  ##
## via an agent-based model. Deme subdivision is assumed in order to model cell  ##
## mixing and spatial contraint.                                                 ##
##                                                                               ##
## Spatial model: pripheral growth                                               ##
## Author: Zheng Hu in Curtis lab@Stanford                                       ##
## Release Date: 2/20/2017                                                       ##
###################################################################################
## Also author: Matt Myers in Raphael Lab @ Princeton

import pickle
import sys,os,math,random
import numpy as np
from collections import Counter
import sets

class deme():
    def __init__(self):
        self.present= 0         ## whether the deme is empty or occupied: 0-empty;1-occupied
        self.background = []    ## the background founder lineage after tumor tranformation
        self.advant = []        ## the advantageous cells

def merge_dicts(dicts):
    a = {}
    for d in dicts:
        a.update(d)
    return a

def createLattice(d):
    """
    Create a 3D cubic lattice with side length of 2d+1 where each site contains a empty deme.
    """
    lattice = {}
    for x in range(0,2*d+1):
        for y in range(0,2*d+1):
            for z in range(0,2*d+1):
                lattice[(x,y,z)] = deme()
    return lattice

def neighbor26((a,b,c)):
    """
    Moore neighbourhood: 26 neighbour sites of (a,b,c).
    """
    neighbor = [(a-1, b-1, c-1),(a-1, b-1, c),(a-1, b-1, c+1),(a-1, b, c-1),(a-1, b, c),(a-1, b, c+1),(a-1, b+1, c-1),(a-1, b+1, c),(a-1, b+1, c+1),(a, b-1, c-1),(a, b-1, c),(a, b-1, c+1),(a, b, c-1),(a, b, c+1),(a, b+1, c-1),(a, b+1, c),(a, b+1, c+1),(a+1, b-1, c-1),(a+1, b-1, c),(a+1, b-1, c+1),(a+1, b, c-1),(a+1, b, c),(a+1, b, c+1),(a+1, b+1, c-1),(a+1, b+1, c),(a+1, b+1, c+1)]
    return neighbor

def neighbor6((a,b,c)):
    """
    von Neumann neighbourhood: 6 neighbour sites of (a,b,c).
    """
    neighbor = [(a-1, b, c),(a+1, b, c),(a, b-1, c),(a, b+1, c),(a, b, c-1),(a, b, c+1)]
    return neighbor

def localNeighbor((a,b,c),r):
    """
    A function to search the local neighbour sites of (a,b,c) within an area of radius r in the 3D cubic lattice.
    """
    neighbor = []
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            for z in range(-r,r+1):
                if pow(x,2)+pow(y,2)+pow(z,2) < pow(r+1,2):
                    neighbor += [(a+x,b+y,c+z)]
    return neighbor

def traceLineage(mlineage,mutid):
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell.
    For example, the input ID (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell

    mlineage - the list that could be used to recover the mutational lineage given the most recent mutation id of a lineage
    mutid - the mutation ID of the most recently occurred mutation in the cell
    """
    recent_muts = mutid.split(',')  # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0]      # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace = recent_muts + trace
        recent_muts = mlineage[first_mut].split(',')
        recent_muts = [int(t) for t in recent_muts]
        first_mut = recent_muts[0]
    return trace

def myTraceLineage(mlineage,mutid):
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell.
    For example, the input ID (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell

    mlineage - the list that could be used to recover the mutational lineage given the most recent mutation id of a lineage
    mutid - the mutation ID of the most recently occurred mutation in the cell
    """
    recent_muts = mutid.split(',')  # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0]      # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace = [",".join([str(x) for x in recent_muts])] + trace
        recent_muts = mlineage[first_mut].split(',')
        recent_muts = [int(t) for t in recent_muts]
        first_mut = recent_muts[0]
    return ["0"] + trace

def lowerORupper(value):
    """
    A function to choose the upper or lower integral value given a non-integral number
    """
    lower_int = int(value)
    upper_int = lower_int+1
    if random.random() < value-lower_int:
        return upper_int
    else:
        return lower_int

def initiateFirstDeme(maxsize,lineage,current_id):
    """
    The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process

    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    current_id - the starting mutation ID
    """
    neu_list = [str(current_id)]
    adv_list = []
    current_deme_size = 1
    while current_deme_size < maxsize:
        n1,n2 = len(neu_list),len(adv_list)                         #n1 and n2 are the current number of neutral founder cells and advantageous cells, respectively
        neu_divcells =  int(n1*birth_rate+1)                        #number of dividing cells of neutral lineage in this generation. The other cells will die in the next generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells = lowerORupper(n2*birth_rate*(1+s_coef))   #number of dividing cells of advantageous lineage in this generation
            adv_list = random.sample(adv_list,adv_divcells)*2
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if n1 > 0:
            new_mut1 = np.random.poisson(mut_rate*n1)               # the total number of mutations occurring in a generation follows Poission distribution with lambda=u*n
            mut_assig1 = Counter(np.random.choice(n1,new_mut1))
            for x1 in mut_assig1.keys():
                nmut = mut_assig1[x1]
                new_mut1 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [neu_list[x1]]
                neu_list[x1] = mut_str
        if n2 > 0:
            new_mut2 = np.random.poisson(mut_rate*n2)
            mut_assig2 = Counter(np.random.choice(n2,new_mut2))
            for x2 in mut_assig2.keys():
                nmut = mut_assig2[x2]
                new_mut2 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv_list[x2]]
                adv_list[x2] = mut_str

        if random.random() < adv_rate*n1:                           # occurence of advantageous mutation on the neutral lineage
            current_id += 1
            current_n1 = len(neu_list)
            lineage += [str(neu_list[current_n1-1])]
            adv_list += [str(current_id)]
            neu_list = neu_list[0:current_n1-1]

    return neu_list,adv_list,current_id,lineage


def demeGrowthFission(neu_list,adv_list,lineage,current_id,current_deme_number):
    """
    A function to simulate deme expansion and fission and keep track of the mutational lineages
    """
    current_deme_size = len(neu_list)+len(adv_list)
    while current_deme_size < 2*deme_size:                          #when the deme size doubles, it will split into two offspring demes
        n1,n2 = len(neu_list),len(adv_list)
        neu_divcells =  lowerORupper(n1*birth_rate)                 #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells =  lowerORupper(n2*birth_rate*(1+s_coef))  #number of dividing cells in this generation
            adv_list = random.sample(adv_list,adv_divcells)*2
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if current_deme_number < 5*pow(10,7)/deme_size:             #stop mutation occurring when the tumor size is larger than 5*10^7 cells. The reason is that late occuring mutations have very small chance to present at detectable frequency even under selection.
            if n1 > 0:
                new_mut1 = np.random.poisson(mut_rate*n1)
                mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    new_mut1 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut1))
                    #if nmut > 1:
                    #    for t in new_mut1:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [neu_list[x1]]
                    neu_list[x1] = mut_str
            if n2 > 0:
                new_mut2 = np.random.poisson(mut_rate*n2)
                mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    new_mut2 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut2))
                    #if nmut > 1:
                    #    for t in new_mut2:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [adv_list[x2]]
                    adv_list[x2] = mut_str

            if random.random() < adv_rate*n1:
                current_id += 1
                current_n1 = len(neu_list)
                lineage += [str(neu_list[current_n1-1])]
                adv_list += [str(current_id)]
                neu_list = neu_list[0:current_n1-1]
            #n1,n2 = len(neu_list),len(adv_list)
    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list),0.5)       # the offpring deme size is determined by a Binomial distribution B(n,0.5)
    else:
        offspring_neu = 0
    neu_list1=neu_list[0:offspring_neu]
    neu_list2=neu_list[offspring_neu:len(neu_list)]
    random.shuffle(adv_list)
    if len(adv_list) > 0:
        offspring_adv = np.random.binomial(len(adv_list),0.5)
    else:
        offspring_adv = 0
    adv_list1=adv_list[0:offspring_adv]
    adv_list2=adv_list[offspring_adv:len(adv_list)]

    return neu_list1,neu_list2,adv_list1,adv_list2,current_id,lineage


def seqProcessing(sp,sample_keys,mlineage,size_par,mean_depth,purity):
    """
    Model the random sampling process in NGS and report the sequencing allele frequencies in a sample of cells

    sp- the lattice space
    sample_keys- the locations for the demes in a bulk sample
    size_par- variance parameter for negative-binomial distribution
    mean_depth- the mean depth of the sequencing
    purity- tumor purity
    """
    all_cur_id = []                                     # all most recently occurred mutations
    all_mut_id = []                                     # all mutations in the sampled cells
    for key in sample_keys:
        smuts = list(sp[key].background + sp[key].advant)
        all_cur_id += smuts
    sample_size = 10000                                 # the number of cells for sequencing analysis
    sample_id = random.sample(all_cur_id,sample_size)
    id_count = Counter(sample_id)
    for x in id_count.keys():
        xlineage = traceLineage(mlineage,x)
        all_mut_id += xlineage*id_count[x]
    mut_count = Counter(all_mut_id)
    prob_par=size_par*1.0/(size_par+mean_depth)
    sampleAF = {}                                       # a dictionary storing the mutation IDs and corresponding depth and allele frequency the seq data
    trueAF = {}
    for x in mut_count.keys():
        true_af = mut_count[x]*0.5*purity/sample_size   # the true allele frequency in the sample
        if true_af > .001:                             # filter mutations with very low frequency that is not detectable by ~1000X sequencing depth
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 15:                        # seq depth cutoff for "calling" a mutation
                var_reads = np.random.binomial(site_depth,true_af)
                seq_af = var_reads*1.0/site_depth
                if var_reads >= 4:                      # variant reads cutof for "calling" a mutation
                    trueAF[str(x)] = (mut_count[x], true_af)
                    sampleAF[str(x)] = (site_depth,seq_af)
    return sampleAF, trueAF

def highMuts(sp,position,mlineage,cutoff):
    """
    Obtain the high-frequency mutations (vaf>cutoff) in a particular deme

    sp - the lattice space
    position - the location of the deme
    mlineage - mutation lineage dictionary
    cutoff - the VAF cutoff for a "high-frequency" mutation, e.g. 0.4
    """
    all_cur_id = sp[position].background + sp[position].advant
    all_mut_id = []
    sample_size = 100
    sample_id = random.sample(all_cur_id,sample_size)
    id_count = Counter(sample_id)
    for y in id_count.keys():
        xlineage = traceLineage(mlineage,y)
        all_mut_id += xlineage*id_count[y]
    mut_count = Counter(all_mut_id)
    highAF_muts = []
    for x in mut_count.keys():
        allele_freq = mut_count[x]*1.0/sample_size
        if allele_freq > cutoff:
            highAF_muts += [x]
    return highAF_muts


def pubMutGenerator(n,size_par,mean_depth,purity):
    """
    A function to generate the public clonal mutations occured during the multi-step tumorigenesis before transformation.

    n- number of clonal mutations
    size_par- variation parameter in the negative binomial distribution
    mean_death- mean seq depth
    purity-sample purity
    """
    prob_par=size_par*1.0/(size_par+mean_depth)
    mean_af = 0.5*purity
    depth_pub = []
    maf_pub = []
    for k in range(0,n):
        correct = 0
        while correct == 0:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 15:
                correct =1
        var_reads = np.random.binomial(site_depth,mean_af)
        site_maf = var_reads*1.0/site_depth
        depth_pub += [site_depth]
        maf_pub += [site_maf]
    return depth_pub,maf_pub


def localSampling(region,sample_number,cutoff):
    """
    A function to sampling the locations of multiple bulk samples in a local region.
    """
    success = 0
    while success == 0:
        locations = random.sample(region,sample_number)
        repeat = sample_number*(sample_number-1)
        minall = 999
        for x in range(0,repeat):
            rs = random.sample(locations,2)
            min_distance = min([abs(rs[0][0]-rs[1][0]),abs(rs[0][1]-rs[1][1]),abs(rs[0][2]-rs[1][2])])
            if min_distance < minall:
                minall = min_distance
        if min_distance > 2*cutoff:
            success = 1
    return locations


def bulkTissueSampling(sp,location,radius):
    """
    A function to sampling a bulk sample in a local region.
    """
    local_region = localNeighbor(location,radius)
    bulk_tissue = []
    for x in local_region:
        if sp[x].present == 1:
            bulk_tissue += [x]
    return bulk_tissue


def lineageDashLink(mlist):
    """
    Transform the mutation lineage from list (e.g [1,3,10,20]) to dash-linked string (e.g. 1-3-10-20)
    """
    if len(mlist) > 0:
        dstring = str(mlist[0])
        for x in mlist[1:len(mlist)]:
            dstring += "-"
            dstring += str(x)
        return dstring
    else:
        return "0"

def missingDepth(vafdata,absent_muts,mean_depth):
    """
    Randomly generate the sequencing depth for the mutation-absent sites across samples
    """
    for x in absent_muts:
        done = 0
        while done == 0:
            missing_depth = np.random.negative_binomial(2,2.0/(2+mean_depth))
            if missing_depth >= 15:
                done = 1
        vafdata[str(x)] = (missing_depth,0)
    return vafdata


def sample(current_keys, mutlineage, seq_depth, filename, longitudinal=False, WGS_muts={}):
    """
    A function to generate a collection of 8 multi-site samples from the surface of the 
    simulated tumor. 
    """
    periphery = [] # the locations of periheral demes on tumor surface
    for key in current_keys:
        neikeys = neighbor26(key)
        for z in neikeys:
            if space[z].present == 0:
                periphery +=[key]
                break


    quadrant1,quadrant2,quadrant3,quadrant4,quadrant5,quadrant6,quadrant7,quadrant8 = [],[],[],[],[],[],[],[] #surface demes in the eight quadrants
    for pky in periphery:
        if pky[0] > rd and pky[1] > rd and pky[2] > rd:
            quadrant1 += [pky]
        if pky[0] < rd and pky[1] < rd and pky[2] < rd:
            quadrant2 += [pky]
        if pky[0] < rd and pky[1] > rd and pky[2] > rd:
            quadrant3 += [pky]
        if pky[0] > rd and pky[1] < rd and pky[2] < rd:
            quadrant4 += [pky]

        if pky[0] > rd and pky[1] > rd and pky[2] < rd:
            quadrant5 += [pky]
        if pky[0] < rd and pky[1] < rd and pky[2] > rd:
            quadrant6 += [pky]
        if pky[0] > rd and pky[1] < rd and pky[2] > rd:
            quadrant7 += [pky]
        if pky[0] < rd and pky[1] > rd and pky[2] < rd:
            quadrant8 += [pky]

    #print "# of demes in the periphery=",len(periphery)
    #print "# of demes in quadrant1=",len(quadrant1)
    #print "# of demes in quadrant2=",len(quadrant2)
    #print "# of demes in quadrant3=",len(quadrant3)
    #print "# of demes in quadrant4=",len(quadrant4)
    #print "# of demes in quadrant5=",len(quadrant5)
    #print "# of demes in quadrant6=",len(quadrant6)
    #print "# of demes in quadrant7=",len(quadrant7)
    #print "# of demes in quadrant8=",len(quadrant8)
    #print

    #p4samples = localSampling(quadrant1,8,1)
    ###multisample == "8samples":
    locat1 = random.choice(quadrant1) # location of bulk tissue1
    locat2 = random.choice(quadrant2)
    locat3 = random.choice(quadrant3)
    locat4 = random.choice(quadrant4)
    locat5 = random.choice(quadrant5)
    locat6 = random.choice(quadrant6)
    locat7 = random.choice(quadrant7)
    locat8 = random.choice(quadrant8)
    sample8 = [locat1,locat2,locat3,locat4,locat5,locat6,locat7,locat8]

    tissue1 = bulkTissueSampling(space,sample8[0],3)
    tissue2 = bulkTissueSampling(space,sample8[1],3)
    tissue3 = bulkTissueSampling(space,sample8[2],3)
    tissue4 = bulkTissueSampling(space,sample8[3],3)
    tissue5 = bulkTissueSampling(space,sample8[4],3)
    tissue6 = bulkTissueSampling(space,sample8[5],3)
    tissue7 = bulkTissueSampling(space,sample8[6],3)
    tissue8 = bulkTissueSampling(space,sample8[7],3)
    print "Average # of demes in the 8 bulks",(len(tissue1)+len(tissue2)+len(tissue3)+len(tissue4)+len(tissue5)+len(tissue6)+len(tissue7)+len(tissue8))/8

    """
    maf1, true1 = seqProcessing(space,tissue1,mutlineage,2,seq_depth,1)
    maf2, true2 = seqProcessing(space,tissue2,mutlineage,2,seq_depth,1)
    maf3, true3 = seqProcessing(space,tissue3,mutlineage,2,seq_depth,1)
    maf4, true4 = seqProcessing(space,tissue4,mutlineage,2,seq_depth,1)
    maf5, true5 = seqProcessing(space,tissue5,mutlineage,2,seq_depth,1)
    maf6, true6 = seqProcessing(space,tissue6,mutlineage,2,seq_depth,1)
    maf7, true7 = seqProcessing(space,tissue7,mutlineage,2,seq_depth,1)
    maf8, true8 = seqProcessing(space,tissue8,mutlineage,2,seq_depth,1)
    """
    mafs, true = seqProcessing(space, tissue1 + tissue2 + tissue3 + tissue4 + tissue5 + tissue6 + tissue7 + tissue8, mutlineage, 2, seq_depth, 1)

    MAF_file = open(filename,"w")
    MAF_file.write("mut_id"+" "+"public"+" "+"depth"+" "+"maf" + "\n")
    #MAF_file.write("mut_id"+" "+"public"+" "+"depth1"+" "+"maf1"+" "+"depth2"+" "+"maf2"+" "+"depth3"+" "+"maf3"+" "+"depth4"+" "+"maf4"+" "+"depth5"+" "+"maf5"+" "+"depth6"+" "+"maf6"+" "+"depth7"+" "+"maf7"+" "+"depth8"+" "+"maf8")
    #MAF_file.write("\n")

    for k in range(0,npub):
        pdepth,pmaf = pubMutGenerator(8,2,seq_depth,1)
        depth = sum(pdepth)
        maf = sum([pdepth[i] * pmaf[i] for i in range(len(pdepth))]) / depth
        MAF_file.write("0"+" "+"1"+" "+str(depth)+" "+str(maf))
        MAF_file.write("\n")

    """
    muts_all = sets.Set(maf1.keys()) | sets.Set(maf2.keys()) | sets.Set(maf3.keys()) | sets.Set(maf4.keys()) |sets.Set(maf5.keys()) | sets.Set(maf6.keys()) |sets.Set(maf7.keys()) | sets.Set(maf8.keys())

    absent1 = muts_all-sets.Set(maf1.keys())
    absent2 = muts_all-sets.Set(maf2.keys())
    absent3 = muts_all-sets.Set(maf3.keys())
    absent4 = muts_all-sets.Set(maf4.keys())
    absent5 = muts_all-sets.Set(maf5.keys())
    absent6 = muts_all-sets.Set(maf6.keys())
    absent7 = muts_all-sets.Set(maf7.keys())
    absent8 = muts_all-sets.Set(maf8.keys())

    maf1 = missingDepth(maf1,absent1,seq_depth)
    maf2 = missingDepth(maf2,absent2,seq_depth)
    maf3 = missingDepth(maf3,absent3,seq_depth)
    maf4 = missingDepth(maf4,absent4,seq_depth)
    maf5 = missingDepth(maf5,absent5,seq_depth)
    maf6 = missingDepth(maf6,absent6,seq_depth)
    maf7 = missingDepth(maf7,absent7,seq_depth)
    maf8 = missingDepth(maf8,absent8,seq_depth)

    samplemix = createMixture([maf1, maf2, maf3, maf4, maf5, maf6, maf7, maf8])
    """
    muts_all = sets.Set(mafs.keys())

    #print 'Found %d mutations at above 0.05' % len(samplemix)
    
    for mt in list(muts_all):
        if longitudinal and not mt in WGS_muts:
            continue

        #MAF_file.write(str(mt)+" "+"0"+" "+str(maf1[mt][0])+" "+str(maf1[mt][1])+" "+str(maf2[mt][0])+" "+str(maf2[mt][1])+" "+str(maf3[mt][0])+" "+str(maf3[mt][1])+" "+str(maf4[mt][0])+" "+str(maf4[mt][1])+" "+str(maf5[mt][0])+" "+str(maf5[mt][1])+" "+str(maf6[mt][0])+" "+str(maf6[mt][1])+" "+str(maf7[mt][0])+" "+str(maf7[mt][1])+" "+str(maf8[mt][0])+" "+str(maf8[mt][1]))
        #MAF_file.write(str(mt) + " 0 " + str(samplemix[mt][0]) + " " + str(samplemix[mt][1]))
        MAF_file.write(str(mt) + " 0 " + str(mafs[mt][0]) + " " + str(mafs[mt][1]))

        MAF_file.write("\n")

    MAF_file.close()

    #truemix = createMixture([true1, true2, true3, true4, true5, true6, true7, true8])
    with open("true_" + filename, "w") as f:
        f.write("mut_id"+" "+"public"+" "+"ncells"+" "+"maf")
        f.write("\n")
        for mt in list(muts_all):
            if longitudinal and not mt in WGS_muts:
                continue
            
            f.write(str(mt)+" "+"0"+" "+str(true[mt][0])+" "+str(true[mt][1]))
            f.write("\n")

    return muts_all

def createMixture(maf_list):
    """ 
    Given a list of dictionaries, where each maps mutation -> (depth, MAF), combines these individual samples into a single bulk sample.
    """
    all_muts = []
    for mafs in maf_list:
        all_muts += list(mafs.keys())
    
    result = {}
    for mut in all_muts:
        total_reads = 0
        alt_reads = 0

        for mafs in maf_list:
            if mut in mafs:
                item = mafs[mut]
                total_reads += item[0]
                alt_reads += item[0] * item[1]
        result[mut] = total_reads, alt_reads / total_reads
    return result


#############main script to simulate a tumor and multi-region sequencing data#########
###parameter intiation###
deme_size = int(sys.argv[1])        # the deme size
#mut_rate = float(sys.argv[2])       # the neutral mutation rate at whole exonic region
#adv_rate = float(sys.argv[3])       # the advantageous mutation rate at each cell generation
#s_coef = float(sys.argv[4])         # the selection coefficient
seed = int(sys.argv[2])             # replication of simulation


# Using default values from repository / recommended parameters from paper
mut_rate = .6
adv_rate = .00001
s_coef = .1

id = seed
WGS_depth = 30
targeted_depth = 3000

np.random.seed(seed=seed)
random.seed(seed)

rd = 60                             # the side length of the 3D space
final_tumor_size = pow(10,7) / 2        # the number of cells in the final tumor
final_deme_number = final_tumor_size/deme_size    # the final number of demes in the tumor
birth_rate = 0.55                   # the birth probability at each cell generation during tumor growth
npub=1                            # the number of public mutation to be generated
seq_depth=80                        # the average sequencing depth
percentage = int(s_coef*100)        # the percentage form of the selection

mut_id = 0
mutlineage = ['0']                  # the lineage tracer
######################################################################################

first_neu,first_adv,mut_id,mutlineage = initiateFirstDeme(deme_size,mutlineage,mut_id)  #the growth of the fisrt deme from single transformed cell

#print(mutlineage[-1])
#print traceLineage(mutlineage, mutlineage[-1])

space = createLattice(rd)
space[(rd,rd,rd)] = deme()                      #initiate the space with a empty deme in the center site (rd,rd,rd)
space[(rd,rd,rd)].present = 1
space[(rd,rd,rd)].background = list(first_neu)
space[(rd,rd,rd)].advant = list(first_adv)
current_keys = [(rd,rd,rd)]
current_deme_number =1                                 #current deme number
surface_keys = [(rd,rd,rd)]
surface_deme_number =1
deme_time_generation = 0

WGS_done = False
deep_seq_iter = 0

# Runs simulation
while current_deme_number < final_deme_number:
    print current_deme_number
    n_cells = current_deme_number * deme_size

    new_keys = []
    for w in range(0,surface_deme_number):             # deme expansion occurs in the surface of a tumor
        ckey = random.choice(current_keys)
        if space[ckey].present == 1:
            rx,ry,rz = ckey[0],ckey[1],ckey[2]
            nei_sites = neighbor26((rx,ry,rz))  # neighbor sites of (rx,ry,rz)
            empty_sites = []                    # the empty neighbor sites
            for key in nei_sites:
                if space[key].present == 0:
                    empty_sites += [key]
            if len(empty_sites) > 0:
                rand_prob = random.random()
                if rand_prob < 1-math.exp(-len(empty_sites)*0.25): # the probability that a deme is chosen for expansion and split in a given step is proportional to the # of empty neighbor sites
                    pre_neu = list(space[(rx,ry,rz)].background)
                    pre_advant = list(space[(rx,ry,rz)].advant)
                    post_neu1,post_neu2,post_adv1,post_adv2,mut_id,mutlineage = demeGrowthFission(pre_neu,pre_advant,mutlineage,mut_id,current_deme_number)
                    nextkey = random.choice(empty_sites)
                    space[ckey].background = list(post_neu1)
                    space[ckey].advant = list(post_adv1)
                    space[nextkey].background = list(post_neu2)
                    space[nextkey].advant = list(post_adv2)
                    space[nextkey].present = 1
                    current_keys += [nextkey]
                    current_deme_number += 1
                    new_keys += [nextkey]
    ###update surface###
    surface_update = list(surface_keys+new_keys)
    surface_keys = []
    for fkey in surface_update:
        neisites = neighbor26(fkey)
        random.shuffle(neisites)
        for key in neisites:
            if space[key].present == 0:
                surface_keys += [fkey]
                break
    surface_deme_number = len(surface_keys)
    current_deme_number = len(current_keys)
    deme_time_generation = 0

    #n_cells = current_deme_number * deme_size

    if current_deme_number > 1000:
        if not WGS_done:
            wgs_muts = sample(current_keys, mutlineage, WGS_depth, "instance%d_WGS.txt" % id)
            WGS_done = True
            print "Initial WGS found %d mutations" % len(wgs_muts)
        else:
            final_muts = sample(current_keys, mutlineage, targeted_depth, "instance%d_longitudinalsample%d.txt" % (id, deep_seq_iter), longitudinal=True, WGS_muts = wgs_muts)
            deep_seq_iter += 1


"""
final_muts = wgs_muts
mut = list(final_muts)[-1]
print mut
print traceLineage(mutlineage, mut)
print mut in final_muts

detected = set([int(x.split(',')[0]) for x in list(final_muts)])
"""
"""
edges = []
to_visit = list(final_muts)
to_skip = {}
for mt in to_visit:
    if mt in to_skip:
        continue
    to_skip[mt] = True
    # trace the lineage of this mutation
    lineage = traceLineage(mutlineage, mt)
    found_lineage = [x for x in lineage if x in detected]
    #print found_lineage
    for i in range(len(found_lineage) - 1):
        to_skip[found_lineage[i]] = True
        edges.append([found_lineage[i], found_lineage[i + 1]])

print len(edges)
"""

with open("instance%d_lineages.txt" % id, "w") as f:
    for mut in wgs_muts:
        f.write(mut + " " + "-".join([str(x) for x in myTraceLineage(mutlineage, mut)]) + "\n")
    

#with open("instance%d_mutlineage.txt" % id, "w") as f:
#    f.write(" ".join(mutlineage))


####visulization of spatial clonal structure in the central plane###
#central_plane = []
#for key in current_keys:
#    if key[2] == rd:
#        central_plane += [key]

#print "# of demes on the central plane=",len(central_plane)

#map_file = open("CloneMap3D_peri_u"+str(mu)+"_birth_rate"+str(birth_rate)+"_s"+str(s_coef)+"_"+str(repl)+".txt","w")
#map_file.write("x"+" "+"y"+" "+"z"+" "+"lineage")
#map_file.write("\n")
#for key in central_plane:
#    cur_muts = highMuts(space,key,mutlineage,0.4)
#    cur_lineage = lineageDashLink(sorted(cur_muts))
#    map_file.write(str(key[0])+" "+str(key[1])+" "+str(key[2])+" "+str(cur_lineage))
#    map_file.write("\n")

#filename = "simulMRS_deme"+str(deme_size)+"_s"+str(percentage)+"percent_8samples_u"+str(mut_rate)+"_"+str(repl)+".txt"
#final_muts = sample(current_keys, mutlineage, seq_depth, filename)