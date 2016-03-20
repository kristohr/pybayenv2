#!/local/bin/python

from locus import *


#Convert form plink (ped) to bayenv format
def ped2bayenv(ped_file):
    
    map_file = ped_file.split(".")[0] + ".map"
    
    ped = open(ped_file, "r")
    map1 = open(map_file, "r")
    ped_lines = ped.readlines()
    map_lines = map1.readlines()

    snp_list = []

    for line in map_lines:
        l = line.split(" ")
        snp_list.append(l[1])

    ped_info = []
    pop_info = []

    for line in ped_lines:
        line = line.strip("\n")
        l = line.split(" ")
        pop_info.append(l[0])
        ped_info.append(l[6:])

    N =  len(ped_info[0])
    
    alleles_list = []
    for i in range(0, len(ped_info)):
        new_line = []
        for j in range(0, N, 2):
            new_line.append(ped_info[i][j] + ped_info[i][j+1])
        new_line = [w.replace('-9', '00') for w in new_line]
        alleles_list.append(new_line)
    pop_names = []
    for pop_name in pop_info:
        if pop_name not in pop_names:
            pop_names.append(pop_name)

    locus_list = [] #List of locus objects
    
    num_pops = len(pop_names)

    for snp in snp_list:
        locus = Locus(snp, int(num_pops))
        locus_list.append(locus)

    for i in range(0, len(alleles_list)):
        pop = pop_names.index(pop_info[i])
        for j in range(0, len(alleles_list[i])):
            al_type = get_allele_type(alleles_list[i][j])
            locus_list[j].update_freqs(al_type, pop)

    return locus_list


#Convert form genepop to bayenv format
def gen_loci(in_file, num_pops):

    locus_list = [] #List of locus names
    
    dataset = open(in_file, 'r')
    lines = dataset.readlines()
    
    new_lines = []
    test = ""
    for line in lines:
        if (len(line) > 2):
            new_lines.append(line)

    file_info = new_lines[0] #File header
    line = new_lines[1] #Locus names
    header = line.replace(" ", "") 
    all_snps = header.split(',')
    for snp in all_snps:
        locus = Locus(snp, int(num_pops))
        locus_list.append(locus)

    new_lines.pop(0) #Removing the two first lines 
    new_lines.pop(0)
    pop = -1 #Population in the
    for line in new_lines:
        if ("Pop" in line or "pop" in line):
            pop += 1
            continue
        data = line.split()
        alleles = data[1::] #The allele data for each population
        for i in range(0, len(alleles)):
            al_type = get_allele_type(alleles[i])
            locus_list[i].update_freqs(al_type, pop)

    return locus_list
