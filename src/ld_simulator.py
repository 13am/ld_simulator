#!/usr/bin/python

from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

import random
import sys
from optparse import OptionParser
from itertools import imap

try:
    import numpy as np
    numpy_imported = True
except:
    numpy_imported = False

def parse_options():
    
    parser = OptionParser()

    parser.add_option("--mutation", 
                      type="float", 
                      action="store", 
                      dest="m",
                      default=5e-8*1000)
    
    parser.add_option("--recombination", 
                      type="float", 
                      action="store", 
                      dest="r",
                      default=0.001)

    parser.add_option("--realistic", 
                      action="store_true", 
                      dest="realistic",
                      default=False)
    
    parser.add_option("--sample-disease", 
                      action="store_true", 
                      dest="sample_disease",
                      default=False)

    parser.add_option("--forever", 
                      action="store_true", 
                      dest="forever",
                      default=False)
    
    parser.add_option("--size", 
                      type="int", 
                      action="store", 
                      dest="size",
                      default=150)

    (options, args) = parser.parse_args()
    
    return options

def pearsonr(x, y):
    assert len(x) == len(y)
    n = len(x)
    sum_x = float(sum(x))
    sum_y = float(sum(y))
    sum_x_sq = sum(map(lambda x: pow(x, 2), x))
    sum_y_sq = sum(map(lambda x: pow(x, 2), y))
    psum = sum(imap(lambda x, y: x * y, x, y))
    num = psum - (sum_x * sum_y/n)
    den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)
    if den == 0: return 0
    return num / den

def make_position_vector(chromosomes, position):
    v = []
    for c in chromosomes:
        if c[position] == "-":
            v.append(0)
        else:
            v.append(1)
    return v

def get_best_tag_position(sample):
    chromosomes = []
    for s in sample:
        chromosomes.append(s[0])
        chromosomes.append(s[1])
    d_pos = 0
    d_vector = []
    for c in chromosomes:
        if 'D' in c:
            d_pos = c.find('D')
            d_vector = make_position_vector(chromosomes, d_pos)
            break
    if d_vector == []:
        return []
    tag_pos = [0]
    tag_r = 0
    for i in range(0, len(chromosomes[0])):
        if i == d_pos:
            continue
        v = make_position_vector(chromosomes, i)
        r = pearsonr(d_vector, v)
        if r > tag_r:
            tag_r = r
            tag_pos = [i]
        elif r == tag_r:
            tag_pos.append(i)
    if tag_r == 0:
        return []
    return tag_pos

def mutate(chromosome, alts):
    chromosome = list(chromosome)
    for pos, base in enumerate(chromosome):
        if random.random() < p_mut:
            chromosome[pos] = random.choice(alts)
    chromosome = "".join(chromosome)
    return chromosome

def numpy_mutate(chromosome, alts):
    chromosome = list(chromosome)
    chromosome_len = len(chromosome)
    n = np.random.binomial(chromosome_len, p_mut)
    while n > 0:
        pos = random.choice(range(chromosome_len))
        chromosome[pos] = random.choice(alts)
        n -= 1
    new_chromosome = "".join(chromosome)
    return new_chromosome

def recombine(a, b):
    cross_overs = []
    for i in range(0, len(a)):
        if random.random() < p_rec:
            cross_overs.append(i)
    for i in cross_overs:
        a_left = a[0:i]
        b_left = a[0:i]
        a_right = a[i:]
        b_right = b[i:]
        a = a_left + b_right
        b = b_left + a_right
    return [a, b]

def numpy_recombine(a, b):
    a_len = len(a)
    n = np.random.binomial(a_len, p_rec)
    cross_overs = random.sample(range(a_len), n)
    for i in cross_overs:
        a_left = a[0:i]
        b_left = a[0:i]
        a_right = a[i:]
        b_right = b[i:]
        a = a_left + b_right
        b = b_left + a_right
    return [a, b]

def count_polymorphic_sites(pool):
    n = 0
    for i in range(0, len(pool[0])):
        for j in pool:
            if pool[0][i] != j[i]:
                n += 1
                break
    return n

def advance_generation(pool, disease_in_pool, n_offspring):
    start_size = len(pool)
    new_pool = []
    disease_mutation = ""
    if numpy_imported:
        mut_fun = numpy_mutate
        rec_fun = numpy_recombine
    else:
        mut_fun = mutate
        rec_fun = recombine
    while len(pool) > 1 and len(new_pool) < start_size:
        a1 = pool.pop()
        b1 = pool.pop()
        for i in range(0, n_offspring):
            recombined = rec_fun(a1,b1)
            a2 = recombined[0]
            b2 = recombined[1]
            if disease_in_pool:
                a2 = mut_fun(a2, ['a', 'c', 'g', 't'])
                b2 = mut_fun(b2, ['a', 'c', 'g', 't'])
            else:
                a2 = mut_fun(a2, ['a', 'c', 'g', 't', 'D'])
                if 'D' in a2:
                    disease_in_pool = True
                    disease_mutation = a2
                    b2 = mut_fun(b2, ['a', 'c', 'g', 't'])
                else:
                    b2 = mut_fun(b2, ['a', 'c', 'g', 't', 'D'])
                    if 'D' in b2:
                        disease_mutation = b2
                        disease_in_pool = True
            new_pool.append(a2)
            new_pool.append(b2)
    random.shuffle(new_pool)
    return [new_pool, disease_mutation]

def make_diploid(sample):
    inds = []
    i = 0
    while i < len(sample) - 1:
        inds.append([sample[i], sample[i+1]])
        i += 2
    return inds

options = parse_options()

l_chr = options.size
n_chr_sample = 48
n_chr_population = 1000*2
n_offspring = 2
p_mut = options.m
p_rec = options.r
burn_in = 500
pool = []

if options.realistic:
    p_mut = 5e-8
    p_rec = 0.01/1000000

#populate the pool with chromosomes
for i in range(0, n_chr_population):
    c = ""
    x = '-'
    for j in range(0, l_chr):
        c = c + x
    pool.append(c)
    
#run simulation
generation = 0

#first run burn-in time
if options.sample_disease:
    while burn_in > 0:
        burn_in -= 1
        generation += 1
        pool = advance_generation(pool, True, n_offspring)
        pool = pool[0]

dah = open('disease_alleles.txt', 'w')
dm = ""
#try:
while True:
    sys.stdout.flush()
    generation += 1
    disease_in_pool = False
    for i in pool:
        if 'D' in i:
            disease_in_pool = True
            dah.write("\n GEN "+ str(generation) + "\t\t" + "5'" + i + "3'")
    if generation > 1:
        pool = advance_generation(pool, disease_in_pool, n_offspring)
        new_dm = pool[1]
        pool = pool[0]
        if new_dm != "":
            dm = new_dm
    sys.stdout.write("\x1b[2J\x1b[H")
    n_poly = count_polymorphic_sites(pool)
    inds = make_diploid(pool)
    sample = []
    #first sample the affecteds
    if options.sample_disease:
        for i in inds:
            if len(sample) >= n_chr_sample/2/2:
                break
            if 'D' in i[0] + i[1]:
                sample.append(i)
    #then sample the unaffecteds
    for i in inds:
        if len(sample) >= n_chr_sample/2:
            break
        if 'D' not in i[0] + i[1]:
            sample.append(i)
    sample_counter = 0
    cases = 0
    
    # join the output lines for the genomes
    genome_lines = ''
    for s in sample:
        sample_counter += 1
        a = s[0]
        b = s[1]
        status = ''
        if options.sample_disease:
            status = 'healthy'
            if 'D' in a + b:
                status = 'sick'
                cases += 1
        genome_lines = genome_lines + '\n' + str(sample_counter) + 'p' + ":\t5'" + a + "3' " + status
        genome_lines = genome_lines + '\n' + str(sample_counter) + 'm' + ":\t5'" + b + "3' " + status
    
    # make the header line
    header_line = "GENERATION %i | %i YEARS PASSED | %i POLYMORPHIC SITES IN THE POPULATION | MUT. RATE=%.2E REC. RATE=%.2E" % (generation, generation*18, n_poly, p_mut, p_rec)
    if options.sample_disease:
        header_line = header_line + " | %i SICK INDS. IN THE SAMPLE" % (cases)
    header_line = header_line + '\n'
        
    # make the disease mutation haplotype line
    # and the asterisk to indicate the tag position
    tag_line = ''
    disease_mutation_line = ''
    if options.sample_disease and cases > 0:
        disease_mutation_line = "dm" + ":\t5'" + dm + "3'\n"
        tag_pos = get_best_tag_position(sample)
        tag_line = [' ' for i in dm]
        for i in tag_pos:
            tag_line[i] = '*'
        tag_line = ''.join(tag_line)
        tag_line = disease_mutation_line.replace(dm, tag_line)
        tag_line = tag_line.replace('dm:', '   ')
        tag_line = tag_line.replace('5\'', '  ')
        tag_line = tag_line.replace('3\'', '  ')
    
    print header_line + disease_mutation_line + tag_line + genome_lines

    try:
        cmd = raw_input()
    except KeyboardInterrupt:
        dah.close()
        sys.exit(0)
    if cmd == "exit":
        sys.exit()
    if cmd == "save":
        ofile = open("ld_simulator_save_" + str(generation) + ".txt", 'w')
        ofile.write(op)
        ofile.close()
    disease_fixed = True
    for i in pool:
        if 'D' not in i:
            disease_fixed = False
            break
    if disease_fixed and options.forever == False:
        print "Disease allele was fixed in population, end of simulation!"
        dah.close()
        sys.exit(0)
    
