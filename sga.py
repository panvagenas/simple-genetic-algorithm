#!/usr/bin/python3
import argparse
from random import random, randrange, randint
from math import ceil, sqrt
from copy import copy
import sys
import time

parser = argparse.ArgumentParser(description='A Simple Genetic Algorithm by Panagiotis Vagenas')

parser.add_argument('-p', '--popsize', dest='popsize', default=50, type=int,
                    help='GA Population Size')
parser.add_argument('-g', '--maxgens', dest='maxgens', default=100, type=int,
                    help='GA Generations')
parser.add_argument('-n', '--nvars', dest='nvars', default=8, type=int,
                    help='GA N Vars')
parser.add_argument('-a', '--antigenpopsize', dest='antigenpopsize', default=1000, type=int,
                    help='Antigens Population Size')
parser.add_argument('-r', '--runtimes', dest='runtimes', default=1, type=int,
                    help='Number of times GA will be called in order to eval Meta-GA chromosome')
parser.add_argument('-x', '--crossover', dest='crossover', default=0.9, type=float,
                    help='Meta-GA Crossover Probability')
parser.add_argument('-m', '--mutation', dest='mutation', default=0.01, type=float,
                    help='Meta-GA Mutation Probability')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', const=True,
                    help='Verbose info to std output')
parser.add_argument('-o', '--output', dest='output', default='', type=str,
                    help='If provided runtime information will be writen to')

args = parser.parse_args()


class Chromosome:
    def __init__(self):
        self.fitness = 0
        self.gene = ''
        self.mate = False
        self.selection_prob = 0
        self.cumulative_prob = 0


class Population:
    def __init__(self, nvars, pop_size, antigene_pop_size):
        self.best = Chromosome()
        self.population = [Chromosome() for x in range(pop_size)]
        self.nvars = nvars
        self.pop_size = pop_size
        self.history = []
        self.mean_history = []
        self.antigene_pop_size = antigene_pop_size
        self.antigene_dna = ''
        self.evaluation_dict = {}

        for chromosome in self.population:
            chromosome.gene = "".join([str(randint(0, 1)) for j in range(0, nvars)])

        self.antigene_dna = "".join([str(randint(0, 1)) for j in range(0, nvars * self.antigene_pop_size)])

    def evaluate(self):
        mean = 0

        for mem in self.population:
            if mem.gene in self.evaluation_dict:
                mem.fitness = self.evaluation_dict.get(mem.gene)
            else:
                mem_dna = mem.gene * self.antigene_pop_size

                xor = int(mem_dna, 2) ^ int(self.antigene_dna, 2)
                mem.fitness = bin(xor).count("1")

                self.evaluation_dict[mem.gene] = mem.fitness

                if mem.fitness > self.best.fitness:
                    self.best = copy(mem)

            mean += mem.fitness

        self.history.append(self.best.fitness)
        self.mean_history.append(mean / self.pop_size)

        return self.best.fitness

    def crossover(self, pxover):
        lovers = 0

        for chromosome in self.population:
            if random() < pxover:
                chromosome.mate = True
                lovers += 1
            else:
                chromosome.mate = False

        if lovers == self.pop_size and lovers % 2 != 0:
            lovers -= 1

        if lovers % 2 != 0:
            while True:
                chromosome = self.population[randrange(0, self.pop_size)]
                if not chromosome.mate:
                    chromosome.mate = True
                    lovers += 1
                    break

        parent_1 = 0
        couples = int(ceil(lovers / 2))
        for c in range(0, couples):
            while not self.population[parent_1].mate:
                parent_1 += 1

            parent_2 = parent_1 + 1
            while not self.population[parent_2].mate:
                parent_2 += 1

            xover_point = randrange(1, self.nvars)

            gene_1 = self.population[parent_1].gene[:]
            gene_2 = self.population[parent_2].gene[:]

            self.population[parent_1].gene = "".join([gene_1[:xover_point], gene_2[xover_point:]])
            self.population[parent_2].gene = "".join([gene_2[:xover_point], gene_1[xover_point:]])

            parent_1 = parent_2 + 1

    def select_new(self):
        buffered_pop = self.population[:]

        total_fitness = 0
        for member in self.population:
            total_fitness += member.fitness

        for member in self.population:
            member.selection_prob = member.fitness / total_fitness

        self.population[0].cumulative_prob = self.population[0].selection_prob

        for i in range(1, self.pop_size):
            self.population[i].cumulative_prob = self.population[i - 1].cumulative_prob + self.population[
                i].selection_prob

        for spin_num in range(0, self.pop_size):
            roulette = random()

            if roulette <= self.population[0].cumulative_prob:
                buffered_pop[spin_num] = copy(self.population[0])
            else:
                for mem in range(1, self.pop_size):
                    if self.population[i - 1].cumulative_prob < roulette <= self.population[i].cumulative_prob:
                        buffered_pop[spin_num] = copy(self.population[i])
                        break

        self.population = buffered_pop

    def mutate(self, pmutation):
        for member in self.population:
            gene = list(member.gene)
            for i in range(0, len(gene)):
                if random() < pmutation:
                    if gene[i] == '0':
                        gene[i] = '1'
                    else:
                        gene[i] = '0'
            member.gene = "".join(gene)


if args.runtimes < 1:
    sys.stderr.write('Runtimes must be an integer above 0. Your input was "{}"\n'.format(args.runtimes))
    exit(1)

print('Running SGA Algorithm...')

if args.output:
    f = open(args.output, 'w')
else:
    f = False

a = [
    '=' * 20, ' Arguments ', '=' * 20, '\n',
    ('Population Size: %d\n'
     'Generations: %d\n'
     'N Vars: %d\n'
     'Antigens Population Size: %d\n'
     'Runtimes: %d\n'
     'Crossover Probability: %f\n'
     'Mutation Probability: %f\n' %
     (args.popsize, args.maxgens, args.nvars, args.popsize, args.runtimes, args.crossover, args.mutation)),
    '=' * 51, '\n'
]
s = "".join(a)

if f:
    f.write(s)

if args.verbose:
    print(s)

mean = 0

pop = False

for r in range(0, args.runtimes):
    pop = Population(args.nvars, args.popsize, args.antigenpopsize)
    best = Chromosome()
    pop.evaluate()

    if args.runtimes > 1 and f:
        s = ''.join([
            '\n\n',
            '=' * 20,
            (' Round %d ' % r),
            '=' * 20,
            '\n\n'
        ])

        if f:
            f.write(s)

    if args.runtimes == 1:
        print("".join(['=' * 20, ' Initial Population ', '=' * 20]))
        for mem in pop.population:
            print(mem.gene)

        print("".join(['=' * 20, ' Antigens ', '=' * 20]))
        for i in range(0, pop.antigene_pop_size, args.nvars):
            print(pop.antigene_dna[i:i + args.nvars])

    if f:
        f.write(("\n%10s  %15s  %15s  %15s  %20s\n"
                 % ("Generation", "Best Value", "Average", "StdDev", "Best Genotype")))

    for generation in range(0, args.maxgens):
        pop.select_new()
        pop.crossover(args.crossover)
        pop.mutate(args.mutation)
        pop.evaluate()

        if args.runtimes > 1:
            sys.stdout.write('Round %d / %d, Generation %d / %d\r'
                             % (r + 1, args.runtimes, generation + 1, args.maxgens))
        else:
            sys.stdout.write('Generation %d / %d\r' % (generation + 1, args.maxgens))
        sys.stdout.flush()

        if f:
            sum_p = 0
            sum_square = 0
            for m in pop.population:
                sum_p += m.fitness
                sum_square += m.fitness ** 2
            square_sum = sum_p * sum_p / pop.pop_size
            dev = sqrt((1.0 / (pop.pop_size - 1)) * (sum_square - square_sum))
            f.write("\n%10d  %15.4f  %15.4f  %15.4f %20s" %
                    (generation + 1, pop.best.fitness, pop.mean_history[generation], dev, pop.best.gene))

    mean += pop.best.fitness

    s = ''.join([
        ('\nBest Member: %s\n' % pop.best.gene),
        ('Fitness: %f\n\n' % pop.best.fitness)
    ])
    if args.verbose:
        print(s)

s = '\n\nSimulation Completed...\n'
print(s)
if f:
    f.write(s)
    if args.runtimes > 1:
        f.write('Mean fitness: %f' % (mean / args.runtimes))
    elif pop:
        f.write('Best Member: %s\n' % pop.best.gene)
        f.write('Fitness: %f' % pop.best.fitness)

if args.runtimes > 1:
    print('Mean fitness: %f' % (mean / args.runtimes))
else:
    print(''.join([
        ('Best Member: %s\n' % pop.best.gene),
        ('Fitness: %f' % pop.best.fitness)
    ]))
pass

