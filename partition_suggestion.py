"""
partition_suggestion.py

purpose: Given Nx, Ny, Nz and a number of processes, suggest a partitioning strategy that would result in
more-or-less cube-shaped partitions

"""
import random
from collections import deque


def random_partition(factors,n_factors):
    """
    factors = list of prime factors of a number
    n_factors = three-tuple [#px, #py, #pz] indicating the number
    of prime factors that should be chosen (at random) for each direction
    
    returns [px,py,pz]
    
    """
    
    l_factors = factors[:] # make a local copy
    
    p_list = [1,1,1]
    
    for d in range(3):
        for i in range(n_factors[d]):
            c = random.choice(l_factors)
            l_factors.remove(c)
            p_list[d]*=c
    return p_list



class Partition:
    def __init__(self,Nx,Ny,Nz,part):
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.px = part[0]
        self.py = part[1]
        self.pz = part[2]
        
        self.set_score()
        
    def set_score(self):
        lx = self.Nx/float(self.px)
        ly = self.Ny/float(self.py)
        lz = self.Nz/float(self.pz)
        
        # obviously, change this measure if it proves crappy.
        
        # estimate surface to volume ratio of the typical partition
        #vol = lx*ly*lz
        #surf = 2.*lx*ly + 2.*lx*lz + 2.*lz*ly
        #self.score = surf/vol
        
        vol = lx*ly*lz
        surf = 2.*lx*ly + 2.*lx*lz + 2.*lz*ly
        interior = vol - surf
        
        self.score = (surf/interior)
        
        
        
        
    def get_score(self):
        return self.score
    
    def get_partition(self):
        return [self.px, self.py, self.pz]



def partitionfunc(n,k,l=1):
    '''n is the integer to partition, k is the length of partitions, l is the min partition element size'''
    if k < 1:
        raise StopIteration
    if k == 1:
        if n >= l:
            yield (n,)
        raise StopIteration
    for i in range(l,n+1):
        for result in partitionfunc(n-i,k-1,i):
            yield (i,)+result


def primes(n):
    """
    return a list containing the prime factorization of positive integer n
    n = positive integer
    
    """
    primfac = []
    d = 2
    while d*d <= n:
        while (n % d) == 0:
            primfac.append(d)  # supposing you want multiple factors repeated
            n //= d
        d += 1
    if n > 1:
       primfac.append(n)
    return primfac

def factors(n):
    """
    n = positive integer
    returns a list of the factors of n in (more-or-less) standard form
    
    """
    return filter(lambda i: n % i == 0, range(1, n + 1))

def part_advisor(Nx,Ny,Nz,num_procs):
    """
    Nx = number of points in the x-direction
    Ny = number of points in the y-direction
    Nz = number of points in the z-direction
    num_procs = the number of partitions to create
    
    returns a suggested px,py,pz
    
    """
        
    p_facts = primes(num_procs)
    p_facts.append(1)
    p_facts.append(1) # to allow effectively 1-D partitioning if that is best....
    
    bestScore = float("inf")
    bestPartition = None
    numTrials = 4000 # not clear to me how big this should be...
    
    for p in partitionfunc(len(p_facts),3):
        #print p
        """
        set up some partitions and keep track of the one with the best score.
        
        """
        p_deque = deque(p);
        for i in range(3):
            p_deque.rotate(1) #shift the groupings
            # take numTrials samples
            for trial in range(numTrials):
                r_part = random_partition(p_facts,p_deque)
                sample_partition = Partition(Nx,Ny,Nz,r_part)
                sample_score = sample_partition.get_score()
                if sample_score < bestScore:
                    bestPartition = Partition(Nx,Ny,Nz,r_part)
            
    
    return bestPartition.get_partition()
        
    """
    partitionfunc will let me generate groupings of the prime factors
    """
    
    """
    if there are fewer than 3 prime factors, then there is no way to solve
    the problem; an error should be returned and the user should be prompted
    to provide a value for num_procs that has more prime factors.
    
    """
    
    if len(p_facts)<3:
        print 'Error!  num_procs is prime and cannot be used for 3D partitioning'
        raise RuntimeError
    
    print p_facts
    
    """
    concept: use the prime factors listed in p_facts
    and put them into 3 groups such that, as close as possible,
    Nx/g1, Ny/g2 and Nz/g3 are nearly equal.
    
    To do this, for each grouping, I will compute the variance of the partition dimensions.
    
    I will then select the grouping that has the lowest variance.
    
    
    1.  Enumerate all of the possible groupings of the prime factors.
    
    2.  Compute the partition dimension variance for each grouping
    
    3.  Pick the smallest one.
    
    """
    
    
if __name__=="__main__":
    """
    write test code here...
    """
    Nx = 150
    Ny = 150
    Nz = 1000
    num_procs = 8
    
    partition = part_advisor(Nx,Ny,Nz,num_procs)
    bestPartition = Partition(Nx,Ny,Nz,partition)
    print 'Best partition found has score = %g \n'%bestPartition.get_score()
    print bestPartition.get_partition()
    
    print 'Block sizes approximately %i x %i x %i'%(Nx/partition[0],Ny/partition[1],Nz/partition[2])
    
    