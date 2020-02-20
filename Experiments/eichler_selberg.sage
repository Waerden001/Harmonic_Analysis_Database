G = DirichletGroup(23)
eps = G.list()[1]

def psi(n):
    return n*prod([1+1/p for p in prime_divisors(n)])

def central(n, weight, level, character):
    for k in range(int(sqrt(n)), int(sqrt(n))+1):
        if k**2 == n:
            return (weight-1)/12*psi(level)*(1/character(k))*n**(weight/2-1)
    return 0

#Trace = central(n=4, weight = 4, level = 25, character = eps)
def symmetric(t,n, k):
    """Compute x^k-x_bar^k/(x-x_bar) for roots of x^2-tx+n
    
    """
    x = t+sqrt(t**2-4*n)
    x_bar = t-sqrt(t**2-4*n)
    
    return (x**k-x_bar**k)/(x-x_bar)


def mu_function(t,m, n, N, character):
    N_m = gcd(N,m)
    c_range = []
    for c in range(N):
        if gcd(c,N) == 1:
            #print(c)
            k=0
            while(c+k*N<N*N_m):
                if (c+k*N)**2-t*(c+k*N)+n % N*N_m == 0:
                    c_range.append(c)
                k +=1
    print(c_range)
                
    return psi(N)/psi(N_m)*sum([1/character(c) for c in c_range])

def weighted_class_number(t,n):
    K.<a> = NumberField(x^2 -t*x + n)
    if K.discriminant() == -4:
        return K.class_number()/2
    if K.discriminant() == -3:
        return K.class_number()/3
    else return K.class_number()

def hyper_unipotent(n, N, k, w):
    ans = 0
    for d in divisors(n):
        tau_range = []
        for tau in divisors(N):
            if N/w.conductor % gcd(N, int(N/tau))and (d-int(n/d)) % gcd(N, N/tau)==0:
                tau_range.append(tau) 
        #TODO: simplify to save time, one loop is already enough
        aux = 0
        for tau in tau_range:
            common_divisor = gcd(tau, int(N/tau))
            y = CRT([d, int(n/d)], [tau, int(N/tau)]) % int(N/common_divisor)
            #TODO: decorate for special case gcd(n, N)>1
            aux += euler_phi(common_divisor)*1/w(y)
        ans +=min(d, n/d)**(k-1)* aux
    return ans

def A_2(k, n, N, w):
    if k == 2 and w.is_trivial():
        return sum([c for c in divisors(n) and gcd(N, int(n/c))==1])
    else return 0

def elliptic(n, k):
    lower_bound = ceil(-2*sqrt(n))
    upper_bound = floor(2*sqrt(n))
    for t in range(lower_bound, upper_bound+1):
        aux = 0
        for m in range(1, t**2-4*n):
        #TODO: sqrt of the range is enough
            if t**2-4*n % m**2 == 0:
                aux += weighted_class_number(int((t**2-4*n)/m**2))*mu_function(t, m, n)
        ans += symmetric(t,n,k)*aux
    return ans
    
def eichler_selberg(n, k, N, w):
    """Compute the trace of T_n---the n-th Hecke operator on S_k(N, w), weight k, level N, cusp forms with Dirichlet character w.
    
    :param: n, int, the n-th Hecke operator
    :param: k, int, the weight of the space of cusp forms, need to be an even number.
    :param: N, int, the level of the space of cusp forms 
    :param: w, DirichletCharacter, one Dirichlet character of (Z/NZ)*-->C* 
    
    """
    reuturn central(n, k, N, w) + elliptic(k,n) + hyper_uniponent() + A_2(k,n,N,w)



#####TEST#########################################################################
