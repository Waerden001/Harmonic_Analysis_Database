def moment(G,n,rho):
    count = 0
    for ele in G:
        count += rho(ele)^n
    return count/G.order()

def abs_moment(G,n,rho):
    count = 0
    for ele in G:
        count += abs(rho(ele))^n
    return count/G.order()

def decompose(rho,G):
    for k, ele in enumerate(G.irreducible_characters()):
        if rho.scalar_product(ele) != 0:
            print(f"({rho.scalar_product(ele)}, c_{k+1})")