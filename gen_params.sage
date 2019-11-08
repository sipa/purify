import sys

if len(sys.argv) != 2:
    print("Usage: sage %s FIELD_SIZE" % __file__)
    sys.exit()

N = int(sys.argv[1])


if not is_prime(N):
    print("Field size is not prime: %i" % N)
    sys.exit()

F = GF(N)

# Find a non-square D in F
for D in range(1, 1000000):
    if not F(D).is_square():
        break
D = F(D)
D2 = D * D
D3 = D * D * D

def embedding_degree(Curve):
    size = Curve.coordinate_ring().base_ring().order()
    order = Curve.order()
    Fg = GF(order)
    x = Fg(size)
    for degree in divisors(order-1):
        if (x**degree == Fg(1)):
            return (order-1)//degree
    return -1

sum_a_b = 1
iter = 0
while True:
    print("Iteration %i..." % iter)
    for a_val in range(0, sum_a_b):
        iter += 1
        a = F(a_val)
        b = F(sum_a_b - a_val)

        # Sanity check for non-singular curve
        if (4 * a * a * a + 27 * b * b) == 0:
            continue

        # Preliminary analysis on E1: y^2 = x^3 + a*x + b
        E1 = EllipticCurve(F, [a, b])
        n1 = E1.order()
        if not is_pseudoprime(n1):
            continue

        # Preliminary analysis on E2: y^2 = x^3 + a*x + b
        E2 = EllipticCurve(F, [a * D2, b * D3])
        n2 = E2.order()
        if not is_pseudoprime(n2):
            continue

        # Full primarily test on both
        if not is_prime(n1) or not is_prime(n2):
            continue

        print("n = %i # Field size" % N)
        print("a = %i # curve equation parameter a" % a)
        print("b = %i # curve equation parameter b" % b)
        print("d = %i # non-square in GF(n)" % D)
        print("n1 = %i # Order of E1: y^2 = x^3 + a*x + b over GF(n)" % n1)
        print("n2 = %i # Order of E2: y^2 = x^3 + a*d^2*x + b*d^3 over GF(n)" % n2)
        print("e1 = (n1 - 1) / %i # Embedding degree of E1" % embedding_degree(E1))
        print("e2 = (n2 - 1) / %i # Embedding degree of E2" % embedding_degree(E2))
        sys.exit()
    sum_a_b += 1

