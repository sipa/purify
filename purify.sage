import hashlib
import sys

# Parameters generated using gen_params.sage for secp256k1
P = 115792089237316195423570985008687907852837564279074904382605163141518161494337
A = 118
B = 339
D = 5
N1 = 115792089237316195423570985008687907853146579067639158218940405176378157516777
N2 = 115792089237316195423570985008687907852528549490510650546269921106658165471899

F.<f> = GF(P)
assert((N1 + N2) == 2*(P + 1))
FF.<SQRT_D> = GF(P*P, modulus=x^2 - D)
EE = EllipticCurve(FF, [A,B])
C1, C2 = crt_basis([N1, N2])

def hash_to_int(data, rang):
    bits = rang.nbits()
    mask = 2 ** bits - 1
    for i in range(256):
        out = "".join(hashlib.sha256(chr(i) + chr(j) + data).digest() for j in range((bits + 255) // 256))
        val = int(out.encode("hex"), 16) & mask
        if val < rang:
            return val

def hash_to_curve(data, p, a, b):
    for i in range(256):
        x = F(hash_to_int(chr(i) + data, p))
        c = x*(x*x + a) + b
        if c.is_square():
            return (x, c.sqrt())

def hash_to_ext_curve(data):
    p1 = hash_to_curve(data + "/1", P, A, B)
    p2 = hash_to_curve(data + "/2", P, A*D*D, B*D*D*D)
    return EE(p1[0], p1[1]) + EE(p2[0] / D, p2[1] / (D * SQRT_D))

def decode_secret(z):
    """Convert a secret from packed integer format to extension curve scalar."""
    z1 = 1 + (z % ((N1 - 1) // 2))
    z2 = 1 + (z // ((N1 - 1) // 2))
    return z1 * C1 + z2 * C2

def encode_public(pp):
    """Convert an extension curve point to a packed x coordinate pair format."""
    p1 = C1 * pp
    p2 = C2 * pp
    x1 = int(p1[0])
    x2 = int(p2[0] * D)
    return x1 + x2 * P

def extract(pp):
    """Extract the constant term of the X coordinate of an extension curve point."""
    return pp[0].polynomial()[0]

GG = hash_to_ext_curve("Generator")

if len(sys.argv) < 2:
    print("Usage: %s gen <seckey>: generate a public key" % __file__)
    print("       %s eval <seckey> <hexmsg>: evaluate the PRF" % __file__)
elif sys.argv[1] == "gen":
    z = decode_secret(int(sys.argv[2], 16))
    print("x=%x # public key" % encode_public(z * GG))
elif sys.argv[1] == "eval":
    z = decode_secret(int(sys.argv[2], 16))
    m = bytearray.fromhex(sys.argv[3])
    print("eval: %x" % extract(z * hash_to_ext_curve("Eval/" + m)))
else:
    print("Unknown command")
