# Purify: A provably-secure PRF with low multiplicative complexity

Purify is a [pseudorandom function](https://en.wikipedia.org/wiki/Pseudorandom_function_family) (PRF) onto the integers
modulo a prime *P*. It does rely on another hash function under the hood, but one that is only applied to the message.
When that hash function is modeled as a random oracle, Purify is provably indistuishable from random assuming the
Elliptic Curve Decision Diffie-Hellman (ECDDH) problem is hard for its chosen parameters.

Furthermore, for any fixed message, Purify can be implemented as an arithmetic circuit that uses less
than *16L/3* multiplication gates, where *L = &lceil;log<sub>2</sub>(P)&rceil;*.

Lastly, the keys used in Purify have a corresponding public key, and Purify remains indistinguishable from random to
an attacker who knows the public key. This means it can function as a VRF, but more generally, it can be used in
zero-knowledge proofs where the output of the function is not known to the verifier. This is specifically what is needed for proving
deterministic nonces for Schnorr signatures over a group of order *P* were computed correctly. In this case,
a single proof that both evaluates Purify for a public message and a secret key, and proves the secret key
corresponds to a known public key, can be written as an arithmetic circuit using no more than *8L - 14* multiplication gates.

## Parameters

* Let *P* be an *L*-bit prime.
* Let *D* be a non-square in *GF(P)*.
* Let *A* and *B* be elements of *GF(P)* such that:
  * The elliptic curve *E<sub>1</sub>* over *GF(P)* with equation *y<sup>2</sup> = x<sup>3</sup> + Ax + B* has prime order *N<sub>1</sub>*.
  * The elliptic curve *E<sub>2</sub>* over *GF(P)* with equation *y<sup>2</sup> = x<sup>3</sup> + AD<sup>2</sup>x + D<sup>3</sup>B* has prime order *N<sub>2</sub>*.
  * The ECDDH problem is assumed to be hard in both *E<sub>1</sub>* and *E<sub>2</sub>*.
* Let *H<sub>1</sub>* be a hash function onto *E<sub>1</sub>* excluding the point at infinity.
* Let *H<sub>2</sub>* be a hash function onto *E<sub>2</sub>* excluding the point at infinity.
* Let *G<sub>1</sub>* be a generator for *E<sub>1</sub>*.
* Let *G<sub>2</sub>* be a generator for *E<sub>2</sub>*.

The tool <code>gen_params.sage</code> can be used to find appropriate values for *A*, *B*, *D*, *N<sub>1</sub>*, and *N<sub>2</sub>* for a given *P*, for example:

    $ sage gen_params.sage 1000000007
    P = 1000000007 # Field size
    A = 17 # curve equation parameter A
    B = 13 # curve equation parameter B
    D = 5 # non-square in GF(P)
    N1 = 999956519 # Order of E1: y^2 = x^3 + A*x + B over GF(P)
    N2 = 1000043497 # Order of E2: y^2 = x^3 + A*D^2*x + B*D^3 over GF(P)
    # E1 = (N1 - 1) / 1 # Embedding degree of E1
    # E2 = (N2 - 1) / 3 # Embedding degree of E2

For values of P that are sufficiently large for cryptographic purposes (256 bits and larger), this
may take several days. For primes equal to common group orders, see the Example Parameters section below.

## Keys

(Private) keys for Purify are tuples *(z<sub>1</sub>, z<sub>2</sub>)*, where *z<sub>1</sub>* is a non-zero integer below *(N<sub>1</sub> + 1) / 2* and
*z<sub>2</sub>* is a non-zero integer below *(N<sub>2</sub> + 1) / 2*. They can be represented as a single integer *z = z<sub>1</sub> - 1 + (z<sub>2</sub> - 1)(N<sub>1</sub> - 1) / 2*.

The corresponding public keys (if desired) are tuples *(x<sub>1</sub>, x<sub>2</sub>)*, where *x<sub>1</sub>* is the X coordinate of *z<sub>1</sub>G<sub>1</sub>* and
*x<sub>2</sub>* is the X coordinate of *z<sub>2</sub>G<sub>2</sub>*. They can be represented as a single integer *x = x<sub>1</sub> + P x<sub>2</sub>*.

The <code>purify.py</code> tool can be used to compute these. **Note that this code is purely for demonstration purposes.**

    $ ./purify.py gen
    z=11427c7268288dddf0cd24af3d30524fd817a91e103e7e02eb28b78db81cb350b3d2562f45fa8ecd711d1becc02fa348cf2187429228e7aac6644a3da2824e93 # private key
    x=9343f981e9c40546061e63f9f4e6f61541c483c8aae8fe27180c490f0faf584d5036a5952b01200d8b0fdb49c83d5f8dcc8ae434e77785c576720d18897bbea5 # public key

## Formula

For a message *m* and key *(z<sub>1</sub>, z<sub>2</sub>)*, *Purify((z<sub>1</sub>, z<sub>2</sub>), m)* can be computed as follows:
* Let *u* be the X coordinate of *z<sub>1</sub>H<sub>1</sub>(m)*.
* Let *v* be the X coordinate of *z<sub>2</sub>H<sub>2</sub>(m)*.
* Return *((u + v / D)(A + u v / D) + 2B) / (u - v / D)<sup>2</sup>*.

The <code>purify.py</code> tool implements this:

    $ ./purify.py eval 11427c7268288dddf0cd24af3d30524fd817a91e103e7e02eb28b78db81cb350b3d2562f45fa8ecd711d1becc02fa348cf2187429228e7aac6644a3da2824e93 01234567
    eval: afae82108c66397451ce376bc95751c398e40eaf8c768d1b18cc9dd4161cee35

## Verification using arithmetic circuits

The <code>purify.py</code> can also construct arithmetic circuits that verify the Purify evaluation as well as correctness of public keys. Specifically:

    $ ./purify.py verifier 01234567 >verifier.py

This generates a Python function <code>verifier(pubkey, output, v)</code> that takes as input the *x* value from above, the output from the evaluation, and
an assignment for all of the circuit's secret variables. It is specific for the message <code>01234567</code> in this case.
The function simply contains a number of assert statements that each verify a relation that must hold between these values.

Input for the verifier can be generated using:

    $ ./purify.py prove 01234567 11427c7268288dddf0cd24af3d30524fd817a91e103e7e02eb28b78db81cb350b3d2562f45fa8ecd711d1becc02fa348cf2187429228e7aac6644a3da2824e93 >proof.py

It invokes the verifier with expected values:

    $ cat proof.py | cut -d , -f 1-2
    verify(0x9343f981e9c40546061e63f9f4e6f61541c483c8aae8fe27180c490f0faf584d5036a5952b01200d8b0fdb49c83d5f8dcc8ae434e77785c576720d18897bbea5, 0xafae82108c66397451ce376bc95751c398e40eaf8c768d1b18cc9dd4161cee35

These are indeed the public key and the evaluation. The third argument to <code>verifier</code> is a long list of secret variables. The proof can be verified too:

    $ cat verifier.py proof.py | python3

**Note that this does not actually implement any zero-knowledge proofs. It only derives the relations that would need to be proven, and the secret values they're over in specific instances.**

## Example parameters

The code in this repository has parameters that correspond to the order of the secp256k1 group:

    P = 115792089237316195423570985008687907852837564279074904382605163141518161494337
    A = 118
    B = 339
    D = 5
    N1 = 115792089237316195423570985008687907853146579067639158218940405176378157516777
    N2 = 115792089237316195423570985008687907852528549490510650546269921106658165471899

Using this 256-bit prime results in verification circuits that have 2030 multiplication gates.

Other target groups can be configured by modifying purify.py. Other parameters are available in comments:

    # Parameters generated using gen_params.sage for Curve25519 (253 bits)
    P = 0x1000000000000000000000000000000014DEF9DEA2F79CD65812631A5CF5D3ED
    A = 95
    B = 78
    D = 2
    N1 = 0x100000000000000000000000000000004E9C306B81CF1C611587B3ED91288DAD
    N2 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDB21C351C4201D4B9A9D124728C31A2F

    # Parameters generated using gen_params.sage for BLS12-381 (255 bits)
    P = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
    A = 245
    B = 46
    D = 5
    N1 = 0x73EDA753299D7D483339D80809A1D804942105BA15136AAC92458EF0CDB43949
    N2 = 0x73EDA753299D7D483339D80809A1D806135A424BEAE94D516DBA710D324BC6BB

    # Parameters generated using gen_params.sage for BN(2,254) (254 bits)
    P = 0x2523648240000001BA344D8000000007FF9F800000000010A10000000000000D
    A = 209
    B = 140
    D = 2
    N1 = 0x2523648240000001BA344D80000000089C9DDF8B4198211E1005BEF4E673BA39
    N2 = 0x2523648240000001BA344D800000000762A12074BE67DF0331FA410B198C45E3

    # Parameters generated using gen_params.sage for Ed448-Goldilocks (446 bits)
    P = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF7CCA23E9C44EDB49AED63690216CC2728DC58F552378C292AB5844F3
    A = 155
    B = 199
    D = 2
    N1 = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF61E19CF8AE93A7F6204DD85972E93B7A4C4733D057799E70F578D05B
    N2 = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF97B2AADADA0A0E9D3D5E94C6CFF0496ACF43EAD9EF77E6B46137B98D
