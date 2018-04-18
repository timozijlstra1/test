q = 7681
n = 256

# define rings, calculate powers of root

Z = IntegerModRing(q)
gen = Z.multiplicative_generator()
root = ZZ(mod(gen^((q-1)/(2*n)),q))
inv_root = root.inverse_mod(q)
inv_n = mod(1/n, q)

omega = ZZ(mod(root^2, q))
inv_omega = omega.inverse_mod(q)

R.<x> = PolynomialRing(IntegerModRing(q))
R.<x> = R.quotient(x^n + 1)






