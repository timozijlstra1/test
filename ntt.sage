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

powers = vector(ZZ,n)
inv_powers = vector(ZZ,n)

rootpowers = vector(ZZ,n)
inv_rootpowers = vector(ZZ,n)

for i in range(n):
	powers[i] = mod(omega^i, q)
	inv_powers[i] = mod(inv_omega^i, q)
	rootpowers[i] = mod(root^i, q)
	inv_rootpowers[i] = mod(inv_root^i, q)








