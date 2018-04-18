q = 7681
n = 256


Z = IntegerModRing(q)
gen = Z.multiplicative_generator()
root = ZZ(mod(gen^((q-1)/(2*n)),q))
inv_root = root.inverse_mod(q)
inv_n = mod(1/n, q)
