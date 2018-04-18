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




p1 = R.random_element()
p2 = R.random_element()

v1 = p1.list()
v1 = vector(v1)
v2 = p2.list()
v2 = vector(v2)

# define vandermonde matrix:

M = matrix(ZZ, n)
for i in range(n):
	for j in range(n):
		M[i,j] = powers[mod(i*j,n)]
invM = matrix(ZZ,n)
for i in range(n):
	for j in range(n):
		invM[i,j] = inv_powers[mod(i*j,n)]


# NTT USING MATRIX
def ntt2(v):
	w= M*v
	for i in range(n):
		w[i] = mod(w[i],q)
	return w

def inv_ntt2(v):
	w= invM*v
	for i in range(n):
		w[i] = mod(w[i],q)
	return w


# NTT USING DEFINITION
def ntt(v):
	v = v.list()
	w = vector(ZZ, n)
	for i in range(n):
		w[i] = 0
		for j in range(n):
			w[i] = mod(w[i] + v[j]*powers[mod(j*i,n)], q)

	return w

def inv_ntt(v):
	v = v.list()
	w = vector(ZZ, n)
	for i in range(n):
		w[i] = 0
		for j in range(n):
			w[i] = mod(w[i] + v[j]*inv_powers[mod(j*i,n)], q)	

	return w

def convolution(v, list):
	v = v.list()
	for i in range(n):
		v[i] = mod(v[i]*list[i], q)
	return vector(v)

# CONVERT TO NTT DOMAIN AND MULTIPLY:
def multiply(p1, p2):
	v = ntt2(convolution(p1, rootpowers))
	w = ntt2(convolution(p2, rootpowers))
	z = vector(ZZ, n)
	for i in range(n):
		z[i] = mod(v[i]*w[i], q)
	n3 = inv_ntt2(z)
	for i in range(n):
		n3[i] = mod(n3[i]*inv_n,q)
	return convolution(n3, inv_rootpowers)


def mult_in_nttdomain(p1,p2):
	z = vector(ZZ, n)
	for i in range(n):
		z[i] = p1[i]*p2[i]
	return z

############# schoolbook


def skbk(a, b):
	a = a.list()
	b = b.list()
	c = vector(ZZ,2*n - 1)
	for i in range(n):
		for j in range(i + 1):
			c[i] = c[i] + a[j]*b[i-j]
	for i in range(n, 2*n):
		
		for j in range(1 + i-n, n):
			c[i] = c[i] + a[j]*b[i - j]
			
	d = vector(ZZ, n)
	for i in range(n-1):
		d[i] = c[i] - c[i+n]
	d[n-1] = c[n-1]	

	for i in range(n):
		d[i] = mod(d[i], q)
	return d

#######  test correctness

def test():
		
	w = multiply(v1, v2)
	p3 = p1*p2
	s = skbk(p1,p2)
	
	l = p3.list()
	for i in range(n):
		if s[i]  - w[i] != 0:
			print 'error'
	return 0

####### speed test


def speedtest():
	print 'ntt:'
	timeit('multiply(p1,p2)')
	print 'schoolbook:'
	timeit('skbk(p1,p2)')
	return 0











############## cooley tukey

def butterfly(a, b, root):
	return mod(a + root*b,q)

def cooley_tukey(a):
	a = a.list()
	d = log(n,2)
	
	for i in range(d):
		for j in range(n/(2^(i+1))):
			a[j] = butterfly(a[2*j], a[2*j]+1, omega^(i))

	return a



def recursive(p, root, n):
	if n == 2:
		return p
	v = vector(ZZ,n)
	print p
	print n
	print v
	for i in range(n/2):
		
		print v 
		print p
		v[i] = butterfly(p[2*i], p[2*i+1], root^(log(n,2)))
	print v
	recursive(v,root, n/2)
			

	return 0






