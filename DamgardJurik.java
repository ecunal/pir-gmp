import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Random;

public class DamgardJurik {

	private BigInteger p, q, d, mu;

	/**
	 * n = p*q, n_s = n^s, n_sp = n^(s+1)
	 */
	public BigInteger n, n_s, n_sp, g;

	private int bitLength;
	public int s;
	
	public long keyGenerationTime;

	// private BigInteger[] invFact; // unused

	public int getS() {
		return s;
	}

	public void setS(int s) {
		this.s = s;
		this.n_s = n.pow(s);
		this.n_sp = n.multiply(n_s);

		if (d != null)
			this.mu = d.modInverse(n_s);
	}

	/**
	 * 
	 * Constructs an instance of the DJ cryptosystem with specified bitLength
	 * 
	 * @param bitLength
	 * @param s
	 */
	public DamgardJurik(int bitLength, int s) {
		keyGeneration(bitLength, s);
	}

	/**
	 * Constructs an instance of the DJ cryptosystem with 512 bits of modulus
	 * 
	 * @param s
	 *            Increasing s in n^s
	 */
	public DamgardJurik(int s) {
		keyGeneration(512, s);
	}

	/**
	 * Initializing a system with known public keys and 512 bit modulus. Only for
	 * encryption, cannot perform decryption.
	 * 
	 * @param s
	 * @param n
	 * @param g
	 */
	public DamgardJurik(int s, BigInteger n, BigInteger g) {
		keyGeneration(bitLength, s, n, g);
	}

	/**
	 * Known public keys with specified modulus length. Only for encryption,
	 * cannot perform decryption.
	 * 
	 * @param bitLength
	 * @param s
	 * @param n
	 * @param g
	 */
	public DamgardJurik(int bitLength, int s, BigInteger n, BigInteger g) {
		keyGeneration(512, s, n, g);
	}

	/**
	 * Sets up the system according to known public keys.
	 * 
	 * @param bitLength
	 * @param s
	 * @param n
	 * @param g
	 */
	public void keyGeneration(int bitLength, int s, BigInteger n, BigInteger g) {

		this.n = n;
		this.g = g;
		this.bitLength = bitLength;
		this.s = s;

		n_s = n.pow(this.s);
		n_sp = n.multiply(n_s);
	}

	/**
	 * Sets up the public key and private key.
	 * 
	 * @param bitLength
	 *            Number of bits of modulus
	 * @param s
	 */
	public void keyGeneration(int bitLength, int s) {

		this.bitLength = bitLength;
		this.s = s;

		long start = System.currentTimeMillis();

		p = BigInteger.valueOf(2).pow(this.bitLength / 2);
		p = p.nextProbablePrime();

		q = p.nextProbablePrime();
		q = q.nextProbablePrime();

/*		System.out.println(p + "\n\n" + q); */
		
		/* new prime generation codes ***/
		
	//	p = BigInteger.probablePrime(bitLength / 2, new SecureRandom());
	//	q = BigInteger.probablePrime(bitLength / 2, new SecureRandom());

		n = p.multiply(q);

		n_s = n.pow(this.s);
		n_sp = n.multiply(n_s);

		g = n.add(BigInteger.ONE);
		d = p.subtract(BigInteger.ONE).multiply(q.subtract(BigInteger.ONE)); // (p-1)*(q-1)
		mu = d.modInverse(n_s); // mu = d^-1 mod n^s

		long end = System.currentTimeMillis();

//		System.out.println("Key generation time for modulus size: " + bitLength
//				+ " and s: " + s + " is: " + (end - start));
		
		keyGenerationTime = end-start;

	}
	
	public static DamgardJurik copy(DamgardJurik o) {
		
		DamgardJurik dj = new DamgardJurik(o.s, o.n, o.g);
		dj.p = o.p;
		dj.q = o.q;
		dj.d = o.d;
		dj.mu = o.mu;
		
		return dj;
	}

	public BigInteger encryption(BigInteger m, BigInteger r) {

		return g.modPow(m, n_sp).multiply(r.modPow(n_s, n_sp)).mod(n_sp);
	}

	public BigInteger encryption(BigInteger m) {

		Random rand = new Random();
		BigInteger r = BigInteger.valueOf(rand.nextInt());

		// long begin = System.currentTimeMillis();
		BigInteger temp1 = g.modPow(m, n_sp);
		// long end = System.currentTimeMillis();

		// System.out.println("g^m mod n^s+1: " + (end - begin));

		// begin = System.currentTimeMillis();
		BigInteger temp2 = r.modPow(n_s, n_sp);
		// end = System.currentTimeMillis();

		// System.out.println("r ussu n^s mod n^s+1: " + (end - begin));

		return temp1.multiply(temp2).mod(n_sp);
	}

	private BigInteger factorial(int k) {

		BigInteger res = BigInteger.valueOf(k);

		k--;

		while (k > 1) {
			res = res.multiply(BigInteger.valueOf(k));
			k--;
		}

		return res;
	}

	private BigInteger find_i(BigInteger c) {

		BigInteger t, t1, t2, nj, f, i = BigInteger.ZERO;

		System.out.println("inside find_i");

		for (int j = 1; j <= s; j++) {

			System.out.println("j: " + j);

			t = c.mod(n.pow(j + 1));
			System.out.println("t: " + t + "\n");

			t1 = t.subtract(BigInteger.ONE).divide(n);
			System.out.println("t1: " + t1 + "\n");

			t2 = i;
			System.out.println("t2: " + t2 + "\n");

			nj = n.pow(j);
			System.out.println("nj: " + nj + "\n");

			for (int k = 2; k <= j; k++) {

				System.out.println("k: " + k);

				i = i.subtract(BigInteger.ONE);
				System.out.println("i: " + i + "\n");

				t2 = t2.multiply(i).mod(nj);
				System.out.println("t2: " + t2 + "\n");

				f = factorial(k);
				f = f.modInverse(nj);
				System.out.println("f inverse: " + f + "\n");

				t1 = t1.subtract(n.pow(k - 1).multiply(f).multiply(t2)).mod(nj);
				System.out.println("t1: " + t1 + "\n");
			}
			i = t1;
		}

		return i;
	}

	public BigInteger decryption(BigInteger c) {

		BigInteger temp = c.modPow(d, n_sp);
		System.out.println("temp: " + temp);
		temp = find_i(temp);

		return temp.multiply(mu).mod(n_s);
	}

	public BigInteger getPublicKeyN() {
		return n;
	}

	public BigInteger getPublicKeyG() {
		return g;
	}
}
