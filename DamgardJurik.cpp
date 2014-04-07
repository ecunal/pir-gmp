#include "DamgardJurik.h"

using namespace std;

DamgardJurik::DamgardJurik(){

	s = 1;
	enc_only = false;
	bit_length = 1024;

	mpz_inits(n, n_s, n_sp, g, p, q, d, mu, NULL);
	keygen();
}

DamgardJurik::DamgardJurik(int bit_length, int s) {

	this->bit_length = bit_length;
	this->s = s;

	mpz_inits(n, n_s, n_sp, g, p, q, d, mu, NULL);

	enc_only = false;

	keygen();
}

DamgardJurik::DamgardJurik(int bit_length, int s, mpz_t n, mpz_t g) {

	mpz_inits(n, n_s, n_sp, g, NULL);

	mpz_set(this->n, n);
	mpz_set(this->g, g);
	this->s = s;
	this->bit_length = bit_length;

	enc_only = true;

	keygen();
}

void DamgardJurik::keygen() {

	init_random();

	if (!enc_only) {

		// Generate p q, random primes of length bitlength/2

		get_random_prime(p);
		get_random_prime(q);

		mpz_mul(n, p, q);
		mpz_add_ui(g, n, 1);
	}

	mpz_pow_ui(n_s, n, s);
	mpz_pow_ui(n_sp, n, (s + 1));

	if (!enc_only) {

		mpz_t p1, q1;
		mpz_inits(p1, q1, NULL);

		mpz_sub_ui(p1, p, 1);
		mpz_sub_ui(q1, q, 1);

		mpz_mul(d, p1, q1);

		mpz_invert(mu, d, n_s);

		mpz_clears(p1, q1, NULL);

	}

}

void DamgardJurik::init_random() {

	int random_dev = open("/dev/urandom", O_RDONLY);
	long random_seed;
	size_t random_len = 0;

	while (random_len < sizeof random_seed) {

		ssize_t result = read(random_dev, ((char*) &random_seed) + random_len,
				(sizeof random_seed) - random_len);
		if (result < 0) {
			cout << "unable to read from urandom" << endl;
		}
		random_len += result;
	}

	close(random_dev);

	gmp_randinit_mt(state);
	gmp_randseed_ui(state, random_seed);

}

void DamgardJurik::get_random_prime(mpz_t result) {

	mpz_urandomb(result, state, (bit_length / 2));
	mpz_nextprime(result, result);

}

void DamgardJurik::encrypt(mpz_t result, const mpz_t m) {

	mpz_t temp1, temp2, r;
	mpz_inits(temp1, temp2, result, r, NULL);
	mpz_urandomb(r, state, bit_length);

	mpz_powm(temp1, g, m, n_sp);
	mpz_powm(temp2, r, n_s, n_sp);

	mpz_mul(result, temp1, temp2);
	mpz_mod(result, result, n_sp);

	mpz_clears(temp1, temp2, r, NULL);

}

void DamgardJurik::find_i(mpz_t result, const mpz_t c) {

	mpz_t t, t1, t2, nj, nj1, f, i;
	mpz_inits(t, t1, t2, nj, nj1, f, i, NULL);
	mpz_set_ui(i, 0);

	for (int j = 1; j <= s; j++) {

		mpz_pow_ui(nj, n, j);
		mpz_pow_ui(nj1, n, j + 1);

		mpz_mod(t, c, nj1);
		mpz_sub_ui(t, t, 1);
		mpz_tdiv_q(t1, t, n);

		mpz_set(t2, i);

		for (int k = 2; k <= j; k++) {

			mpz_sub_ui(i, i, 1);

			mpz_mul(t2, t2, i);
			mpz_mod(t2, t2, nj);

			mpz_fac_ui(f, k);
			mpz_invert(f, f, nj);

			mpz_t nkm1; // n^(k-1)
			mpz_init(nkm1);
			mpz_pow_ui(nkm1, n, k-1);

			mpz_mul(f, f, nkm1);
			mpz_mod(f, f, nj);

			mpz_mul(f, f, t2);
			mpz_mod(f, f, nj);

			mpz_sub(t1, t1, f);
			mpz_mod(t1, t1, nj);

		}

		mpz_set(i, t1);

	}

}

void DamgardJurik::decrypt(mpz_t result, const mpz_t c) {

	mpz_t temp1;
	mpz_init(temp1);

	mpz_powm(temp1, c, d, n_sp);
	find_i(temp1, temp1);

	mpz_mul(result, temp1, mu);
	mpz_mod(result, result, n_s);

}

void DamgardJurik::set_s(int s) {

	this->s = s;
	mpz_pow_ui(n_s, n, s);
	mpz_mul(n_sp, n, n_s);

	if(!enc_only) {
		mpz_invert(mu, d, n_s);
	}

}

int main() {

	DamgardJurik dj(1024, 1);

	mpz_t plain, cipher;
	mpz_inits(plain, cipher, NULL);

	mpz_set_ui(plain, 15);

	cout << plain << endl;

	dj.encrypt(cipher, plain);

	cout << cipher << endl;

	dj.decrypt(plain, cipher);

	cout << "decrypted: " << plain << endl;

	return 0;
}
