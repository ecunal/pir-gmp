#include "DamgardJurik.h"

DamgardJurik::DamgardJurik(int bit_length, int s) {
	
	this->bit_length = bit_length;
	this->s = s;
	keygen(false);
}

DamgardJurik::DamgardJurik(int bit_length, int s, mpz_t n, mpz_t g) {
	
	mpz_set(this->n, n);
	mpz_set(this->g, g);
	this->s = s;
	this->bit_length = bit_length;
	
	keygen(true);
}

void DamgardJurik::keygen(bool pub) {
	
	
	
}