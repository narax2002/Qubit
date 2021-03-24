#ifndef Qubit
#include <iostream>
#include <complex>
#include <math.h>
using namespace std;

const double pi = 4.0 * atan(1.0);

class Qubit
{
private:
//public:
	int n;
	int size;
	complex<double>* q;
public:
	Qubit();
	Qubit(const int n);
	void Initial();
	friend ostream& operator<<(ostream& os, const Qubit& Q);
	double* Qnorm();
	void PrintQnorm();
	void PrintError(const int nv);

	// One-qubit gates
	void H(const int idx); // Hadamard (H) gate
	void X(const int idx); // Pauli-X gate
	void SRX(const int idx); // Square Root of X
	void Y(const int idx); // Pauli-Y gate
	void Z(const int idx); // Pauli-Z gate
	void S(const int idx);
	void Sd(const int idx); // S dagger
	void T(const int idx);
	void Td(const int idx); // T dagger
	void R(const int idx, const double phi); // Phase shift gate
	

	// Multi-qubit gates
	void SWAP(const int idx_a, const int idx_b); // swap gate
	void CX(const int idx_a, const int idx_b); // Controlled X or Controlled NOT
	void CZ(const int idx_a, const int idx_b); // Controlled Z
	void CCX(const int idx_a, const int idx_b, const int idx_c); // Controlled Controlled X or Controlled Controlled NOT

	void Toffoli(const int idx); // Toffoli

	void CR(const int idx_a, const int idx_b, const double phi); // Contrulled Phase shift gate


	// Algorithms
	void Grover(const int k); // Grover
	void Bernstein_Vazirani(const int s); // Bernstein-Vazirani
	void Shor(); // Shor

	void QFT();
	void FFT();
	void IFFT();
	void FFT_iter(); // Not use
	void Multi(Qubit& q_a, Qubit& q_b);
};
#endif
