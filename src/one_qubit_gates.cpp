#include "qubit.h"

/// <summary>
/// one-qubit gate 정의
/// H: Hadamard, X: Pauli X, Y: Pauli Y, Z: Pauli X
/// S, T, R: Phase shift
/// </summary>

void Qubit::H(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double>* temp = new complex<double>[len];
	//complex<double> temp[len];
	for (int i = 0; i < len; ++i) {
		if ((i & ref) == 0) {
			int j = i | ref;
			temp[i] = (this->q[i] + this->q[j]) / sqrt(2);
			temp[j] = (this->q[i] - this->q[j]) / sqrt(2);
		}
	}
	delete this->q;
	this->q = temp;
}

void Qubit::X(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	for (int i = 0; i < len; ++i) {
		if ((i & ref) == 0) {
			int j = i | ref;
			complex<double> temp = this->q[i];
			this->q[i] = this->q[j];
			this->q[j] = temp;
		}
	}
}

void Qubit::SRX(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double> c_a(0.5, 0.5);
	complex<double> c_b(0.5, -0.5);
	for (int i = 0; i < len; ++i) {
		if ((i & ref) == 0) {
			int j = i | ref;
			complex<double> temp_a = this->q[i];
			complex<double> temp_b = this->q[j];
			this->q[i] = c_a * temp_a + c_a * temp_b;
			this->q[j] = c_b * temp_a + c_b * temp_b;
		}
	}
}

void Qubit::Y(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double> c(0.0, 1.0);
	for (int i = 0; i < len; ++i) {
		if ( !(i & ref)) {
			int j = i | ref;
			complex<double> temp = this->q[i]; // idx번째 큐빗이 |0>
			this->q[i] = c * this->q[j];
			this->q[j] = - c * temp;
		}
	}
}

void Qubit::Z(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	for (int i = 0; i < len; ++i) {
		if (i & ref) {
			this->q[i] *= -1.0;
		}
	}
}

void Qubit::S(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double> c(0.0, 1.0);
	for (int i = 0; i < len; ++i) {
		if (i & ref) {
			this->q[i] *= c;
		}
	}
}

void Qubit::Sd(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double> c(0.0, -1.0);
	for (int i = 0; i < len; ++i) {
		if (i & ref) {
			this->q[i] *= c;
		}
	}
}

void Qubit::T(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double> c(1.0/sqrt(2), 1.0/sqrt(2));
	for (int i = 0; i < len; ++i) {
		if (i & ref) {
			this->q[i] *= c;
		}
	}
}

void Qubit::Td(const int idx)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double> c(1.0 / sqrt(2), -1.0 / sqrt(2));
	for (int i = 0; i < len; ++i) {
		if (i & ref) {
			this->q[i] *= c;
		}
	}
}

void Qubit::R(const int idx, double phi)
{
	if (idx >= this->n) this->PrintError(1);
	int len = this->size;
	int ref = 1 << idx;
	complex<double> c = polar(1.0, phi);
	for (int i = 0; i < len; ++i) {
		if (i & ref) {
			this->q[i] *= c;
		}
	}
}