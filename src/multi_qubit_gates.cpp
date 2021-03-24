#include "qubit.h"

/// <summary>
/// multi-qubit gate ����
/// SWAP, CX: Controlled X, CZ: Controlled Z
/// CCX: Controlled Controlled X
/// </summary>

void Qubit::SWAP(const int idx_a, const int idx_b)
{
	if (this->n < 2) this->PrintError(2);
	if (idx_a >= this->n || idx_b >= this->n) this->PrintError(1);

	int len = this->size;
	int ref_a = 1 << idx_a;
	int ref_b = 1 << idx_b;
	for (int i = 0; i < len; ++i) {
		if ((i & ref_a) && ((i & ref_b) == 0)) {
			//int j = i - ref_a + ref_b;
			int j = i ^ (ref_a + ref_b);
			complex<double> temp = this->q[i];
			this->q[i] = this->q[j];
			this->q[j] = temp;
		}
	}
}

void Qubit::CX(const int idx_a, const int idx_b)
{
	if (this->n < 2) this->PrintError(2);
	if (idx_a >= this->n || idx_b >= this->n) this->PrintError(1);

	int len = this->size;
	int ref_a = 1 << idx_a;
	int ref_b = 1 << idx_b;
	for (int i = 0; i < len; ++i) {
		if (((i & ref_a) == 0) && (i & ref_b)) {
			int j = i | ref_a;
			complex<double> temp = this->q[i];
			this->q[i] = this->q[j];
			this->q[j] = temp;
		}
	}
}

void Qubit::CZ(const int idx_a, const int idx_b)
{
	if (this->n < 2) this->PrintError(2);
	if (idx_a >= this->n || idx_b >= this->n) this->PrintError(1);

	int len = this->size;
	int ref_a = 1 << idx_a;
	int ref_b = 1 << idx_b;
	complex<double> c(0.0, 1.0);
	for (int i = 0; i < len; ++i) {
		if (((i & ref_a) == 0) && (i & ref_b)) {
			int j = i | ref_a;
			complex<double> temp = this->q[i];
			this->q[i] = c * this->q[j];
			this->q[j] = -c * temp;
		}
	}
}

void Qubit::CCX(const int idx_a, const int idx_b, const int idx_c)
{
	if (this->n < 3) this->PrintError(2);
	if (idx_a >= this->n || idx_b >= this->n || idx_c >= this->n) this->PrintError(1);

	int len = this->size;
	int ref_a = 1 << idx_a;
	int ref_b = 1 << idx_b;
	int ref_c = 1 << idx_c;
	for (int i = 0; i < len; ++i) {
		if (((i & ref_a) == 0) && (i & ref_b) && (i & ref_c)) {
			int j = i | ref_a;
			complex<double> temp = this->q[i];
			this->q[i] = this->q[j];
			this->q[j] = temp;
		}
	}
}

void Qubit::Toffoli(const int idx)
{
	if (this->n < 2) this->PrintError(2);
	int len = this->size;
	int ref = 1 << idx;
	
	int i = len - 1;
	int j = i ^ ref;
	complex<double> temp = this->q[i];
	this->q[i] = this->q[j];
	this->q[j] = temp;

}

void Qubit::CR(const int idx_a, const int idx_b, const double phi)
{
	if (this->n < 2) this->PrintError(2);
	int len = this->size;
	int ref_a = 1 << idx_a;
	int ref_b = 1 << idx_b;
	complex<double> c = polar(1.0, phi);
	for (int i = 0; i < len; ++i) {
		if ((i & ref_a) && (i & ref_b)) {
			this->q[i] *= c;
		}
	}
}