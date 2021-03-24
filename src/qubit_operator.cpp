#include "qubit.h"

/// <summary>
/// qubit operator
/// </summary>

Qubit::Qubit()
{
	n = 1;
	size = 2;
	q = new complex<double>[size];

	q[0] = complex<double>(1.0, 0.0);
	q[1] = complex<double>(0.0, 0.0);
}

Qubit::Qubit(const int nv) {
	n = nv;
	size = 1 << nv;
	q = new complex<double>[size];

	q[0] = complex<double>(1.0, 0.0);
	for (int i = 1; i < size; ++i) {
		q[i] = complex<double>(0.0, 0.0);
	}
}

void Qubit::Initial() {
	this->q[0] = complex<double>(1.0, 0.0);
	for (int i = 1; i < this->size; ++i) {
		q[i] = complex<double>(0.0, 0.0);
	}
}

ostream& operator<<(ostream& os, const Qubit& Q)
{
	//os << "(";
	os << Q.q[0];
	for (int i = 1; i < Q.size; ++i) {
		os << "," << Q.q[i];
	}
	//os << ")";

	return os;
}

double* Qubit::Qnorm()
{
	int len = this->size;
	double* result = new double[len];
	for (int i = 0; i < len; ++i) {
		result[i] = norm(this->q[i]);
	}

	return result;
}

void Qubit::PrintQnorm()
{
	int len = this->size;
	double* temp = this->Qnorm();
	cout << "(";
	for (int i = 0; i < len - 1; ++i) {
		cout << temp[i] << ", ";
	}
	cout << temp[len - 1] << ")" << endl;
}

void Qubit::PrintError(const int n)
{
	switch(n) {
	case 1:
		cout << "Array index is out of bound" << endl;
		exit(1);
	case 2:
		cout << "The number of qubits is insufficient" << endl;
		exit(1);
	case 3:
		cout << "Can not implement this algorithm" << endl;
		exit(1);
	}
}
