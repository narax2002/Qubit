#include "qubit.h"

/// <summary>
/// Grover algorithm
/// </summary>

void Qubit::Grover(const int k)
{
	int dim = this->n - 1;	// 마지막 큐빗은 ancilla
	int len = 1 << dim;

	if (len <= k) this->PrintError(3);
	this->Initial();

	int MaxIter = (int)ceil(pi * pow(2.0, (double)dim/2) / 4);
	cout << "MaxIter = " << MaxIter << endl;

	/// <summary>
	/// State Preparation
	for (int i = 0; i < dim; ++i) this->H(i);
	this->X(dim);
	this->H(dim);
	/// </summary>


	for (int iter = 0; iter < MaxIter; ++iter) {
		/// <summary>
		/// Quantum Oracle
		int yk = k + (1 << dim);
		complex<double> temp = this->q[k];
		this->q[k] = this->q[yk];
		this->q[yk] = temp;
		//this->Toffoli(dim);
		/// </summary>
		
		/// <summary>
		/// Action identity - 2 * |psi><psi|
		for (int i = 0; i < dim; ++i) this->H(i);

		/// <summary>
		/// Action identity - 2 * |0><0|
		/// |0> --> -|0> and |x> --> |x> otherwise
		for (int i = 0; i < dim; ++i) this->X(i); // |0> --> |2^dim-1>
		this->q[(2 << dim) - 1] *= -1.0;
		this->q[(1 << dim) - 1] *= -1.0;
		/*
		this->H(dim); // |0> --> |+> and |1> --> |->
		this->Toffoli(dim); // if x == 2^dim-1 then |-> --> -|->
		this->Z(dim); // |+> --> |-> and |-> --> |+>
		this->Toffoli(dim); // if x == 2^dim -1 then |-> --> -|->
		this->Z(dim); // |+> --> |-> and |-> --> |+>
		this->H(dim); // |+> --> |0> and -|-> --> -|1>
		/**/

		for (int i = 0; i < dim; ++i) this->X(i); // |2^dim-1> --> |0>
		/// </summary>

		for (int i = 0; i < dim; ++i) this->H(i);
		/// </summary>

		
		
		/*
		double* a_temp = this->Qnorm();
		cout << "(" << a_temp[0] + a_temp[len];
		for (int i = 1; i < (1 << dim); ++i) {
			cout << "," << a_temp[i] + a_temp[i + len];
		}
		cout << ")" << endl;
		/**/
	}

	/// <summary>
	/// Mesuarement
	/// </summary>
	double* a_temp = this->Qnorm();
	int max_idx = 0;
	for (int i = 1; i < len; ++i) {
		if (a_temp[i] + a_temp[i + len] > a_temp[max_idx] + a_temp[max_idx + len]) {
			max_idx = i;
		}
	}
	cout << "k=" << max_idx << endl;
}
