#include "qubit.h"

/// <summary>
/// Bernstein-Vazirani
/// Find Hidden string Boolean f_s(x) = <s,x>
/// </summary>



void Qubit::Bernstein_Vazirani(const int s)
{
	int dim = this->n - 1;
	int len = 1 << dim;

	if (len <= s) this->PrintError(3);
	this->Initial();

	/// <summary>
	/// State Preparation
	this->Initial();
	this->X(dim);
	for (int i = 0; i <= dim; ++i) this->H(i);
	/// </summary>
	
	/// <summary>
	/// Classical Hiddem String
	/// For ancilla qubit (dim) |y> --> |y+f_s(x)>
	for (int i = 0; i < len; ++i) {
		int vec_dot = s & i;
		int f = 0;
		while (vec_dot) {
			f = f ^ (vec_dot & 1);
			vec_dot = vec_dot >> 1;
		}

		if (f) {
			complex<double> temp = this->q[i];
			this->q[i] = this->q[i + len];
			this->q[i + len] = temp;
		}
	}
	/// </summary>

	/// <summary>
	/// Basis Change
	/// </summary>
	for (int i = 0; i < dim; ++i) this->H(i);

	/// <summary>
	/// Measurement
	/// </summary>
	
	//cout << *this << endl;
	/**/
	double* a_temp = this->Qnorm();
	int max_idx = 0;
	for (int i = 1; i < len; ++i) {
		if (a_temp[i] + a_temp[i + len] > a_temp[max_idx] + a_temp[max_idx + len]) {
			max_idx = i;
		}
	}
	cout << "s=" << max_idx << endl;
	/**/
}
