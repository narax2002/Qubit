#include "qubit.h"

void Qubit::QFT()
{
	int nv = this->n;
	for (int i = 0; i < nv / 2; ++i) {
		this->SWAP(i, nv - 1 - i);
	}

	for (int i = 0; i < nv; ++i) {
		int k = 1 << i;
		for (int j = 0; j < i; ++j) {
			this->CR(j, i, pi / k);
			k = k >> 1;
		}
		this->H(i);
	}
}

/**/
void Qubit::FFT()
{
	int n = this->n;
	int d1 = 1 << (n - 1);
	for (int m = 1; m <= n; ++m) {
		int dm = 1 << (n - m);
		int D = 1 << m;
		complex<double> xi = polar(1.0, -2.0 * pi / D);
		Qubit temp(n);

		D = 1 << (m - 1);

		for (int i = 0; i < dm; ++i) {

			complex<double> eta(1.0, 0.0);
				for (int k = 0; k < D; ++k) {
					// temp[k1] = q[k2] + eta * q[k2+d]
					// temp[k1+kd] = q[k2] - eta * q[k2+d]
					int k1 = i + k * dm;
					int k2 = i + 2 * k * dm;
					temp.q[k1] = this->q[k2] + eta * this->q[k2 + dm];
					temp.q[k1 + d1] = this->q[k2] - eta * this->q[k2 + dm];

					eta *= xi;
				}
			
		}
		for (int i = 0; i < this->size; ++i) {
			this->q[i] = temp.q[i];
		}
	}
}
/**/

/**/
void Qubit::IFFT()
{
	int n = this->n;
	int d1 = 1 << (n - 1);
	for (int m = 1; m <= n; ++m) {
		int dm = 1 << (n - m);
		int D = 1 << m;
		complex<double> xi = polar(1.0, 2.0 * pi / D);
		Qubit temp(n);

		D = 1 << (m - 1);
		for (int i = 0; i < dm; ++i) {
			complex<double> eta(1.0, 0.0);
			for (int k = 0; k < D; ++k) {
				// temp[k1] = q[k2] + eta * q[k2+d]
				// temp[k1+kd] = q[k2] - eta * q[k2+d]
				int k1 = i + k * dm;
				int k2 = i + 2 * k * dm;
				temp.q[k1] = this->q[k2] + eta * this->q[k2 + dm];
				temp.q[k1 + d1] = this->q[k2] - eta * this->q[k2 + dm];

				eta *= xi;
			}
		}
		for (int i = 0; i < this->size; ++i) {
			this->q[i] = temp.q[i];
		}
	}
	for (int i = 0; i < this->size; ++i) {
		this->q[i] /= this->size;
	}
}
/**/

/**/
void Qubit::Multi(Qubit& q_a, Qubit& q_b)
{
	for (int i = 0; i < this->size; ++i) {
		this->q[i] = q_a.q[i] * q_b.q[i];
	}
}
/**/

/**/
void Qubit::FFT_iter()
{
	complex<double> w = polar(1.0, 2.0 * pi / this->size);
	if (this->n == 1) {
		complex<double> temp_a = this->q[0];
		complex<double> temp_b = this->q[1];
		this->q[0] = (temp_a + temp_b) / sqrt(2.0);
		this->q[1] = (temp_a + w * temp_b) / sqrt(2.0);
		return;
	}
	
	int nv = this->n - 1;
	Qubit q_even(nv);
	Qubit q_odd(nv);

	int len = q_even.size;
	for (int i = 0; i < len; ++i) {
		q_even.q[i] = this->q[2 * i];
		q_odd.q[i] = this->q[2 * i + 1];
	}

	q_even.FFT_iter();
	q_odd.FFT_iter();

	for (int i = 0; i < len; ++i) {
		this->q[i] = (q_even.q[i] + pow(w, i) * q_odd.q[i]) / sqrt(2.0);
		this->q[len + i] = (q_even.q[i] - pow(w, i) * q_odd.q[i]) / sqrt(2.0);
	}
	return;
}
/**/
