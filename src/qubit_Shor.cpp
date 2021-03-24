#include "qubit.h"

/// <summary>
/// Shor algorithm
/// Reference
/// Experimental realization of Shor's quantum factoring algorithm using nuclear magnetic resonance, Nature 2001.
/// The case N=15, a=7
/// </summary>


void Qubit::Shor()
{
	this->H(0);
	this->H(1);
	this->H(2);

	this->CX(4, 2);
	this->CX(5, 2);
	this->CX(5, 3);
	this->CCX(3, 1, 5);
	this->CX(5, 3);
	this->CX(4, 6);
	this->CCX(6, 1, 4);
	this->CX(4, 6);

	this->H(0);
	this->CR(0, 1, - pi / 2);
	this->H(1);
	this->CR(0, 2, - pi / 4);
	this->CR(1, 2, - pi / 2);
	this->H(2);
	
	double* a_temp = this->Qnorm();
	double temp = 0;
	for (int j = 0; j < 16; ++j) temp += a_temp[8 * j];
	cout << "(" << temp;
	for (int i = 1;i < 8;++i) {
		temp = 0;
		for (int j = 0; j < 16; ++j) temp += a_temp[8 * j + i];
		cout << ", " << temp;
	}
	cout << ")" << endl;

	cout << *this << endl;
}