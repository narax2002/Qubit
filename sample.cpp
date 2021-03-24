#include "qubit.h"
#include <omp.h>
#include <chrono>

/// <summary>
/// Personal study  (Qubit algorithm implementation personal Study)
/// Refference
/// "Quantum algorithm implementations for beginners", 10 Apr 2018
/// </summary>

chrono::high_resolution_clock::time_point t1, t2;
#define START() t1 = chrono::high_resolution_clock::now();
#define END() t2 = chrono::high_resolution_clock::now();
#define PRINTTIME(msg, repeat) cout << "* " << msg << " time = " << (double)chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / repeat << " ms" << endl;



int main()
{
	printf("There are %d max threads\n", omp_get_max_threads());
	printf("There are %d threads currently\n", omp_get_num_threads());

	#pragma omp parallel
	{
		if(omp_in_parallel()) {
			printf("In parallel: thread number = %d / %d\n", omp_get_thread_num(), omp_get_num_threads());
		} else {
			printf("Out parallel: thread number = %d / %d\n", omp_get_thread_num(), omp_get_num_threads());
		}
	}
	
	// std::cout << fixed;
	// std::cout.precision(4);
	
	// int N = 18;

	// START();
	// Qubit q1(N), q2(N), q3(N);
	// END();
	// PRINTTIME("Gnerate qubit", 3);

	// START();
	// q1.X(1);
	// q2.X(2);
	// END();
	// PRINTTIME("X gate", 2);
	
	// START();
	// q1.FFT();
	// q2.FFT();
	// END();
	// PRINTTIME("FFT", 2);

	// START();
	// q1.IFFT();
	// q2.IFFT();
	// END();
	// PRINTTIME("IFFT", 2);

	// START();
	// q1.QFT();
	// q2.QFT();
	// END();
	// PRINTTIME("QFT", 2);

	// START();
	// q1.IFFT();
	// q2.IFFT();
	// END();
	// PRINTTIME("IFFT", 2);

	// START();
	// q1.FFT_iter();
	// q2.FFT_iter();
	// END();
	// PRINTTIME("FFT_iter", 2);

	// START();
	// q1.IFFT();
	// q2.IFFT();
	// END();
	// PRINTTIME("IFFT", 2);

	return 0;
}
