#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
#include <mutex>
#include <cmath>

class Timer {
	public:
		Timer() {
			auto m_StartTimepoint = std::chrono::high_resolution_clock::now();
		}
		~Timer() {
			Stop();
		}
		void Stop() {
			auto endTimepoint = std::chrono::high_resolution_clock::now();

			auto start = std::chrono::time_point_cast<std::chrono::seconds>(m_StartTimepoint).time_since_epoch().count();
			auto end = std::chrono::time_point_cast<std::chrono::seconds>(endTimepoint).time_since_epoch().count();

			auto ms = (end - start) / 1000;
			double s = ms / 1000;

			std::cout << ms << "ms (" << s << "s)\n";
			
		}
private:
	std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
};

bool isPerfectSquare(long double x) {
	// Find floating point value of
	// square root of x.
	if (x >= 0) {

		long long sr = sqrt(x);

		// if product of square root
		//is equal, then
		// return T/F
		return (sr * sr == x);
	}
	// else return false if n<0
	return false;
}

void thr_bench_nines(long start, long offset, long threadcount) {
	long a, b, c, d, e, f, g, h, i;
	long cycles = 0; long matches = 0; long best = 0;
	e = start * start;
	long nmlimit = e;
	long bestn = 0, bestm = 0;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	for (long n = 1; n < nmlimit; n++) {

		for (long m = 1 + offset; m < nmlimit; m += threadcount) {
			if (n == m) { goto end; }
			if (n + m >= e) { break; }
			
			a = e + n;
			if (!isPerfectSquare(a)) { goto end; }
			matches++;
			b = e - n - m;
			if (!isPerfectSquare(b)) { goto end; }
			matches++;
			c = e + n;
			if (!isPerfectSquare(c)) { goto end; }
			matches++;
			d = e - n + m;
			if (!isPerfectSquare(d)) { goto end; }
			matches++;
			f = e + n - m;
			if (!isPerfectSquare(f)) { goto end; }
			matches++;
			g = e - n;
			if (!isPerfectSquare(g)) { goto end; }
			matches++;
			h = e + n + m;
			if (!isPerfectSquare(h)) { goto end; }
			matches++;
			i = e - n;
			if (!isPerfectSquare(i)) { goto end; }
			matches++;
		end:
			cycles++;
			

			if (matches > best) {
				best = matches;
				bestn = n; bestm = m;
			}
			matches = 1;
		}
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	double cps = (double)cycles / t_time.count();

	std::cout << std::fixed;
	std::cout << "best:" << best;
	std::cout << " n:" << bestn;
	std::cout << " m:" << bestm;
	std::cout << " cycles:" << cycles / 1000000 << "m";
	std::cout << " time:" << t_time.count();
	std::cout << " cps:" << cps / 1000000 << "m\n";
}
void thr_all(long start, long offset, long threadcount) {
	long a, b, c, d, e, f, g, h, i;
	long cycles = 0; long matches = 0; long best = 0;
	e = start * start;
	long nmlimit = e;
	long bestn = 0, bestm = 0;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	for (long n = 1; n < nmlimit; n++) {

		for (long m = 1 + offset; m < nmlimit; m += threadcount) {
			if (n == m) { goto end; }
			if (n + m >= e) { break; }

			a = e + n;
			if (!isPerfectSquare(a)) { goto end; }
			matches++;
			b = e - n - m;
			if (!isPerfectSquare(b)) { goto end; }
			matches++;
			c = e + n;
			if (!isPerfectSquare(c)) { goto end; }
			matches++;
			d = e - n + m;
			if (!isPerfectSquare(d)) { goto end; }
			matches++;
			f = e + n - m;
			if (!isPerfectSquare(f)) { goto end; }
			matches++;
			g = e - n;
			if (!isPerfectSquare(g)) { goto end; }
			matches++;
			h = e + n + m;
			if (!isPerfectSquare(h)) { goto end; }
			matches++;
			i = e - n;
			if (!isPerfectSquare(i)) { goto end; }
			matches++;
		end:
			cycles++;


			if (matches > best) {
				best = matches;
				bestn = n; bestm = m;
			}
			matches = 1;
		}
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	double cps = (double)cycles / t_time.count();

	std::cout << std::fixed;
	std::cout << "best:" << best;
	std::cout << " n:" << bestn;
	std::cout << " m:" << bestm;
	std::cout << " cycles:" << cycles / 1000000 << "m";
	std::cout << " time:" << t_time.count();
	std::cout << " cps:" << cps / 1000000 << "m\n";
}
//long double sqrt(long double n) {
//	long double x = n / 2;
//	while (1) {
//		x = 0.5 * (x + (n / x));
//		if (x * x >= n - 0.0000001 && x * x <= n + 0.0000001) {
//			break;
//		}
//	}
//	return x;
//}

int main()
{
	
	long numthreads = 15;
	long e = 425;

	std::vector<std::thread> thr(numthreads);

	thr.resize(numthreads);

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	for (long i = 0; i < numthreads; i++) {
		thr[i] = std::thread(thr_all, e, i, numthreads);
		thr[i].join();
	}
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	std::cout << "time:" << t_time.count();

	std::cout << " done" << "\n";
}