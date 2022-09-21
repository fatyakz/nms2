#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
#include <mutex>
#include <math.h>
#include <iomanip>

void pSquare(long e, long n, long m) {

	

	std::vector<long> squares(9);
	std::vector<long double> roots(9);

	int longest = 0;

	//build squares
	squares[4] = e * e;

	squares[0] = e + n;
	squares[1] = e - n - m;
	squares[2] = e + n;
	squares[3] = e - n + m;
	
	squares[5] = e + n - m;
	squares[6] = e - n;
	squares[7] = e + n + m;
	squares[8] = e - n;

	//copy sqrt of squares into roots
	for (int i = 0; i < squares.size(); i++) {
		
		roots[i] = sqrt(squares[i]);

		//find longest digits
		if (sizeof(squares[i] > longest)) {
			longest = sizeof(squares[i]);
		}
		if (sizeof(roots[i] > longest)) {
			longest = sizeof(roots[i]);
		}
	}

	std::cout << std::setprecision(8);
	//std::cout << std::setw(longest);
	std::cout << std::setw(12);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);

	std::cout << std::setw(12)  << roots[0] << ":" << roots[1] << ":" << roots[2] << " -> "
		<< squares[0] << ":" << squares[1] << ":" << squares[2] << "\n"
		<< roots[3] << ":" << roots[4] << ":" << roots[5] << " -> "
		<< squares[3] << ":" << squares[4] << ":" << squares[5] << "\n"
		<< roots[6] << ":" << roots[7] << ":" << roots[8] << " -> "
		<< squares[6] << ":" << squares[7] << ":" << squares[8] << "\n";



}

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

bool square(int n) {

	// If ceil and floor are equal
	// the number is a perfect
	// square
	if (ceil((double)sqrt(n)) == floor((double)sqrt(n))) {
		return true;
	}
	else {
		return false;
	}
}

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

uint_fast64_t g_cycles = 0;
uint_fast32_t e = 425;
uint_fast32_t gbest;
uint_fast32_t gbestn;
uint_fast32_t gbestm;

std::vector<uint_fast64_t> jobs;
std::vector<bool> complete;


void thr_jobs(uint_fast32_t start, uint_fast32_t offset, uint_fast32_t threadcount, uint_fast32_t job) {
	uint_fast32_t a, b, c, d, e, f, g, h, i;
	uint_fast64_t cycles = 0; uint_fast32_t matches = 0; uint_fast32_t best = 0;
	e = start * start;
	uint_fast32_t nmlimit = e;
	uint_fast32_t bestn = 0, bestm = 0;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	uint_fast64_t n = 0;

	start:
	for (uint_fast64_t j = 0; j < complete.size(); j++) {
		if (!complete[j]) { uint_fast64_t n = jobs[j]; break; }
	}
	
	for (uint_fast32_t m = 1; m < nmlimit; m ++) {
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
		if (n < nmlimit) { goto start; }
	}
	

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	double cps = (double)cycles / t_time.count();

	std::mutex tex;

	tex.lock();
	complete[job] = true;
	if (best > gbest) { 
		gbest = best; 
		gbestn = bestn;
		gbestm = bestm;
	}
	g_cycles += cycles;
	/*
	std::cout << std::fixed;
	std::cout << "best:" << best;
	std::cout << " n:" << bestn;
	std::cout << " m:" << bestm;
	std::cout << " cycles:" << cycles / 1000000 << "m";
	std::cout << " time:" << t_time.count();
	std::cout << " cps:" << cps / 1000000 << "m\n";
	*/
	//std::cout << "[" << job << ":" << best << "]";

	tex.unlock();
}

void thr_nines(uint_fast32_t start, uint_fast32_t offset, uint_fast32_t threadcount) {
	std::mutex tex;
	uint_fast32_t a, b, c, d, e, f, g, h, i;
	uint_fast64_t cycles = 0; uint_fast32_t matches = 0; uint_fast32_t best = 0;
	e = start * start;
	uint_fast32_t nmlimit = e;
	uint_fast32_t bestn = 0, bestm = 0;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	for (uint_fast32_t n = 1; n < nmlimit; n++) {
		for (uint_fast32_t m = 1 + offset; m < nmlimit; m += threadcount) {
			//if (n == m) { goto end; }
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
	

	tex.lock();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	cycles = (e * e * e * e) / 2;

	double cps = (double)cycles / t_time.count();

	std::cout << std::fixed;
	std::cout << "best:" << best;
	std::cout << " n:" << bestn;
	std::cout << " m:" << bestm;
	std::cout << " cycles:" << cycles / 1000000 << "m";
	std::cout << " time:" << t_time.count();
	std::cout << " cps:" << cps / 1000000 << "m\n";
	tex.unlock();
}

void thr_vector(uint_fast32_t start, uint_fast32_t offset, uint_fast32_t threadcount) {
	std::mutex tex;
	uint_fast32_t a, b, c, d, e, f, g, h, i;
	uint_fast64_t cycles = 0; uint_fast32_t matches = 0; uint_fast32_t best = 0;
	e = start * start;
	uint_fast32_t nmlimit = e;
	uint_fast32_t bestn = 0, bestm = 0;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::vector<uint_fast32_t> squares;
	squares.resize(nmlimit * nmlimit);

	for (uint_fast32_t n = 1; n < nmlimit; n++) {
		for (uint_fast32_t m = 1 + offset; m < nmlimit; m += threadcount) {
			//if (n == m) { goto end; }
			if (n + m >= e) { break; }

			a = e + n;
			b = e - n - m;
			c = e + n;
			d = e - n + m;
			f = e + n - m;
			g = e - n;
			h = e + n + m;
			i = e - n;

			squares[n + m + 0] = a;
			squares[n + m + 1] = b;
			squares[n + m + 2] = c;
			squares[n + m + 3] = d;
			squares[n + m + 4] = e;
			squares[n + m + 5] = f;
			squares[n + m + 6] = g;
			squares[n + m + 7] = h;
			squares[n + m + 8] = i;
			squares[n + m + 9] = n;
			squares[n + m + 10] = m;
		
		}
	}

	for (uint_fast32_t n = 0; n < squares.size(); n += 11) {

		if (isPerfectSquare(squares[n + 0])) { matches++; }
		if (isPerfectSquare(squares[n + 1])) { matches++; }
		if (isPerfectSquare(squares[n + 2])) { matches++; }
		if (isPerfectSquare(squares[n + 3])) { matches++; }
		if (isPerfectSquare(squares[n + 5])) { matches++; }
		if (isPerfectSquare(squares[n + 6])) { matches++; }
		if (isPerfectSquare(squares[n + 7])) { matches++; }
		if (isPerfectSquare(squares[n + 8])) { matches++; }

	end:
		if (matches > best) {
			best = matches;
			bestn = squares[n + 10]; bestm = squares[n + 9];
		}
		matches = 1;
	}

	tex.lock();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	double cps = (double)cycles / t_time.count();

	g_cycles += cycles;

	std::cout << std::fixed;
	std::cout << "best:" << best;
	std::cout << " n:" << bestn;
	std::cout << " m:" << bestm;
	std::cout << " cycles:" << cycles / 1000000 << "m";
	std::cout << " time:" << t_time.count();
	std::cout << " cps:" << cps / 1000000 << "m\n";
	tex.unlock();
}

void thr_ifless(uint_fast32_t start, uint_fast32_t offset, uint_fast32_t threadcount) {
	std::mutex tex;
	uint_fast32_t a, b, c, d, e, f, g, h, i;
	uint_fast64_t cycles = 0; uint_fast32_t matches = 0; uint_fast32_t best = 0;
	e = start * start;
	uint_fast32_t nmlimit = e;
	uint_fast32_t bestn = 0, bestm = 0;
	long long sr = 0;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	for (uint_fast32_t n = 4; n < nmlimit; n++) {
		for (uint_fast32_t m = 1 + offset; m < nmlimit; m += threadcount) {

			if (n + m >= e) { break; }

			a = e + n;
			b = e - n - m;
			c = e + n;
			d = e - n + m;
			f = e + n - m;
			g = e - n;
			h = e + n + m;
			i = e - n;
			
			sr = (long double)sqrt(a);
			matches += (sr * sr == a);

			sr = (long double)sqrt(b);
			matches += (sr * sr == b);

			sr = (long double)sqrt(c);
			matches += (sr * sr == c);

			sr = (long double)sqrt(d);
			matches += (sr * sr == d);

			sr = (long double)sqrt(f);
			matches += (sr * sr == f);

			sr = (long double)sqrt(g);
			matches += (sr * sr == g);

			sr = (long double)sqrt(h);
			matches += (sr * sr == h);

			sr = (long double)sqrt(i);
			matches += (sr * sr == i);

			best += (matches - best) * (matches >= best);
			bestn += (n - bestn) * (matches >= best);
			bestm += (m - bestm) * (matches >= best);

			matches = 0;
			
		}
	}


	
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	double cps = (double)cycles / t_time.count();

	tex.lock();
	g_cycles += cycles;

	std::cout << std::fixed;
	std::cout << "best:" << best;
	std::cout << " n:" << bestn;
	std::cout << " m:" << bestm;
	std::cout << " cycles:" << cycles / 1000000 << "m";
	std::cout << " time:" << t_time.count();
	std::cout << " cps:" << cps / 1000000 << "m\n";

	if (best >= gbest) { gbest = best; gbestn = bestn; gbestm = bestm; }

	tex.unlock();
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
			b = e - n - m;
			c = e + n;
			d = e - n + m;
			f = e + n - m;
			g = e - n;
			h = e + n + m;
			i = e - n;
			
			if (isPerfectSquare(i)) { matches++; }
			if (isPerfectSquare(h)) { matches++; }
			if (isPerfectSquare(g)) { matches++; }
			if (isPerfectSquare(f)) { matches++; }
			if (isPerfectSquare(d)) { matches++; }
			if (isPerfectSquare(c)) { matches++; }
			if (isPerfectSquare(b)) { matches++; }
			if (isPerfectSquare(a)) { matches++; }
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

	if (best >= gbest) { gbest = best; gbestn = bestn; gbestm = bestm; }

	std::cout << std::fixed;
	std::cout << "best:" << best;
	std::cout << " n:" << bestn;
	std::cout << " m:" << bestm;
	std::cout << " cycles:" << cycles / 1000000 << "m";
	std::cout << " time:" << t_time.count();
	std::cout << " cps:" << cps / 1000000 << "m\n";
}

int main()
{

	const uint_fast32_t numthreads = 15;
	

	std::vector<std::thread> thr(numthreads);

	std::cout << "e:";
	std::cin >> e;

	thr.resize(numthreads);

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	
	for (uint_fast32_t t = 0; t < numthreads; t++) {
		thr[t] = std::thread(thr_nines, e, t, numthreads);
	}


	for (uint_fast32_t i = 0; i < numthreads; i++) {
		thr[i].join();
	}
		

	
	g_cycles = (e * e * e * e) / 2;

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	std::cout << "time:" << t_time.count() << "\n";
	std::cout << "cycles:" << g_cycles << "\n";
	std::cout << "cps:" << (g_cycles / t_time.count()) / 1000000 << "m\n";
	std::cout << "e:" << e << "m\n";
	std::cout << "best:" << gbest << " n:" << gbestn << " m:" << gbestm << "\n";

	pSquare(e, gbestn, gbestm);
}