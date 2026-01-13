#include "crystal.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
	const char* path = (argc > 1) ? argv[1] : "tests/sample.dump";
	int iterations = (argc > 2) ? std::atoi(argv[2]) : 1000;

	LOG logdt{};
	std::vector<DATA> vt;
	std::vector<std::vector<DATA2>> bt;
	std::vector<int> typejudge(256, 0);

	auto start = std::chrono::steady_clock::now();
	long long total_atoms = 0;
	int frames = 0;
	for (int i = 0; i < iterations; ++i) {
		std::ifstream in(path);
		if (!in.is_open()) {
			std::cerr << "failed to open benchmark input\n";
			return 1;
		}
		while (input(in, vt, bt, logdt, typejudge, true)) {
			total_atoms += logdt.Ncen;
			frames++;
		}
	}
	auto end = std::chrono::steady_clock::now();
	double seconds = std::chrono::duration<double>(end - start).count();
	double time_per_atom = (total_atoms > 0) ? (seconds / total_atoms) : 0.0;
	std::cout << "frames=" << frames << " atoms=" << total_atoms
	          << " seconds=" << seconds
	          << " time_per_atom=" << time_per_atom << "\n";
	return 0;
}
