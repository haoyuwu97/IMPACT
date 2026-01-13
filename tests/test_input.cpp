#include "crystal.h"

#include <fstream>
#include <iostream>
#include <vector>

int main() {
	LOG logdt{};
	std::vector<DATA> vt;
	std::vector<std::vector<DATA2>> bt;
	std::vector<int> typejudge(256, 0);

	{
		std::ifstream in("tests/sample.dump");
		if (!in.is_open()) {
			std::cerr << "failed to open sample.dump\n";
			return 1;
		}
		if (!input(in, vt, bt, logdt, typejudge, true)) {
			std::cerr << "fast input failed on frame 1\n";
			return 1;
		}
		if (logdt.Ncen <= 0 || logdt.Mol < 0) {
			std::cerr << "unexpected parsed counts\n";
			return 1;
		}
		if (!input(in, vt, bt, logdt, typejudge, true)) {
			std::cerr << "fast input failed on frame 2\n";
			return 1;
		}
	}

	{
		std::ifstream in("tests/sample.dump");
		if (!in.is_open()) {
			std::cerr << "failed to open sample.dump for slowio\n";
			return 1;
		}
		if (!input(in, vt, bt, logdt, typejudge, false)) {
			std::cerr << "slow input failed on frame 1\n";
			return 1;
		}
	}

	{
		std::ifstream in("tests/partial.dump");
		if (!in.is_open()) {
			std::cerr << "failed to open partial.dump\n";
			return 1;
		}
		if (input(in, vt, bt, logdt, typejudge, true)) {
			std::cerr << "fast input unexpectedly succeeded on partial dump\n";
			return 1;
		}
	}

	{
		const char* follow_path = "tests/follow_runtime.dump";
		const char* first_frame =
			"ITEM: TIMESTEP\n"
			"0\n"
			"ITEM: NUMBER OF ATOMS\n"
			"4\n"
			"ITEM: BOX BOUNDS xy xz yz pp pp pp\n"
			"0 10 0\n"
			"0 10 0\n"
			"0 10 0\n"
			"ITEM: ATOMS id mol type x y z ix iy iz\n"
			"1 1 1 1 1 1 0 0 0\n"
			"2 1 1 2 1 1 0 0 0\n"
			"3 1 1 3 1 1 0 0 0\n"
			"4 1 1 4 1 1 0 0 0\n"
			"ITEM: TIMESTEP\n"
			"1\n"
			"ITEM: NUMBER OF ATOMS\n"
			"4\n";
		const char* second_tail =
			"ITEM: BOX BOUNDS xy xz yz pp pp pp\n"
			"0 10 0\n"
			"0 10 0\n"
			"0 10 0\n"
			"ITEM: ATOMS id mol type x y z ix iy iz\n"
			"1 1 1 1.1 1 1 0 0 0\n"
			"2 1 1 2.1 1 1 0 0 0\n"
			"3 1 1 3.1 1 1 0 0 0\n"
			"4 1 1 4.1 1 1 0 0 0\n";
		{
			std::ofstream out(follow_path);
			out << first_frame;
		}
		std::ifstream in(follow_path);
		if (!in.is_open()) {
			std::cerr << "failed to open follow_runtime.dump\n";
			return 1;
		}
		if (!input(in, vt, bt, logdt, typejudge, true)) {
			std::cerr << "follow test failed to read first frame\n";
			return 1;
		}
		std::streampos pos = in.tellg();
		if (input(in, vt, bt, logdt, typejudge, true)) {
			std::cerr << "follow test unexpectedly read partial frame\n";
			return 1;
		}
		{
			std::ofstream out(follow_path, std::ios::app);
			out << second_tail;
		}
		in.clear();
		in.seekg(pos);
		if (!input(in, vt, bt, logdt, typejudge, true)) {
			std::cerr << "follow test failed after append\n";
			return 1;
		}
	}

	return 0;
}
