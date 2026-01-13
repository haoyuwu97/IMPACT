#include "crystal.h"

#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdio.h>
#include <vector>

using namespace std;

namespace {
	bool read_non_empty_line(istream& in, string& line) {
		while (std::getline(in, line)) {
			if (!line.empty()) {
				return true;
			}
		}
		return false;
	}

	bool parse_two_doubles(const string& line, double& a, double& b) {
		const char* ptr = line.c_str();
		char* end = nullptr;
		a = std::strtod(ptr, &end);
		if (end == ptr) return false;
		ptr = end;
		b = std::strtod(ptr, &end);
		return end != ptr;
	}

	bool parse_atom_line(const string& line, int& id, int& molid, int& ty,
						 double& x, double& y, double& z,
						 int& ix, int& iy, int& iz) {
		const char* ptr = line.c_str();
		char* end = nullptr;
		id = std::strtol(ptr, &end, 10);
		if (end == ptr) return false;
		ptr = end;
		molid = std::strtol(ptr, &end, 10);
		if (end == ptr) return false;
		ptr = end;
		ty = std::strtol(ptr, &end, 10);
		if (end == ptr) return false;
		ptr = end;
		x = std::strtod(ptr, &end);
		if (end == ptr) return false;
		ptr = end;
		y = std::strtod(ptr, &end);
		if (end == ptr) return false;
		ptr = end;
		z = std::strtod(ptr, &end);
		if (end == ptr) return false;
		ptr = end;
		ix = std::strtol(ptr, &end, 10);
		if (end == ptr) return false;
		ptr = end;
		iy = std::strtol(ptr, &end, 10);
		if (end == ptr) return false;
		ptr = end;
		iz = std::strtol(ptr, &end, 10);
		return end != ptr;
	}

	bool parse_bounds_tokens(istream& in, double& lo, double& hi) {
		if (!(in >> lo >> hi)) {
			return false;
		}
		string rest;
		std::getline(in, rest);
		return true;
	}

	bool build_atoms(int N,
					 vector<DATA>& vt,
					 vector<vector<DATA2>>& bt,
					 LOG& logdt,
					 const vector<int>& typejudge,
					 const vector<int>& id,
					 const vector<int>& molid,
					 const vector<int>& ty,
					 const vector<double>& x,
					 const vector<double>& y,
					 const vector<double>& z,
					 const vector<int>& X,
					 const vector<int>& Y,
					 const vector<int>& Z) {
		const int molid_sentinel = -1;
		const int ty_sentinel = -1;
		logdt.Ncen=0;
		logdt.Mol=-1;
		logdt.Nmol.clear();
		logdt.Nmol.reserve(N);
		bt.clear();
		bt.reserve(N);
		vt.assign(N + 1, DATA());

		for(int i=0;i<N;++i){
			const int prev_molid = (i > 0) ? molid[i - 1] : molid_sentinel;
			const int prev2_molid = (i > 1) ? molid[i - 2] : molid_sentinel;
			const int prev_ty = (i > 0) ? ty[i - 1] : ty_sentinel;
			bool include_curr = (ty[i] < 0 || ty[i] >= static_cast<int>(typejudge.size())) ? true : (typejudge[ty[i]] == 0);
			bool include_prev = (prev_ty < 0 || prev_ty >= static_cast<int>(typejudge.size())) ? true : (typejudge[prev_ty] == 0);
			if(include_curr || (i > 0 && include_prev)){
				if(prev_molid==molid[i]){
					if(i > 1 && prev2_molid==molid[i]){
						vt[logdt.Ncen].cenid=id[i-1];
						vt[logdt.Ncen].ty=ty[i-1];
						vt[logdt.Ncen].cenx=x[i-1];
						vt[logdt.Ncen].ceny=y[i-1];
						vt[logdt.Ncen].cenz=z[i-1];
						vt[logdt.Ncen].cenix=X[i-1];
						vt[logdt.Ncen].ceniy=Y[i-1];
						vt[logdt.Ncen].ceniz=Z[i-1];
						vt[logdt.Ncen].x=x[i]+X[i]*logdt.box_x-x[i-2]-X[i-2]*logdt.box_x;
						vt[logdt.Ncen].y=y[i]+Y[i]*logdt.box_y-y[i-2]-Y[i-2]*logdt.box_y;
						vt[logdt.Ncen].z=z[i]+Z[i]*logdt.box_z-z[i-2]-Z[i-2]*logdt.box_z;
						vt[logdt.Ncen].mol=molid[i-1];
						
						logdt.Ncen++;
					}if(i > 1 && prev2_molid!=molid[i]){
						vt[logdt.Ncen].cenid=id[i-1];
						vt[logdt.Ncen].ty=ty[i-1];
						vt[logdt.Ncen].cenx=x[i-1];
						vt[logdt.Ncen].ceny=y[i-1];
						vt[logdt.Ncen].cenz=z[i-1];
						vt[logdt.Ncen].cenix=X[i-1];
						vt[logdt.Ncen].ceniy=Y[i-1];
						vt[logdt.Ncen].ceniz=Z[i-1];
						vt[logdt.Ncen].x=x[i-1]+X[i-1]*logdt.box_x-x[i-2]-X[i-2]*logdt.box_x;
						vt[logdt.Ncen].y=y[i-1]+Y[i-1]*logdt.box_y-y[i-2]-Y[i-2]*logdt.box_y;
						vt[logdt.Ncen].z=z[i-1]+Z[i-1]*logdt.box_z-z[i-2]-Z[i-2]*logdt.box_z;
						vt[logdt.Ncen].mol=molid[i-1];
						
						logdt.Ncen++;
					}
					
					if(i==N-1){
						vt[logdt.Ncen].cenid=id[i];
						vt[logdt.Ncen].ty=ty[i];
						vt[logdt.Ncen].cenx=x[i];
						vt[logdt.Ncen].ceny=y[i];
						vt[logdt.Ncen].cenz=z[i];
						vt[logdt.Ncen].cenix=X[i];
						vt[logdt.Ncen].ceniy=Y[i];
						vt[logdt.Ncen].ceniz=Z[i];
						vt[logdt.Ncen].x=x[i]+X[i]*logdt.box_x-x[i-1]-X[i-1]*logdt.box_x;
						vt[logdt.Ncen].y=y[i]+Y[i]*logdt.box_y-y[i-1]-Y[i-1]*logdt.box_y;
						vt[logdt.Ncen].z=z[i]+Z[i]*logdt.box_z-z[i-1]-Z[i-1]*logdt.box_z;
						vt[logdt.Ncen].mol=molid[i];
						
						logdt.Ncen++;
					}
				}
				if(prev_molid!=molid[i]){
					logdt.Mol++;
					logdt.Nmol.push_back(-1);
					bt.push_back(vector<DATA2>());
					
					if(i>1){
						vt[logdt.Ncen].cenid=id[i-1];
						vt[logdt.Ncen].ty=ty[i-1];
						vt[logdt.Ncen].cenx=x[i-1];
						vt[logdt.Ncen].ceny=y[i-1];
						vt[logdt.Ncen].cenz=z[i-1];
						vt[logdt.Ncen].x=x[i-1]+X[i-1]*logdt.box_x-x[i-2]-X[i-2]*logdt.box_x;
						vt[logdt.Ncen].y=y[i-1]+Y[i-1]*logdt.box_y-y[i-2]-Y[i-2]*logdt.box_y;
						vt[logdt.Ncen].z=z[i-1]+Z[i-1]*logdt.box_z-z[i-2]-Z[i-2]*logdt.box_z;
						vt[logdt.Ncen].mol=molid[i-1];
						
						logdt.Ncen++;
					}
					
				}
				logdt.Nmol[logdt.Mol]++;
				DATA2 atom;
				atom.x=x[i];
				atom.y=y[i];
				atom.z=z[i];
				atom.ix=X[i];
				atom.iy=Y[i];
				atom.iz=Z[i];
				atom.realx=x[i]+X[i]*logdt.box_x;
				atom.realy=y[i]+Y[i]*logdt.box_y;
				atom.realz=z[i]+Z[i]*logdt.box_z;
				atom.ty=ty[i];
				atom.dttid=0;
				bt[logdt.Mol].push_back(atom);
			}
		}
		vt.resize(logdt.Ncen + 1);
		vt[logdt.Ncen].mol=0;

		return true;
	}
}

static bool input_fast(istream &in, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, const vector<int>& typejudge){
	string item;
	string label;
	if (!(in >> item >> label)) {
		return false;
	}
	if (item != "ITEM:" || label != "TIMESTEP") {
		return false;
	}
	string rest;
	std::getline(in, rest);
	int timestep = 0;
	if (!(in >> timestep)) {
		return false;
	}
	std::getline(in, rest);
	string label2;
	string label3;
	if (!(in >> item >> label2 >> label3)) {
		return false;
	}
	if (item != "ITEM:" || label2 != "NUMBER" || label3 != "OF") {
		return false;
	}
	string atoms_word;
	if (!(in >> atoms_word)) {
		return false;
	}
	if (atoms_word != "ATOMS") {
		return false;
	}
	int N = 0;
	if (!(in >> N)) {
		return false;
	}
	if (N <= 0) {
		return false;
	}
	std::getline(in, rest);
	if (!(in >> item >> label2 >> label3)) {
		return false;
	}
	if (item != "ITEM:" || label2 != "BOX" || label3 != "BOUNDS") {
		return false;
	}
	std::getline(in, rest);
	if (!parse_bounds_tokens(in, logdt.box_xlo, logdt.box_xhi)) return false;
	if (!parse_bounds_tokens(in, logdt.box_ylo, logdt.box_yhi)) return false;
	if (!parse_bounds_tokens(in, logdt.box_zlo, logdt.box_zhi)) return false;
	logdt.box_x=logdt.box_xhi-logdt.box_xlo;
	logdt.box_y=logdt.box_yhi-logdt.box_ylo;
	logdt.box_z=logdt.box_zhi-logdt.box_zlo;
	if (!(in >> item >> label2)) return false;
	if (item != "ITEM:" || label2 != "ATOMS") {
		return false;
	}
	std::getline(in, rest);

	vector<int> ty(N + 2);
	vector<int> id(N + 2);
	vector<int> molid(N + 5);
	vector<double> x(N + 2);
	vector<double> y(N + 2);
	vector<double> z(N + 2);
	vector<int> X(N + 2);
	vector<int> Y(N + 2);
	vector<int> Z(N + 2);
	for (int i = 0; i < N; ++i) {
		if (!(in >> id[i] >> molid[i] >> ty[i] >> x[i] >> y[i] >> z[i] >> X[i] >> Y[i] >> Z[i])) {
			return false;
		}
		std::getline(in, rest);
	}
	return build_atoms(N, vt, bt, logdt, typejudge, id, molid, ty, x, y, z, X, Y, Z);
}

static bool input_line(istream &in, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, const vector<int>& typejudge){
	string line;
	if (!read_non_empty_line(in, line)) {
		return false;
	}
	if (line.find("ITEM: TIMESTEP") != 0) {
		return false;
	}
	if (!read_non_empty_line(in, line)) return false;
	if (!read_non_empty_line(in, line)) return false;
	if (line.find("ITEM: NUMBER OF ATOMS") != 0) return false;
	if (!read_non_empty_line(in, line)) return false;
	int N = std::atoi(line.c_str());
	if (N <= 0) return false;

	if (!read_non_empty_line(in, line)) return false;
	if (line.find("ITEM: BOX BOUNDS") != 0) return false;
	if (!read_non_empty_line(in, line)) return false;
	double w = 0.0, w1 = 0.0, m1 = 0.0, n1 = 0.0, m2 = 0.0, n2 = 0.0;
	if (!parse_two_doubles(line, w, w1)) return false;
	if (!read_non_empty_line(in, line)) return false;
	if (!parse_two_doubles(line, m1, n1)) return false;
	if (!read_non_empty_line(in, line)) return false;
	if (!parse_two_doubles(line, m2, n2)) return false;
	if (!read_non_empty_line(in, line)) return false;
	if (line.find("ITEM: ATOMS") != 0) return false;

	vector<int> ty(N + 2);
	vector<int> id(N + 2);
	vector<int> molid(N + 5);
	vector<double> x(N + 2);
	vector<double> y(N + 2);
	vector<double> z(N + 2);
	vector<int> X(N + 2);
	vector<int> Y(N + 2);
	vector<int> Z(N + 2);

	logdt.box_xlo=w;
	logdt.box_xhi=w1;
	logdt.box_ylo=m1;
	logdt.box_yhi=n1;
	logdt.box_zlo=m2;
	logdt.box_zhi=n2;
	logdt.box_x=w1-w;
	logdt.box_y=n1-m1;
	logdt.box_z=n2-m2;

	for(int i=0;i<N;++i){
		if (!read_non_empty_line(in, line)) return false;
		if (!parse_atom_line(line, id[i], molid[i], ty[i], x[i], y[i], z[i], X[i], Y[i], Z[i])) {
			return false;
		}
	}
	return build_atoms(N, vt, bt, logdt, typejudge, id, molid, ty, x, y, z, X, Y, Z);
}

bool input(istream &in, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, const vector<int>& typejudge, bool fast_io){
	if (fast_io) {
		return input_fast(in, vt, bt, logdt, typejudge);
	}
	return input_line(in, vt, bt, logdt, typejudge);
}
