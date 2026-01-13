#include <stdio.h>
#include <fstream>
#include <vector>

using namespace std;

struct DATA
{
	double cenx,ceny,cenz,x,y,z,sop;
	int cenid,mol,ty,site,site_size,cenix,ceniy,ceniz,moltype_sop,moltype_dtt;
	int btid_k,btid_s;
};

struct DATA2
{
	double x,y,z,realx,realy,realz,vx,vy,vz;
	int ty,dttid,site,site_size,ix,iy,iz;
};

struct LOG
{
	double box_x,box_y,box_z,sop_avg,dtt_avg,dttc_avg;
	double box_xlo,box_xhi,box_ylo,box_yhi,box_zlo,box_zhi;
	int Ncen;
	
	int Ncount;
	
	int Ntail_s,Ntie_s,Nloop_s;
	double Ltail_s,Ltie_s,Lloop_s;
	int Ntail_d,Ntie_d,Nloop_d;
	double Ltail_d,Ltie_d,Lloop_d;
	
	double sop_crystal_pro,sop_crystal_pro_v;	
	double sop_cluster_gy_max,sop_cluster_gy_avg,sop_cluster_gy_x_max,sop_cluster_gy_x_avg,sop_cluster_gy_y_max,sop_cluster_gy_y_avg,sop_cluster_gy_z_max,sop_cluster_gy_z_avg,sop_cluster_size_avg;	
	double sop_cluster_size_avgW;
	int sop_crystal_num,sop_cluster_num,sop_cluster_size_max;
	int sop_cluster_num2,Ncl_all_sop;
	vector<int> sop_cluster_size;
	vector<double> sop_cluster_gy;
	vector<double> sop_cluster_gyx;
	vector<double> sop_cluster_gyy;
	vector<double> sop_cluster_gyz;
	vector<double> sop_cluster_cx;
	vector<double> sop_cluster_cy;
	vector<double> sop_cluster_cz;
	vector<double> sop_cluster_vx;
	vector<double> sop_cluster_vy;
	vector<double> sop_cluster_vz;
	
	vector<int> ncluster;
	
	double crystal_pro,crystal_pro_v;	
	double cluster_gy_max,cluster_gy_avg,cluster_gy_x_max,cluster_gy_x_avg,cluster_gy_y_max,cluster_gy_y_avg,cluster_gy_z_max,cluster_gy_z_avg,cluster_size_avg;	
	double cluster_size_avgW;
	int crystal_num,cluster_num,cluster_size_max;
	int cluster_num2,Ncl_all;
	vector<int> cluster_size;
	vector<double> cluster_gy;
	vector<double> cluster_gyx;
	vector<double> cluster_gyy;
	vector<double> cluster_gyz;
	vector<double> cluster_cx;
	vector<double> cluster_cy;
	vector<double> cluster_cz;
	vector<double> cluster_vx;
	vector<double> cluster_vy;
	vector<double> cluster_vz;
	vector<int> dtt_num;
	int Mol,Ndtt;
	vector<int> Nmol;
	double probe;

	void reset_sop_clusters(size_t n) {
		sop_cluster_size.assign(n, 0);
		sop_cluster_gy.assign(n, 0.0);
		sop_cluster_gyx.assign(n, 0.0);
		sop_cluster_gyy.assign(n, 0.0);
		sop_cluster_gyz.assign(n, 0.0);
		sop_cluster_cx.assign(n, 0.0);
		sop_cluster_cy.assign(n, 0.0);
		sop_cluster_cz.assign(n, 0.0);
		sop_cluster_vx.assign(n, 0.0);
		sop_cluster_vy.assign(n, 0.0);
		sop_cluster_vz.assign(n, 0.0);
	}

	void reset_dtt_clusters(size_t n) {
		cluster_size.assign(n, 0);
		cluster_gy.assign(n, 0.0);
		cluster_gyx.assign(n, 0.0);
		cluster_gyy.assign(n, 0.0);
		cluster_gyz.assign(n, 0.0);
		cluster_cx.assign(n, 0.0);
		cluster_cy.assign(n, 0.0);
		cluster_cz.assign(n, 0.0);
		cluster_vx.assign(n, 0.0);
		cluster_vy.assign(n, 0.0);
		cluster_vz.assign(n, 0.0);
	}

	void reset_clusters(size_t n) {
		ncluster.assign(n, 0);
	}
};

bool input(istream &in, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, const vector<int>& typejudge, bool fast_io);
void calculate_1(double sop_r, double sop_j, int mode_c, double cluster_p[], ofstream& out4, vector<DATA>& vt, LOG& logdt, int Nt, int mode_v, int v_dump, ofstream& out7, int mode_chain_sop);
void calculate_2(double dtt_j, double dtt_sj, int mode_c, double cluster_p[], ofstream& out6, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, int Nt, int mode_v, int v_dump, ofstream& out8, int mode_chain_dtt);
void volume_cry1(vector<DATA>& vt, LOG& logdt, double sop_j, int v_dump, ofstream& out7, int Nt);
void volume_cry2(vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, double dtt_sj, int v_dump, ofstream& out8, int Nt);
void cluster_1(double sop_j, double cluster_p[], vector<DATA>& vt, LOG& logdt);
void cluster_2(double dtt_j, double dtt_sj, double cluster_p[], vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt);
int findproclust(int sn, LOG& logdt);
void output(int Ntemp, int mode_sop, int mode_dtt, int mode_c_sop, int mode_c_dtt, ofstream& out1, ofstream& out2, ofstream& out3, ofstream& out5, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, int mode_v, int mode_chain_sop, int mode_chain_dtt);
void conformation_sop(vector<DATA>& vt, LOG& logdt, double sop_j);
void conformation_dtt(vector<vector<DATA2>>& bt, vector<DATA>& vt, LOG& logdt, double dtt_sj);
