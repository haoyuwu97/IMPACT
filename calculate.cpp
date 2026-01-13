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
void calculate_1(double sop_r, double sop_j, int mode_c, double cluster_p[], ofstream& out4, vector<DATA>& vt, LOG& logdt, int Nt, int mode_v, int v_dump, ofstream& out7, int mode_chain_sop){
	
	if(mode_c!=1&&mode_c!=0){
		cout << "ERROR: Wrong cluster mode set" << endl;
		exit(0);
	}
	
	double rx=0.0,ry=0.0,rz=0.0,vec_AB=0.0,vec_A=0.0,vec_B=0.0,CosD=0.0;
	logdt.sop_crystal_num=0;
	logdt.sop_avg=0.0;
	int Nsop;

	int ncen = logdt.Ncen;
	static std::vector<int> cell_index;
	static std::vector<int> head;
	static std::vector<int> next;
	static std::vector<double> prev_cenx;
	static std::vector<double> prev_ceny;
	static std::vector<double> prev_cenz;
	static double prev_box_x = 0.0;
	static double prev_box_y = 0.0;
	static double prev_box_z = 0.0;
	static double prev_sop_r = 0.0;
	bool reuse_cells = false;
	int nx = (sop_r > 0.0) ? int(logdt.box_x / sop_r) : 1;
	int ny = (sop_r > 0.0) ? int(logdt.box_y / sop_r) : 1;
	int nz = (sop_r > 0.0) ? int(logdt.box_z / sop_r) : 1;
	if (nx < 1) nx = 1;
	if (ny < 1) ny = 1;
	if (nz < 1) nz = 1;
	int ncell = nx * ny * nz;
	if (!cell_index.empty() && prev_cenx.size() == static_cast<size_t>(ncen) && prev_box_x == logdt.box_x && prev_box_y == logdt.box_y && prev_box_z == logdt.box_z && prev_sop_r == sop_r) {
		double threshold = (sop_r > 0.0) ? 0.25 * sop_r : 0.0;
		if (threshold > 0.0) {
			double max_disp2 = 0.0;
			for (int k = 0; k < ncen; ++k) {
				double dx = vt[k].cenx - prev_cenx[k];
				double dy = vt[k].ceny - prev_ceny[k];
				double dz = vt[k].cenz - prev_cenz[k];
				double d2 = dx * dx + dy * dy + dz * dz;
				if (d2 > max_disp2) max_disp2 = d2;
				prev_cenx[k] = vt[k].cenx;
				prev_ceny[k] = vt[k].ceny;
				prev_cenz[k] = vt[k].cenz;
			}
			reuse_cells = max_disp2 <= (threshold * threshold);
		}
	}
	if (!reuse_cells) {
		head.assign(ncell, -1);
		next.assign(ncen, -1);
		cell_index.assign(ncen, 0);
		prev_cenx.assign(ncen, 0.0);
		prev_ceny.assign(ncen, 0.0);
		prev_cenz.assign(ncen, 0.0);
		prev_box_x = logdt.box_x;
		prev_box_y = logdt.box_y;
		prev_box_z = logdt.box_z;
		prev_sop_r = sop_r;
		for (int k = 0; k < ncen; ++k) {
			prev_cenx[k] = vt[k].cenx;
			prev_ceny[k] = vt[k].ceny;
			prev_cenz[k] = vt[k].cenz;
			double fx = (vt[k].cenx - logdt.box_xlo) / logdt.box_x;
			double fy = (vt[k].ceny - logdt.box_ylo) / logdt.box_y;
			double fz = (vt[k].cenz - logdt.box_zlo) / logdt.box_z;
			fx -= floor(fx);
			fy -= floor(fy);
			fz -= floor(fz);
			int ix = int(fx * nx);
			int iy = int(fy * ny);
			int iz = int(fz * nz);
			if (ix == nx) ix = nx - 1;
			if (iy == ny) iy = ny - 1;
			if (iz == nz) iz = nz - 1;
			int cell = (ix * ny + iy) * nz + iz;
			cell_index[k] = cell;
			next[k] = head[cell];
			head[cell] = k;
		}
	}

	for(int k=0;k<logdt.Ncen;++k){
		Nsop=0;
		vt[k].sop=0.0;
		vec_A=vt[k].x*vt[k].x + vt[k].y*vt[k].y + vt[k].z*vt[k].z;
		vec_A=sqrt(vec_A);
		int cell = cell_index[k];
		int cx = cell / (ny * nz);
		int cy = (cell / nz) % ny;
		int cz = cell % nz;
		for (int dx = -1; dx <= 1; ++dx) {
			int nx_i = (cx + dx + nx) % nx;
			for (int dy = -1; dy <= 1; ++dy) {
				int ny_i = (cy + dy + ny) % ny;
				for (int dz = -1; dz <= 1; ++dz) {
					int nz_i = (cz + dz + nz) % nz;
					int neighbor = (nx_i * ny + ny_i) * nz + nz_i;
					for (int s = head[neighbor]; s != -1; s = next[s]) {
						if (k == s) continue;
						
						rx=vt[k].cenx-vt[s].cenx;
						ry=vt[k].ceny-vt[s].ceny;
						rz=vt[k].cenz-vt[s].cenz;
						if(rx>logdt.box_x/2){rx=rx-logdt.box_x;}
						if(rx<-logdt.box_x/2){rx=rx+logdt.box_x;}
						if(ry>logdt.box_y/2){ry=ry-logdt.box_y;}
						if(ry<-logdt.box_y/2){ry=ry+logdt.box_y;}
						if(rz>logdt.box_z/2){rz=rz-logdt.box_z;}
						if(rz<-logdt.box_z/2){rz=rz+logdt.box_z;}
						rx=rx*rx+ry*ry+rz*rz;
						rx=sqrt(rx);
						if(rx<=sop_r){
							Nsop++;
							vec_AB=vt[k].x*vt[s].x + vt[k].y*vt[s].y + vt[k].z*vt[s].z;
							vec_B=vt[s].x*vt[s].x + vt[s].y*vt[s].y + vt[s].z*vt[s].z;
							vec_B=sqrt(vec_B);
							CosD=vec_AB/(vec_A*vec_B);
							vt[k].sop+=(3*CosD*CosD-1)/2;
						}
					}
				}
			}
		}
		if(Nsop!=0){
			vt[k].sop=vt[k].sop/Nsop;
		}else{
			vt[k].sop=vt[k].sop;
		}
		logdt.sop_avg+=vt[k].sop;
		if(vt[k].sop>=sop_j){
			logdt.sop_crystal_num++;
		}
	}
	logdt.sop_avg=logdt.sop_avg/logdt.Ncen;
	logdt.sop_crystal_pro=double(logdt.sop_crystal_num)/logdt.Ncen;

	if(mode_v==1){
		volume_cry1(vt, logdt, sop_j, v_dump, out7, Nt);
	}
	
	if(mode_chain_sop==1){
		conformation_sop(vt, logdt, sop_j);
	}
	
	if(mode_c==1){
		cluster_1(sop_j, cluster_p, vt, logdt);	

		out4<<"ITEM: TIMESTEP"<<endl;
		out4<<Nt<<endl;
		out4<<"ITEM: NUMBER OF ATOMS"<<endl;
		out4<<logdt.sop_cluster_num<<endl;
		out4<<"ITEM: BOX BOUNDS xy xz yz pp pp pp"<<endl;
		out4<<logdt.box_xlo<<" "<<logdt.box_xhi<<" "<<0.0<<endl;
		out4<<logdt.box_ylo<<" "<<logdt.box_yhi<<" "<<0.0<<endl;
		out4<<logdt.box_zlo<<" "<<logdt.box_zhi<<" "<<0.0<<endl;
		out4<<"ITEM: ATOMS id x y z vx vy vz rg rgx rgy rgz size"<<endl;
		logdt.sop_cluster_size_max=0;
		logdt.sop_cluster_gy_max=0;
		logdt.sop_cluster_gy_x_max=0;
		logdt.sop_cluster_gy_y_max=0;
		logdt.sop_cluster_gy_z_max=0;
		logdt.sop_cluster_gy_avg=0;
		logdt.sop_cluster_gy_x_avg=0;
		logdt.sop_cluster_gy_y_avg=0;
		logdt.sop_cluster_gy_z_avg=0;
		logdt.sop_cluster_size_avg=0;
		logdt.sop_cluster_size_avgW=0;
		for(int q=0;q<logdt.sop_cluster_num;++q){
			out4 << q+1 << " " << logdt.sop_cluster_cx[q] << " " << logdt.sop_cluster_cy[q] << " " << logdt.sop_cluster_cz[q] << " " << logdt.sop_cluster_vx[q] << " " << logdt.sop_cluster_vy[q] << " " << logdt.sop_cluster_vz[q] << " " << logdt.sop_cluster_gy[q] << " " << logdt.sop_cluster_gyx[q] << " " << logdt.sop_cluster_gyy[q] << " " << logdt.sop_cluster_gyz[q] << " " << logdt.sop_cluster_size[q]<<endl;
			/*logdt.sop_cluster_gy_avg+=logdt.sop_cluster_gy[q]/double(logdt.sop_cluster_num);
			logdt.sop_cluster_gy_x_avg+=logdt.sop_cluster_gyx[q]/double(logdt.sop_cluster_num);
			logdt.sop_cluster_gy_y_avg+=logdt.sop_cluster_gyy[q]/double(logdt.sop_cluster_num);
			logdt.sop_cluster_gy_z_avg+=logdt.sop_cluster_gyz[q]/double(logdt.sop_cluster_num);*/
			
			logdt.sop_cluster_gy_avg+=logdt.sop_cluster_gy[q]*logdt.sop_cluster_size[q]/double(logdt.Ncl_all_sop);
			logdt.sop_cluster_gy_x_avg+=logdt.sop_cluster_gyx[q]*logdt.sop_cluster_size[q]/double(logdt.Ncl_all_sop);
			logdt.sop_cluster_gy_y_avg+=logdt.sop_cluster_gyy[q]*logdt.sop_cluster_size[q]/double(logdt.Ncl_all_sop);
			logdt.sop_cluster_gy_z_avg+=logdt.sop_cluster_gyz[q]*logdt.sop_cluster_size[q]/double(logdt.Ncl_all_sop);
			logdt.sop_cluster_size_avg+=double(logdt.sop_cluster_size[q])/double(logdt.sop_cluster_num);
			logdt.sop_cluster_size_avgW+=double(logdt.sop_cluster_size[q]*logdt.sop_cluster_size[q])/double(logdt.Ncl_all_sop);
			
			if(logdt.sop_cluster_size[q]>logdt.sop_cluster_size_max) logdt.sop_cluster_size_max=logdt.sop_cluster_size[q];
			if(logdt.sop_cluster_gy[q]>logdt.sop_cluster_gy_max) logdt.sop_cluster_gy_max=logdt.sop_cluster_gy[q];
			if(logdt.sop_cluster_gyx[q]>logdt.sop_cluster_gy_x_max) logdt.sop_cluster_gy_x_max=logdt.sop_cluster_gyx[q];
			if(logdt.sop_cluster_gyy[q]>logdt.sop_cluster_gy_y_max) logdt.sop_cluster_gy_y_max=logdt.sop_cluster_gyy[q];
			if(logdt.sop_cluster_gyz[q]>logdt.sop_cluster_gy_z_max) logdt.sop_cluster_gy_z_max=logdt.sop_cluster_gyz[q];
		}
	}

}

void calculate_2(double dtt_j, double dtt_sj, int mode_c, double cluster_p[], ofstream& out6, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, int Nt, int mode_v, int v_dump, ofstream& out8, int mode_chain_dtt){
	
	if(mode_c!=1&&mode_c!=0){
		cout << "ERROR: Wrong cluster mode set" << endl;
		exit(0);
	}
	
	double bond1x=0.0,bond1y=0.0,bond1z=0.0,bond2x=0.0,bond2y=0.0,bond2z=0.0;
	double vec_AB=0.0,vec_A=0.0,vec_B=0.0,CosD=0.0,p2=0.0;
	logdt.crystal_num=0;
	logdt.crystal_pro=0.0;
	logdt.dtt_avg=0;
	logdt.dttc_avg=0;
	logdt.Ndtt=0;
	logdt.dtt_num.assign(logdt.Ncen + logdt.Mol + 2, 0);
	auto ensure_dtt = [&](int id) {
		if (id >= static_cast<int>(logdt.dtt_num.size())) {
			logdt.dtt_num.resize(id + 1, 0);
		}
	};
	int Ntrans=0;

	for(int k=0;k<logdt.Mol+1;++k){
		logdt.Ndtt++;
		bt[k][0].dttid=logdt.Ndtt;
		ensure_dtt(logdt.Ndtt);
			logdt.dtt_num[logdt.Ndtt]++;
			
		vt[Ntrans].btid_k=k;
		vt[Ntrans].btid_s=0;
		Ntrans++;
			
		for(int s=1;s<logdt.Nmol[k];++s){
			vt[Ntrans].btid_k=k;
			vt[Ntrans].btid_s=s;
			Ntrans++;
			if(bt[k][s].dttid!=0){
				continue;
			}
			
			bond1x=bt[k][s].realx-bt[k][s-1].realx;
			bond1y=bt[k][s].realy-bt[k][s-1].realy;
			bond1z=bt[k][s].realz-bt[k][s-1].realz;
			vec_A=bond1x*bond1x+bond1y*bond1y+bond1z*bond1z;
			vec_A=sqrt(vec_A);

			bond2x=bt[k][s+1].realx-bt[k][s].realx;
			bond2y=bt[k][s+1].realy-bt[k][s].realy;
			bond2z=bt[k][s+1].realz-bt[k][s].realz;
			vec_B=bond2x*bond2x+bond2y*bond2y+bond2z*bond2z;
			vec_B=sqrt(vec_B);
				
			vec_AB=bond1x*bond2x+bond1y*bond2y+bond1z*bond2z;
			CosD=vec_AB/(vec_A*vec_B);
			p2=(3*CosD*CosD-1)/2;
			if(p2>=dtt_j){
				bt[k][s].dttid=bt[k][s-1].dttid;
				ensure_dtt(logdt.Ndtt);
				logdt.dtt_num[logdt.Ndtt]++;
				if(s==logdt.Nmol[k]-1){
					bt[k][s+1].dttid=bt[k][s-1].dttid;
					ensure_dtt(logdt.Ndtt);
					logdt.dtt_num[logdt.Ndtt]++;
				}
			}else{
				logdt.Ndtt++;
				bt[k][s].dttid=logdt.Ndtt;
				ensure_dtt(logdt.Ndtt);
				logdt.dtt_num[logdt.Ndtt]++;
			}
		}
		if(bt[k][logdt.Nmol[k]].dttid==0){
			logdt.Ndtt++;
			bt[k][logdt.Nmol[k]].dttid=logdt.Ndtt;
			ensure_dtt(logdt.Ndtt);
			logdt.dtt_num[logdt.Ndtt]++;
		}
		vt[Ntrans].btid_k=k;
		vt[Ntrans].btid_s=logdt.Nmol[k];
		Ntrans++;
	}

	int Ncal1=0;
	int Ncal2=0;
	for(int ss=1;ss<=logdt.Ndtt;ss++){
		if(double(logdt.dtt_num[ss])>1){
			Ncal1++;
			logdt.dtt_avg+=logdt.dtt_num[ss];
			if(double(logdt.dtt_num[ss])>dtt_sj){
				Ncal2++;
				logdt.crystal_num+=logdt.dtt_num[ss];
				logdt.crystal_pro+=double(logdt.dtt_num[ss])/logdt.Ncen;
				logdt.dttc_avg+=logdt.dtt_num[ss];
			}
		}
	}
	logdt.dtt_avg=logdt.dtt_avg/Ncal1;
	if(Ncal2!=0){
		logdt.dttc_avg=logdt.dttc_avg/Ncal2;
	}
	
	if(mode_v==1){
		volume_cry2(vt, bt, logdt, dtt_sj, v_dump, out8, Nt);
	}
	
	if(mode_chain_dtt==1){
		conformation_dtt(bt, vt, logdt, dtt_sj);
	}

	if(mode_c==1){
		cluster_2(dtt_j, dtt_sj, cluster_p, vt, bt, logdt);

		out6<<"ITEM: TIMESTEP"<<endl;
		out6<<Nt<<endl;
		out6<<"ITEM: NUMBER OF ATOMS"<<endl;
		out6<<logdt.cluster_num<<endl;
		out6<<"ITEM: BOX BOUNDS xy xz yz pp pp pp"<<endl;
		out6<<logdt.box_xlo<<" "<<logdt.box_xhi<<" "<<0.0<<endl;
		out6<<logdt.box_ylo<<" "<<logdt.box_yhi<<" "<<0.0<<endl;
		out6<<logdt.box_zlo<<" "<<logdt.box_zhi<<" "<<0.0<<endl;
		out6<<"ITEM: ATOMS id x y z vx vy vz rg rgx rgy rgz size"<<endl;
		logdt.cluster_size_max=0;
		logdt.cluster_gy_max=0;
		logdt.cluster_gy_x_max=0;
		logdt.cluster_gy_y_max=0;
		logdt.cluster_gy_z_max=0;
		logdt.cluster_gy_avg=0;
		logdt.cluster_gy_x_avg=0;
		logdt.cluster_gy_y_avg=0;
		logdt.cluster_gy_z_avg=0;
		logdt.cluster_size_avg=0;
		logdt.cluster_size_avgW=0;
		for(int q=0;q<logdt.cluster_num;++q){
			out6 << q+1 << " " << logdt.cluster_cx[q] << " " << logdt.cluster_cy[q] << " " << logdt.cluster_cz[q] << " " << logdt.cluster_vx[q] << " " << logdt.cluster_vy[q] << " " << logdt.cluster_vz[q] << " " << logdt.cluster_gy[q] << " " << logdt.cluster_gyx[q] << " " << logdt.cluster_gyy[q] << " " << logdt.cluster_gyz[q] << " " << logdt.cluster_size[q]<<endl;
			/*logdt.cluster_gy_avg+=logdt.cluster_gy[q]/double(logdt.cluster_num);
			logdt.cluster_gy_x_avg+=logdt.cluster_gyx[q]/double(logdt.cluster_num);
			logdt.cluster_gy_y_avg+=logdt.cluster_gyy[q]/double(logdt.cluster_num);
			logdt.cluster_gy_z_avg+=logdt.cluster_gyz[q]/double(logdt.cluster_num);*/
			
			logdt.cluster_gy_avg+=logdt.cluster_gy[q]*logdt.cluster_size[q]/double(logdt.Ncl_all);
			logdt.cluster_gy_x_avg+=logdt.cluster_gyx[q]*logdt.cluster_size[q]/double(logdt.Ncl_all);
			logdt.cluster_gy_y_avg+=logdt.cluster_gyy[q]*logdt.cluster_size[q]/double(logdt.Ncl_all);
			logdt.cluster_gy_z_avg+=logdt.cluster_gyz[q]*logdt.cluster_size[q]/double(logdt.Ncl_all);
			logdt.cluster_size_avg+=double(logdt.cluster_size[q])/double(logdt.cluster_num);
			logdt.cluster_size_avgW+=double(logdt.cluster_size[q]*logdt.cluster_size[q])/double(logdt.Ncl_all);
			
			if(logdt.cluster_size[q]>logdt.cluster_size_max) logdt.cluster_size_max=logdt.cluster_size[q];
			if(logdt.cluster_gy[q]>logdt.cluster_gy_max) logdt.cluster_gy_max=logdt.cluster_gy[q];
			if(logdt.cluster_gyx[q]>logdt.cluster_gy_x_max) logdt.cluster_gy_x_max=logdt.cluster_gyx[q];
			if(logdt.cluster_gyy[q]>logdt.cluster_gy_y_max) logdt.cluster_gy_y_max=logdt.cluster_gyy[q];
			if(logdt.cluster_gyz[q]>logdt.cluster_gy_z_max) logdt.cluster_gy_z_max=logdt.cluster_gyz[q];
		}
	}

}