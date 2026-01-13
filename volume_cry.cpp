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
void volume_cry1(vector<DATA>& vt, LOG& logdt, double sop_j, int v_dump, ofstream& out7, int Nt){
	
	double lxx=logdt.box_x/logdt.probe;
	double lyy=logdt.box_y/logdt.probe;
	double lzz=logdt.box_z/logdt.probe;
	int lx=int(lxx);
	int ly=int(lyy);
	int lz=int(lzz);
	
	int mx = lx + 1;
	int my = ly + 1;
	int mz = lz + 1;
	vector<int> matrix(mx * my * mz, 0);
	auto idx = [&](int i, int j, int k) { return (i * my + j) * mz + k; };
	
	double ixx,iyy,izz;
	int Nc=0,ix,iy,iz;
	for(int k=0;k<logdt.Ncen;++k){
		if(vt[k].sop>=sop_j){
			ixx=(vt[k].cenx-logdt.box_xlo)/logdt.probe;
			iyy=(vt[k].ceny-logdt.box_ylo)/logdt.probe;
			izz=(vt[k].cenz-logdt.box_zlo)/logdt.probe;
			ix=int(ixx);
			iy=int(iyy);
			iz=int(izz);
			if(ix<0) ix=0;
			if(ix>lx) ix=lx;
			if(iy<0) iy=0;
			if(iy>ly) iy=ly;
			if(iz<0) iz=0;
			if(iz>lz) iz=lz;
			matrix[idx(ix, iy, iz)]=1;
		}
	}
	
	if(v_dump==1){
		out7<<"ITEM: TIMESTEP"<<endl;
		out7<<Nt<<endl;
		out7<<"ITEM: NUMBER OF ATOMS"<<endl;
		out7<<lx*ly*lz<<endl;
		out7<<"ITEM: BOX BOUNDS xy xz yz pp pp pp"<<endl;
		out7<<logdt.box_xlo<<" "<<logdt.box_xhi<<" "<<0.0<<endl;
		out7<<logdt.box_ylo<<" "<<logdt.box_yhi<<" "<<0.0<<endl;
		out7<<logdt.box_zlo<<" "<<logdt.box_zhi<<" "<<0.0<<endl;
		out7<<"ITEM: ATOMS id type x y z"<<endl;
	}
	
	for(int i=0;i<lx;i++){
		for(int j=0;j<ly;j++){
			for(int s=0;s<lz;s++){
				if(matrix[idx(i, j, s)]==1){
					Nc++;
				}
				if(v_dump==1){
					out7<<i*ly*lz+j*lz+s+1<<" "<<matrix[idx(i, j, s)]<<" "<<(i+0.5)*logdt.probe+logdt.box_xlo<<" "<<(j+0.5)*logdt.probe+logdt.box_ylo<<" "<<(s+0.5)*logdt.probe+logdt.box_zlo<<endl;
				}
			}
		}
	}

	logdt.sop_crystal_pro_v=Nc/double(lx*ly*lz);

}

void volume_cry2(vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, double dtt_sj, int v_dump, ofstream& out8, int Nt){
	
	double lxx=logdt.box_x/logdt.probe;
	double lyy=logdt.box_y/logdt.probe;
	double lzz=logdt.box_z/logdt.probe;
	int lx=int(lxx);
	int ly=int(lyy);
	int lz=int(lzz);
	
	int mx = lx + 1;
	int my = ly + 1;
	int mz = lz + 1;
	vector<int> matrix(mx * my * mz, 0);
	auto idx = [&](int i, int j, int k) { return (i * my + j) * mz + k; };
	
	double ixx,iyy,izz;
	int Nc=0,kk,ss,id,ix,iy,iz;
	for(int k=0;k<logdt.Ncen;++k){
		kk=vt[k].btid_k;
		ss=vt[k].btid_s;
		id=bt[kk][ss].dttid;
		if(double(logdt.dtt_num[id])>dtt_sj){
			ixx=(vt[k].cenx-logdt.box_xlo)/logdt.probe;
			iyy=(vt[k].ceny-logdt.box_ylo)/logdt.probe;
			izz=(vt[k].cenz-logdt.box_zlo)/logdt.probe;
			ix=int(ixx);
			iy=int(iyy);
			iz=int(izz);
			if(ix<0) ix=0;
			if(ix>lx) ix=lx;
			if(iy<0) iy=0;
			if(iy>ly) iy=ly;
			if(iz<0) iz=0;
			if(iz>lz) iz=lz;
			matrix[idx(ix, iy, iz)]=1;
		}
	}
	
	if(v_dump==1){
		out8<<"ITEM: TIMESTEP"<<endl;
		out8<<Nt<<endl;
		out8<<"ITEM: NUMBER OF ATOMS"<<endl;
		out8<<lx*ly*lz<<endl;
		out8<<"ITEM: BOX BOUNDS xy xz yz pp pp pp"<<endl;
		out8<<logdt.box_xlo<<" "<<logdt.box_xhi<<" "<<0.0<<endl;
		out8<<logdt.box_ylo<<" "<<logdt.box_yhi<<" "<<0.0<<endl;
		out8<<logdt.box_zlo<<" "<<logdt.box_zhi<<" "<<0.0<<endl;
		out8<<"ITEM: ATOMS id type x y z"<<endl;
	}
	
	for(int i=0;i<lx;i++){
		for(int j=0;j<ly;j++){
			for(int s=0;s<lz;s++){
				if(matrix[idx(i, j, s)]==1){
					Nc++;
				}
				if(v_dump==1){
					out8<<i*ly*lz+j*lz+s+1<<" "<<matrix[idx(i, j, s)]<<" "<<(i+0.5)*logdt.probe+logdt.box_xlo<<" "<<(j+0.5)*logdt.probe+logdt.box_ylo<<" "<<(s+0.5)*logdt.probe+logdt.box_zlo<<endl;
				}
			}
		}
	}

	logdt.crystal_pro_v=Nc/double(lx*ly*lz);
	
}
