#include "crystal.h"

#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdio.h>
#include <vector>

using namespace std;
void cluster_1(double sop_j, double cluster_p[], vector<DATA>& vt, LOG& logdt){
	int nmtotal=0,temp;
	int totcluster=1;
	int kcluster,clustertrue;
	logdt.reset_clusters(logdt.Ncen + 1);
	logdt.ncluster[0]=1;
	vector<int> clustnb(logdt.Ncen, -1);
	vector<int> cluster2(logdt.Ncen, 0);
	cluster2[0]=1;
	double rx,ry,rz,vec_AB,vec_A,vec_B,CosD,p2;
	
	double Rdf=cluster_p[1];
	double cluster_judge=cluster_p[2];
	logdt.reset_sop_clusters(logdt.Ncen + 1);
	int ncen = logdt.Ncen;
	static vector<int> head;
	static vector<int> next;
	static vector<int> cell_index;
	static vector<double> prev_cenx;
	static vector<double> prev_ceny;
	static vector<double> prev_cenz;
	static double prev_box_x = 0.0;
	static double prev_box_y = 0.0;
	static double prev_box_z = 0.0;
	static double prev_Rdf = 0.0;
	bool reuse_cells = false;
	int nx = (Rdf > 0.0) ? int(logdt.box_x / Rdf) : 1;
	int ny = (Rdf > 0.0) ? int(logdt.box_y / Rdf) : 1;
	int nz = (Rdf > 0.0) ? int(logdt.box_z / Rdf) : 1;
	if (nx < 1) nx = 1;
	if (ny < 1) ny = 1;
	if (nz < 1) nz = 1;
	int ncell = nx * ny * nz;
	if (!cell_index.empty() && prev_cenx.size() == static_cast<size_t>(ncen) && prev_box_x == logdt.box_x && prev_box_y == logdt.box_y && prev_box_z == logdt.box_z && prev_Rdf == Rdf) {
		double threshold = (Rdf > 0.0) ? 0.25 * Rdf : 0.0;
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
		prev_Rdf = Rdf;
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

	for(int k=1;k<logdt.Ncen;++k){
		nmtotal++;
		vec_A=vt[nmtotal].x*vt[nmtotal].x + vt[nmtotal].y*vt[nmtotal].y + vt[nmtotal].z*vt[nmtotal].z;
		vec_A=sqrt(vec_A);
		int nbtotal=-1;
		int clustmin=10000000;
		std::fill(clustnb.begin(), clustnb.end(), -1);
		int cell = cell_index[nmtotal];
		int cx = cell / (ny * nz);
		int cy = (cell / nz) % ny;
		int cz = cell % nz;
		for(int dx=-1; dx<=1; ++dx){
			int nx_i = (cx + dx + nx) % nx;
			for(int dy=-1; dy<=1; ++dy){
				int ny_i = (cy + dy + ny) % ny;
				for(int dz=-1; dz<=1; ++dz){
					int nz_i = (cz + dz + nz) % nz;
					int neighbor = (nx_i * ny + ny_i) * nz + nz_i;
					for(int i=head[neighbor]; i!=-1; i=next[i]){
						if(i>=nmtotal){
							continue;
						}
						rx=vt[nmtotal].cenx-vt[i].cenx;
						ry=vt[nmtotal].ceny-vt[i].ceny;
						rz=vt[nmtotal].cenz-vt[i].cenz;
						if(rx>logdt.box_x/2){rx=rx-logdt.box_x;}
						if(rx<-logdt.box_x/2){rx=rx+logdt.box_x;}
						if(ry>logdt.box_y/2){ry=ry-logdt.box_y;}
						if(ry<-logdt.box_y/2){ry=ry+logdt.box_y;}
						if(rz>logdt.box_z/2){rz=rz-logdt.box_z;}
						if(rz<-logdt.box_z/2){rz=rz+logdt.box_z;}
						rx=rx*rx+ry*ry+rz*rz;
						rx=sqrt(rx);
						if(rx<=Rdf && vt[nmtotal].sop>=sop_j && vt[i].sop>=sop_j){
							vec_AB=vt[nmtotal].x*vt[i].x + vt[nmtotal].y*vt[i].y + vt[nmtotal].z*vt[i].z;
							vec_B=vt[i].x*vt[i].x + vt[i].y*vt[i].y + vt[i].z*vt[i].z;
							vec_B=sqrt(vec_B);
							CosD=vec_AB/(vec_A*vec_B);
							p2=(3*CosD*CosD-1)/2;
							if(p2>=cluster_judge){
								kcluster=cluster2[i];
								clustertrue=findproclust(kcluster,logdt);
								bool seen=false;
								for(int s=0;s<=nbtotal;++s){
									if(clustertrue==clustnb[s]){
										seen=true;
										break;
									}
								}
								if(seen){
									continue;
								}
								nbtotal++;
								clustnb[nbtotal]=clustertrue;
								clustmin=min(clustertrue, clustmin);
							}
						}
					}
				}
			}
		}
 
		
		if(nbtotal==-1){
			totcluster++;
			cluster2[nmtotal]=totcluster;
			logdt.ncluster[totcluster]=1;
			goto part2;
		}else{
			cluster2[nmtotal]=clustmin;
			temp=0;
			for(int m=0;m<=nbtotal;++m){
				clustertrue=clustnb[m];
				temp+=logdt.ncluster[clustertrue];
			}
			logdt.ncluster[clustmin]=temp;
			logdt.ncluster[clustmin]+=1;
			for(int m=0;m<=nbtotal;++m){
				if(clustnb[m]!=clustmin){
					int kk=clustnb[m];
					logdt.ncluster[kk]=-clustmin;
				}
			}
		}
		part2:continue;
	}
	
	vector<int> nl(nmtotal, 0);
	vector<int> site(nmtotal, 0);

	vector<double> id_centroidx0(nmtotal, 0.0);
	vector<double> id_centroidy0(nmtotal, 0.0);
	vector<double> id_centroidz0(nmtotal, 0.0);
	vector<double> id_centroidx1(nmtotal, 0.0);
	vector<double> id_centroidy1(nmtotal, 0.0);
	vector<double> id_centroidz1(nmtotal, 0.0);
	vector<double> id_vx(nmtotal, 0.0);
	vector<double> id_vy(nmtotal, 0.0);
	vector<double> id_vz(nmtotal, 0.0);
	vector<double> id_gyration(nmtotal, 0.0);
	vector<double> id_gyration_gyx2(nmtotal, 0.0);
	vector<double> id_gyration_gyy2(nmtotal, 0.0);
	vector<double> id_gyration_gyz2(nmtotal, 0.0);
	for(int i=0;i<nmtotal;++i){
		double xx=vt[i].cenx;
		double yy=vt[i].ceny;
		double zz=vt[i].cenz;
		
		kcluster=cluster2[i];
		int idc=findproclust(kcluster,logdt);
		site[i]=idc;
		nl[idc]++;
		if(id_centroidx0[idc]==0.0) id_centroidx0[idc]=xx;
		if(id_centroidy0[idc]==0.0) id_centroidy0[idc]=yy;
		if(id_centroidz0[idc]==0.0) id_centroidz0[idc]=zz;
		
		if((xx-id_centroidx0[idc])>0.5*logdt.box_x){
			xx=xx-logdt.box_x;
		}
		if((xx-id_centroidx0[idc])<-0.5*logdt.box_x){
			xx=xx+logdt.box_x;
		}
		if((yy-id_centroidy0[idc])>0.5*logdt.box_y){
			yy=yy-logdt.box_y;
		}
		if((yy-id_centroidy0[idc])<-0.5*logdt.box_y){
			yy=yy+logdt.box_y;
		}
		if((zz-id_centroidz0[idc])>0.5*logdt.box_z){
			zz=zz-logdt.box_z;
		}
		if((zz-id_centroidz0[idc])<-0.5*logdt.box_z){
			zz=zz+logdt.box_z;
		}
		id_centroidx1[idc]+=xx;
		id_centroidy1[idc]+=yy;
		id_centroidz1[idc]+=zz;
		
		double rr=sqrt(vt[i].x*vt[i].x+vt[i].y*vt[i].y+vt[i].z*vt[i].z);
		id_vx[idc]=vt[i].x/rr;
		id_vy[idc]=vt[i].y/rr;
		id_vz[idc]=vt[i].z/rr;
	}
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>0){
			
			id_centroidx1[i]=id_centroidx1[i]/nl[i];
			id_centroidy1[i]=id_centroidy1[i]/nl[i];
			id_centroidz1[i]=id_centroidz1[i]/nl[i];

			if(id_centroidx1[i]<logdt.box_xlo) id_centroidx1[i]+=logdt.box_x;
			if(id_centroidx1[i]>logdt.box_xhi) id_centroidx1[i]-=logdt.box_x;
			
			if(id_centroidy1[i]<logdt.box_ylo) id_centroidy1[i]+=logdt.box_y;
			if(id_centroidy1[i]>logdt.box_yhi) id_centroidy1[i]-=logdt.box_y;
			
			if(id_centroidz1[i]<logdt.box_zlo) id_centroidz1[i]+=logdt.box_z;
			if(id_centroidz1[i]>logdt.box_zhi) id_centroidz1[i]-=logdt.box_z;
			
		}
	}
	double gy_x, gy_y, gy_z;
	for(int i=0;i<nmtotal;++i){	
		
		gy_x=vt[i].cenx-id_centroidx1[site[i]];
		gy_y=vt[i].ceny-id_centroidy1[site[i]];
		gy_z=vt[i].cenz-id_centroidz1[site[i]];
		if(gy_x>logdt.box_x/2){gy_x=gy_x-logdt.box_x;}
		if(gy_x<-logdt.box_x/2){gy_x=gy_x+logdt.box_x;}
		if(gy_y>logdt.box_y/2){gy_y=gy_y-logdt.box_y;}
		if(gy_y<-logdt.box_y/2){gy_y=gy_y+logdt.box_y;}
		if(gy_z>logdt.box_z/2){gy_z=gy_z-logdt.box_z;}
		if(gy_z<-logdt.box_z/2){gy_z=gy_z+logdt.box_z;}
		id_gyration_gyx2[site[i]]=id_gyration_gyx2[site[i]]+(gy_x*gy_x)/nl[site[i]];
		id_gyration_gyy2[site[i]]=id_gyration_gyy2[site[i]]+(gy_y*gy_y)/nl[site[i]];
		id_gyration_gyz2[site[i]]=id_gyration_gyz2[site[i]]+(gy_z*gy_z)/nl[site[i]];
		id_gyration[site[i]]=id_gyration[site[i]]+(gy_x*gy_x+gy_y*gy_y+gy_z*gy_z)/nl[site[i]];
		if(nl[site[i]]>1){
			vt[i].site=site[i]+1;
		}else{
			vt[i].site=0;
		}
		
		vt[i].site_size=nl[site[i]];
	}
	int cl_last=0;
	int cl_last2=0;
	int cl_all=0;
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>1){
			logdt.sop_cluster_gy[cl_last]=sqrt(id_gyration[i]);
			logdt.sop_cluster_gyx[cl_last]=sqrt(id_gyration_gyx2[i]);
			logdt.sop_cluster_gyy[cl_last]=sqrt(id_gyration_gyy2[i]);
			logdt.sop_cluster_gyz[cl_last]=sqrt(id_gyration_gyz2[i]);
			logdt.sop_cluster_size[cl_last]=nl[i];
			logdt.sop_cluster_cx[cl_last]=id_centroidx1[i];
			logdt.sop_cluster_cy[cl_last]=id_centroidy1[i];
			logdt.sop_cluster_cz[cl_last]=id_centroidz1[i];
			logdt.sop_cluster_vx[cl_last]=id_vx[i];
			logdt.sop_cluster_vy[cl_last]=id_vy[i];
			logdt.sop_cluster_vz[cl_last]=id_vz[i];
			cl_last++;
			cl_all+=nl[i];
		}
		if(nl[i]>10){
			cl_last2++;
		}
	}
	logdt.sop_cluster_num=cl_last;
	logdt.sop_cluster_num2=cl_last2;
	logdt.Ncl_all_sop=cl_all;
}

void cluster_2(double dtt_j, double dtt_sj, double cluster_p[], vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt){
	int nmtotal=0,temp;
	int totcluster=0;
	int kcluster,clustertrue;
	logdt.reset_clusters(logdt.Ncen + 1);
	logdt.ncluster[0]=1;
	vector<int> clustnb(logdt.Ncen + 1, -1);
	vector<int> cluster2(logdt.Ncen + 1, 0);
	cluster2[0]=1;
	double rx,ry,rz,vec_AB,vec_A,vec_B,CosD,p2;
	
	double Rdf=cluster_p[1];
	double cluster_judge=cluster_p[2];
	logdt.reset_dtt_clusters(logdt.Ncen + 1);
	int ncen = logdt.Ncen;
	static vector<int> head;
	static vector<int> next;
	static vector<int> cell_index;
	static vector<double> prev_cenx;
	static vector<double> prev_ceny;
	static vector<double> prev_cenz;
	static double prev_box_x = 0.0;
	static double prev_box_y = 0.0;
	static double prev_box_z = 0.0;
	static double prev_Rdf = 0.0;
	bool reuse_cells = false;
	int nx = (Rdf > 0.0) ? int(logdt.box_x / Rdf) : 1;
	int ny = (Rdf > 0.0) ? int(logdt.box_y / Rdf) : 1;
	int nz = (Rdf > 0.0) ? int(logdt.box_z / Rdf) : 1;
	if (nx < 1) nx = 1;
	if (ny < 1) ny = 1;
	if (nz < 1) nz = 1;
	int ncell = nx * ny * nz;
	if (!cell_index.empty() && prev_cenx.size() == static_cast<size_t>(ncen) && prev_box_x == logdt.box_x && prev_box_y == logdt.box_y && prev_box_z == logdt.box_z && prev_Rdf == Rdf) {
		double threshold = (Rdf > 0.0) ? 0.25 * Rdf : 0.0;
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
		prev_Rdf = Rdf;
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

	int k1,k2,s1,s2,id1,id2;
	
	for(int k=1;k<logdt.Ncen;++k){
		nmtotal++;
		k1=vt[nmtotal].btid_k;
		s1=vt[nmtotal].btid_s;
		id1=bt[k1][s1].dttid;
		vec_A=vt[nmtotal].x*vt[nmtotal].x + vt[nmtotal].y*vt[nmtotal].y + vt[nmtotal].z*vt[nmtotal].z;
		vec_A=sqrt(vec_A);
		int nbtotal=-1;
		int clustmin=10000000;
		std::fill(clustnb.begin(), clustnb.end(), -1);
		int cell = cell_index[nmtotal];
		int cx = cell / (ny * nz);
		int cy = (cell / nz) % ny;
		int cz = cell % nz;
		for(int dx=-1; dx<=1; ++dx){
			int nx_i = (cx + dx + nx) % nx;
			for(int dy=-1; dy<=1; ++dy){
				int ny_i = (cy + dy + ny) % ny;
				for(int dz=-1; dz<=1; ++dz){
					int nz_i = (cz + dz + nz) % nz;
					int neighbor = (nx_i * ny + ny_i) * nz + nz_i;
					for(int i=head[neighbor]; i!=-1; i=next[i]){
						if(i>=nmtotal){
							continue;
						}
						k2=vt[i].btid_k;
						s2=vt[i].btid_s;
						id2=bt[k2][s2].dttid;
						rx=vt[nmtotal].cenx-vt[i].cenx;
						ry=vt[nmtotal].ceny-vt[i].ceny;
						rz=vt[nmtotal].cenz-vt[i].cenz;
						if(rx>logdt.box_x/2){rx=rx-logdt.box_x;}
						if(rx<-logdt.box_x/2){rx=rx+logdt.box_x;}
						if(ry>logdt.box_y/2){ry=ry-logdt.box_y;}
						if(ry<-logdt.box_y/2){ry=ry+logdt.box_y;}
						if(rz>logdt.box_z/2){rz=rz-logdt.box_z;}
						if(rz<-logdt.box_z/2){rz=rz+logdt.box_z;}
						rx=rx*rx+ry*ry+rz*rz;
						rx=sqrt(rx);
						if(rx<=Rdf && double(logdt.dtt_num[id1])>dtt_sj && double(logdt.dtt_num[id2])>dtt_sj){
							vec_AB=vt[nmtotal].x*vt[i].x + vt[nmtotal].y*vt[i].y + vt[nmtotal].z*vt[i].z;
							vec_B=vt[i].x*vt[i].x + vt[i].y*vt[i].y + vt[i].z*vt[i].z;
							vec_B=sqrt(vec_B);
							CosD=vec_AB/(vec_A*vec_B);
							p2=(3*CosD*CosD-1)/2;
							if(p2>=cluster_judge){
								kcluster=cluster2[i];
								clustertrue=findproclust(kcluster,logdt);
								bool seen=false;
								for(int s=0;s<=nbtotal;++s){
									if(clustertrue==clustnb[s]){
										seen=true;
										break;
									}
								}
								if(seen){
									continue;
								}
								nbtotal++;
								clustnb[nbtotal]=clustertrue;
								clustmin=min(clustertrue, clustmin);
							}
						}
					}
				}
			}
		}

	
		if(nbtotal==-1){
			totcluster++;
			cluster2[nmtotal]=totcluster;
			logdt.ncluster[totcluster]=1;
			goto part2;
		}else{
			cluster2[nmtotal]=clustmin;
			temp=0;
			for(int m=0;m<=nbtotal;++m){
				clustertrue=clustnb[m];
				temp+=logdt.ncluster[clustertrue];
			}
			logdt.ncluster[clustmin]=temp;
			logdt.ncluster[clustmin]+=1;
			for(int m=0;m<=nbtotal;++m){
				if(clustnb[m]!=clustmin){
					int kk=clustnb[m];
					logdt.ncluster[kk]=-clustmin;
				}
			}
		}
		part2:continue;
	}
	
	vector<int> nl(nmtotal, 0);
	vector<int> site(nmtotal, 0);

	vector<double> id_centroidx0(nmtotal, 0.0);
	vector<double> id_centroidy0(nmtotal, 0.0);
	vector<double> id_centroidz0(nmtotal, 0.0);
	vector<double> id_centroidx1(nmtotal, 0.0);
	vector<double> id_centroidy1(nmtotal, 0.0);
	vector<double> id_centroidz1(nmtotal, 0.0);
	vector<double> id_vx(nmtotal, 0.0);
	vector<double> id_vy(nmtotal, 0.0);
	vector<double> id_vz(nmtotal, 0.0);
	vector<double> id_gyration(nmtotal, 0.0);
	vector<double> id_gyration_gyx2(nmtotal, 0.0);
	vector<double> id_gyration_gyy2(nmtotal, 0.0);
	vector<double> id_gyration_gyz2(nmtotal, 0.0);
	for(int i=0;i<nmtotal;++i){
		double xx=vt[i].cenx;
		double yy=vt[i].ceny;
		double zz=vt[i].cenz;
		
		kcluster=cluster2[i];
		int idc=findproclust(kcluster,logdt);
		site[i]=idc;
		nl[idc]++;
		if(id_centroidx0[idc]==0.0) id_centroidx0[idc]=xx;
		if(id_centroidy0[idc]==0.0) id_centroidy0[idc]=yy;
		if(id_centroidz0[idc]==0.0) id_centroidz0[idc]=zz;
		
		if((xx-id_centroidx0[idc])>0.5*logdt.box_x){
			xx=xx-logdt.box_x;
		}
		if((xx-id_centroidx0[idc])<-0.5*logdt.box_x){
			xx=xx+logdt.box_x;
		}
		if((yy-id_centroidy0[idc])>0.5*logdt.box_y){
			yy=yy-logdt.box_y;
		}
		if((yy-id_centroidy0[idc])<-0.5*logdt.box_y){
			yy=yy+logdt.box_y;
		}
		if((zz-id_centroidz0[idc])>0.5*logdt.box_z){
			zz=zz-logdt.box_z;
		}
		if((zz-id_centroidz0[idc])<-0.5*logdt.box_z){
			zz=zz+logdt.box_z;
		}
		id_centroidx1[idc]+=xx;
		id_centroidy1[idc]+=yy;
		id_centroidz1[idc]+=zz;
		
		double rr=sqrt(vt[i].x*vt[i].x+vt[i].y*vt[i].y+vt[i].z*vt[i].z);
		id_vx[idc]=vt[i].x/rr;
		id_vy[idc]=vt[i].y/rr;
		id_vz[idc]=vt[i].z/rr;
	}
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>0){
			
			id_centroidx1[i]=id_centroidx1[i]/nl[i];
			id_centroidy1[i]=id_centroidy1[i]/nl[i];
			id_centroidz1[i]=id_centroidz1[i]/nl[i];

			if(id_centroidx1[i]<logdt.box_xlo) id_centroidx1[i]+=logdt.box_x;
			if(id_centroidx1[i]>logdt.box_xhi) id_centroidx1[i]-=logdt.box_x;
			
			if(id_centroidy1[i]<logdt.box_ylo) id_centroidy1[i]+=logdt.box_y;
			if(id_centroidy1[i]>logdt.box_yhi) id_centroidy1[i]-=logdt.box_y;
			
			if(id_centroidz1[i]<logdt.box_zlo) id_centroidz1[i]+=logdt.box_z;
			if(id_centroidz1[i]>logdt.box_zhi) id_centroidz1[i]-=logdt.box_z;
			
		}
	}
	double gy_x, gy_y, gy_z;
	for(int i=0;i<nmtotal;++i){	
		
		gy_x=vt[i].cenx-id_centroidx1[site[i]];
		gy_y=vt[i].ceny-id_centroidy1[site[i]];
		gy_z=vt[i].cenz-id_centroidz1[site[i]];
		if(gy_x>logdt.box_x/2){gy_x=gy_x-logdt.box_x;}
		if(gy_x<-logdt.box_x/2){gy_x=gy_x+logdt.box_x;}
		if(gy_y>logdt.box_y/2){gy_y=gy_y-logdt.box_y;}
		if(gy_y<-logdt.box_y/2){gy_y=gy_y+logdt.box_y;}
		if(gy_z>logdt.box_z/2){gy_z=gy_z-logdt.box_z;}
		if(gy_z<-logdt.box_z/2){gy_z=gy_z+logdt.box_z;}
		id_gyration_gyx2[site[i]]=id_gyration_gyx2[site[i]]+(gy_x*gy_x)/nl[site[i]];
		id_gyration_gyy2[site[i]]=id_gyration_gyy2[site[i]]+(gy_y*gy_y)/nl[site[i]];
		id_gyration_gyz2[site[i]]=id_gyration_gyz2[site[i]]+(gy_z*gy_z)/nl[site[i]];
		id_gyration[site[i]]=id_gyration[site[i]]+(gy_x*gy_x+gy_y*gy_y+gy_z*gy_z)/nl[site[i]];
		if(nl[site[i]]>1){
			bt[vt[i].btid_k][vt[i].btid_s].site=site[i]+1;
		}else{
			bt[vt[i].btid_k][vt[i].btid_s].site=0;
		}
		
		bt[vt[i].btid_k][vt[i].btid_s].site_size=nl[site[i]];
		
	}
	int cl_last=0;
	int cl_last2=0;
	int cl_all=0;
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>1){
			logdt.cluster_gy[cl_last]=sqrt(id_gyration[i]);
			logdt.cluster_gyx[cl_last]=sqrt(id_gyration_gyx2[i]);
			logdt.cluster_gyy[cl_last]=sqrt(id_gyration_gyy2[i]);
			logdt.cluster_gyz[cl_last]=sqrt(id_gyration_gyz2[i]);
			logdt.cluster_size[cl_last]=nl[i];
			logdt.cluster_cx[cl_last]=id_centroidx1[i];
			logdt.cluster_cy[cl_last]=id_centroidy1[i];
			logdt.cluster_cz[cl_last]=id_centroidz1[i];
			logdt.cluster_vx[cl_last]=id_vx[i];
			logdt.cluster_vy[cl_last]=id_vy[i];
			logdt.cluster_vz[cl_last]=id_vz[i];
			cl_last++;
			cl_all+=nl[i];
		}
		if(nl[i]>10){
			cl_last2++;
		}
	}
	logdt.cluster_num=cl_last;
	logdt.cluster_num2=cl_last2;
	logdt.Ncl_all=cl_all;

}

int findproclust(int sn, LOG& logdt){
	int r,t,k;
	
	r=sn;
	t=-logdt.ncluster[r];
	if(t<0) goto find1;
	
	find2:
	r=t;
	k=t;
	t=-logdt.ncluster[k];
	if(t<0){
		logdt.ncluster[sn]=-r;
	}else{
		goto find2;
	}
	
	find1:
		return r;

}
