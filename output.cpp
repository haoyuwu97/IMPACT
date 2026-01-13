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
void output(int Ntemp, int mode_sop, int mode_dtt, int mode_c_sop, int mode_c_dtt, ofstream& out1, ofstream& out2, ofstream& out3, ofstream& out5, vector<DATA>& vt, vector<vector<DATA2>>& bt, LOG& logdt, int mode_v, int mode_chain_sop, int mode_chain_dtt){
	
	out2<<"ITEM: TIMESTEP"<<endl;
	out2<<Ntemp<<endl;
	out2<<"ITEM: NUMBER OF ATOMS"<<endl;
	out2<<logdt.Ncen<<endl;
	out2<<"ITEM: BOX BOUNDS xy xz yz pp pp pp"<<endl;
	out2<<logdt.box_xlo<<" "<<logdt.box_xhi<<" "<<0.0<<endl;
	out2<<logdt.box_ylo<<" "<<logdt.box_yhi<<" "<<0.0<<endl;
	out2<<logdt.box_zlo<<" "<<logdt.box_zhi<<" "<<0.0<<endl;
	
	out2<<"ITEM: ATOMS id mol type x y z ix iy iz ";
	if(mode_sop==1){
		out2<<"sop ";
		if(mode_c_sop==1){
			out2<<"cluster(sop) cluster_size(sop) ";
		}
		if(mode_chain_sop==1){
			out2<<"mol_sop ";
		}
	}
	if(mode_dtt==1){
		out2<<"dtt ";
		if(mode_c_dtt==1){
			out2<<"cluster(dtt) cluster_size(dtt) ";
		}
		if(mode_chain_dtt==1){
			out2<<"mol_dtt ";
		}
	}
	out2<<endl;
	int kk,ss,id;
	for(int k=0;k<logdt.Ncen;++k){
		out2<< vt[k].cenid << " " << vt[k].mol << " " << vt[k].ty << " " << vt[k].cenx << " " << vt[k].ceny << " " << vt[k].cenz << " " << vt[k].cenix << " " << vt[k].ceniy << " " << vt[k].ceniz << " ";
		if(mode_sop==1){
			out2<< vt[k].sop << " ";
			if(mode_c_sop==1){
				out2<< vt[k].site << " " << vt[k].site_size<<" ";
			}
			if(mode_chain_sop==1){
				out2<< vt[k].moltype_sop<<" ";
			}
		}
		if(mode_dtt==1){
			kk=vt[k].btid_k;
			ss=vt[k].btid_s;
			id=bt[kk][ss].dttid;
			out2<< logdt.dtt_num[id] << " ";
			if(mode_c_dtt==1){
				out2<< bt[kk][ss].site << " " << bt[kk][ss].site_size<<" ";
			}
			if(mode_chain_dtt==1){
				out2<< vt[k].moltype_dtt<<" ";
			}
		}
		out2<<endl;
	}

	
	if(logdt.Ncount==1){
		out1<<"Time ";
		if(mode_sop==1){
			out1<<"Ncry(sop) Cry(sop)% ";
			if(mode_v==1){
				out1<<"Cry_vol(sop)% ";
			}
			out1<<"<sop> ";
			if(mode_chain_sop==1){
				out1<<"Ntail_s Ntie_s Nloop_s <Ltail>_s <Ltie>_s <Lloop>_s ";
			}
		}
		if(mode_dtt==1){
			out1<<"Ncry(dtt) Cry(dtt)% ";
			if(mode_v==1){
				out1<<"Cry_vol(dtt)% ";
			}
			out1<<"<dtt> <dttc> ";
			if(mode_chain_dtt==1){
				out1<<"Ntail_d Ntie_d Nloop_d <Ltail>_d <Ltie>_d <Lloop>_d ";
			}
		}
		out1<<endl;
	}
	out1<< Ntemp <<" ";	
	if(mode_sop==1){
		out1<< logdt.sop_crystal_num <<" "<< logdt.sop_crystal_pro<<" ";
		if(mode_v==1){
			out1<< logdt.sop_crystal_pro_v <<" ";
		}
		out1<< logdt.sop_avg <<" ";
		if(mode_chain_sop==1){
			out1<< logdt.Ntail_s <<" "<< logdt.Ntie_s <<" "<< logdt.Nloop_s <<" "<< logdt.Ltail_s <<" "<< logdt.Ltie_s <<" "<< logdt.Lloop_s <<" ";
		}
	}
	if(mode_dtt==1){
		out1<< logdt.crystal_num <<" "<< logdt.crystal_pro<<" ";
		if(mode_v==1){
			out1<< logdt.crystal_pro_v <<" ";
		}
		out1<< logdt.dtt_avg <<" "<< logdt.dttc_avg<<" ";
		if(mode_chain_dtt==1){
			out1<< logdt.Ntail_d <<" "<< logdt.Ntie_d <<" "<< logdt.Nloop_d <<" "<< logdt.Ltail_d <<" "<< logdt.Ltie_d <<" "<< logdt.Lloop_d <<" ";
		}
	}
	out1<<endl;
	
	if(mode_c_sop==1){
		if(logdt.Ncount==1){
			out3 << "Nc" <<" "<< "Nc(>10)" <<" "<<"<Cn>" <<" "<< "<Cn>_w" <<" "<< "Cn_max" <<" "<< "<Rg>" <<" "<< "<Rgx>" <<" "<< "<Rgy>" <<" "<< "<Rgz>" <<" "<< "Rg_max" <<" "<< "Rgx_max" <<" "<< "Rgy_max" <<" "<< "Rgz_max" << endl;
		}
		out3 << logdt.sop_cluster_num <<" "<< logdt.sop_cluster_num2 <<" "<< logdt.sop_cluster_size_avg <<" "<< logdt.sop_cluster_size_avgW <<" "<< logdt.sop_cluster_size_max <<" "<< logdt.sop_cluster_gy_avg <<" "<< logdt.sop_cluster_gy_x_avg <<" "<< logdt.sop_cluster_gy_y_avg <<" "<< logdt.sop_cluster_gy_z_avg <<" "<< logdt.sop_cluster_gy_max <<" "<< logdt.sop_cluster_gy_x_max <<" "<< logdt.sop_cluster_gy_y_max <<" "<< logdt.sop_cluster_gy_z_max << endl;
	}
	if(mode_c_dtt==1){
		if(logdt.Ncount==1){
			out5 << "Nc" <<" "<< "Nc(>10)" <<" "<<"<Cn>" <<" "<< "<Cn>_w" <<" "<< "Cn_max" <<" "<< "<Rg>" <<" "<< "<Rgx>" <<" "<< "<Rgy>" <<" "<< "<Rgz>" <<" "<< "Rg_max" <<" "<< "Rgx_max" <<" "<< "Rgy_max" <<" "<< "Rgz_max" << endl;
		}
		out5 << logdt.cluster_num <<" "<< logdt.cluster_num2 <<" "<< logdt.cluster_size_avg <<" "<< logdt.cluster_size_avgW <<" "<< logdt.cluster_size_max <<" "<< logdt.cluster_gy_avg <<" "<< logdt.cluster_gy_x_avg <<" "<< logdt.cluster_gy_y_avg <<" "<< logdt.cluster_gy_z_avg <<" "<< logdt.cluster_gy_max <<" "<< logdt.cluster_gy_x_max <<" "<< logdt.cluster_gy_y_max <<" "<< logdt.cluster_gy_z_max << endl;
	}
	
}