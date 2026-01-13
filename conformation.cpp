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
void conformation_sop(vector<DATA>& vt, LOG& logdt, double sop_j){
	
	vector<int> segment(logdt.Ncen + 10, 0);
	vector<int> segment_type(logdt.Ncen, 0);
	vector<int> head(logdt.Ncen, 0);
	int tail;
	int Ns=-1;
	double vec_A,vec_B,vec_AB,CosD;

	for(int k=0;k<logdt.Ncen;++k){
		if(vt[k].sop>=sop_j){
			segment[k]=-1;
		}else{
			if(k==0||(k!=0&&vt[k].mol!=vt[k-1].mol)){
				Ns++;
				segment[k]=Ns;
				segment_type[Ns]=2;
			}
			if(k!=0&&vt[k].mol==vt[k-1].mol){
				if(vt[k].mol!=vt[k+1].mol){
					segment[k]=Ns;
					segment_type[Ns]=2;
				}
				if(vt[k].mol==vt[k+1].mol){
					if(vt[k-1].sop>=sop_j){
						Ns++;
						segment[k]=Ns;
						segment_type[Ns]=3;
						head[Ns]=k-1;
					}if(vt[k+1].sop>=sop_j){
						segment[k]=Ns;
						if(segment_type[Ns]==3){
  						tail=k+1;			
  						vec_A=vt[head[Ns]].x*vt[head[Ns]].x + vt[head[Ns]].y*vt[head[Ns]].y + vt[head[Ns]].z*vt[head[Ns]].z;
  						vec_A=sqrt(vec_A);
  						vec_AB=vt[head[Ns]].x*vt[tail].x + vt[head[Ns]].y*vt[tail].y + vt[head[Ns]].z*vt[tail].z;
  						vec_B=vt[tail].x*vt[tail].x + vt[tail].y*vt[tail].y + vt[tail].z*vt[tail].z;
  						vec_B=sqrt(vec_B);
  						CosD=vec_AB/(vec_A*vec_B);
  						if(CosD>=0){
  							segment_type[Ns]=4;
  						}
						}
					}if(vt[k-1].sop<sop_j&&vt[k+1].sop<sop_j){
						segment[k]=segment[k-1];
					}
				}
			}
		}
	}
	int Nss=0;
	logdt.Ntail_s=0;
	logdt.Ntie_s=0;
	logdt.Nloop_s=0;

	for(int k=0;k<logdt.Ncen;++k){
		if(segment[k]==-1){
			vt[k].moltype_sop=1;
		}else if(segment[k]!=-1){
			Nss=segment[k];
			vt[k].moltype_sop=segment_type[Nss];
			if(segment_type[Nss]==2){
				logdt.Ntail_s++;
			}
			if(segment_type[Nss]==3){
				logdt.Nloop_s++;
			}
			if(segment_type[Nss]==4){
				logdt.Ntie_s++;
			}
		}
	}
	
	int n2=0,n3=0,n4=0;
	for(int i=0;i<Ns+1;i++){
		if(segment_type[i]==2) n2++;
		if(segment_type[i]==3) n3++;
		if(segment_type[i]==4) n4++;
	}

	if(n2!=0){
		logdt.Ltail_s=double(logdt.Ntail_s)/n2;
	}else{
		logdt.Ltail_s=0.0;
	}

	if(n3!=0){
		logdt.Lloop_s=double(logdt.Nloop_s)/n3;
	}else{
		logdt.Lloop_s=0.0;
	}

	if(n4!=0){
		logdt.Ltie_s=double(logdt.Ntie_s)/n4;
	}else{
		logdt.Ltie_s=0.0;
	}
	
}

void conformation_dtt(vector<vector<DATA2>>& bt, vector<DATA>& vt, LOG& logdt, double dtt_sj){
	
	vector<int> segment(logdt.Ncen + 10, 0);
	vector<int> segment_type(logdt.Ncen, 0);
	vector<int> head(logdt.Ncen, 0);
	int tail;
	int Ns=-1;
	double vec_A,vec_B,vec_AB,CosD;
	
	int k0,s0,id0,k1,s1,id1,k2,s2,id2;

	for(int k=0;k<logdt.Ncen;++k){
	
		k0=vt[k].btid_k;
		s0=vt[k].btid_s;
		id0=bt[k0][s0].dttid;
		
		if(k!=0){
  		k1=vt[k-1].btid_k;
  		s1=vt[k-1].btid_s;
  		id1=bt[k1][s1].dttid;
		}
			
		k2=vt[k+1].btid_k;
		s2=vt[k+1].btid_s;
		id2=bt[k2][s2].dttid;
			
		if(double(logdt.dtt_num[id0])>dtt_sj){
			segment[k]=-1;
		}else{
			if(k==0||(k!=0&&vt[k].mol!=vt[k-1].mol)){
				Ns++;
				segment[k]=Ns;
				segment_type[Ns]=2;
			}else if(vt[k].mol==vt[k-1].mol){
				if(vt[k].mol!=vt[k+1].mol){
					segment[k]=Ns;
					segment_type[Ns]=2;
				}
				if(vt[k].mol==vt[k+1].mol){
					if(double(logdt.dtt_num[id1])>dtt_sj){
						Ns++;
						segment[k]=Ns;
						segment_type[Ns]=3;
						head[Ns]=k-1;

					}if(double(logdt.dtt_num[id2])>dtt_sj){
						segment[k]=Ns;
						if(segment_type[Ns]==3){                                                
  						tail=k+1;                            
              vec_A=vt[head[Ns]].x*vt[head[Ns]].x + vt[head[Ns]].y*vt[head[Ns]].y + vt[head[Ns]].z*vt[head[Ns]].z;
  						vec_A=sqrt(vec_A);
  						vec_AB=vt[head[Ns]].x*vt[tail].x + vt[head[Ns]].y*vt[tail].y + vt[head[Ns]].z*vt[tail].z;
  						vec_B=vt[tail].x*vt[tail].x + vt[tail].y*vt[tail].y + vt[tail].z*vt[tail].z;
  						vec_B=sqrt(vec_B);
  						CosD=vec_AB/(vec_A*vec_B);
  						if(CosD>=0){
  							segment_type[Ns]=4;
  						}
            }					
					}if(double(logdt.dtt_num[id1])<=dtt_sj&&double(logdt.dtt_num[id2])<=dtt_sj){
						segment[k]=segment[k-1];
					}
				}
			}
		}
	}

	int Nss=0;
	logdt.Ntail_d=0;
	logdt.Ntie_d=0;
	logdt.Nloop_d=0;
	for(int k=0;k<logdt.Ncen;++k){
		if(segment[k]==-1){
			vt[k].moltype_dtt=1;
		}else if(segment[k]!=-1){
			Nss=segment[k];
			vt[k].moltype_dtt=segment_type[Nss];
			if(segment_type[Nss]==2){
				logdt.Ntail_d++;
			}
			if(segment_type[Nss]==3){
				logdt.Nloop_d++;
			}
			if(segment_type[Nss]==4){
				logdt.Ntie_d++;
			}
		}
	}
	
	int n2=0,n3=0,n4=0;
	for(int i=0;i<Ns+1;i++){
		if(segment_type[i]==2) n2++;
		if(segment_type[i]==3) n3++;
		if(segment_type[i]==4) n4++;
	}

	if(n2!=0){
		logdt.Ltail_d=double(logdt.Ntail_d)/n2;
	}else{
		logdt.Ltail_d=0.0;
	}

	if(n3!=0){
		logdt.Lloop_d=double(logdt.Nloop_d)/n3;
	}else{
		logdt.Lloop_d=0.0;
	}

	if(n4!=0){
		logdt.Ltie_d=double(logdt.Ntie_d)/n4;
	}else{
		logdt.Ltie_d=0.0;
	}
	
}