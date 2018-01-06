// gasNozzle.cpp : 定义控制台应用程序的入口点。
// copyright by daiyuqiang@dlut.edu.cn
#include <stdio.h>
#include <math.h>
#include "idealGasSteady.h"
//=================================================
int main()
{
	printf("Copyright by daiyuqiang@dlut.edu.cn\n");

	//收缩喷嘴=========================================================================
	double pin=1.01325e5;//上游压力	
	double Tin=293;//上游温度
	double pamb=8e3;//背压
	
	
	//临界压力
    double Pcr=criticalPress(pin);
	printf("Pcr/Pin=%f\t [Pa]\n",Pcr/pin);

	//单位面积的通流能力
	double vo=0.0;
	double To=0.0;
	double mass_dA=massflux_dA_convergingNozzle(pin,Tin,pamb,vo,To);
	printf("mass_dA=%f\t[kg/m^2]\tVout=%f \t[m/s] Tout=%f\t[K]\n",mass_dA,vo,To);
	
	//给定流量求候补面积
	double massflowRate=0.18;
	double Area=massflowRate/mass_dA;
	printf("A=%f\t[m2] d=%f\t[m]\n",Area,sqrt(4.*Area/PI));


	//=========================================================================
	//激波前后关系
	double pt=1.52e6;
	double Tt=500;
	double ma=3.9376;//激波前mach
	double denst=pt/Rg/Tt;
	//激波前静参数
	double p1=pt* ma_to_pst(ma);
	double T1=Tt* ma_to_Tst(ma);
	double dens1=denst* ma_to_Denst(ma);
	double u1=ma*sqrt(Cpv*Rg*T1);
	//激波后静参数
	double p2=p1*p21_shock(ma);
	double T2=T1*T21_shock(ma);
	double dens2=dens1*dens21_shock(ma);
	double deltas=s21_shock(ma);
	double u2=velo21_shock(ma)*u1;

	//R-H关系验证
	//活塞运动产生的激波，波前p1，T1和波后p2已知
	p1=1.013e5;
	T1=288;
	p2=1.1143e5;
	double p21=p2/p1;
	dens1=p1/Rg/T1;
	dens2=dens1*D21_from_p21(p21);
	T2=T1*T21_from_p21(p21);
	double c1=sqrt(Cpv*Rg*T1);
	double c2=sqrt(Cpv*Rg*T2);
	double vs=c1*Vsc1_from_p21(p21);
	double vg=c1*Vgc1_from_p21(p21);
	
	//验证面积比与ma和lam关系
	double ar=0.2/0.01972;
	double lam1,lam2,ma1,ma2;
	Ar_to_ma(ar,0.5,4,ma1,ma2);
	lam1=ma_to_lamda(ma1);
	lam2=ma_to_lamda(ma2);

	Ar_to_lamda(ar,0,2.,lam1,lam2);
	ma1=lamda_to_ma(lam1);
	ma2=lamda_to_ma(lam2);

	//double Artest1=ma_to_Ar(ma1);
	//double Artest2=ma_to_Ar(ma2);
	double Artest1=lamda_to_Ar(lam1);
	double Artest2=lamda_to_Ar(lam2);
	
	//验证缩放喷管
	pt=20.265e5;
	Tt=500;
	double Acr=0.01792;
	double Ao=0.2;
	double Ax=0;
	double mx=0;
	double lamx=0;
	pamb=4.053e5;

	double Pi=pIt_convDivNoz(ma2)*pt;
	double Pii=pIIt_convDivNoz(ma2)*pt;
	double Piii=pIIIt_convDivNoz(ma1)*pt;	
	double qcr;
	double pto=0;
	if(pamb<Piii&&pamb>Pii){
		//pii<p<piii,在渐扩段产生激波
		//候部达到音速,求临界流量
		 qcr=qcr_dA(pt,pt/Rg/Tt)*Acr;
		 pto=pt_from_qApTt(qcr,Ao,pamb,Tt,1e5);
		 mx=pt21_to_ma1(pto/pt,10);//激波前Ma1
		 Ax=ma_to_Ar(mx)*Acr;//激波出现在Ax截面上
		 lamx=ma_to_lamda(mx);//激波前lamda 

	}



	


	return 0;
}

