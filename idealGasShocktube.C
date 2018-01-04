// gasNozzle.cpp : �������̨Ӧ�ó������ڵ㡣
// copyright by daiyuqiang@dlut.edu.cn
//ref[1] ����������������������ը���嶯��ѧ,1987,6
//ref[2] ������  ������ ��ΰ�⣬���嶯��ѧ������������ѧ�����磬2006,7
//����˵��
/*
 I,II    ����˶������Ĳ�ǰ�Ͳ���
 1 2     �����˶��Ĳ�ǰ�Ͳ���
 U       ������Բ�ǰ���������ٶȣ�����Vs����
 D       ���������ٶ�
 c       ����
 v       �ٶ�
 
 vI =v1-D =        U;
 vII=v2-D =(v2-v1)+U;
 Ms      ������������Ms= U/c1 = D - v1(rightRuning)> 1.0
            or     Ms=-U/c1 =v1 -  D( leftRuning)<-1.0 
 */
#include <iostream>
#include <math.h>
#include "idealGasShocktube.h"

//================================================
int main()
{
	std::cout<< ("Copyright by daiyuqiang@dlut.edu.cn\n")<< std::endl;	
//��֤v2-v1/c1 ��Ms,ref[1] p320 ex(8-1)
    //����ֹ״̬�������ٵ�״̬1
    double p0=7.83e5;
    double T0=288;
    double p1=1e5;
    double ma1=pst_to_ma(p1/p0);
    double v1,v2=0.;
    std::cout<< "Ma1=   "<< ma1  <<std::endl;
    double c10=ma_to_Sndst(ma1);
    std::cout<< "c10=c1/c0=   "  <<c10  <<std::endl;
    double T10=ma_to_Tst(ma1);
    std::cout<< "T10=T1/T0=   "  <<T10  <<std::endl;
    double c1=sqrt(Rg*Cpv*T0*ma_to_Tst(ma1));
    std::cout<< "c1=   "  <<c1  <<std::endl;
    
    double v21c1=(v2-ma1*c1)/c1;
    double ms=v21c1_to_Ms(v21c1,Cpv);
    std::cout<< "shockcmach Ms=   "  <<ms  <<std::endl;
    if (ms<-1) std::cout << "leftRuning shockwave found" << std::endl;
    double shockVa=ms*c1+ma1*c1;
    std::cout<<  "sbsolute veocity of the moving shock is : "<< shockVa << std:: endl;
    double p2=p1*p21_from_Ms(ms,Cpv);
    std::cout<<  "pressure behind the moving shock is : "<< p2 << std:: endl;
    
 // ref[2] p278 ex16
    v1=150;
    double T1=300;
    p1=1.5e5;
    v2=0;
    c1=sqrt(Rg*Cpv*T1);
    ma1=v1/c1;
    ms=v21c1_to_Ms((v2-v1)/c1,Cpv);
    std::cout<< "shockcmach Ms=   "  <<ms  <<std::endl;
    std::cout<< "shockcmach absolue velocity=   "  <<ms*c1+v1  <<std::endl;
    if (ms<-1) std::cout << "leftRuning shockwave found" << std::endl;
    std::cout<<   "pressure behind the moving shock is : p2=   [v21c_to_p21ͨ�ó��������]     "<<   p1* v21c1_to_p21((v2-v1)/c1,Cpv) << std::endl;
    std::cout<<   "pressure behind the moving shock is : p2=   [v21c_to_Ms,Ms_to_p21]   "          <<  p1* p21_from_Ms( ms, Cpv)   << std::endl;
//�˶�����ǰ���ϵ
	double pt=1.52e6;
	double Tt=500;
	double ma=3.9376;//����ǰmach
	double denst=pt/Rg/Tt;
	//�˶�����ǰ������
	            p1=pt* ma_to_pst(ma);
	            T1=Tt* ma_to_Tst(ma);
	double dens1=denst* ma_to_Denst(ma);
	double u1=ma*sqrt(Cpv*Rg*T1);
//R-H��ϵ��֤
//�����˶������ļ�������ǰp1��T1�Ͳ���p2��֪
	p1=1.013e5;
	T1=288;
	p2=1.1143e5;
    double p21=p2/p1;
	std::cout<< "p21= "<< p21 << std::endl;
    dens1=p1/Rg/T1;
	std::cout<< "des1= "<< dens1 << std::endl;
    double dens2=dens1*D21_from_p21(p21,Cpv);
    std::cout<< "dens2= "<< dens2 << std::endl;
	double T2=T1*T21_from_p21(p21,Cpv);
    std::cout<< "T2=    "<< T2    << std::endl;
	       c1=sqrt(Cpv*Rg*T1);
	double c2=sqrt(Cpv*Rg*T2);
	std::cout<< "c1=    "<< c1 << std::endl;
    std::cout<< "c2=    "<< c2 << std::endl;
 
//ref[2] p226-227 ex15
//(1)
    std::cout<< "case 1 p21=2.5 �ھ�ֹ�����д�����" << std::endl;
    p21=2.5;
    p21=1.656;
    ma1=0;
    v21c1=v21c1_from_p21(p21,Cpv,1);
    //std::cout<< "!!!v2c1=     "<< v21c1+ma1<< std::endl;
  
    double ma2=0;
    //ma2=v2/c2=(v21c1*c1+v1)/c2= (v21c1*c1+v1)/c1*(c1/c2= (v21c1+ma1)/c21=(v21c1+ma1)/sqrt(T21)
    ma2=(v21c1+ma1)/sqrt(T21_from_p21(p21,Cpv));
    std::cout<< "v2/c2=Ma2    "<< ma2 << std::endl;
    std::cout<< "Ma2^2    "<< ma2*ma2 << std::endl;
	
//(2)
    std::cout<< "case 2 p21=5 �ھ�ֹ�����д�����" << std::endl;
    ma1=0;
    p21=5;
    v21c1=v21c1_from_p21(p21,Cpv,1);
    ma2=(v21c1+ma1)/sqrt(T21_from_p21(p21,Cpv));
    std::cout<< "v2/c2=Ma2    "<< ma2 << std::endl;
    std::cout<< "Ma2^2    "<< ma2*ma2 << std::endl;
//(3-1)
    std::cout<< "case 3-1 p21=2.5 ��ma1=-2.0�������������򴫲�(rightRuning)��" << std::endl;
    ma1=-2.0;
    p21=2.5;
    v21c1=v21c1_from_p21(p21,Cpv,1);//1-rightRuning
    ma2=(v21c1+ma1)/sqrt(T21_from_p21(p21,Cpv));
    std::cout<< "v2/c2=Ma2    "<< ma2 << std::endl;
    std::cout<< "Ma2^2    "<< ma2*ma2 << std::endl;
    std::cout<< "Ms=      "<< Ms_from_p21(p21,Cpv,1)<< std::endl; 
	
//(3-2)
    std::cout<< "case 3-2 p21=2.5 ��ma1=2.0�������������򴫲���leftRuning��" << std::endl;
    ma1=2.0;
    p21=2.5;
    v21c1=v21c1_from_p21(p21,Cpv,0);//0-leftRuning
    ma2=(v21c1+ma1)/sqrt(T21_from_p21(p21,Cpv));
    std::cout<< "v2/c2=Ma2    "<< ma2 << std::endl;
    std::cout<< "Ma2^2    "<< ma2*ma2 << std::endl;
    std::cout<< "Ms=      "<< Ms_from_p21(p21,Cpv,0)<< std::endl; 
    
//(4-1)
    std::cout<< "case 4-1 p21= 5 ��ma1=-2.0�������������򴫲�(rightRuning)��" << std::endl;
    ma1=-2.0;
    p21=5;
    v21c1=v21c1_from_p21(p21,Cpv,1);//1-rightRuning
    ma2=(v21c1+ma1)/sqrt(T21_from_p21(p21,Cpv));
    std::cout<< "v2/c2=Ma2    "<< ma2 << std::endl;
    std::cout<< "Ma2^2    "<< ma2*ma2 << std::endl;
    std::cout<< "Ms=      "<< Ms_from_p21(p21,Cpv,1)<< std::endl; 
	
//(4-2)
    std::cout<< "case 4-2 p21= 5 ��ma1=2.0�������������򴫲���leftRuning��" << std::endl;
    ma1=2.0;
    p21=5;
    v21c1=v21c1_from_p21(p21,Cpv,0);//0-leftRuning
    ma2=(v21c1+ma1)/sqrt(T21_from_p21(p21,Cpv));
    std::cout<< "v2/c2=Ma2    "<< ma2 << std::endl;
    std::cout<< "Ma2^2    "<< ma2*ma2 << std::endl;
    std::cout<< "Ms=      "<< Ms_from_p21(p21,Cpv,0)<< std::endl; 

// �����ھ�ֹ����ķ�����֤  ref[1] p342 ex 8-4
    std::cout<< " �˶������ھ�ֹ����ķ�����֤ "<<  std::endl;

    double absVi=655;
                    v1=0;
    double snd1=340;
                    p1=1e5;
              ms=(absVi-v1)/snd1;
    std::cout<< " incidence shock mach: Ms =  "<< ms << std::endl;
    double Mr =Mr_from_Ms_shockRef_sWall(ms, Cpv,1);
    std::cout<< " reflect shock mach: Mr =  "<< Mr<< std::endl;
    p21=p21_from_Ms(ms,Cpv);
    c2=c1*sqrt(T21_from_p21(p21,Cpv));
    v2 =v21c1_from_p21(p21,Cpv,1)*c1+v1;//���䲨������,���伤��Ϊ����
    std::cout<< "  zone2 gas Velocity: v2=  "<<  v2  <<std::endl;

    double c3 =c1*c31_from_p21_shockRef_sWall(p21,Cpv);
   
    std::cout<< "  zone2,33  Snd Velocity: c2=  "<<  c2 << "   c3=   " <<  c3 <<std::endl;
    double Dr=Dri_from_Ms_shockRef_sWall(ms, Cpv);
    std::cout<< " reflect shock Velocity: absVr =  "<< Mr*c2 + v2<< std::endl;  //Dr =Ur+v2 
    //���伤������ٶ�Ur Ϊ���伤����� ���伤���Ĳ�ǰ�ٶȶ��ԣ�Ҳ�������伤���Ĳ������٣�����v2,����  Ur = Dr -v2
    //����ѹǿp3    
    std::cout<< " reflect Press:  P3 =  "<< p31_from_p21_shockRef_sWall(p21, Cpv) *p1<< std::endl;  //Dr =Ur+v2 
    
    //���伤���ڹ̱��Ϸ��䣬��v2-v1���루v3-v2��Ӧ�ô�С��ȷ����෴�����䲨�ͷ��䲨һ������һ�����У���Ӧ��v1=v3
    //��֤v3=0, ���Լ�����Ҫ����v21c1����ȷ��
    double p32=p32_from_p21_shockRef_sWall( p21, Cpv);
    double v32c2=v21c1_from_p21(p32,Cpv,0); //���з��伤��
    
    v21c1=v21c1_from_p21(p21,Cpv,1);//�������伤��
    std::cout<< " v2-v1=   "<< v21c1 *c1<< "     v3-v2=   " << v32c2*c2  << "    v3=   "<< v32c2*c2+v2<< std::endl;  
    
   //�������˶�����ķ�����֤
    std::cout<< " �˶��������˶�����ķ�����֤ "<<  std::endl;
    absVi= -655;
    double vw1=-100;  //  
    double vw2=-100;//389.125 make the v2=vw2
    double vw3=-100;
    v1=vw1;
    double v3=vw3;
    
    std::cout<< " Moving Wall velocity: Vw1 =  "<< vw1<< "     Vw2=  "<<vw2 << "    Vw3= "<< vw3 << std::endl;
    snd1=340;
    p1=1.013e5;
//
    ms=(absVi-v1)/snd1;  // v1!=0
    std::cout<< " incidence shock mach: Ms =  "<< ms << std::endl;
    p21=p21_from_Ms(ms,Cpv);
    c2=c1*sqrt(T21_from_p21(p21,Cpv));
    v2 =v21c1_from_p21(p21,Cpv,1)*c1+v1;//���䲨������,���伤��Ϊ����
    std::cout<< "  zone2 gas Velocity: v2=  "<<  v2  <<std::endl;
    //
    Mr  =   Mr_from_p21_shockRefAsShock_mWall( p21,0, Cpv, vw1/snd1,vw3/snd1);
    std::cout<< " reflecive shock mach: Mr =  "<< Mr<< std::endl;   
    p32= p32_from_p21_shockRefAsShock_mWall( p21,0, Cpv, vw1/snd1,vw3/snd1);
    std::cout<< " pressure 3 / pressure2: p32 =  "<< p32<< std::endl;
    double  p31=p32*p21;
    std::cout<< " pressure 3 / pressure1: p31 =  "<< p31<< std::endl;
    std::cout<< " pressure 2 / pressure1: p21 =  "<< p21<< std::endl;

//�˶������ڿ��ڶ˵ķ�����֤
     std::cout<< " �˶������ڿ��ڶ˵ķ�����֤,case 1 "<<  std::endl;
    ms=1.25;
    p21=p21_from_Ms(ms,Cpv);
    std::cout<<" p21=    "<< p21<< std::endl;
    v1=0;
    c1=340;
    p1=1e5;
    double pa=p1;
    bool rRun=1;
    
    double dxdtc1_lead;
    double dxdtc1_trail;
    if(case_from_p21_shockRef_opening(p21,pa/p1,Cpv, v1/c1,rRun )==1){
        double ma3=ma3_from_p21_shockRefOpen_case1(p21, pa/p1, Cpv,  v1/c1, rRun);
        std::cout<<" Ma3=    "<< ma3<< std::endl;
        dxdtc1_from_p21_shockRefOpen_case12(p21, pa/p1,  Cpv,v1/c1,rRun, dxdtc1_lead, dxdtc1_trail);
         std::cout<<" dxdtc1_lead=    "<<dxdtc1_lead <<" dxdtc1_trail=    " << dxdtc1_trail << std::endl;
    }
//case 2
std::cout<< " �˶������ڿ��ڶ˵ķ�����֤,case 2 "<<  std::endl;
    ms=1.5;
    p21=p21_from_Ms(ms,Cpv);
    std::cout<<" p21=    "<< p21<< std::endl;
    v1=0;
    c1=340;
    p1=1e5;
    pa=p1;
    rRun=1;
    double c21=sqrt(T21_from_p21(p21,Cpv));
    double pStarp1,cStarc1;
    if(case_from_p21_shockRef_opening(p21,pa/p1,Cpv, v1/c1,rRun )==2){
        double ma3=ma3_from_p21_shockRefOpen_case1(p21, pa/p1, Cpv,  v1/c1, rRun);
        std::cout<<" Ma3=    "<< ma3<< std::endl;
        dxdtc1_from_p21_shockRefOpen_case12(p21, pa/p1,  Cpv,v1/c1,rRun, dxdtc1_lead, dxdtc1_trail);
         std::cout<<" dxdtc1_lead=    "<<dxdtc1_lead <<" dxdtc1_trail=    " << dxdtc1_trail << std::endl;
         cStar_pStar_from_p21_shockRefOpen_case2(p21, pa/p1, Cpv,v1/c1,rRun, cStarc1, pStarp1);
         std::cout<<" c*/c1=    "<<cStarc1<<" c*/c2=    "<< cStarc1/c21 <<" p*/p1=    " <<  pStarp1<< std::endl;
    }
//case 3
std::cout<< " �˶������ڿ��ڶ˵ķ�����֤,case 3 "<<  std::endl;
    ms=2.5;
    p21=p21_from_Ms(ms,Cpv);
    std::cout<<" p21=    "<< p21<< std::endl;
    v1=0;
    c1=340;
    p1=1e5;
    pa=p1;
    rRun=1;
   c21=sqrt(T21_from_p21(p21,Cpv));
   double v2c1=v21c1_from_p21(p21,Cpv,rRun);
    if(case_from_p21_shockRef_opening(p21,pa/p1,Cpv, v1/c1,rRun )==3)
                 std::cout<<" v2/c1=    "<<v2c1 << std::endl;
   //
    
 std::cout<< " �˶�����+�˶�������������֤,case L+R "<<  std::endl;   
  //ref[1] p347 ex8-7  
 p1=1.e5;
 T1=300.;
 c1=sqrt(Cpv*Rg*T1);
 v1=153.;
 double v2r=0.;
 double p2l1=5.0;
 double p2r1=v21c1_to_p21((v2r-v1)/c1, Cpv);
 std::cout<<" p2r1=    "<<p2r1 << std::endl;
 double msa ,  msb, msc, msd,    
               v2lc1, v2rc1, v3c1,
               m2l, m2r, m3l, m3r,
               /*p31,  */      T3l1,    T3r1, 
              c2l1,        c2r1,         c3l1,    c3r1,
              p32l,       p32r,
              Dac1,           Dbc1,            Dcc1,       Ddc1;

 lrRunShocksIntersection_from_pr(    p2l1, p2r1, Cpv, v1/c1,
                                                                    msa ,        msb,        msc,      msd,    
                                                                    v2lc1, v2rc1, v3c1, 
                                                                    m2l, m2r, m3l, m3r,
                                                                    p31,        T3l1,    T3r1, 
                                                                    c2l1,        c2r1,         c3l1,    c3r1,
                                                                    p32l,        p32r,
                                                                    Dac1,           Dbc1,            Dcc1,       Ddc1   );
    
    std::cout<< "Msa=   "<< msa <<" Msb=    "<< msb<<"   Msc="<< msc<<"   Msd=" << msd << std::endl;
    std::cout<< "V2l/c1=   "<< v2lc1 <<" V2r/c1=    "<< v2rc1<<"   V3/c1="<< v3c1 << std::endl;
    std::cout<< "Ma[2l]=   "<< m2l <<" Ma[2r]=    "<< m2r<<"   Ma[3l]="<< m3l <<"   Ma[3r]=" << m3r << std::endl;
    std::cout<< "p3/p1=   "<< p31 <<" T3l/T1=    "<< T3l1<<"   T3r1="<< T3r1  << std::endl;
    std::cout<< "c2l/c1=   "<< c2l1 <<"c2r/c1=    "<< c2r1<<"  c3l/c1="<< c3l1<<"  c3r/c1=" << c3r1 << std::endl;
    std::cout<< "p3l/p2l=   "<< p32l <<"             p3r/p2r=    "<< p32r << std::endl;
     std::cout<< "AbsV_Sa/c1=   "<< Dac1<<"   AbsV_Sb/c1=    "<< Dbc1<<"   AbsV_Sc/c1=   "<< Dcc1<<"   AbsV_Sd/c1=  " << Ddc1 << std::endl;
     //
    //ref[2] p301 ex---------------------------------------------------------------------------------------------------------------------
     //================
   std::cout<< "............ANOTHER EXAMPLE.................."<<std::endl;  
  p1=0.689e5;
 T1=294.;
 c1=sqrt(Cpv*Rg*T1);
 v1=0.;
  v2r=0.;
 p2l1=2.0;
 p2r1=4.0;
 lrRunShocksIntersection_from_pr(    p2l1, p2r1, Cpv, v1/c1,
                                                                    msa ,        msb,        msc,      msd,    
                                                                    v2lc1, v2rc1, v3c1, 
                                                                    m2l, m2r, m3l, m3r,
                                                                    p31,        T3l1,    T3r1, 
                                                                    c2l1,        c2r1,         c3l1,    c3r1,
                                                                    p32l,        p32r,
                                                                    Dac1,           Dbc1,            Dcc1,       Ddc1   );
    
    std::cout<< "Msa=   "<< msa <<" Msb=    "<< msb<<"   Msc="<< msc<<"   Msd=" << msd << std::endl;
    std::cout<< "V2l/c1=   "<< v2lc1 <<" V2r/c1=    "<< v2rc1<<"   V3/c1="<< v3c1 << std::endl;
    std::cout<< "Ma[2l]=   "<< m2l <<" Ma[2r]=    "<< m2r<<"   Ma[3l]="<< m3l <<"   Ma[3r]=" << m3r << std::endl;
    std::cout<< "p3/p1=   "<< p31 <<" T3l/T1=    "<< T3l1<<"   T3r1="<< T3r1  << std::endl;
    std::cout<< "c2l/c1=   "<< c2l1 <<"c2r/c1=    "<< c2r1<<"  c3l/c1="<< c3l1<<"  c3r/c1=" << c3r1 << std::endl;
    std::cout<< "p3l/p2l=   "<< p32l <<"             p3r/p2r=    "<< p32r << std::endl;
     std::cout<< "AbsV_Sa/c1=   "<< Dac1<<"   AbsV_Sb/c1=    "<< Dbc1<<"   AbsV_Sc/c1=   "<< Dcc1<<"   AbsV_Sd/c1=  " << Ddc1 << std::endl;   
     
    //
//ͬ�򼤲�׷�Ϸ���Ϊ���Ͳ�����֤
std::cout<< "............ͬ�򼤲�׷�Ϸ���Ϊ���Ͳ�����֤.................."<<std::endl;  
 p1=1.e5;
 T1=300.;
 c1=sqrt(Cpv*Rg*T1);
 v1=0.;
 p21=2.0;
 p32=2.0;
 rRun=1;
 double v4c1,      v5c1,    p41,      p51;  
 double Msa, Msb, Msc;
 SbChaseSa_reflectExp_from_Pr( p21,       p32,     Cpv,      v1/c1, rRun,
                                                                         Msa ,        Msb,        Msc,      
                                                                         v4c1,       v5c1,      p41,        p51  );
 std::cout<< "Msa=  "<< Msa <<"    Msb=  "<< Msb <<"   Msc=  "<< Msc <<  std::endl;
 std::cout<<" p41=p51=   " << p41  << "   p31=" << p21*p32 << std::endl;
 
 
 //
 //����+�Ӵ��淴��Ϊ���Ͳ�����֤
 p1=1.e5;
 p2=p1;
 T1=300.;
 c1=sqrt(Cpv1*Rg1*T1);
 c21=1.0;
 c2=c21*c1;
 T2=pow(c21*c1,2)/Rg2/Cpv2;
 v2=0;
 v1=v2;
 p32=2.0;
 rRun=0;
 double rho21=T1/T2;
   if(case_SaInterJ_from_P32( p32,   Cpv1, Cpv2,    c21) ==0)    {
   // std::cout<< "...........����+�Ӵ��淴��Ϊ���Ͳ�����֤.................."<<std::endl;  
    std::cout<< "The shock meets interface reflects as a RAREFACTION."<< std::endl;
    SaInterJ_reflectExp_from_P32(     p32,           v2/c2,        Cpv1, Cpv2,    c21,   rho21, rRun,
                                     Msa , Msb, v4c1,p41, p51  ) ; 
 
    std::cout<< "Msa=  "<< Msa <<"    Msb=  "<< Msb  <<  std::endl;
    std::cout<<" p41=p51=   " << p41  << "   p31=" <<p32 << std::endl;
    std::cout<<" v4/c1=   " << v4c1  << std::endl;}
 //
 //����+�Ӵ��淴��Ϊ��������֤
 
  if(case_SaInterJ_from_P32( p32,   Cpv1, Cpv2,    c21) ==1)   { 
     std::cout<< "The shock meets interface reflects as a SHOCK."<< std::endl;
     // std::cout<< "...........����+�Ӵ��淴��Ϊ��������֤.................."<<std::endl;  
    SaInterJ_reflectShock_from_P32(  p32,  v2/c2,  Cpv1, Cpv2,  c21,  rho21, rRun,
                                     Msa ,Mr,Msb, v4c1,p41, p51  ) ;  
    std::cout<< "Msa=  "<< Msa <<"    Msb=  "<< Msb  << "  Mr =   "<< Mr << std::endl;
    std::cout<<" p41=p51=   " << p41  << "   p31=" <<p32 << std::endl;
    std::cout<<" v4/c1=   " << v4c1  << std::endl; }

    
    
    //��������֤
    //ref[1] p361 ex8-8
    std::cout<< "...........��������֤1.................."<<std::endl;   
    double pH=1.e7;
    double Th=300;
    double pL=1.e5;
    double pHL=pH/pL;
    double Tl=300;
    double CpvH=1.4;
    double CpvL=1.4;
    double RgHL=1.;
    double Thl=Th/Tl; 
    rRun=1;
    double p61,p71,p81;
    double T21,T31,T51,T61,T71,T81;
    double Mra,Mrb;
    iGasShockTube(pHL, CpvH,CpvL, RgHL,  Thl,rRun,
            Msa, Mra,Mrb,
            p21, p51,p71,p61,p81,v2c1,
            T21, T31,T51,T71,T81,T61);
    std::cout<< "Msa=  "<< Msa <<"  Mra=  "<< Mra <<"        Mrb=  " <<  Mrb   <<  std::endl;
    std::cout<<" p31=p21=   " << p21  << std::endl;
    std::cout<<" v2/c1=   " << v2c1  << std::endl; 
    
    std::cout<< "T2=  "<< T21*Tl<<"       T3=  "<< T31*Tl  << std::endl;
    std::cout<< "p2=  "<< p21*pL<<"      Ma2=  "<<v2c1/ sqrt( T21_from_p21(p21,CpvL))  <<std::endl;
    std::cout<< "p5= "<< p51*pL<<"       T5= "<< T51*Tl<< std::endl;
    std::cout<<" p6= "<<p61*pL<<" p7= "<<p71*pL<<" p8= "<<p81*pL <<std::endl;
    std::cout<<" T6= "<<T61*Tl <<" T7= "<<T71*Tl <<" T8= "<<T81*Tl <<std::endl;
    //test 2 
    //ref[1] p361 ex8-8
  std::cout<< "...........��������֤2.................."<<std::endl;        
     CpvH=1.4;
     CpvL=1.4;
     RgHL=1.;
     pH=2.e5;
    double  cH=344;
    double  cL=343.5;
    pL=1.e5;
    pHL=pH/pL;
    double RgH=287.1;
    double RgL=287.1;
    Th=cH*cH/CpvH/RgH;
    Tl= cL*cL /CpvL /RgL;
    RgHL=RgH/RgL;
    Thl=Th/Tl; 
    rRun=1;
    p61,p71,p81, T21,T31,T51,T61,T71,T81;
    iGasShockTube(pHL, CpvH,CpvL, RgHL,  Thl,rRun,
            Msa, Mra,Mrb,
            p21, p51,p71,p61,p81,v2c1,
            T21, T31,T51,T71,T81,T61);
    std::cout<<"rRun=0-----Left Incident Shockwave" << std::endl;
    std::cout<< "Msa=  "<< Msa <<"  Mra=  "<< Mra <<"        Mrb=  " <<  Mrb   <<  std::endl;
    std::cout<<" p31=p21=   " << p21  << std::endl;
    std::cout<<" v2/c1=   " << v2c1  << std::endl; 
    
    std::cout<< "T2=  "<< T21*Tl<<"       T3=  "<< T31*Tl   << std::endl;
    std::cout<< "p2=  "<< p21*pL<<"      Ma2=  "<<  v2c1/sqrt(T21_from_p21(p21,CpvL))  <<std::endl;
    std::cout<< "p5= "<< p51*pL<<"       T5= "<< T51*Tl<< std::endl;
    std::cout<<" p6= "<<p61*pL<<" p7= "<<p71*pL<<" p8= "<<p81*pL <<std::endl;
    std::cout<<" T6= "<<T61*Tl <<" T7= "<<T71*Tl <<" T8= "<<T81*Tl <<std::endl;
    
    //     �����ܵ���֤3
    std::cout<< "  �����ܵ���֤3"<<std::endl;
    
    pH=20.e5;
    Th=300;
    pL=1.e5;
    pHL=pH/pL;
    Tl=300;
    CpvH=1.4;
    CpvL=1.4;
    RgHL=1.;
     Thl=Th/Tl; 
    rRun=1;
    p61,p71,p81;
    T21,T31,T51,T61,T71,T81;
    Mra,Mrb;
    iGasShockTube(pHL, CpvH,CpvL, RgHL,  Thl,rRun,
            Msa, Mra,Mrb,
            p21, p51,p71,p61,p81,v2c1,
            T21, T31,T51,T71,T81,T61);
    std::cout<< "Msa=  "<< Msa <<"  Mra=  "<< Mra <<"        Mrb=  " <<  Mrb   <<  std::endl;
    std::cout<<" p31=p21=   " << p21  << std::endl;
    std::cout<<" v2/c1=   " << v2c1  << std::endl; 
    
    std::cout<< "T2=  "<< T21*Tl<<"       T3=  "<< T31*Tl  << std::endl;
    std::cout<< "p2=  "<< p21*pL<<"      Ma2=  "<<v2c1/ sqrt( T21_from_p21(p21,CpvL))  <<std::endl;
    std::cout<< "p5= "<< p51*pL<<"       T5= "<< T51*Tl<< std::endl;
    std::cout<<" p6= "<<p61*pL<<" p7= "<<p71*pL<<" p8= "<<p81*pL <<std::endl;
    std::cout<<" T6= "<<T61*Tl <<" T7= "<<T71*Tl <<" T8= "<<T81*Tl <<std::endl;
    
    return 0;}