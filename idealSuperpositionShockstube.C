// 激波对撞管程序v1.0
// copyright by daiyuqiang@dlut.edu.cn
//符号说明
/*
 +______________________________________________________________________+
 |      11L               *                              /|                     |\                                         *        11R    |
 |                     Se                                 /||         ||           ||  \                                   \      Se             |
 |                *     /                               /  ||          ||          ||       \                                 \         *         |
 |        *           /                              /     ||          ||          ||            \                               \  10R   *   |
   *   10L     /                             /        ||            ||          ||              \                              \       Sd    |       
 |    Sd       /                              /          ||            ||          ||                  \                            \ *          |
 |        *   /                            /               ||           ||          ||                        \                    *  \ 9R     |
 |   9L   / *                         /                 JbL       ||       JbR                            \          *          \        |
 |       /          *              /                        ||          ||          ||                                   \*                  \     |
 |    /                    *   /                            ||          ||          ||                               *       \                  \  |                            
 | /                        /  *           7L            ||    8L  ||    8R ||         7R            *               \               \ |
 | \                     /         Sc                      ||  \       ||        / ||                     *                        \            /|
 |   \                /                *                    || SRb  || SRb ||                *                                \        / |
 |     \           /                        *               ||      \   ||   /      ||             Sc                                     \   /  |
 |        \    /                                  *         ||  6L / || \ 6R ||         *                                             /\    |
 |         \/                                        *      ||   SRa  SRa   ||     *                                               /   \   |
 |         /\                                           *   ||  /       ||       \ || *                                                 /        \|
 |   /        \                                              //* 5L  || 5R *\\                                               /           \|
 |  \           \                                        //        *    ||   Sb    \\                                         /              /|
 |     \            \                                //                 *                \\                                    Rl            /  |
 |        \           \           3L           //      2L     *      *     2R     \\             3R            /             /      |
 |            \         \                      //             *                  *           \\                       /           Rt          |
 |               \         \                 Ja           *                          *          \\                 /          /                 |
 |                   \       \           //        Sa                                  Sa      Ja            /        /                      |
 |      4L            \     \      //       *                 (1)                        *      \\     /     /                   4R    |
 |                           \  \ //    *                                                         *    \\/   /                                 |
 +_____________\\+________________________________* \\/__________ ______  + 
 zone 4L,4R     -  high p with CpvH, RgH
 zone 1             -  low p  with CpvL, RgL
 */
 //=========================================================================
#include <iostream>
#include <math.h>
#include "idealGasShocktube.h"
int main()
{
	std::cout<< "Copyright by daiyuqiang@dlut.edu.cn\n"<< std::endl;	
    std::cout<< " Supersition of two moving shocks for ideal gases."<< std::endl;
    
    double p4L  =     1.e6;
    double p4R  =     1.e6;
    double T4L  =     300.;
    double T4R  =     300.;
    double p1   =      2e5;
    double T1   =      300;
    double v1   =      0;
    
    double Cpv4L    =   1.4;
    double Cpv4R    =   1.4;
    double Cpv1     =   1.4;
    double Rg4L     =   287.1;
    double Rg4R     =   287.1;
    double Rg1       =   287.1;
    double c1          =    sqrt(Cpv1*Rg1*T1);
    bool rRun          =    1; //incident shock running?
    
    double p2L1,p3L1,p4L1,p5L1,p6L1,p7L1,p8L1,p9L1,p10L1,p11L1;
    double T2L1,T3L1,T4L1,T5L1,T6L1,T7L1,T8L1,T9L1,T10L1,T11L1;
    double p2R1,p3R1,p4R1,p5R1,p6R1,p7R1,p8R1,p9R1,p10R1,p11R1;
    double T2R1,T3R1,T4R1,T5R1,T6R1,T7R1,T8R1,T9R1,T10R1,T11R1;
    double MaL,MbL,McL, MdL,MeL,MaR,MbR,McR, MdR,MeR;
    
    double M2L,M3L,M4L,M5L,M6L,M7L,M8L,M9L,M10L,M11L;
    double M2R,M3R,M4R,M5R,M6R,M7R,M8R,M9R,M10R,M11R;
    
    double n_Msb,n_Msc,n_p51,n_p71,n_p81,n_T71,n_T51,n_T81;
    double v2Lc1,v2Rc1,v5c1;
    
    //调用激波管函数，求1L,2L,3L,4L,9L参数
    /*
    void iGasShockTube(double pHL,double CpvH, double CpvL,double RgHL, double Thl,
                   double& Msa, double& Msb,double& Msc,
                   double& p21, double& p51,double& p71,double& p61,double& p81,double& v2c1,
                   double& T21, double& T31,double& T51,double& T71,double& T81,double& T61){
     */ 
    //左侧
   iGasShockTube(p4L/p1, Cpv4L,Cpv1, Rg4L/Rg1,  T4L/T1,rRun,//incident shock rRun
            MaL, n_Msb,n_Msc,
            p2L1, n_p51,n_p71,p9L1,n_p81,v2Lc1,
            T2L1, T3L1,n_T51,n_T71,n_T81,T9L1);
    std::cout<< "   MaL=  "<< MaL <<  std::endl;
    std::cout<<"    p3L1=p2L1=   " << p2L1  << std::endl;
    std::cout<<"    v2L/c1=   " << v2Lc1  << std::endl;     
    std::cout<< "   T2L=  "<< T2L1*T1<<"       T3L=  "<< T3L1*T1   << std::endl;
    std::cout<< "   p2L=  "<< p2L1*p1<<"      Ma2L=  "<<v2Lc1/ sqrt( T21_from_p21(p2L1,Cpv1))  <<std::endl;
    std::cout<<"    p9L= "<<p9L1*p1<<std::endl;
    std::cout<<"    T9L= "<<T9L1*T1<<std::endl;
    //右侧
    iGasShockTube(p4R/p1, Cpv4L,Cpv1, Rg4L/Rg1,  T4R/T1,!rRun,//incident shock lRun
            MaR, n_Msb,n_Msc,
            p2R1, n_p51,n_p71,p9R1,n_p81,v2Rc1,
           T2R1, T3R1,n_T51,n_T71,n_T81,T9R1);
    std::cout<< "       MaR=  "<< MaR <<  std::endl;
    std::cout<<"        p3R1=p2R1=   " << p2R1  << std::endl;
    std::cout<<"        v2R/c1=   " << v2Rc1  << std::endl;     
    std::cout<< "       T2R=  "<< T2R1*T1<<"       T3R=  "<< T3R1*T1   << std::endl;
    std::cout<< "       p2R=  "<< p2R1*p1<<"      Ma2R=  "<<v2Rc1/ sqrt( T21_from_p21(p2R1,Cpv1))  <<std::endl;
    std::cout<<"        p9R= "<<p9R1*p1<<std::endl;
    std::cout<<"        T9R= "<<T9R1*T1<<std::endl;
    //
    // 反向碰撞Sa(rRun)+Sb(lRun)，zone 1---波前,zone2l --Sa波后；zone2r-Sb波后,zone3l  --Sc波后；zone3r-Sd 波后
    // Sa--rRunning incident shock; Sb--lRunning incident shock;Sc--Sb‘s refraction shock, lRunning
    // Sd--Sa's refraction shock,rRunning;  J-contact interface  
    double MaL_test, MaR_test; 
    double p52L,p52R;
    double Dac1,Dbc1,Dcc1,Ddc1;
   double  c2L1,    c2R1,  c5L1,   c5R1;
    // 
    lrRunShocksIntersection_from_pr( p2L1,  p2R1,  Cpv1,  v1/c1,
                                    MaL_test ,  MaR_test,  MbR,   MbL,    
                                    v2Lc1,      v2Rc1,      v5c1,   
                                    M2L,        M2R,       M5L,  M5R,
                                    p5L1,       T5L1,       T5R1, 
                                    c2L1,     c2R1,      c5L1,   c5R1,
                                    p52L,    p52R,       
                                    Dac1,  Dbc1, Dcc1, Ddc1)  ;                                    
  std::cout<< "     MaL_=    " <<MaL_test<< ",      MaR_="  << MaR_test <<",    MbR_= "<< MbR<< ",       MbL_=" <<MbL<< std::endl;
  std::cout<< "     v2Lc1=  " << v2Lc1<< ",  v2Rc1="  << v2Rc1 <<", v5c1= "<< v5c1 <<std::endl;
  std::cout<< "     M2L=    " <<M2L<< ",        M2R="  << M2R <<",  M5L= "<< M5L<< ",       M5R=" <<M5R << std::endl;           
  std::cout<< "     p5L=p5R=    " << p5L1*p1<< ",     T5L1=T5R1= "  << T5L1 <<", T5L=T5R= "<< T5R1*T1 <<std::endl;
  //std::cout<< "     c2L1=" << c2L1<< ",     c2R1="  << c2R1<<"  ,c5L1= "<< c5L1 <<",      c5R1= "<< c5R1 <<std::endl;
  std::cout<< "     p5L2L=p52L=  "<< p52L << " ,     p5R2R=p52R= "<< p52R <<std::endl;
  std::cout<< "     Dac1=   " <<Dac1<< ",      Dbc1="  <<  Dbc1<<", Dcc1=  "<<Dcc1<< ",    Ddc1=" <<Ddc1 << std::endl;          
  //
  //透射激波SbL+JaL  以及 SbR + JaR接触面相遇 
  /*利用S+J的结果注意各区的对应关系:
       // 介质 :              I           II          II          II          I               II          I                                           II              
       // 形参 :              1           2           3          4           5             Sa         Sb          Ja           Jb          Rr(Sr） 
       //---------      (Ja后) (Ja前 )(入波Sa后)(Rr后)(透波Sb后)
       // 实参R侧:      3R       2R          5R         6R        7R         SbR       ScR          JaR          JbR           RaR 
       // 实参L侧：    3L        2L          5L         6L        7R          SbL        Sc  L        JaL          JbL            RaL
       // 介       质:       H           L           L            L          H           L             H                           L                 H 
       //
        void SaInterJ_reflectExp_from_P32( double  p32,  double  v2c2,    double  Cpv1,double Cpv2,  double  c21,   double rho21, bool rRunning,
                                     double& Msa ,   double& Msb,       double& v4c1,    double& p41,     double& p51  ) 
      */
       //右侧rRun=1
      double c23R=sqrt(Cpv1/Cpv4R)*sqrt( Rg1/Rg4R )*sqrt(T2R1/T3R1);
      double rho23R=(T3R1/T2R1)*(Rg4R/Rg1);
      double v6Rc3R,p63R,p73R;     
       if(case_SaInterJ_from_P32( p52R,   Cpv4R, Cpv1,    c23R) ==0)    {
            std::cout<< "    The RIGHT side  shock meets interface reflects as a RAREFACTION."<< std::endl;    
            SaInterJ_reflectExp_from_P32(     p52R,           M2R,        Cpv4R, Cpv1,   c23R,   rho23R, rRun ,
                                            MbR , McR, v6Rc3R,p63R, p73R  ) ;  
            std::cout<< "     MSbR=" << MbR<< ",      MaScR= "  << McR  <<",   v6R/c3R= "<< v6Rc3R << std::endl;
            std::cout<< "     p8R/p3R = " << p63R << ",  p7R/3R="  << p73R  <<std::endl;  
              }
     else { 
            std::cout<< "     The RIGHT  side shock meets interface reflects as a SHOCK."<< std::endl;
            double MrR;//reflective shock _R
            SaInterJ_reflectShock_from_P32(  p52R,           M2R,        Cpv4R, Cpv1,    c23R,   rho23R, rRun ,
                                                    MbR ,MrR,  McR, v6Rc3R,p63R, p73R  ) ;  
            std::cout<< "     MSbR=" << MbR << " MrR=    "<< MrR  << ",      MScR= "  << McR  <<",   v6R/c3R= "<< v6Rc3R << std::endl;
            std::cout<< "     p6R/p3R = " << p63R << ",  p7R/3R="  << p73R  <<std::endl;      
            std::cout<< "     p6R = " << p63R *p2R1*p1<< ",  p7R ="  << p73R  *p2R1*p1<<std::endl;       }

    //left lRun
      double c23L=sqrt(Cpv1/Cpv4L)*sqrt( Rg1/Rg4L )*sqrt(T2L1/T3L1);
      double rho23L=(T3L1/T2L1)*(Rg4L/Rg1);
      double v6Lc3L,p63L,p73L;     
       if(case_SaInterJ_from_P32( p52L,   Cpv4L, Cpv1,    c23L) ==0)    {
            std::cout<< "    The LEFT side  shock meets interface reflects as a RAREFACTION."<< std::endl;    
            SaInterJ_reflectExp_from_P32(     p52L,           M2L,        Cpv4L, Cpv1,    c23L,   rho23L,!rRun ,
                                            MbL , McL, v6Lc3L,p63L, p73L  ) ;  
            std::cout<< "     MSbL=" << MbL<< ",      MaScL= "  << McL  <<",   v6L/c3L= "<< v6Lc3L << std::endl;
            std::cout<< "     p6L/p3L = " << p63L << ",  p7L/3L="  << p73L  <<std::endl;  
            std::cout<< "     p6L = " << p63L *p2L1*p1<< ",  p7 ="  << p73L  *p2L1*p1<<std::endl;       }
      else { 
            std::cout<< "        The LEFT  side shock meets interface reflects as a SHOCK."<< std::endl;
            double MrL;//reflective shock _L
            SaInterJ_reflectShock_from_P32(  p52L,           M2L,        Cpv4L, Cpv1,    c23L,   rho23L, !rRun ,
                                                    MbL ,MrL,  McL, v6Lc3L,p63L, p73L  ) ;  
            std::cout<< "     MSbL=" << MbL << " MrL=    "<< MrL  << ",      MScL= "  << McL  <<",   v6L/c3L= "<< v6Lc3L << std::endl;
            std::cout<< "     p6L/p3L = " << p63L << ",  p7L/p3L ="  << p73L  <<std::endl;     
            std::cout<< "     p6L = " << p63L *p2L1*p1<< ",  p7L ="  << p73L *p2L1*p1<<std::endl;    }

  //
   // 反向碰撞Sa(rRun)+Sb(lRun)，zone 1---波前,zone2l --Sa波后；zone2r-Sb波后,zone3l  --Sc波后；zone3r-Sd 波后
    // Sa(SRaL)--rRunning incident shock; Sb(SRaR)--lRunning incident shock;Sc(SRbL)--Sb‘s refraction shock, lRunning
    // Sd(SRbR)--Sa's refraction shock,rRunning;  J-contact interface  
    double MsraL_test, MsraR_test, MsrbL, MsrbR; 
    double p65L= p63L *p2L1*p1/(p5L1*p1);
    double p65R=p63R *p2R1*p1 /(p5L1*p1);//p5R1=p5L1;
    std::cout<< "p65L= " <<p65L << "p65R= " <<p65R << std::endl;
    double  DarLc5L,DarRc5R,DbrLc5L,DbrRc5R;                                   
    double v6Lc5L,      v6Rc5R,      v8Lc5L;    
    double  c2L1_n,    c2R1_n,  c5L1_n,   c5R1_n;
    // 求8区参数
    lrRunShocksIntersection_from_pr(p65L,  p65R,  Cpv1,  M5L,
                                    MsraL_test ,  MsraR_test,  MsrbL,   MsrbR,    
                                    v6Lc5L,      v6Rc5R,      v8Lc5L,   
                                    M6L,        M6R,       M8L,  M8R,
                                    p5L1,       T5L1,       T5R1, 
                                    c2L1_n,     c2R1_n,      c5L1_n,   c5R1_n,
                                    p65L,    p65R,       
                                    DarLc5L,DarRc5R,DbrLc5L,DbrRc5R)  ;                                    
  std::cout<< "     MsraL_=    " <<MsraL_test<< ",      MsraR_="  << MsraR_test <<",    MsrbR= "<< MsrbR<< ",       MsrbL=" <<MsrbL<< std::endl;
  std::cout<< "     v6L/c5L=  " << v6Lc5L<< ",  v6R/c5R="  << v6Rc5R <<", v8Lc5L= "<< v8Lc5L <<std::endl;
  std::cout<< "     M6L=    " <<M6L<< ",        M6R="  << M6R <<",  M8L= "<< M8L<< ",       M8R=" <<M8R << std::endl;           
  std::cout<< "     p5L=p5R=    " << p5L1*p1<< ",     T5L1=T5R1= "  << T5L1 <<", T5L=T5R= "<< T5R1*T1 <<std::endl;
  std::cout<< "     p6L5L=p65L=  "<< p65L << " ,     p6R5R=p65R= "<< p65R <<std::endl;
  std::cout<< "     DarLc5L=   " <<DarLc5L<< ",      DarRc5R="  <<  DarRc5R<<", DbrLc5L=  "<<DbrLc5L<< ",    DbrRc5R=" <<DbrRc5R << std::endl;          
  double p8L=p21_from_Ms(MsrbL,Cpv1)*( p63L *p2L1*p1);//p6L
  double p8R=p21_from_Ms(MsrbR,Cpv1)*( p63R *p2R1*p1);//p6R
  std::cout<< "     p8L=    " << p8L << ",     p8R= "  << p8R <<std::endl;

  
  
    return 0;}
