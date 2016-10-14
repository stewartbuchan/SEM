#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <random>
#include <cmath>
#include <string.h>
#include <new>
#include <fstream>
#include <list>
#include <iterator>
#include <ctime>

using namespace std;

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

/*cosmological parameters*/
double omega_m=0.25;
double omega_l=1.-omega_m;
double omega_k=0;
double h=0.73;
double spectral_index=1.0;
double omega_baryon=0.045;
double sigma_8=0.9;
double omega_0=1.;
double domega=0.05;
double dc0=1.69;

//  //random number generator

std::default_random_engine generator(time(0));
std::normal_distribution<double> my_normal_distribution(0,1.0); //normal with mean=0 standard dev=1
std::uniform_real_distribution<double> my_uniform_distribution(0,1);

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


struct MatZ {
double logM;
double logZ;
MatZ(double M, double Z){logM=M;logZ=Z;}
};


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_deltavir(double z){
double Esq=omega_m*pow((1.+z),3)+omega_l;
double omegaMz=omega_m*pow((1.+z),3)/Esq;
double delta=(18.*M_PI*M_PI+82.*(1-omegaMz)-39*(1-omegaMz)*(1-omegaMz))/omegaMz;
return delta;
}

double get_Hz(double z){
double Hz=sqrt(omega_m*pow((1+z),3)+omega_l);
return Hz;
}

double Mmin(double z){
//Neistein+2006 eqn 20
return log10(0.5*1.52*pow(10.0,9.0)*pow((get_deltavir(z)/101.0),-0.5)/get_Hz(z));
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


std::list<double> deriv(std::list<double> x, std::list<double> y){
std::list<double> dydx;
std::list<double>::iterator itx = x.begin();
std::list<double>::iterator ity = y.begin();
double si,x0,x1,x2,y0,y1,y2;
int counter=0;
while (itx!=x.end()&&ity!=y.end()){
//cout <<counter << "  " << *ity << "  " << *itx << endl;
if(counter>2){x0=x1;x1=x2;y0=y1;y1=y2;}
if(counter==0){x0=*itx;y0=*ity;itx++;ity++;}
else if(counter==1){x1=*itx;y1=*ity;si=(y1-y0)/(x1-x0);itx++;ity++;dydx.push_back(si);}
else{x2=*itx;y2=*ity; si=((y1-y0)/(x1-x0)+(y2-y1)/(x2-x1))/2.;itx++;ity++;dydx.push_back(si);}
counter++;
}
si=(y2-y1)/(x2-x1);
dydx.push_back(si);
return dydx;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double integrate(double deltax, std::list<double> y){
// //intergrates y(x) for evenly spaced x
// //I=sum(0.5*deltax(y_i+y_i+1)) = deltax*(sum(y)-0.5*(y1+y_n))
std::list<double>::iterator ity=y.begin();
double I;
while(ity!=y.end()){I+=*ity;ity++;}
I-=0.5*(y.front()+y.back());
I*=deltax;

return I;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


std::list<double> interpol(std::list<double> x1, std::list<double> x2, std::list<double> y1)
{
// //linear interpolation; x1 y1 is library and x2 is input.
std::list<double> y2;
std::list<double>::iterator itx1=x1.begin(),itx2=x2.begin(),ity1=y1.begin();
double prex,prey;
while(itx2!=x2.end()){
//cout << "WSDFSF  " << *itx1 << "  " << *ity1 << "  " << *itx2 << endl;
while(*itx1 > *itx2&&itx2!=x2.end()){/*cout << "SDGAGF  " << *itx1 << "  " << *ity1 << "  " << *itx2 << endl;*/y2.push_back(((*itx2-prex)/(*itx1-prex)*(*ity1-prey))+prey);itx2++;}
prex=*itx1;prey=*ity1;
itx1++;ity1++;
}
return y2;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



double M_to_sig(double M)
{
double gammam=omega_m*h*exp(-omega_baryon * (1.+sqrt(2.*h)/omega_m));
double u=3.804*0.0001*gammam*pow((1/omega_m),(1./3.))*pow(10,M/3);

double g=64.087*pow((1.+(1.074*(pow(u,(0.3)))) - (1.581*(pow(u,(0.4)))) + (0.954*(pow(u,(0.5)))) - (0.185*(pow(u,(0.6))))),(-10.));
u=32.0*gammam;
double gs=64.087* pow(( 1.+ (1.074*(pow(u,(0.3)))) - (1.581*(pow(u,(0.4)))) + (0.954*(pow(u,(0.5)))) - (0.185*(pow(u,(0.6))))),(-10.));
double f=(g*g)/(gs*gs);
double s=(f * sigma_8 * sigma_8);
// //returns sigma
return sqrt(s);
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double D_z_white(double z)
{
// //the linear growth factor D(log(1+Z))
double Ez=sqrt(omega_l+omega_k*pow(10,2*z)+omega_m*pow(10,3*z));
double omega_m_z=omega_m*pow(10,3*z)/(Ez*Ez);
double omega_l_z=omega_l/(Ez*Ez);
double gz=2.5*omega_m_z/(pow(omega_m_z,(4./7.))-omega_l_z+(1.+omega_m_z/2.)*(1.+omega_l_z/70.));
double gz0=2.5*omega_m/(pow(omega_m,(4./7.))-omega_l+(1.+omega_m/2.)*(1.+omega_l/70.));

double Dz=(gz/gz0)/pow(10,z);

return Dz;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_z(std::list<double> Zlib, std::list<double> Dlib, int n3, double zp, double domega)
{
// //get redshift by interpolating the linear growth factor using domega=dc(zp)-dc(z)
double Dp=D_z_white(zp);
std::list<double> D,Z;
D.push_back(1/((domega/dc0)+(1./Dp)));
Z=interpol(Dlib, D, Zlib);

return Z.back();
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



MatZ splitting_P08(std::list<double> siglib, std::list<double> Mlib, double Mmin, double logM, double logZ, list<MatZ> *prolist){
// //full algarthm derived in Parkinson 2008 (MNRAS 282,557-564) Appendix A
double G0=0.57,gamma1=0.38,gamma2=-0.01, epsilon1=0.1,epsilon2=0.1;
double tmp,qres=pow(10,(Mmin-logM)), alphah, Vres, Vh, beta, B, mu, eta, I, deltaZ, F, ures, c,J, r1,r2,r3, qtry, sigmatry,R, logZp;
std::list<double> q,logsiglib,alphalib,tmplist;
bool sucessfulsplit;

std::list<double>::iterator itM=Mlib.begin(),itsig=siglib.begin();
while(itM!=Mlib.end()){q.push_back(pow(10,(*itM-logM)));itM++;}
double sigma2=M_to_sig(logM), sigmah=M_to_sig(logM-log10(2)), sigmares=M_to_sig(Mmin);
while(itsig!=siglib.end()){logsiglib.push_back(log10(*itsig));itsig++;}
alphalib=deriv(Mlib,logsiglib);
Vres=sigmares*sigmares/pow(((sigmares*sigmares)-(sigma2*sigma2)),1.5);
Vh=sigmah*sigmah/(pow(((sigmah*sigmah)-(sigma2*sigma2)),1.5));
beta=log10(Vres/Vh)/log10(2*qres);
B=Vh*pow(2,beta);
tmplist.push_back(logM-log10(2));
mu=alphah=-(interpol(Mlib,tmplist,alphalib).back());
tmplist.clear();
eta=beta-1-gamma1*mu;

double ddelta_dz=dc0*((1/D_z_white(log10(pow(10,logZ)+0.000001)))-(1/D_z_white(log10(pow(10,logZ)-0.000001))))/0.000002;

// //eqn A5 P08 just integrating S(q)
I=sqrt(2./M_PI)*B*alphah*G0*pow(2.,(-mu*gamma1))*(pow(((dc0/D_z_white(logZ))/sigma2),gamma2))*(pow((sigmah/sigma2),gamma1))*ddelta_dz*(1/eta)*(pow(0.5,eta)-pow(qres,eta));
// //double check!!

cout << "all the constants!!  " << I << "  " << B << "  " << beta << "  " << alphah << "  " << sigma2 << "  " << sigmah << "  " << eta<< "  " << (pow(0.5,eta)-pow(qres,eta)) << "  " << qres << "  " << Vres << "  " << Vh << endl;

// //calculate redshift step
deltaZ=epsilon2/I;
tmp=epsilon1*sqrt(2.*(sigmah*sigmah-sigma2*sigma2))/ddelta_dz;
if(deltaZ>tmp){cout << "swapped\n"; deltaZ=tmp;}

// //calcualte J (eqn A7) 
ures=sigma2/(sqrt(sigmares*sigmares-sigma2*sigma2));

// //integrate to get j(ures) 
std::list<double> u,ju;
tmp=ures/1000.; // //check?
c=tmp;
while(c<=ures){u.push_back(c);ju.push_back(pow((1.+(1./(c*c))),(gamma1/2.)));c+=tmp;}
J=integrate(tmp,ju);
//cout << "J(u)  " << J  << "  " << ures << "  " <<  gamma1 << endl;

// //calculate the fraction of mass below resolution limit (eqn A6)
F=sqrt(2./M_PI)*J*(G0/sigma2)*pow(((dc0/D_z_white(logZ))/sigma2),gamma2)*ddelta_dz*deltaZ;

sucessfulsplit=false;
r1=my_uniform_distribution(generator);
if(r1<I*deltaZ){
// //try for splitting
r2=my_uniform_distribution(generator);
qtry=pow((pow(qres,eta)+(pow(2,-eta)-pow(qres,eta))*r2),(1/eta)); // //choose random mass to split from power law
sigmatry=M_to_sig(log10(qtry)+logM);
tmplist.push_back(log10(qtry)+logM);
double alphaq=-interpol(Mlib,tmplist,alphalib).back();

// //eqn A3 to compinsate for 
R=(alphaq/alphah)*(sigmatry*sigmatry/pow((sigmatry*sigmatry-sigma2*sigma2),1.5))/(B*pow(qtry,beta))*pow((pow(2*qtry,mu)*sigmatry/sigmah),gamma1);
tmplist.clear();
r3=my_uniform_distribution(generator);

if(r3<R){sucessfulsplit=true;}
}
//cout << "F  " <<F <<endl;
logZp=log10(pow(10,logZ)+deltaZ); // //calculate new redshift
MatZ Macc(log10(F)+logM,logZp);  // //Msmooth=F*M0
if(sucessfulsplit){prolist->push_back(MatZ(logM+log10(1-F-qtry),logZp));prolist->push_back(MatZ(logM+log10(qtry),logZp));} // //if the split was sucessful; 2 progenitors, one with mass qtry*M0 and the other (1-qtry-F)*M0
else{prolist->push_back(MatZ(logM+log10(1-F),logZp));} //else; just subtract Msmooth
//cout << "deltaZ  " << deltaZ<<endl;
return Macc;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


int main(int argc,char *argv[]){
if(argc < 2) {
	printf("You must provide at least two arguments: the initial mass and the file to save the mergre tree to\n");
	exit(0);
}

double Initial_mass=strtod(argv[1],NULL);  //starting mass
//string filename=strtod(argv[1],NULL);
//double minM=Initial_mass-2;  //smallest mass/mass resolution
double minM=9.0;  //smallest mass/mass resolution


// report settings
cout << "cosmology   " << Initial_mass << "  " << minM << "  "<< omega_m << "  " << h << "  " << omega_baryon << "  " << sigma_8 << "  " << dc0 << "  " <<endl;

int no_its = 1; // number of iterations
int i,c,nrun;
double Z, tmpZ;

std::list<double> Mlib,siglib;
int p=100, n2=p*(Initial_mass+1)-p*(Initial_mass-6);
for(i=0;i<n2;i++) {Mlib.push_back(Initial_mass-6+(double)i/p);siglib.push_back(M_to_sig(Mlib.back()));}

int n3=3*p;
p=100;
std::list<double> Zlib,Dlib,deltalib,alphalib;
for(i=0;i<n3;i++){Zlib.push_front((double)i/p); Dlib.push_front(D_z_white(Zlib.front()));deltalib.push_front(dc0/Dlib.front());}


//main code
std::list<MatZ> prolist;
for(nrun=0;nrun<no_its;nrun++){
Z=log10(1.0);
std::list<MatZ> all_pros_list;
std::list<MatZ> all_pros_tmp;
all_pros_list.push_back(MatZ(Initial_mass,Z));
cout << "START" << endl;
do{ //while main pro > 2* m min
cout<< "mpro " << all_pros_list.begin()->logZ << "  " << all_pros_list.begin()->logM << endl;
tmpZ=Z;
c=0;

for(std::list<MatZ>::iterator pro = all_pros_list.begin(); pro != all_pros_list.end(); pro++){ //for each subhalo
std::list<MatZ> prolist;
MatZ Macc(0,0);
//cout << "TESTQ£  " << pro->logM << 
if(pro->logM>Mmin(pow(10,pro->logZ)-1.0)+(log10(2))){Macc=splitting_P08(siglib, Mlib, Mmin(pow(10,pro->logZ)-1.0), pro->logM, pro->logZ, &prolist);}
else{Macc=*pro;}
c+=1;

std::list<MatZ>::iterator list_iter=prolist.begin();
while(list_iter!=prolist.end()){cout << c << "  " <<  list_iter->logZ << "  " << list_iter->logM << endl;list_iter++;}
cout << c << " acc " << Macc.logZ << "  " << Macc.logM << "  " << 1 << endl;
all_pros_tmp.splice(all_pros_tmp.end(),prolist);
}
//add subhalos to a temporary list
all_pros_list.clear();
all_pros_list.splice(all_pros_list.end(),all_pros_tmp);
cout << "MINM  " << Mmin(pow(10,all_pros_list.begin()->logZ)-1.0) << "  " << (pow(10,all_pros_list.begin()->logZ)-1.0) << "  " << all_pros_list.begin()->logM << endl;
}while(all_pros_list.begin()->logM>Mmin(pow(10,all_pros_list.begin()->logZ)-1.0)+log10(2));
}


return 0;
}
