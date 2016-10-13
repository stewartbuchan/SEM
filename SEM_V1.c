#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stack>
#include <list>
#include <vector>
#include <iterator>
#include <random>
#include <cmath>


using namespace std;


/*cosmological parameters*/
double omega_m=0.25;
double omega_l=1.-omega_m;
double h=0.73;
double spectral_index=1.0;
double omega_baryon=0.045;
double sigma_8=0.9;
double omega_0=1.;
double domega=0.1;
double dc_0=1.69;
double G=4.32; //  km^2 * kpc^3 * Mpc^-2 * s^-2 * Msun^-1
double forb=0.;
bool includemerger=true;
bool includeshutdown=false;
bool includediskregrowth=true;
bool includeclumpy_accretion=false;
bool includebathtub_model_accretion=false;
bool includegasdissipation=false;
bool includeconcentration=false;
bool includeexpansion=true;
//bool includeDMsubstructure=false;
bool reinitialiseDM=false;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0,1.0);

class satillite_halo;
typedef std::list<satillite_halo> listgals;  //list of satilite halos that are merging onto the main branch
typedef std::list<listgals> listlistgals;  //list of those lists 



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


std::list<double> interpol(std::list<double> x1, std::list<double> x2, std::list<double> y1)
	{
		//linear interpolation; x1 y1 is library and x2 is input.

		cout << "asdfasf  " <<  x1.size() << "  " << x2.size() << "  " << y1.size() << endl; 

		std::list<double> y2;
		std::list<double>::iterator itx1=x1.begin(),itx2=x2.begin(),ity1=y1.begin();
		double prex,prey;
		while(itx2!=x2.end())
		{
			cout <<"blub_blub  " << *itx1 << "  " << *ity1 << "  " << *itx2 << endl;
			while(*itx1 > *itx2&&itx2!=x2.end())
				{
					cout <<  *itx1 << "  " << *ity1 << "  " << *itx2 << endl;y2.push_back(((*itx2-prex)/(*itx1-prex)*(*ity1-prey))+prey);itx2++;
				}
			prex=*itx1;prey=*ity1;
			itx1++;ity1++;
		}
		cout << "sdfsfg  " << y2.front()<<endl;
		return y2;
	}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double Mdm_to_Mstellar(double logM, double logZ)
	{
		// Z = 1+z
		double Z=pow(10,logZ);
		float m10=11.59, m11=1.195, n10=0.0351, n11=-0.0257, beta10=1.376, beta11=-0.826, gamma10=0.608, gamma11=0.329;
		double logM1=m10+m11*(Z-1)/Z;
		double n=n10+n11*(Z-1)/Z;
		double beta=beta10+beta11*(Z-1)/Z;
		double gamma=gamma10+gamma11*(Z-1)/Z;
		double logm=logM+log10(2*n)-log10(pow(10,gamma*(logM-logM1))+pow(10,-beta*(logM-logM1)));
		//cout <<logm<<endl;
		return logm;
	}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double get_Mgas(double logMstar, double logZ,double maxgas)
	{
		//Stewart et al 2009 apj 702:307
		double Z=pow(10,logZ);
		double alpha=0.59*pow(Z,0.45);
		double logMgas=logMstar+log10(0.04)-alpha*(logMstar-log10(4.5*pow(10,11)));
		//cout << "GAS!  " << logMgas << "  " << maxgas<<endl;
		if(logMgas>maxgas)
		{
			/*cout <<"ops"<<endl;*/ logMgas=maxgas;
		}
		return logMgas;
	}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



/*
double get_Rvir(double Mhalo,double Z)
	{
		//Guo et al 2010
		//Z=log(1+Z) Mhalo=log(Mhalo)
		double H_Z=(h*h*10000)*(omega_m*pow(10.,3.*Z)+omega_l);
		return pow(pow(10,Mhalo)*G/100./H_Z,(1./3.));
	}
*/


double get_Rvir(double Mhalo,double Z)
	{
		double omegaZ=omega_m*pow(10,3.*Z)/(omega_m*pow(10,3*Z)+omega_l);
		double dc=10.*M_PI*M_PI+82*(1.-omegaZ)-39.*(1.-omegaZ)*(1.-omegaZ);
		double alz=31.*pow(omega_m/omegaZ*dc/(18.*M_PI*M_PI),(-1./3.))*(7./pow(10,Z));
		return (0.7/h)*alz*pow(10,(Mhalo-12.)/3.)/1000;
	}


double get_Rdisk(double Mdisk,double Z)
	{
		return (-1)-0.4*Z+0.14*Mdisk+(0.39-0.14)*log10(1+pow(10,(Mdisk-10.59988)));
	}

double get_Rbulge(double M1,double R1,double M2,double R2, double fc,double gas) //after Major merger work out the new radius of the bulge by magic (lets have another look at the Literature)
	{
		//double c=0.5;
		//double f=forb;
		double f0=0.25;
		if(M2<11)
			{
				R2*=2;
			}
		double AA,BB;
		if(R1<-3)
			{
				AA=0.;
			}
		else
			{
				AA=pow(10.,(2.*M1-R1));
			}
		if(R2<-3)
			{
				BB=0.;
			}
		else
			{
				BB=pow(10.,(2.*M2-R2));
			}
		//cout << "testR  " << AA << "  " << BB << "  " << endl;
		double p=AA+BB+fc*pow(10.,(M1+M2))/(pow(10.,R1)+pow(10.,R2));
		double r=(pow(10.,M1)+pow(10.,M2))*(pow(10.,M1)+pow(10.,M2))/p;
		if(includegasdissipation)
			{
				double fgas=pow(10.,gas)/(pow(10.,M1)+pow(10.,M2));
				r=r/(1.+fgas/f0);
			}
		return log10(r);
	}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_T_dynamical_friction(double Mh,double Ms,double Z)
	{
		//using formula in Boylan & kolchin 2008
		double A=0.9,B=1.,C=0.6,D=0.1; //revised values from McCavana et al 2012
		double T_dynamical=(0.1*9.78/(h*sqrt(omega_m*pow(10,3.*Z)+omega_l)));  // UNITS!!!!
		cout << "Td  " << Mh << "  " << Ms << "  " << Z << "  " << T_dynamical<<endl;
		double M_ratio=pow(10,(Mh-Ms));
		double coeff=A*T_dynamical*pow(M_ratio,B)/log(1+(M_ratio));
		double eta=0.23*distribution(generator)+0.5;
		//cout << eta<<endl;
		if(eta>=1.)
			{
				eta=0.999;
			}
		if(eta<=0.)
			{
				eta=0.001;
			}
		double eps=sqrt(1-eta*eta);
		double rvir=get_Rvir(Mh,Z)/1000;
		double rper=rvir*pow(eta,2.17);
		double rc=rper/(1.-eps);
		//cout << "TDF" << "  " <<T_dynamical <<endl;
		double T_dynamical_friction=coeff*exp(C*eta)*pow(rc/rvir,D);
		cout << "TDF" << "  " << Mh << "  " << Ms << "  " << pow(10,Z)-1. << "  " <<  rc/rvir << "  " << T_dynamical_friction << "  " << T_dynamical << "  " << T_dynamical_friction/T_dynamical <<endl;

		return T_dynamical_friction;
	}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_time_from_redshift(double Z) //returns a time dependent on the cosmology as omega_m
	{
		//returns time in Gyrs from log(1+Z)
		if(omega_m==1)
			{
				return 9.777505969*(2./3./h)*(pow(10.,(-Z*3./2.)));
			}
		else
			{
				return 9.777505969*(2./3/h/sqrt(1.-omega_m))*asinh(sqrt((1.-omega_m)/omega_m)/pow(10,(Z*(3./2))));
			}
	}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_M_halo_shutdown(double logZ)
	{

		//see Cattaneo+06 for details

		double logZc=log10(3.+1);
		double logMshock=12.;
		double K;

		if(logZ<logZc)
			{
				K=1.3*(pow(10,logZ)-pow(10,logZc));
			}
		else
			{
				K=0;
			}

		return logMshock+K;
	}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


// sub halo class
class satillite_halo 
	{
		double mass; //mass of DM in halo
		double redshift; //redshift of impact
		double Tdy; //dynamical time - from impact to merger
		double stellar_bulge; //mass of stars in bulge
		double stellar_disk; //mass of stars in disk
		double gas; //mass of gas in disk - assuming no gas in bulge
		double Rhl_disk; //half light radius of disk
		double Rhl_bulge; //hald light redius of bulge
		bool had_major_merger;
		bool shutdown;
		bool hadextememerger;
		bool cgt4;
		float test;
		double mu; // mass ratio of merger


		public:
		satillite_halo(double M, double Z,double new_stellar_bulge, double new_stellar_disk, double new_gas)
			{
				mass=M;
				gas=gas;
				redshift=Z;
				stellar_bulge=new_stellar_bulge;
				stellar_disk=new_stellar_disk;
				gas=new_gas;
				shutdown=false;
				had_major_merger=false;
				hadextememerger=false;
				cgt4=false;
				test=0.1;
			}

		double return_mass(){return mass;}
		void put_DM(double new_DM){mass=new_DM;}
		void add_DM(double new_DM){mass=log10(pow(10,mass)+pow(10,new_DM));}
		void subtract_DM(double new_DM){mass=log10(pow(10,mass)-pow(10,new_DM));}
		double return_gas(){return gas;}
		void put_gas(double new_gas){gas=new_gas;}
		void add_gas(double add_mass){gas=log10(pow(10,gas)+pow(10,add_mass));}
		double return_stellar_bulge() {return stellar_bulge;}
		void add_stellar_bulge(double add_mass){stellar_bulge=log10(pow(10,stellar_bulge)+pow(10,add_mass));}
		double return_stellar_disk() {return stellar_disk;}
		void add_stellar_disk(double new_stellar){stellar_disk=log10(pow(10,stellar_disk)+pow(10,new_stellar));}
		double return_stellar_total(){return(log10(pow(10,stellar_disk)+pow(10,stellar_bulge)));}
		double return_baryons() {return log10(pow(10,stellar_bulge)+pow(10,stellar_disk)+pow(10,gas));}

		double return_Rhl_bulge() {return Rhl_bulge;}
		double return_Rhl_disk() {return Rhl_disk;}
		void put_Rhl_bulge(double new_Rhl_bulge){Rhl_bulge=new_Rhl_bulge;}
		void put_Rhl_disk(double new_Rhl_disk){Rhl_disk=new_Rhl_disk;}
		double return_Rhl(){return (Rhl_bulge*pow(10,stellar_bulge)+Rhl_disk*(pow(10,stellar_disk)+pow(10,gas)))/(pow(10,stellar_bulge)+pow(10,stellar_disk)+pow(10,gas));}


		double return_redshift() {return redshift;}
		void put_Z(double Z){redshift=Z;}
		//double return_Tdy() {return 10000000000000000.;}
		//double return_Tdy() {return 0.;}
		double return_Tdy() {return Tdy;}
		double put_Tdy(double nTdy){Tdy=nTdy;}
		bool return_had_major_merger(){return had_major_merger;}
		bool return_shutdown(){return shutdown;}
		bool return_extrememerger(){return hadextememerger;}
		bool put_had_major_merger(){had_major_merger=true;}
		void put_shutdown(){shutdown=true;}
		void put_extrememerger(){hadextememerger=true;}
		void put_mu(double new_mu){mu=new_mu;}
		double return_mu(){return mu;}
		bool return_cgt4(){return cgt4;}
		void put_cgt4(){cgt4=true;}
		void output_vars(){cout<< mass << "  " << redshift << "  " << gas << "  " << stellar_bulge << "  " << stellar_disk << "  " << return_stellar_total() << "  " << Rhl_disk << "  " << Rhl_bulge << "  " << return_Rhl() << "  " << had_major_merger << "  " << shutdown << "  " << hadextememerger << "  " << cgt4 << endl;}
	};

// accreted mass class class
class accreted_mass 
	{
		double mass;
		double redshift;
		bool addmass;
		public:
		accreted_mass(double M, double Z,bool read_add)
			{
				mass=M;
				redshift=Z;
				addmass=read_add;
			}
		double return_mass() {return mass;}
		double return_redshift() {return redshift;}
		double return_addmass(){return addmass;}
	};

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


satillite_halo clumpyaccretion(satillite_halo cent,double Z1, double Z2,double Z,double gas_frac_crit) //mechanisim to grow bulge in situ that possibly happens at high redshift(Z>~2)
	{
		//cosmic flow puts too much gas in disk disck gets clumpy efficently torques to collapse disk to make a bulge
		double k=1;

		if(pow(10,(cent.return_gas()-cent.return_baryons()))>gas_frac_crit)
			{
				double dgasdt=k*25.*pow(10,cent.return_stellar_disk()-11.)*pow(10,1.5*(Z-log10(3))); //rate of change of gas mass in solar mass/yr

				/////check this!!! which frame is dt in? and does it matter? 
				double dt=get_time_from_redshift(Z2)-get_time_from_redshift(Z1); //time piriod where gas is accreting in Gyr

				double dgas=log10(dgasdt*dt)+9; //log(change in gas mass) in log(solar mass)
				if(dgas>cent.return_gas())
					{
						dgas=cent.return_gas();
					} //if there isn't enough gas in disk; use it all up!!

				//cent.put_Rhl_bulge(log10(pow(10,(cent.return_Rhl_bulge()+cent.return_stellar_bulge()))/(pow(10,cent.return_stellar_bulge())+pow(10,dgas))));

				cout << "CLUMPS  " << pow(10,cent.return_redshift())-1. << "  " << cent.return_mass() << "  " <<  cent.return_stellar_disk() << "  " << cent.return_gas() << "  " << dgas << "  " << dgas/cent.return_stellar_disk() <<endl;

				//double  rb_old=cent.return_Rhl_bulge();

				double rclump=cent.return_Rhl_disk()*M_PI/4.*cent.return_gas()/cent.return_baryons();

				if(cent.return_stellar_bulge()<3.)
					{
						cent.put_Rhl_bulge(rclump);
					}
				else
					{
						cent.put_Rhl_bulge(get_Rbulge(cent.return_stellar_bulge(),cent.return_Rhl_bulge(),dgas,rclump, 2.0/0.5,cent.return_gas()));
					}

				/*
				double logaa=2.*cent.return_stellar_total();
				double bb=pow(10,(2.*cent.return_stellar_bulge()-cent.return_Rhl_bulge()));
				double cc=pow(10,(2.*cent.return_stellar_disk()-cent.return_Rhl_disk()));
				double dd=pow(10,cent.return_stellar_bulge()+cent.return_stellar_disk())/(pow(10,cent.return_Rhl_bulge())+pow(10,cent.return_Rhl_disk()));
				cent.put_Rhl_bulge(logaa-log10(bb+1.*cc+(2./0.5)*dd));
				*/

				cent.add_stellar_bulge(dgas);
				cent.add_gas(-dgas);
				//cout << "clumpyaccretion  " << dgasdt << "  " << dgas << "  " << rb << "  " << cent.return_Rhl_bulge() << "  " <<  cent.return_mass() << "  " << Z << endl;
			}

		return cent;

	}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


void get_ssfr(satillite_halo cent, double Z1,double Z2)
	{
		double logA11=log10(0.0324)+11, alpha=3.45, beta=-0.35; 
		double dmdt=pow(10,logA11+(beta+1)*(cent.return_stellar_disk()-11)+alpha*Z1);
		double dt=get_time_from_redshift(Z2)-get_time_from_redshift(Z1); //time piriod in Gyr //again check!!
		double dm=log10(dmdt*dt);
		cent.add_gas(-dm);
		cent.add_stellar_disk(dm);
		cent.put_Rhl_disk(get_Rdisk(cent.return_stellar_total(),Z1));
	}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

void bathtub_model_accretion(satillite_halo cent, double Z1,double Z2)
	{
		double beta=0.14,alpha=2.4,mu=0.54,epsilon=0.02;
		double fga=1., eta=1.;
		double dt=get_time_from_redshift(Z2)-get_time_from_redshift(Z1); //time piriod where gas is accreting in Gyr
		double td=0.071*get_time_from_redshift(Z1);

		double Mdotb=80+pow(10,(1+beta)*cent.return_mass()-12)*pow(10,mu*Z1)/3*omega_baryon/omega_m;
		double Mdotsf=cent.return_gas()*epsilon/td;
		double Mdotg=fga*Mdotb-(mu+eta)*Mdotsf;
		double Mdots=(1.-fga)*Mdotb+mu*Mdotsf;

		cent.add_gas(log10(dt*pow(10,9)*Mdotg));
		cent.add_stellar_bulge(log10(dt*pow(10,9)*Mdots));
		cout << "bathtub_model_accretion  " << Mdots << "  "  << cent.return_mass() << "  " << Z1 <<endl;
	}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double starburst(satillite_halo central, satillite_halo satillite) //from Shankar 2014a
	{
		if(satillite.return_baryons()>central.return_baryons())
			{
				return log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()))+log10(0.54)+0.7*(central.return_baryons()-satillite.return_baryons());
			}
		else
			{
				return log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()))+log10(0.54)+0.7*(satillite.return_baryons()-central.return_baryons());
			}
	}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

void disk_regowth(satillite_halo central, satillite_halo satillite, double *gas, double *stellar_bulge,double *r_bulge, double *r_disk)
	{

		double newstars=starburst(central,satillite);
		cout << "gas  " << *gas << endl;
		if(newstars>=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas())))
			{
				newstars=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()));*gas=-300;
			}
		else
			{
				*gas=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas())-pow(10,newstars));
			}
		*stellar_bulge=log10(pow(10,central.return_baryons())+pow(10,satillite.return_baryons())-pow(10,*gas));
		*r_bulge=get_Rbulge(central.return_baryons(),central.return_Rhl(),satillite.return_baryons(),satillite.return_Rhl(),forb/0.5,*gas);  // // need to think of a better way of doing this!!!
		r_disk=r_bulge; // // and this line!!
		cout << "gas  " << *gas << "  " << newstars  << endl;

	}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double get_concentration(double logM, double logZ)
	{

		double z=pow(10,logZ)-1.;
		double a=0.520+(0.905-0.520)*exp(-0.617*pow(z,1.21));
		double b=-0.101+0.026*z;

		double logc200=a+b*(logM-12.);
		return pow(10,logc200);
	}



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


// // // this bit does the initializing either as disks/bulges or finds a galaxy from the catallogue \\  \\  \\

void find_galaxy(double DM, double Z, listlistgals *mergercat, std::stack<satillite_halo> *stack_Shalos)
	{
		
		bool onlydisk=true;
		bool onlybulge=false;
		bool findfromcatalogue=false; //these three set what galaxy type will be assigened to the DM halo
		
		listlistgals::iterator itZ = mergercat->begin(); //iterator for listlistgals points to the first item in listlistgals mergercat
		listgals::iterator itH; //generic iterator for listgals
			std::list<double> tmpDM,tmpG,tmpB,tmpD,newDM,tmpRD,tmpRB; //tempory lists
		double stars; //number of stars 
		bool foundhalo=false; //has matched it to a preexisting halo
		cout << "new_satillite  " << DM << "  " << Z<<endl;
		
		std::advance(itZ,(int)(Z*100)); // iterates from z = whatever to the redshift of the halo
		cout << "ok1\n";
	//	cout << DM << "  " << Z << "  " << mergercat->size() << "  " << itZ->size() << "  " << itZ->back().return_mass() << "  " << itZ->front().return_mass() << endl;
		
		if(findfromcatalogue)  //hunts through everything already simulated and finds a halo within a givin redhift range and simmilar mass and 'copys and pasts' to the halo to get a mtach
		{
			if(itZ->size()>6)
				{ 
					//maybe getrid of this condition?
					if(itZ->back().return_mass()<DM)
						{
							cout << "ok2  " << itZ->back().return_mass() << "\n";
							if(itZ->front().return_mass()>DM)
								{
									cout << "ok3  " << itZ->front().return_mass() << "\n";	
									foundhalo=true;
									itH=itZ->begin();
									tmpDM.clear();tmpG.clear();tmpB.clear();tmpD.clear();tmpRB.clear();tmpRD.clear();newDM.clear();
									newDM.push_front(DM);
									cout << itZ->size() << "  " << mergercat->size() << endl;
									while(itH!=itZ->end())
										{
											//cout<<"FOUND IT !!!!!!!!!!!!!!!!!!!\n";
											tmpDM.push_back(itH->return_mass());tmpG.push_back(itH->return_gas());tmpD.push_back(itH->return_stellar_disk());tmpB.push_back(itH->return_stellar_bulge());tmpRB.push_back(itH->return_Rhl_bulge());tmpRD.push_back(itH->return_Rhl_disk());
											//	cout <<  "ok4\n";			
											itH++;
										}
									cout <<"ok!\n";
									stack_Shalos->push(satillite_halo(DM,Z,interpol(tmpDM,newDM,tmpB).back(),interpol(tmpDM,newDM,tmpD).back(),interpol(tmpDM,newDM,tmpG).back()));
									stack_Shalos->top().put_Rhl_bulge(interpol(tmpDM,newDM,tmpRB).back()); stack_Shalos->top().put_Rhl_disk(interpol(tmpDM,newDM,tmpRD).back());

									cout << "we chose one!!  "; stack_Shalos->top().output_vars();
								}
							else if(DM-5<itZ->front().return_mass())
								{
									cout << "LARGER  " << DM -itZ->front().return_mass() <<endl; stack_Shalos->push(itZ->front());stack_Shalos->top().put_DM(DM);stack_Shalos->top().put_Z(Z);
								}
							else
								{
									cout<<"AHHHHHHHHH\n";itH=itZ->begin();while(itH!=itZ->end()){itH->output_vars();itH++;} exit(EXIT_FAILURE);
								}
						}
				}
		}
				
		if(!foundhalo || onlydisk) //this uses abundance matching (analytic fits we will get from somewhere) to slam a disk galaxy into the halo if a galaxy has not been found or you just want a disk
			{
				stars = Mdm_to_Mstellar(DM,Z);
				stack_Shalos->push(satillite_halo(DM,Z,-300,stars,get_Mgas(stars,Z,DM+log10(omega_baryon/omega_m))));
				stack_Shalos->top().put_Rhl_disk(get_Rdisk(stars,Z));
				stack_Shalos->top().put_Rhl_bulge(-300);
			}


	/*	
		if(oblybulge)
			{
				stars=Mdm_to_Mstellar(DM,Z);
				stack_Shalos->push(satillite_halo(DM,Z,stars,-300,get_Mgas(stars,Z,DM+log10(omega_baryon/omega_m))));
				stack_Shalos->top().put_Rhl_disk(get_Rdisk(stars,Z));
			}

	*/

		if(includeshutdown) //quenching 
			{
				if(get_M_halo_shutdown(Z) < DM) //if the halo mass is above a threshold that shock heats the infall gas quenching star formation then shut the galaxy down (stops star formation)
					{
						stack_Shalos->top().put_shutdown();
					}
			}

		cout << "added  ";
		stack_Shalos->top().output_vars(); //prints the data for the halo we just appended to the stack
	}



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



satillite_halo merge(satillite_halo central, satillite_halo satillite,double Z)
	{

		double DM=log10(pow(10,central.return_mass())+pow(10,satillite.return_mass()));
		double stellar_bulge=central.return_stellar_bulge(),stellar_disk=central.return_stellar_disk(),gas=central.return_gas();
		double mu=satillite.return_mu();
		double r_bulge=-300,r_disk=-300;
		bool major_merger=false;

		if(includemerger) // allows the code to iterate wiout all those messy mergers
			{
				if(mu<0.3)
					{				
						// //minor merger//
						cout << "minor  "<<mu << "  " <<endl;
						satillite.output_vars(); //prints stuff about the satilite
						central.output_vars(); //prints stuff about the central
						gas = log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()));//add gass mass together
						stellar_bulge = log10(pow(10,central.return_stellar_bulge())+pow(10,satillite.return_stellar_bulge())+pow(10,satillite.return_stellar_disk())); //total stellar mass of satilite + bulge mass of satilite (assumption good bad? maybe)
						stellar_disk = central.return_stellar_disk(); //disk remains unchanged
						double newstars = starburst(central,satillite); //making new stars in the disk
						gas = log10(pow(10,gas)-newstars); //deduct new stars from gas res
						stellar_disk = log10(pow(10,stellar_disk)+newstars); //add new stars to disk mass
						r_bulge = get_Rbulge(central.return_stellar_bulge(), central.return_Rhl_bulge(), satillite.return_stellar_total(),satillite.return_Rhl(),forb/0.5,gas); //from Shankar 14a conserving orbital energies therefore getting an expansion in galaxy size
						r_disk = get_Rdisk(log10(pow(10,stellar_disk)+pow(10,stellar_bulge)),Z); //reinitilize the disk size Shen et al year? (famous paper on sizes from SDSS)
					}
				else
					{
						// //major merger//
						major_merger=true;
						cout << "MAGOR!!!!!  "<<mu << "  " <<satillite.return_redshift() << "  " << satillite.return_mass() << "  " << central.return_mass() << "  " << satillite.return_stellar_total() << "  " << central.return_stellar_total() << "  " << satillite.return_gas() << "  " << central.return_gas() << "  " << satillite.return_Rhl() << "  " << central.return_Rhl() <<endl; //really long print thing
						cout<< "major  " << central.return_mass() << "  " << central.return_baryons() << "  " << satillite.return_mass() << "  " << satillite.return_baryons() << "  " << satillite.return_redshift() << endl; //really long print thing
						double newstars=starburst(central,satillite); //making new stars
						cout << "gas  " << gas << endl;
						if(newstars>=log10(pow(10,central.return_gas())+pow(10,satillite.return_gas()))) //if all the gas is used use all the gas avaliable
							{
								newstars = log10(pow(10,central.return_gas())+pow(10,satillite.return_gas())); //everything is converted to stars
								gas = -300; //makes gas really small to avoig crapping out a logarithm
							}
						else
							{
								gas = log10(pow(10,central.return_gas())+pow(10,satillite.return_gas())-pow(10,newstars)); //almost everything is converted to stars
							}
						stellar_bulge=log10(pow(10,central.return_baryons())+pow(10,satillite.return_baryons())-pow(10,gas)); // put disk stars in bulge
						stellar_disk=-300; //makes disk mass really small to avoig crapping out a logarithm
						r_bulge=get_Rbulge(central.return_baryons(),central.return_Rhl(),satillite.return_baryons(),satillite.return_Rhl(),forb/0.5,gas);  // // need to think of a better way of doing this!!!
						r_disk=-300; //logarithm safe
						cout << "gas  " << gas << "  " << newstars  << endl; //prints 
					}
			}
		satillite_halo new_cent(DM,Z,stellar_bulge,stellar_disk,gas); //creats a new instance of the class to hold the new galaxy that has been created in the Major merger
		new_cent.put_Rhl_bulge(r_bulge);
		new_cent.put_Rhl_disk(r_disk);
		//sets some bools about how the galaxy is going to act in the future
		if(satillite.return_extrememerger() || central.return_extrememerger())
			{
				new_cent.put_extrememerger();
			}
		if(central.return_shutdown())
			{
				new_cent.put_shutdown();
			}
		if(central.return_cgt4())
			{
				new_cent.put_cgt4();
			}
		if(major_merger || central.return_had_major_merger())
			{
				new_cent.put_had_major_merger();
			}
		cout << "new cent  "; new_cent.output_vars();
		return new_cent; //returns the new galaxy into the next timestep
	}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



//mergeing --- handles all of the data for merging and outputs to merge() to do the actual merging
std::stack<satillite_halo> merge_all_galaxies(std::stack<satillite_halo> *stack_Shalo, std::stack<accreted_mass> *stack_acc, std::stack<accreted_mass> *stack_mpro) //sent stacks in order full of the information we need
	{

		std::stack<satillite_halo> merged_halos; //where we put things that have been merged
		satillite_halo cent(-300,-300,-300,-300,-300);
		double top_redshift,next_redshift, MDM=0.0, Mstellar=0.0,T;
		std::list<satillite_halo> to_be_merged; //tempory stack to hold all of the things at redshift whatever that can then be neatly passed to merge()

		cent = stack_Shalo->top(); //class of the main progenotor
		MDM = cent.return_mass(); //DM mass of the main pro
		Mstellar = cent.return_stellar_total(); // Stellar mass of disk+bulge of main pro
		stack_Shalo->pop(); //pops the main pro out of the stack to avoid self merge

		while(! stack_acc->empty()) //while there is still stuff in the acc mass stack (proxy for redshif counter as there is allways mass to accret at every redshift)
			{
				////

				//cout << "SDFDSF  " << stack_Shalo->top().return_mass() <<endl;// << "  "  << stack_Shalo->top().return_redshift() << "  " << stack_acc->top().return_mass() << "   "<< stack_acc->top().return_redshift()<<endl;

				top_redshift=(*stack_acc).top().return_redshift(); //works out the zmax by looing at the top item of the acc mass stack
				accreted_mass tmp(stack_acc->top().return_mass(),stack_acc->top().return_redshift(),stack_acc->top().return_addmass()); //creats a tempory instance of the accreted mass class so we dont loose current top whilst looking at the next redshift
				stack_acc->pop(); //pops that thing we just saved
				//looks at the next redshift step 
				if(! stack_acc->empty()) //is the next stack isnt empty retruns its redshift value as the next redshift value
					{
						next_redshift=(*stack_acc).top().return_redshift(); 
					}
				else
					{
						next_redshift=0; //if there wasnt a next redshit then we are at presant day hurrah
					}
				stack_acc->push(tmp); // put top element that we saved back on top the stack

				T = get_time_from_redshift(next_redshift); //converts redshift into time

				if(!stack_Shalo->empty()) //if the halo stack isnt empty ie there are still mergers to do
					{
						while(! stack_Shalo->empty() && stack_Shalo->top().return_redshift() >= next_redshift) // while the stack isnt empty and we still have steps in redshift do stuff
							{
								to_be_merged.push_front(stack_Shalo->top()); // flipping the satalite halo stack into the to be merged
								if(MDM < to_be_merged.begin()->return_mass()-1.) //error handling such that if the mergiee large compared to main prog then handle carfully
									{
										to_be_merged.begin()->put_extrememerger();
									}
								to_be_merged.begin()->put_Tdy(T+get_T_dynamical_friction(MDM, to_be_merged.begin()->return_mass(),top_redshift)); //setting the time of the merger using the dynamical friction timescale + time of the merger of the halo
								to_be_merged.begin()->put_mu(pow(10,(to_be_merged.begin()->return_baryons()-cent.return_baryons()))); //define the ratio: Msat/Mprog this is where we define the size of the merger
								cout << "mudfgdfg  " << to_be_merged.begin()->return_mu() << "  " << to_be_merged.begin()->return_baryons() << "  " <<cent.return_baryons() <<endl; //outputting to terminal mud.... is for searching purpouses
								stack_Shalo->pop(); //popping the satalite off the top of the satilites
							}
					}


				//using a linked list we cycle through to be merged at every step such that whilst we are merging haloes if we come across a timestep where a galactic merger happens we do it
				std::list<satillite_halo>::iterator it = to_be_merged.begin(); //iteriator to to be merged so we can cycle through list of to be merged 
				while (it!= to_be_merged.end()) //while there are still elements in the list to look at
					{
						if(it->return_Tdy()<T)//if current timestep is closer to the present day then the time I defined it to be merged
							{
								cent = merge(cent,*it,top_redshift); //cent is the new main progentior 
								MDM = log10(pow(10.,MDM)+pow(10.,it->return_mass())); //adding the DM mass of central and satalite together done now to avoid over bulking galaxys when reinitilising
								//cout << "M stellar  " << it->return_baryons() << endl;
								Mstellar = log10(pow(10.,Mstellar)+pow(10.,it->return_stellar_total())); //stellar mass of new Main progentior
								it = to_be_merged.erase(it); //erases the merge
							}
						else
							{
								if(! it->return_shutdown()) // if it hasnt been shutdown
									{
										if(includeshutdown && get_M_halo_shutdown(next_redshift)<it->return_mass()) //if shutdown is a thing and the threshold has been met
											{
												it->put_shutdown(); //shutdown is threshold met
											}
										else
											{
												get_ssfr(*it,top_redshift,next_redshift); //have some star formation in the disk (somewhere mass is added to disk) 
											}
									}
								++it;
							}
					}
				cout << "MDMtest  " << MDM << "  " << stack_acc->top().return_mass() << "  " << stack_acc->top().return_addmass() <<endl; //prints

				if(reinitialiseDM) //if we need to reinitilise the DM halo (if its messy between timesteps)
					{
						MDM=stack_mpro->top().return_mass();
						cent.put_DM(MDM);
					}
				else //we just add the new dark matter onto the existing halo 
					{
						if(stack_acc->top().return_addmass())
							{
								cent.add_DM(stack_acc->top().return_mass()); //adding the DM
								MDM=log10(pow(10.0,MDM)+pow(10.0,stack_acc->top().return_mass())); //setting the variable
							}
						else //error handling for loosing mass (in messy mergers)
							{
								cent.subtract_DM(stack_acc->top().return_mass()); //
								MDM=log10(pow(10.0,MDM)-pow(10.0,stack_acc->top().return_mass())); //
							}
					}

				cent.put_Z(top_redshift); //putting variables into the class
				cout << "MDMtest2  " << MDM << "  " << stack_acc->top().return_mass() << "  " << stack_acc->top().return_addmass() <<endl; //prints
				if(includeshutdown) //if we are using shutdown
					{
						if(! cent.return_shutdown()) //if not shutdown
							{
								if(get_M_halo_shutdown(top_redshift)<cent.return_mass()) //if above shutdown threshold
									{
										cout <<"shutdown!!\n";//print
										cent.put_shutdown(); //shutdown
									}
							}
					}

				if(! cent.return_had_major_merger() && ! cent.return_shutdown()) //if it hasnt shutdown and hasnt had a mjor merger
					{
						cout << "reinitialise\n"; //then reinitilise everything
						if(cent.return_stellar_total() < Mdm_to_Mstellar(cent.return_mass(),top_redshift)) //if stellar mass is too low
							{
								cent.add_stellar_disk(log10(pow(10,Mdm_to_Mstellar(cent.return_mass(),top_redshift))-pow(10,cent.return_stellar_total()))); //accrete some stellar mass
							}
						if(cent.return_gas() < get_Mgas(cent.return_stellar_total(),top_redshift,cent.return_mass()+log10(omega_baryon/omega_m))) //if gas is too low 
							{
								cent.put_gas(get_Mgas(cent.return_stellar_total(),top_redshift,cent.return_mass()+log10(omega_baryon/omega_m))); //add some cold gas
							}
						if(includeclumpy_accretion)
							{
								cout << "clumpy_before  " << cent.return_Rhl_bulge()<<endl;
								cent=clumpyaccretion(cent,top_redshift,next_redshift,top_redshift,0.); // deals with clumpy accretion
								cout << "clumpy_after  " << cent.return_Rhl_bulge()<<endl; //prints
							}
						if(includebathtub_model_accretion)
							{
								bathtub_model_accretion(cent,top_redshift,next_redshift); //deals with bathtub
							}
					}
				
				if(includeexpansion) //if expansion is being included (dont worry too much) adaidbatic expansion after a violent accretion event
					{
						if(! cent.return_cgt4())
							{
								if(get_concentration(cent.return_mass(), top_redshift)>4)
									{
										cent.put_cgt4();cent.put_Rhl_bulge(2*cent.return_Rhl_bulge());cout << "EXPANSION  " << top_redshift << "  " << cent.return_mass() << "  "  << cent.return_baryons() << endl;
									}
							}
					} //include mass loss?!

				//double Mdisk=log10(pow(10,cent.return_gas())+pow(10,cent.return_stellar_disk()));
				double rdisk = get_Rdisk(cent.return_stellar_total(),top_redshift); //re initilise the disk
				cent.put_Rhl_disk(rdisk); //put the disk in the class
				cout << "merge_tree  " << cent.return_shutdown()  << "  "; //print
				cent.output_vars(); //print
				merged_halos.push(cent); puts this instance of the class on the stack
				stack_acc->pop(); //removes the top of the accretion stack 
				stack_mpro->pop(); //removes the top of the proginator stack these two things essentially move to the next redshift step
			} ///////// end of loop all stacks should be empty except to be merged as some dynamical fiction timescales go beond Z=0 


		cout << "end  " << stack_Shalo->empty() << "  " <<  to_be_merged.empty() << "\n"; //prints

		cout << "central ";
		cent.output_vars();
		//cout << "ok\n";
		if(! stack_Shalo->empty())
			{
				while (! stack_Shalo->empty())
				{
					cout << "satillite  ";
					stack_Shalo->top().output_vars(); //prints to inform of error
					stack_Shalo->pop();
				}
			}
		if(! stack_acc->empty())
			{
				cout << "acc "; 
				while(! stack_acc->empty())
					{
						cout << stack_acc->top().return_mass() << "  ";
						stack_acc->pop();//prints to inform of error
					} 
				cout << endl;
			}
		if(! to_be_merged.empty())
			{
				while(! to_be_merged.empty())
					{
						cout<<"merge   ";
						to_be_merged.front().output_vars();//prints to inform of error
						to_be_merged.pop_front();
					}
			} //these three things empty all stacks to avoid mem leak if error

		cout << "length  " << merged_halos.size() <<endl;
		return merged_halos;
	} 



/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


 
// main program



int main(int argc,char *argv[]) 
	{
		double fileDM,fileZ,stars;
		double endtreeDM,endtreeZ;

		int i;

		cout << "cosmology   " << omega_m << "  " << h << "  " << omega_baryon << "  " << sigma_8 << "  " <<endl;
		cout << "mergers?  " << includemerger << endl;
		cout << "halo shutdown? " << includeshutdown <<endl;
		cout << "disk regrowth?  " << includediskregrowth << endl;
		cout << "clumpy accretion?  " << includeclumpy_accretion <<endl;
		cout << "bathtub model accretion?  " << includebathtub_model_accretion << endl;
		cout << "gas dissipation?  " << includegasdissipation << endl;
		cout << "forb  " << forb << endl;

		listgals emp; //emp is an empty list or one step in redshift
		listlistgals mergercat; //the whole merger tree stored as sepreate lists where each list is a step in redshift
		while(i<=6000) //adds 6000 empty lists to the merger catalogue
			{
				mergercat.push_back(emp);
				i++;
			}
		cout <<"mergetcat!!!  " << mergercat.size()<<endl;
		listlistgals::iterator itZ; //iterator over redshift
		listgals::iterator itH; //iterates over haloes to be merged
		std::stack<satillite_halo> tmpmergedstack; //tempory stack...

		//stack of satilite halos
		std::stack<satillite_halo> stack_Shalos;
		//stack of accreted DM
		std::stack<accreted_mass> stack_acc; //storing non resolved DM mass
		//stack of mass of the main progenitor
		std::stack<accreted_mass> stack_mpro; 


		// get merger tree filename fron argument //
		// decode arguments
		if(argc <2) //error handling
			{
				printf("You must provide at least one argument\n");
				exit(0);
			}
		// report settings
		string filename=argv[1]; //passing the filename of the merger catalogue
		cout << filename << endl;
		int k=-1;
		//while (k<600){
		cout<<"RESET\n";
		k+=1;
		ifstream myfile; //mergertree stored as text file 
		myfile.open (filename);
		string line; //defines a vsting varible to read lines from the file
		int j;
		i=0;
		string tmp;



		while ( getline (myfile,line)) //loops over all the lines in the file
			{
				std::istringstream ss(line);
				std::istream_iterator<std::string> begin(ss), end;
				std::vector<std::string> words(begin, end); //seperates line into space seperated stuff (eg a word)

				if(words[0]=="START") //defines start of merger tree
					{
						i+=1;
					}
				else if(words[0]=="mpro") //defines a main progonator (backbone of the merger tree)
					{
						stack_mpro.push(accreted_mass((double)atof(words[1].c_str()), (double)atof(words[2].c_str()), 1)); //sends theaccretedmass function Mass, Redshift & A true value stores that instance of that class in a stack called m pro
						j=0;
					}
				else if(words[0]=="1") //looking at main progeniter branch
					{
						if(words[1]!="acc") //looking at halo (ie not the accreted mass)
						{
							if(j==0)
								{
									endtreeDM=(double)atof(words[2].c_str()); //defines the mass of the halo at the end of the tree 
									endtreeZ=(double)atof(words[1].c_str()); //defines the redshift of the halo at the end of the tree
								}
							else
								{
									fileDM=(double)atof(words[2].c_str()); //saves the mass of the halo at...
									fileZ=(double)atof(words[1].c_str()); //a given redshift
									find_galaxy(fileDM,fileZ,&mergercat,&stack_Shalos); //looping through listlistgal (mergercat) to find a simmilar halo and match to it
								}
							j+=1;
						}
						else if(words[1]=="acc") //looking at the DM halo accretion
						{
							if(reinitialiseDM) //if the DM is being re initilised 
								{
									stack_acc.push(accreted_mass(-300,(double)atof(words[2].c_str()),1));
								}
							else // if the DM is being accreated
								{
									stack_acc.push(accreted_mass((double)atof(words[3].c_str()),(double)atof(words[2].c_str()),(bool)atof(words[4].c_str())));
								}
						}
					}


				if(words[0]=="START") //top of merger treeat z = 0
					{
						find_galaxy(endtreeDM,endtreeZ,&mergercat,&stack_Shalos); // populate the original main progeniter halo at zmax
						
						cout << "NEW\n";
						tmpmergedstack = merge_all_galaxies(&stack_Shalos,&stack_acc,&stack_mpro); //original branch and all merging branches are populated such that we can 'pop' them off the stack in order and merge
						cout <<"DONE\n";
						
						while(!tmpmergedstack.empty()) //if its not empty empty and print
							{
								if(tmpmergedstack.top().return_mass()==tmpmergedstack.top().return_mass()){
									itZ=mergercat.begin();
									std::advance(itZ, (int)(tmpmergedstack.top().return_redshift()*100.));
									itH=itZ->begin();
									while(itH->return_mass()>tmpmergedstack.top().return_mass())
										{
											itH++;
										}
									itZ->insert(itH,tmpmergedstack.top());
								}
								tmpmergedstack.pop();
							}
						while(!tmpmergedstack.empty())
							{
								cout<<tmpmergedstack.top().return_redshift()<<endl;tmpmergedstack.pop();
							}
						cout << "done merging"<< endl;

						while(! stack_Shalos.empty())
							{
								cout << "Shalo sdaf  " << stack_Shalos.top().return_mass()<<"  " << stack_Shalos.top().return_redshift() <<endl; stack_Shalos.pop();
							}
						while(! stack_acc.empty())
							{
								cout << "acc sdaf " << stack_acc.top().return_mass()<<"  " << stack_acc.top().return_redshift() <<endl; stack_acc.pop();
							}
					}
			}
		myfile.close();

		return(0);
	}

