/******************************************************************************
*****	Langevin dynamics calculation of a one-dimensional harmonic	*******
*****	oscillator.	  	      	   		   		*******
*****									*******
*****   Random number generator from Numerical recipes in C++ (3rd ed.) *******
*****   Normal deviate generation adapted from Numerical recipes in C   *******
*****                                                                   *******
*****   Input file is a column of numbers called tcf.in with the        *******
*****   following entries in this order:                                *******
*****   timestep, total time, initial position, mass, friction coeff.   *******
*****   temperature, force constant, verlet   (0 for on, 1 for off),    *******
*****   random force (1 for on, 0 for off)                              *******
******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "nr3.h"
#include "ran.h"
using namespace std;

// Function prototypes
Ran myran(23464362);
void bivargasdev(double&, double&, double, double, double, double, Ullong);
void force(double, double&);
void integrator(const double[], double, int, double&);
Ullong seed;
Ran myran2(seed);

main()
{
  //Declare variables
  const double k_b=1.;
  double dt, total_time, mass, gamma, temp, force_const;
  double initial_posn, posn, posn_tally, posn2_tally;
  double vel, initial_vel, vel_tally, vel2_tally;
  double c0, c1, c2;
  double time = 0.;
  int count = 0, i = 0, j = 0, k=0;                     
  int lang_on;                                                     // flag that controls langevin vs. velocity verlet
  int random_force;                                                // flag that controls random force
  double f;                                                        // force
 static double posn_posn_corfun[1000000], vel_vel_corfun[1000000];      // Position position, velocity-velocity correlation functions
  int paths = 100;                                                      // Number of paths requested
  int steps_between_initials = 100;                                     // Number of steps between initial conditions
  Ullong seed;                                                    // This is the seed for the "dynamic" uniform deviate generator.
  int nsteps;                                                  // The number of steps to propagate the initial condition
  double integral = 0.;
  int freeparticle;

  // Open input and output files
  ifstream input;
  ofstream output;
  ofstream trajectory;
  ofstream trajectoryv;
  ofstream energy;
  ofstream tcf;
  input.open("tcf.in", ios::in);
  output.open("tcf.out", ios::out);
  trajectory.open("trajectory.out", ios::out);
  trajectoryv.open("velocities.out", ios::out);
  energy.open("energies.out", ios::out);

  // Read in data from langevin.in
  input >> dt >> total_time >> initial_posn >> initial_vel >> mass >> gamma >> temp >> force_const >> lang_on >> seed >> freeparticle;
  
  // Initialize useful varaibles
  double gdt = gamma * dt;
  double dt2 = dt*dt;
  int n = total_time/dt;
  double beta = 1. / (k_b * temp);

  // Set up the parameters if a velocity verlet algorithm is preferred intstead...
  if (lang_on == 0)
    {
      c0 = 1.;
      c1 = 1.;
      c2 = 0.5;
    }
  else
    {
      c0 = exp ( - gdt );
      c1 = (1. - c0) / gdt;
      c2 = (1. - c1) / gdt;
    }
  double c02 = exp ( -2 * gdt );

  double stdev_r = sqrt( ( dt / (mass * gamma * beta) ) * (2. - (3. - 4.*c0 + c02)/gdt));
  double stdev_v = sqrt((1. - c02)/( mass * beta));
  double correlation_rv = pow(1. - c0, 2) / ( beta * stdev_r * stdev_v * mass * gamma);

  static double parameter1 = sqrt(1- correlation_rv*correlation_rv);
  static double parameter2 = correlation_rv * stdev_v / stdev_r;

  // Print the headers for the output file and the trajectory files
  output << "******************************************************************************" << endl;
  output << "*****     Langevin dynamics for a single particle in a one dimensional  ******" << endl;
  output << "*****     harmonic potential (report file)                              ******" << endl;
  output << "*****                                                                   ******" << endl;
  output << "*****     Nick Preketes, 22 July 2008                                   ******" << endl;
  output << "******************************************************************************" << endl << endl << endl;

  output << "SPECIFICS OF SIMULATION" << endl;
  output << "Timestep:              " << dt << endl;
  output << "Timestep squared:      " << dt2 << endl;
  output << "Total simulation time: " << total_time << endl;
  output << "Initial position:      " << initial_posn << endl;
  output << "Mass:                  " << mass << endl;
  output << "Friction constant:     " << gamma << endl;
  output << "Temperature:           " << temp << endl;
  output << "Force constant:        " << force_const << endl;
  output << "Coefficient 0:         " << c0 << endl;
  output << "Coefficient 1:         " << c1 << endl;
  output << "Coefficient 2:         " << c2 << endl;
  output << "Std. Deviation r:      " << stdev_r << endl;
  output << "Std. Deviation v:      " << stdev_v << endl;
  output << "Correlation rv:        " << correlation_rv << endl;
  output << "Paramter 1:            " << parameter1 << endl;
  output << "Parameter 2:           " << parameter2 << endl;
  output << "P1 * sig_v:            " << parameter1 * stdev_v << endl;
  output << "Number of paths:       " << paths << endl;
  output << "Langevin? (1 for yes, 0 for no) : " << lang_on << endl;
  output << "Steps between initial conditions: " << steps_between_initials << endl;

  trajectory << "******************************************************************************" << endl;
  trajectory << "*****     Langevin dynamics for a single particle in a one dimensional  ******" << endl;
  trajectory << "*****     harmonic potential (trajectory file)                          ******" << endl;
  trajectory << "*****                                                                   ******" << endl;
  trajectory << "*****     Nick Preketes, 22 July 2008                                   ******" << endl;
  trajectory << "******************************************************************************" << endl << endl << endl;

  trajectory << "   Time   " << "     " << "    Position     " << "      " << "     Velocity    " << endl;
  trajectory << "----------------------------------------------------" << endl;


  for (i=0; i < paths; i++)
    {
      // Set the number of steps to propagate the initial condition
      // If this is the first trajectory, run it for a very long time to allow equilibration with the bath.
      // This length is the specified length of trajectory t_total.
      if (i == 0)
	nsteps = n;
      else
	nsteps = steps_between_initials;

      for (j = 0; j < nsteps; j++)
	{
	  double rand_r, rand_v;
	  if (lang_on == 0 || random_force == 0)
	    {
	      rand_r = rand_v = 0;
	    }
	  if (lang_on == 1)
	    bivargasdev(rand_r, rand_v, parameter1, parameter2, stdev_v, stdev_r,0);

	  if (rand_r != rand_r || rand_v != rand_v)
	    {
	    cout << "Random variables are not numbers due to negative square root!" << endl;
	    break;
	    }
	  

	  if (freeparticle == 0)
	    {
	      // Call the force function and advance initial positions and half advance initial velocities...
	      force(initial_posn, f);
	      initial_posn += c1 * dt * initial_vel + c2 * dt2 * f / mass + rand_r;
	      initial_vel = c0 * initial_vel + (c1- c2) * dt * f /mass;

	      // Call the force function again and fully advance velocities...
	      force(initial_posn, f);
	      initial_vel += rand_v + c2 * dt * f / mass;
	    }
	  if (freeparticle == 1)
	    {
	      initial_posn += c1 * dt * initial_vel  + rand_r;
	      initial_vel = c0 * initial_vel;

	      // fully advance velocities...
	      initial_vel += rand_v;
	    }
	}
    
  // Initialize position and velocity ( and time).
  posn = initial_posn;
  vel = initial_vel;
  time = 0;

  // Update the seed before entering propagation...
  /*seed += 1;*/

  for (k = 0; k < n; k++)
    {
      double rand_r, rand_v;
      if (lang_on == 0 || random_force == 0)
	{
	  rand_r = rand_v = 0;
	}

      if (lang_on == 1)
	bivargasdev(rand_r, rand_v, parameter1, parameter2, stdev_v, stdev_r,seed);

      if (rand_r != rand_r || rand_v != rand_v)
	{
	  cout << "Random variables are not numbers due to negative square root!" << endl;
	  break;
	}

      if (freeparticle == 0)
	{
      	  // Call the force function and advance  positions and half advance  velocities...
	  force( posn, f);
	  posn += c1 * dt * vel + c2 * dt2 * f / mass + rand_r;
	  vel = c0 * vel + (c1- c2) * dt * f /mass;

	  // Call the force function again and fully advance velocities...
	  force(posn, f);
	  vel += rand_v + c2 * dt * f / mass;
	  time += dt;
	}
      if (freeparticle == 1)
	{
	  posn += c1 * dt * vel + rand_r;
	  vel = c0 * vel;

	  // fully advance velocities...
	  vel += rand_v;
	  time += dt;
	}

	if (i == paths-1)
	  {
	   // Print trajectory for the last run through...
	      trajectory << time << "    " << setw(12) << posn << endl;
	      posn_tally += posn;
	      posn2_tally += posn*posn;
	      vel_tally += vel;
	      vel2_tally += vel * vel;
	  }
    }    
}

  //Print the results to the output file.
  output << "******************************************************************************" << endl;
  output << "************                    RESULTS                           ************" << endl;
  output << "************ Average position:         " << setprecision(16) << posn_tally/n  << endl;
  output << "************ Average position squared: " << setprecision(16) << posn2_tally/n  << endl;
  output << "************ Average velocity:         " << setprecision(16) << vel_tally/n  << endl;
  output << "************ Average velocity squared: " << setprecision(16) << vel2_tally/n  << endl;
  output << "*******************************************************************************" << endl << endl << endl;
  output << "JOB COMPLETE." << endl;
  //End main function
  return 0;
}

/* <<<<<<<<<<<<<<<<<< Acceleration function definition >>>>>>>>>>>>>>>>>>>>>>*/

void force(/* in */ double posn, 
		/* out */ double& f)
{
  f = -4*posn*posn*posn+10*posn ;
}


/*<<<<<<< Normal deviate generator for bivariate normal distribution >>>>>>>>*/

void bivargasdev(/* out */ double& rand_r, 
		 /* out */ double& rand_v, 
		 /* in  */ double parameter1, 
		 /* in  */ double parameter2,
		 /* in  */ double stdev_v,
		 /* in  */ double stdev_r,
		 /* in  */ Ullong seed)
{
  double fac,rsq,v1,v2;

  if (seed == 0)
    {
	do
		{
		v1=2.0*myran.doub()-1.0;
		v2=2.0*myran.doub()-1.0;
		rsq=v1*v1+v2*v2;
		}
	while (rsq >= 1.0 || rsq == 0.0);
	fac=sqrt(-2.0*log(rsq)/rsq);
	rand_r = v2*fac*stdev_r;
	rand_v = v1 * fac * stdev_v * parameter1 + rand_r * parameter2;
    }
  else
    {
    	do
		{
		v1=2.0*myran2.doub()-1.0;
		v2=2.0*myran2.doub()-1.0;
		rsq=v1*v1+v2*v2;
		}
	while (rsq >= 1.0 || rsq == 0.0);
	fac=sqrt(-2.0*log(rsq)/rsq);
	rand_r = v2*fac*stdev_r;
	rand_v = v1 * fac * stdev_v * parameter1 + rand_r * parameter2;
    }
}



