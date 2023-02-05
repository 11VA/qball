////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// VectorPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>
#include <math/d3vector.h>
#include <qball/Basis.h>
#include "UnitCell.h"
#include "Messages.h"
#include "Context.h"

#ifndef VECTORPOTENTIAL_H
#define VECTORPOTENTIAL_H

using namespace std;

class VectorPotential {

public:

  enum class Dynamics {
    NONE,
    POLARIZATION
  };

VectorPotential(Dynamics dyn, const D3vector & initial_external, const D3vector & initial_induced,const D3vector & initial_velocity, const D3vector & initial_accel, double laser_freq, D3vector laser_amp, string envelope_type, double envelope_center, double envelope_width, double lrc_alpha):
    dynamics_(dyn),
    initial_external_(initial_external),
    initial_induced_(initial_induced),
    initial_velocity_(initial_velocity),
    initial_accel_(initial_accel),
    laser_freq_(laser_freq),
    laser_amp_(laser_amp),
    envelope_type_(envelope_type),
    envelope_center_(envelope_center),
    envelope_width_(envelope_width),
    lrc_alpha_(lrc_alpha) {

    external_ = initial_external_;
    induced_ =  initial_induced_;
    value_ = induced_ + external_;
    value2_ = norm(value_);
    velocity_ = initial_velocity_;
    accel_ = initial_accel_;

    if(norm(external_) > 1e-15 && norm(laser_amp) > 1e-15) {
        Messages::fatal("Cannot specify a vector potential and a laser at the same time.");
    }

    if(fabs(laser_freq) < -1e-15 && norm(laser_amp) > 1e-15) {
        Messages::fatal("The laser_freq cannot be zero. Zero is specifically for static electrical field");
    }


}

  double * get_kpgpa(const Basis & basis) const {
    const double * kpgpa2 = get_kpgpa2(basis);
    double * kpgpa = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      kpgpa[ig] = sqrt(kpgpa2[ig]);
    }
    delete [] kpgpa2;
    return kpgpa;
  }
  
  double * get_kpgpa2(const Basis & basis) const {
    double * kpgpa2 = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      kpgpa2[ig] = basis.kpg2_ptr()[ig] + value2();
      kpgpa2[ig] -= 2 * value_[0]*basis.kpgx_ptr(0)[ig];
      kpgpa2[ig] -= 2 * value_[1]*basis.kpgx_ptr(1)[ig];
      kpgpa2[ig] -= 2 * value_[2]*basis.kpgx_ptr(2)[ig];
    }
    return kpgpa2;
  }

  double * get_kpgpai(const Basis & basis) const {
    const double * kpgpa2 = get_kpgpa2(basis);
    double * kpgpai = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      if(kpgpa2[ig] > 0.0){
	kpgpai[ig] = 1.0/sqrt(kpgpa2[ig]);
      } else {
	kpgpai[ig] = 0.0;
      }
    }
    delete [] kpgpa2;
    return kpgpai;
  }
  
  double * get_kpgpax(const Basis & basis, int j) const {
    if(j == 0){ 
      double * kpgpax = new double[3*basis.localsize()];
      for(int j2 = 0; j2 < 3; j2++){
	for(int ig = 0; ig < basis.localsize(); ig++){
	  kpgpax[j2*basis.localsize() + ig] = basis.kpgx_ptr(j2)[ig] - value_[j2];
	}
      }
      return kpgpax;
    } else {
      double * kpgpax = new double[basis.localsize()];
      for(int ig = 0; ig < basis.localsize(); ig++){
	kpgpax[ig] = basis.kpgx_ptr(j)[ig] - value_[j];
      }
      return kpgpax;
    }
  }
  
  const D3vector & value() const {
    return value_;
  }

  const double & value2() const {
    return value2_;
  }

  void calculate_acceleration(const double & dt, const D3vector& total_current, const UnitCell & cell){
    //update the velocity to time t - dt/2
    induced_ += 0.5*dt*velocity_ ;
    velocity_ += 0.5*dt*accel_;

    if(dynamics_ == Dynamics::POLARIZATION){
      accel_ = -4.0*M_PI*total_current/cell.volume();
        if(lrc_alpha_ > 1.0e-4) {
            accel_ = lrc_alpha_*total_current/cell.volume();
        }
    } else {
      accel_ = D3vector(0.0, 0.0, 0.0);
    }

    //update the velocity to time t
    velocity_ += 0.5*dt*accel_;
    induced_ += 0.5*dt*velocity_;
  }
  
  void propagate(double time, const double & dt){
     
    // evaluation of analytic form
    if(norm(laser_amp_) > 0.0 && envelope_type_ == "constant" && laser_freq_>0.0){ 
         double laser_periods=envelope_center_;
       if ( laser_periods == 0) external_ = -sin(laser_freq_*time)*laser_amp_/laser_freq_;
       else if ( laser_periods >0){
           double Tend=laser_periods*(4.136/(laser_freq_*27.2114))*41.341374575751; // atomic time
           if (time<Tend){
               external_ = -sin(laser_freq_*time)*laser_amp_/laser_freq_;
               if ((time+dt)>Tend){
                   eop_=1;
               }
           }
           else{
               eop_=0;
               external_=D3vector (0,0,0);
           }
       }
    }
    if(norm(laser_amp_) > 0.0 && envelope_type_ == "cos") external_ = (cos(laser_freq_*time)-1)*laser_amp_/laser_freq_;
    // Static external field
    if(norm(laser_amp_) > 0.0 && envelope_type_ == "constant" && laser_freq_ == 0.0) external_ = -laser_amp_ * time;
    // Static external field
	if(norm(laser_amp_) > 0.0 && envelope_type_ == "delta" && laser_freq_ == 0.0) {
		if (time > envelope_center_) external_ = -laser_amp_;
		else external_ = D3vector (0,0,0);
	}
    // Delta field
    if(norm(laser_amp_) > 0.0 && envelope_type_ == "delta" && laser_freq_ == 0.0) external_ = laser_amp_;
    // Numerical integration for guassian d(A/c)/dt = - E = - Amp * Normalization_factor * cos(wt) * Gaussing(center,width)
    if(norm(laser_amp_) > 0.0 && envelope_type_ == "gaussian") external_ +=  - dt*(laser_amp_/(envelope_width_ * sqrt(M_PI*2.0))) * cos(laser_freq_*time) * exp(-((time-envelope_center_)*(time-envelope_center_))/(2.0*envelope_width_*envelope_width_)) ;

// A = A0*F(t) = A0 * cos^2(wt/2N)* sin(wt)
    if(norm(laser_amp_) > 0.0 && envelope_type_ == "Asin2" ) {
        double N_sin = envelope_width_;
        double omega_Aext = laser_freq_ /27.21138505;  // eV2hartree
        double envelope_NT_  = N_sin*2.0*M_PI/omega_Aext;
        D3vector A0_ext;
        for(int i=0; i<3 ; i++) {
            A0_ext[i]= 5.338e-9*sqrt(laser_amp_[i]);
        }

        bool test_pulse=false;
        if(test_pulse) {
            double fac1 = N_sin/(N_sin-1.0);
            double fac1_inv  = 1.0/fac1;
            double fac2 = N_sin/(N_sin+1.0);
            double fac2_inv  = 1.0/fac2;
            if(time < envelope_NT_) {
                if(N_sin < 1.1) { // actually N=1

                    external_=A0_ext*(1.0/omega_Aext*pow(sin(omega_Aext*time/2.0),2)- 4.0/ omega_Aext*pow(sin(omega_Aext*time),2)  );
                }
                else {

                    external_=A0_ext*(1.0/omega_Aext*pow(sin(omega_Aext*time/2.0),2) -fac1/2.0/omega_Aext*pow(sin(omega_Aext*time/2.0*fac1_inv),2) -fac2/2.0/omega_Aext*pow(sin(omega_Aext*time/2.0*fac2_inv),2));
                }

            }
            else {
                external_ = D3vector(0.0, 0.0, 0.0);
            }
        }
        else {
            double tt =time - 0.5*envelope_NT_;
            if(fabs(tt) < 0.5*envelope_NT_) {
                external_=-A0_ext/omega_Aext*pow(cos(M_PI*tt/envelope_NT_),2)* sin(  omega_Aext*tt);
            }
            else {
                external_ = D3vector(0.0, 0.0, 0.0);
            }
        }
    }

    value_ = induced_ + external_;
    value2_ = norm(value_);
  }
  int get_eop(){
      return eop_;
  }
  const D3vector & vp_velocity() const {
    return velocity_;
  }

  const D3vector & vp_induced() const {
    return induced_;
  }  

  const D3vector & vp_accel() const {
    return accel_;
  }  
  
private:
  Dynamics dynamics_;

  D3vector initial_external_;
  D3vector external_;
  D3vector initial_induced_;  
  D3vector induced_;
  D3vector value_;
  D3vector initial_velocity_;  
  D3vector velocity_;
  D3vector initial_accel_;
  D3vector accel_;
  double value2_;

  double laser_freq_;
  D3vector laser_amp_;
  string envelope_type_;
  double envelope_width_;
  double envelope_center_;
  double lrc_alpha_;
  int eop_=0;
};
#endif

// Local Variables:
// mode: c++
// End:

