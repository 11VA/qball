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
// MRRKWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MRRKWavefunctionStepper.C,v 1.8 2011-06-02 15:56:19 schleife Exp $

#include <config.h>

#include "MRRKWavefunctionStepper.h"
#include "SelfConsistentPotential.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
MRRKWavefunctionStepper::MRRKWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s)
    : tddt_(tddt), WavefunctionStepper(wf,tmap), ef_(ef), s_(s),k1(s.wf),k2(s.wf),k3(s.wf),Vk2(s.wf),Vk3(s.wf) {
}


////////////////////////////////////////////////////////////////////////////////
void MRRKWavefunctionStepper::update(Wavefunction& dwf) {
    //k1=phi0
    copy(k1,wf_);
    k2.clear();

    //ef_.KE(k2); //k2=T|k1>
    ////k2=k1-idt/2*T|k1>
    //for ( int ispin = 0; ispin < k2.nspin(); ispin++) {
    //    for ( int ikp = 0; ikp < k2.nkp(); ikp++ ) {
    //        k2.sd(ispin, ikp)->c() *= -0.5*complex<double>(0,1)*tddt_; //k2=-idt/2*T|k1>
    //        k2.sd(ispin, ikp)->c().axpy(1,k1.sd(ispin,ikp)->c()); //k2=k1-i*dt/2*T|k1>
    //    }
    //}

    //evaluate the kinetic energy analytically
    ef_.expKE(k2,tddt_/2);


//k3=k1-i*dt/2*T|k2>-i*dt*V(k2)
    copy(k2); //set wf to k2
    updateV();
    Vk2.clear();
    ef_.Vnl(Vk2);
    ef_.Vr(Vk2); //get V(k2)

    k3.clear();
    ef_.KE(k3); //k3=T|k2>

    for ( int ispin = 0; ispin < k2.nspin(); ispin++) {
        for ( int ikp = 0; ikp < k2.nkp(); ikp++ ) {
            k3.sd(ispin, ikp)->c() *= -0.5*complex<double>(0,1)*tddt_; //k3=-idt/2*T|k2>
            k3.sd(ispin, ikp)->c().axpy(-complex<double>(0,1)*tddt_,Vk2.sd(ispin,ikp)->c()); //k3=-i*dt/2*T|k2>-i*dt*V|k2>
            k3.sd(ispin, ikp)->c().axpy(1,k1.sd(ispin,ikp)->c()); //k3=k1-i*dt/2*T|k2>-i*dt*V|k2>
        }
    }
    copy(k3);//set wf to k3
    updateV();

    Vk3.clear();
    ef_.Vnl(Vk3);
    ef_.Vr(Vk3); //get V(k3)

    copy(wf_,k2); //phi1=k2
    for ( int ispin = 0; ispin < k2.nspin(); ispin++) {
        for ( int ikp = 0; ikp < k2.nkp(); ikp++ ) {
            wf_.sd(ispin, ikp)->c().axpy(1,k3.sd(ispin,ikp)->c()); //phi1=k2+k3
        }
    }

    ef_.KE_in_place(wf_); //phi1=T|k2+k3>
    for ( int ispin = 0; ispin < k2.nspin(); ispin++) {
        for ( int ikp = 0; ikp < k2.nkp(); ikp++ ) {
            wf_.sd(ispin, ikp)->c().axpy(1,Vk2.sd(ispin,ikp)->c()); //phi1=T|k2+k3>+Vk2
            wf_.sd(ispin, ikp)->c().axpy(1,Vk3.sd(ispin,ikp)->c()); //phi1=T|k2+k3>+V(k2)+V(k3)
            wf_.sd(ispin,ikp)->c()*=-complex<double>(0,1)*0.5*tddt_;  //phi1=-i*dt/2[T|k2+k3>+V(k2)+V(k3)]
            wf_.sd(ispin, ikp)->c().axpy(1,k1.sd(ispin,ikp)->c()); //phi1=k1-i*dt/2[T|k2+k3>+V(k2)+V(k3)]
            s_.hamil_wf->sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c(); //copy to hamil_wf
        }
    }

}

void MRRKWavefunctionStepper::copy(Wavefunction& newwf, const Wavefunction& oldwf) {
    tmap_["copy"].start();
    for ( int ispin = 0; ispin < oldwf.nspin(); ispin++) {
        for ( int ikp = 0; ikp < oldwf.nkp(); ikp++ ) {
            newwf.sd(ispin, ikp)->c() = oldwf.sd(ispin, ikp)->c();
        }
    }
    tmap_["copy"].stop();
}
void MRRKWavefunctionStepper::copy(const Wavefunction& oldwf) {
    tmap_["copy"].start();
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++) {
        for ( int ikp = 0; ikp < wf_.nkp(); ikp++ ) {
            wf_.sd(ispin, ikp)->c() = oldwf.sd(ispin, ikp)->c();
            s_.hamil_wf->sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
        }
    }
    tmap_["copy"].stop();
}
void MRRKWavefunctionStepper::updateV() {
    tmap_["updateV"].start();
    ef_.hamil_cd()->update_density();
    ef_.update_hamiltonian();
    ef_.update_vhxc();
    tmap_["updateV"].stop();
}
