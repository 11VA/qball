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
// PETScUnit.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef PETSCUNIT_H
#define PETSCUNIT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/SlaterDet.h>

class PETScUnit : public Var {
    Sample *s;

public:

    char const*name ( void ) const {
        return "petsc_unit";
    };

    int set ( int argc, char **argv ) {
        if ( argc < 3) {
            if ( ui->oncoutpe() )
                cout << " <ERROR> petsc_unit takes 2 or more values: eps unit_k savefreq </ERROR>" << endl;
            return 1;
        }

        float v = atof(argv[1]);
        s->ctrl.petsc_unit_eps = v;
        s->ctrl.petsc_unit_k = atoi(argv[2]);
        if (argc>3){
            s->ctrl.petsc_savefreq = atoi(argv[3]);
            if(s->ctrl.petsc_savefreq==0){
                s->ctrl.petsc_savefreq=1;
            }
            if (argc>4){
                if (strcmp(argv[4],"KE")==0){
                    s->ctrl.petsc_KE=true;
                    s->ctrl.petsc_Vr=false;
                }
                else if(strcmp(argv[4],"KEVr")==0){
                    s->ctrl.petsc_KE=true;
                    s->ctrl.petsc_Vr=true;
                }
                else if(strcmp(argv[4],"Vr")==0){
                    s->ctrl.petsc_KE=false;
                    s->ctrl.petsc_Vr=true;
                }
            }
        }
        s->ctrl.petsc_unit = true;

        return 0;
    }

    string print (void) const {
        ostringstream st;
        st.setf(ios::left,ios::adjustfield);
        st << setw(10) << name() << "eps = ";
        st.setf(ios::right,ios::adjustfield);
        st << setw(10) << s->ctrl.petsc_unit_eps<<" unit_k = ";
        st << setw(10) << s->ctrl.petsc_unit_k;
        st << setw(10) << s->ctrl.petsc_savefreq;
        st << setw(10) <<"compute KE="<<s->ctrl.petsc_KE;
        st << setw(10) <<"compute Vr="<<s->ctrl.petsc_Vr;
        return st.str();
    }

    PETScUnit(Sample *sample) : s(sample) {
        s->ctrl.petsc_unit_eps = 0;
        s->ctrl.petsc_unit_k = 0;
        s->ctrl.petsc_savefreq = 1;
        s->ctrl.petsc_unit = false;
        s->ctrl.petsc_Vr=true;
        s->ctrl.petsc_KE=true;
    };
};
#endif

// Local Variables:
// mode: c++
// End:
