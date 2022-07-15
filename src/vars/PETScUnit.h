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
                cout << " <ERROR> petsc_unit takes 2-3 values: eps unit_k [usehdf5? 1:0] </ERROR>" << endl;
            return 1;
        }

        float v = atof(argv[1]);
        s->ctrl.petsc_unit_eps = v;
        s->ctrl.petsc_unit_k = atoi(argv[2]);
        if (argc==4)
            s->ctrl.petsc_hdf5 = atoi(argv[3]);
        //if (argc==4) {
        //    if (argv[3]=="hdf5")
        //        s->ctrl.petsc_hdf5=true;
        //    else if (argv[3]=="binary")
        //        s->ctrl.petsc_hdf5=false;
        //    else {
        //        if ( ui->oncoutpe() ){
        //            cout << " <ERROR> please specify whether output in hdf5 or binary format </ERROR>" << endl;
        //        }
        //        return 1;
        //    }
        //}

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
        st << setw(10) << s->ctrl.petsc_hdf5;
        return st.str();
    }

    PETScUnit(Sample *sample) : s(sample) {
        s->ctrl.petsc_unit_eps = 0;
        s->ctrl.petsc_unit_k = 0;
        s->ctrl.petsc_hdf5 = false;
        s->ctrl.petsc_unit = false;
    };
};
#endif

// Local Variables:
// mode: c++
// End:
