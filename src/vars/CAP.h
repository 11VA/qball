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
// CAP.h Set the complex absorbing potential for absorbing boundary condition
//
////////////////////////////////////////////////////////////////////////////////
#include <config.h>
#ifndef CAP_H
#define CAP_H
#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include <qball/Sample.h>
using namespace std;
class CAP : public Var{
    Sample *s;
    public:
    char const* name (void) const {return "cap";};
    int set (int argc, char ** argv){
         const string shape=argv[1];
         s->ctrl.cap_shape=shape;
         s->ctrl.cap_axis=atoi(argv[2]);
         s->ctrl.cap_start=atof(argv[3]);
         s->ctrl.cap_center=atof(argv[4]);
         if (argc>4){
             for(int i=5;i<argc;i++){
                 s->ctrl.cap_params.push_back(atof(argv[i]));
             }
         }
         const string error = "CAP should be used as set cap [type] [axis] [shape] [start] [center] [cap parameters..] (start and center of the cap should be in lattice unit)";
         if (argc<4){
             ui->error(error);
             return 1;
         }
         else if (shape!="linear" && shape!="power" && shape!="sin2"){
             ui->error(error);
             ui->error("CAP shape must be in linear, power or sin2, not "+shape);
             return 1;
         }
         else if (s->ctrl.cap_axis>2){
             ui->error(error);
             ui->error("axis must be 0, 1 or 2");
             return 1;
         }
         if (argv[1]=="linear" && s->ctrl.cap_params.size()!=1){
             ui->error(error);
             ui->error("if the shape is linear, only W_0 is needed. V_cap = -i*W_0*|r-R|/dR, for R<r<R+2dR");
             return 1;
         }
         else if (argv[1]=="sin2" && s->ctrl.cap_params.size()!=1){
             ui->error(error);
             ui->error("if the shape is sin2, only eta is needed. V_cap = -i*eta*sin^2((r-R)*pi/(2*dR)), for R<r<R+2dR");
             return 1;
         }
         else if  (argv[1]=="sin2" && s->ctrl.cap_params.size()!=1){
             ui->error(error);
             ui->error("if the shape is delta, only W is needed. V_cap = -i*W, for R<r<R+2dR");
             return 1;
         }
         s->ctrl.has_cap=true;
         return 0;
    }
    string print (void) const {
        ostringstream st;
        st.setf(ios::left,ios::adjustfield);
        st << setw(10) << name() << " = ";
        st.setf(ios::right,ios::adjustfield);
        st << setw(10) << s->ctrl.cap_shape;
        st << " " << setw(10) <<  s->ctrl.cap_axis << " " << setw(10) << s->ctrl.cap_start<<" "<<setw(10)<<s->ctrl.cap_center;
        return st.str();
    }
    CAP(Sample *sample): s(sample) {s->ctrl.has_cap=false;}
};

#endif
