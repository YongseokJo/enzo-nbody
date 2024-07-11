#include "../global.h"
#include "../defs.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

//void getBlockTimeStep(double dt, int& TimeLevel, double &TimeStep);
void direct_sum(double *x, double *v, double r2, double vx,
                double mass, double mdot, double (&a)[3], double (&adot)[3]);


void generate_Matrix(double a[3], double (&A)[3][4]) {
 
    A[0][0] =  a[0];
    A[0][1] = -a[1];
    A[0][2] =  a[2];
    A[0][3] =  a[3];

    A[1][0] =  a[1];
    A[1][1] =  a[0];
    A[1][2] = -a[3];
    A[1][3] = -a[2];

    A[2][0] =  a[2];
    A[2][1] =  a[3];
    A[2][2] =  a[0];
    A[2][3] =  a[1];

}


void Binary::getStumpffCoefficients(double z){


    double z4; // 4*z, used for calculation of time polynomial approximations 

    
    // note that the coefficient converges to 1
    // so for higher orders, we can assume that cn ~ 1

    // use the recursion relation to obtain the coefficients
    // the recursion relation is c(n) = 1 - z*c(n+2)/((n+1)(n+2))
    // also note that the index starts from 0. 


    cn[stumpffN]   = 1-z*1/((stumpffN+1)*(stumpffN+2));
    cn[stumpffN-1] = 1-z*1/((stumpffN)*(stumpffN+1));

    for (int i=(stumpffN-2); i>2; i--) {
        cn[i] = 1-z*cn[i+2]/((i+1)*(i+2));
    }

    cn[2] = 1.0; // is this right? needs checking
    cn[1] = 1.0;
    cn[0] = 1.0;
    

    // also obtain the corresponding coefficients for 4z

    z4 = 4*z;

    cn_4z[stumpffN]   = 1-z*1/((stumpffN+1)*(stumpffN+2));
    cn_4z[stumpffN-1] = 1-z*1/((stumpffN)*(stumpffN+1));

    for (int i=(stumpffN-2); i>2; i--) {
        cn_4z[i] = 1-z*cn[i+2]/((i+1)*(i+2));
    }

    cn_4z[2] = 1.0; // is this right? needs checking
    cn_4z[1] = 1.0;
    cn_4z[0] = 1.0;


}




void Binary::InitializeBinary(double current_time) {

    
    double x[Dim],xdot[Dim];
    double r2;

    double dxi[Dim],dvi[Dim],dxj[Dim],dvj[Dim];
    double dr2i, _dr3i, dr2j, _dr3j, dxdvi, dxdvj;

    double Q2dot[4];
    double L[3][4], Ldot[3][4], L2dot[3][4];  // Levi-civita transformation matrix

    double P[Dim], Pdot[Dim]; // perturbed acceleration in cartesian coordinates
    double Ptot; // total perturbed force

    int TimeLevelTmp;
    double TimeStepTmp;


    // variables related with time step and stumpff coefficients

    double dtau_temp, dtau; // time in u-coordinate system
    double dtau2, dtau3, dtau4, dtau5, dtau6;
    double dt_ks; // ks time in physical units
    double z;

    double dt, dt2, dt3;
    double rinv, rinv2, rinv3, rinv4, rinv5;

    double mdot = 0.0;

    Particle* ptclI;
    Particle* ptclJ;

    
    ptclI = ptclCM->BinaryParticleI;
    ptclJ = ptclCM->BinaryParticleJ;

    // just in case
    ptclI->predictParticleSecondOrderIrr(current_time);
    ptclJ->predictParticleSecondOrderIrr(current_time);


    r2 = 0;
    
    // define the relative coordinates for binaries

    for (int dim=0; dim<Dim; dim++) {
        x[dim] = ptclI->PredPosition[dim] - ptclJ->PredPosition[dim];
        xdot[dim] = ptclI->PredVelocity[dim] - ptclJ->PredVelocity[dim];
        r2 += x[dim]*x[dim];
    }

    r = sqrt(r2);


    // TRANSFORMATION OF COORDINATES

    // calculate the initial values for new coordinate system

    if (x[0]<0) {

        u[2] = 0.0;
        u[1] = sqrt(0.5*(r-x[0]));
        u[0] = 0.5*x[1]/u[1];
        u[3] = 0.5*x[2]/u[1];

    } else {

        u[3] = 0.0;
        u[0] = sqrt(0.5*(r+x[0]));
        u[1] = 0.5*x[1]/u[0];
        u[2] = 0.5*x[2]/u[0];

    }

    fprintf(binout, "\n >>>>>>>>>>>>from BinaryRoutines.cpp<<<<<<<<<< \n");

    fprintf(binout, "KS coordinates - u1:%e, u2:%e, u3:%e, u4:%e \n", u[0],u[1],u[2],u[3]);
    //fflush(binout);

    // form the Levi-civita transformation matrix based on initial values

    generate_Matrix(u,L);


    // calculate the velocity in new coordinates using Levi-civita matrix
    // and also the binding energy per unit reduced mass, represented as h

    h = 0.0; //binding energy per unit reduced mass

    for (int dim=0; dim<4; dim++) {

        udot[dim] = 0.5*(L[0][dim]*xdot[0] + L[1][dim]*xdot[1] + L[2][dim]*xdot[2]);
        u_pred[dim] = u[dim];
        udot_pred[dim] = udot[dim];

        h += 2*udot[dim]*udot[dim];
    }

    h -= ptclCM->Mass/r;
    a = -0.5*ptclCM->Mass/h;


    // form the list of perturbers : refer to kslist.f
    // save the list of perturbers in the neighbor list of ith particle. search inside COM particle list

    // change the step of the 2nd particle so it won't be calculated in the list.
    // need to add exclusion of particles CHECK!!

    // Initialize the perturbing force

    for (int dim=0; dim<Dim; dim++) {
        P[dim] = 0.0;
        Pdot[dim] = 0.0;
    }


    // CALCULATION OF PERTURBER FORCE AND FORCE POLYNOMIALS

    fprintf(binout, "# of binary neighbors = %d\n",ptclCM->NumberOfAC);

    // predict the perturber position and calculate the perturbing force

    for (Particle* ptcl: ptclCM->ACList) {

        ptcl->predictParticleSecondOrderIrr(current_time);

        dr2i = 0.0;
        dr2j = 0.0;

        dxdvi = 0.0;
        dxdvj = 0.0;

        for (int dim=0; dim<Dim; dim++) {

            dxi[dim] = ptcl->PredPosition[dim] - ptclI->PredPosition[dim];
            dvi[dim] = ptcl->PredVelocity[dim] - ptclI->PredVelocity[dim];
            dr2i += dxi[dim]*dxi[dim];
            dxdvi += dxi[dim]*dvi[dim];

            dxj[dim] = ptcl->PredPosition[dim] - ptclJ->PredPosition[dim];
            dvj[dim] = ptcl->PredVelocity[dim] - ptclJ->PredVelocity[dim];
            dr2j += dxj[dim]*dxj[dim];
            dxdvj += dxj[dim]*dvj[dim];

        }

        _dr3i = 1/(dr2i)/sqrt(dr2i);
        _dr3j = 1/(dr2j)/sqrt(dr2j);

	//fprintf(binout, "dxi = %e %e %e  dvi = %e %e %e \n",dxi[0],dxi[1],dxi[2],dvi[0],dvi[1],dvi[2]);
    	//fprintf(binout, "dxj = %e %e %e  dvj = %e %e %e \n \n",dxj[0],dxj[1],dxj[2],dvj[0],dvj[1],dvj[2]);


        for (int dim=0; dim<Dim; dim++) {

            // ith particle
            P[dim]    += ptcl->Mass*_dr3i*dxi[dim];
            Pdot[dim] += ptcl->Mass*_dr3i*(dvi[dim] - 3*dxi[dim]*dxdvi/dr2i);//+mdot*_dr3i*dxi[dim];

            // jth partilce
            P[dim]    -= ptcl->Mass*_dr3j*dxj[dim];
            Pdot[dim] -= ptcl->Mass*_dr3j*(dvj[dim] - 3*dxj[dim]*dxdvj/dr2j);//+mdot*_dr3j*dxj[dim];
        }

        //direct_sum(dxi ,dvi, dr2i, dxdvi, ptcl->Mass, mdot, P, Pdot);
        //direct_sum(dxj ,dvj, dr2j, dxdvj, -ptcl->Mass, mdot, P, Pdot);
    }

    // multiply the seperation r to dadot using relation t' = R
    // because da and dadot should be converted to KS coordinates

    for (int dim=0; dim<Dim; dim++) {
        Pdot[dim] *= r;
    }

    Ptot = sqrt(P[0]*P[0] + P[1]*P[1] + P[2]*P[2]);
    fprintf(binout, "Perturbing force - Px:%e, Py:%e, Pz:%e \n", P[0],P[1],P[2]);

    gamma = Ptot*r*r/ptclCM->Mass;


    // scale the perturbing force by modification factor...?
    // CHECK


    // calculate derivatives for KS motion (refer to kspoly.f)
    // relative variables are already initialized on



    // caculation of lower order derivatives

    generate_Matrix(u,L);
    generate_Matrix(udot,Ldot);
    
    for (int dim=0; dim<4; dim++) {

        Q[dim]     = L[0][dim]*P[0] + L[1][dim]*P[1] + L[2][dim]*P[2];
        Qdot[dim]  = Ldot[0][dim]*P[0] + Ldot[1][dim]*P[1] + Ldot[2][dim]*P[2] \
                   + L[0][dim]*Pdot[0] + L[1][dim]*Pdot[1] + L[2][dim]*Pdot[2];

        u2dot[dim] = 0.5*h*u[dim] + 0.5*r*Q[dim];

        rdot      += 2*u[dim]*udot[dim];
        r2dot     += 2*(udot[dim]*udot[dim] + u[dim]*u2dot[dim]);

        hdot      += 2*udot[dim]*Q[dim];
        h2dot     += 2*(u2dot[dim]*Q[dim] + udot[dim]*Qdot[dim]);
        
    }


    // calculation of more higher derivatives
    // use the lower order derivation calculated above for calculation


    generate_Matrix(u2dot,L2dot);

    for (int dim=0; dim<4; dim++) {

        Q2dot[dim] =     (L2dot[0][dim]*P[0]   + L2dot[1][dim]*P[1]   + L2dot[2][dim]*P[2]) \
                   + 2.0*(Ldot[0][dim]*Pdot[0] + Ldot[1][dim]*Pdot[1] + Ldot[2][dim]*Pdot[2]);
                // + 2.0*(L[0][dim]*P2dot[0]   + L[1][dim]*P2dot[1]   + L[2][dim]*P2dot[2])

        u3dot[dim] = 0.5*(h*udot[dim] + hdot*u[dim]) + 0.5*(rdot*Q[dim] + r*Qdot[dim]);

        r3dot     += 2*(3*u2dot[dim]*udot[dim] + u3dot[dim]*u[dim]);

        h3dot     += 2*(u3dot[dim]*Q[dim] + 2*u2dot[dim]*Qdot[dim]+ 2*udot[dim]*Q2dot[dim]);
    
    }

    // calculate the highest derivatives


    for (int dim=0; dim<4; dim++) {

        u4dot[dim] = 0.5*(h2dot*u[dim] + hdot*udot[dim] + h*u2dot[dim] \
                        + r2dot*Q[dim] + rdot*Qdot[dim] + r*Q2dot[dim]);
        u5dot[dim] = 0.5*(h3dot*u[dim] + 3*h2dot*udot[dim] + 3*hdot*u2dot[dim] + h*u3dot[dim] \
                   +      r3dot*Q[dim] + 3*r2dot*Qdot[dim] + 3*rdot*Q2dot[dim]);

        h4dot     += 2*(u4dot[dim]*Q[dim] + 3*u3dot[dim]*Qdot[dim] + 3*u2dot[dim]*Q2dot[dim]);

        r4dot     += 2*(u4dot[dim]*u[dim] + 4*u3dot[dim]*udot[dim] + 3*u2dot[dim]*u2dot[dim]);
        r5dot     += 2*(u5dot[dim]*u[dim] + 5*u4dot[dim]*udot[dim] + 10*u3dot[dim]*u2dot[dim]);

    }


    // change the radius derivatives to time derivatives
    // do this for convinient expression of relative variables

    tdot = r;
    t2dot = rdot;
    t3dot = r2dot;
    t4dot = r3dot;
    t5dot = r4dot;
    t6dot = r5dot;


    // obtain the apropriate time step for binary

    dtau_temp = std::min(r/ptclCM->Mass,0.5*abs(h));
    dtau = 0.8*eta*sqrt(dtau_temp)/pow((1 + 1000.0 * gamma), 1.0/3.0);

    dtau2 = dtau*dtau;
    dtau3 = dtau2*dtau;
    dtau4 = dtau3*dtau;
    dtau5 = dtau4*dtau;
    dtau6 = dtau5*dtau;


    // obtain the Stumpff coefficients

    z = -0.5*h*dtau;
    getStumpffCoefficients(z);


    // obtain the physical interval corresponding to dtau using stumpff coefficients
    // note that r = t', r'=t'' and so on...

    TimeStep = (tdot*dtau + t2dot*dtau2/2 + t3dot*dtau3/6 + t4dot*dtau4/24 \
             + t5dot*cn_4z[5]*dtau5/120 + t6dot*cn_4z[6]*dtau6/720)/EnzoTimeStep;

    // also, convert the time step into block steps. 

//    getBlockTimeStep(TimeStep, TimeLevelTmp, TimeStepTmp);

 //   TimeLevel = TimeLevelTmp;
    dTau = dtau;

    // save the initial seperation
    r0 = r;

    fprintf(binout, "Perturbing term Q - Q0:%e, Q1:%e, Q2:%e, Q3: %e \n", Q[0], Q[1], Q[2], Q[3]);
    fprintf(binout, "Perturbing term Qdot - Qdot0:%e, Qdot1:%e, Qdot2:%e, Qdot3: %e \n", Qdot[0], Qdot[1], Qdot[2], Qdot[3]);
    //fflush(binout);

    fprintf(binout, "derivatives of r: r0 = %e r1 = %e r2 = %e r3 = %e r4 = %e r5 = %e \n", r, rdot, r2dot, r3dot, r4dot, r5dot);
    fprintf(binout, "derivatives of h: h0 = %e h1 = %e h2 = %e h3 = %e h4 = %e \n", h, hdot, h2dot, h3dot, h4dot);
    //fflush(binout);
    
    fprintf(binout, "derivatives of u[0]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[0], udot[0], u2dot[0], u3dot[0], u4dot[0], u5dot[0]);
    fprintf(binout, "derivatives of u[1]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[1], udot[1], u2dot[1], u3dot[1], u4dot[1], u5dot[1]);
    fprintf(binout, "derivatives of u[2]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[2], udot[2], u2dot[2], u3dot[2], u4dot[2], u5dot[2]);
    fprintf(binout, "derivatives of u[3]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[3], udot[3], u2dot[3], u3dot[3], u4dot[3], u5dot[3]);
    //fflush(binout);

    fprintf(binout,"Stumpff Coefficinets : c1 = %e c2 = %e c3 = %e c4 = %e c5 = %e \n \n ", cn[1], cn[2], cn[3], cn[4], cn[5]);
    //fflush(binout);

    fprintf(binout, "\n >>>>>>>>>>>>from BinaryRoutines.cpp<<<<<<<<<< \n");
    //fflush(binout);


}
