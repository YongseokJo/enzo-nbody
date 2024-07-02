#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"


void generate_Matrix(double a[3], double (&A)[3][4]);
void direct_sum(double *x, double *v, double r2, double vx,
                double mass, double mdot, double (&a)[3], double (&adot)[3]);
//void getBlockTimeStep(double dt, int& TimeLevel, double &TimeStep);


// refer to ksint.f

void Binary::KSIntegration(double next_time, int &calnum){

    double binaryCalTime;

    calnum = 0;

    while ((CurrentTime + TimeStep) <= next_time) {

        // first predict the future position of the binary particle. 

        binaryCalTime = CurrentTime + TimeStep;
        calnum += 1;

        //fprintf(binout,"-----------------------------------------------------------------------------------------------\n");
        //fprintf(binout, "In KSRegularlizationIntegration.cpp KSIntegration, after position prediction for %dth time\n",calnum);
	//fflush(binout);


        predictBinary(binaryCalTime);

        // if there are zero neighbors for binary particles, calculate by unperturbed equations

        IntegrateBinary(binaryCalTime);
    }

    if ((r>2*r0)||(TimeStep>2*KSTime)) {
        isTerminate = true;
    }

//    fprintf(binout,"KSREGULARLIZATIONINTEGRATION");
//    fflush(binout);
//    fprintf(binout,"total number of Calculation = %d \n", calnum);
//    fflush(binout);
//    fprintf(binout,"binary time after calculation is %e \n", binaryCalTime); 
//    fflush(binout);

}



void Binary::predictBinary(double next_time) {

    double dt, dt2, dt3;
    double dtCM, dt2CM, dt3CM;
    double dtau, dtau2, dtau3, dtau4;
    double rinv, rinv2, rinv3, rinv4, rinv5;

    double R[Dim], Rdot[Dim];
    double Rinv;
    double ratioM;

    double L[3][4];

    Particle* ptclI;
    Particle* ptclJ;

    ptclI = ptclCM->BinaryParticleI;
    ptclJ = ptclCM->BinaryParticleJ;


    // time interval for prediction

    dtCM = (next_time - ptclCM->CurrentTimeIrr)*EnzoTimeStep;
    dt2CM = dtCM*dtCM;
    dt3CM = dt2CM*dtCM;

    dt = (next_time - CurrentTime)*EnzoTimeStep;
    dt2 = dt*dt;
    dt3 = dt2*dt;


    // inverse of r

    rinv = 1/r;
    rinv2 = rinv/r;
    rinv3 = rinv2/r;
    rinv4 = rinv3/r;
    rinv5 = rinv4/r;


    // the binary particle prediction consists of two parts
    // first, we need to predict the future position of the center of mass

    for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] = ptclCM->Position[dim] + ptclCM->Velocity[dim]*dtCM + ptclCM->a_tot[dim][0]*dt2CM/2 + ptclCM->a_tot[dim][1]*dt3CM/6;
		ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim] + ptclCM->a_tot[dim][0]*dtCM + ptclCM->a_tot[dim][1]*dt2CM/2;
	}


    // then we need to convey this information to the single particles as well
    // originally, the prediction order differs depending on conditions
    // but for now I would only use the 3rd order prediction

    // dtau = dtau0 - 

    dtau = dt*rinv - t2dot*dt2*rinv3/2 + t2dot*t2dot*dt3*rinv5/2 - t3dot*dt3*rinv4/6;

    // if (abs(dtau)>dTau) {
    //     dtau = 0.8*dTau;
    // }

    dtau2 = dtau*dtau;
    dtau3 = dtau2*dtau;
    dtau4 = dtau3*dtau;


    for (int dim=0; dim<4; dim++) {
        u_pred[dim]    = u[dim]    + udot[dim]*dtau  + u2dot[dim]*dtau2/2 + u3dot[dim]*dtau3/6 + u4dot[dim]*dtau4/24;
        udot_pred[dim] = udot[dim] + u2dot[dim]*dtau + u3dot[dim]*dtau2/2 + u4dot[dim]*dtau3/6;
    }

    
    // convert things back to global coordinates
    // r is the distance vector between the two binaries, and R is the distance between them
    // using the predicted R and predicted cm values, we obtain the position of each binary component at a future time


    R[0]   = u_pred[0]*u_pred[0] - u_pred[1]*u_pred[1] - u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3];
    R[1]   = 2*(u_pred[0]*u_pred[1] - u_pred[2]*u_pred[3]);
    R[2]   = 2*(u_pred[0]*u_pred[2] + u_pred[1]*u_pred[3]);
    ratioM = ptclJ->Mass/ptclCM->Mass;


    for (int dim=0; dim<Dim; dim++) {
        ptclI->PredPosition[dim] = ptclCM->PredPosition[dim] + ratioM*R[dim];
        ptclJ->PredPosition[dim] = ptclI->PredPosition[dim] - R[dim];
    }


    // do the same thing for velocity components


    generate_Matrix(u_pred,L);

    Rinv = 1/(u_pred[0]*u_pred[0] + u_pred[1]*u_pred[1] + u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3]) ;


    for (int dim=0; dim<Dim; dim++) {

        Rdot[dim] = 0.0;

        for (int dimu=0; dimu<4; dimu++) {
            Rdot[dim] += 2*L[dim][dimu]*udot_pred[dim]*Rinv;
        }
    }


    for (int dim=0; dim<Dim; dim++) {
        ptclI->PredVelocity[dim] = ptclCM->PredVelocity[dim] + ratioM*Rdot[dim];
        ptclJ->PredVelocity[dim] = ptclI->PredVelocity[dim] - Rdot[dim];
    }

    //fprintf(binout, "CM particle: x = %e, y = %e, z = %e", ptclCM->PredPosition[0], ptclCM->PredPosition[1], ptclCM->PredPosition[2]);
    //fprintf(binout, "CM particle: vx = %e, vy = %e, vz = %e", ptclCM->PredVelocity[0], ptclCM->PredVelocity[1], ptclCM->PredVelocity[2]);
    //fflush(binout);

    //fprintf(binout, "Ith particle: x = %e, y = %e, z = %e", ptclI->PredPosition[0], ptclI->PredPosition[1], ptclI->PredPosition[2]);
    //fprintf(binout, "Ith particle: vx = %e, vy = %e, vz = %e", ptclI->PredVelocity[0], ptclI->PredVelocity[1], ptclI->PredVelocity[2]);
    //fflush(binout);

    //fprintf(binout, "Jth particle: x = %e, y = %e, z = %e", ptclJ->PredPosition[0], ptclJ->PredPosition[1], ptclJ->PredPosition[2]);
    //fprintf(binout, "Jth particle: vx = %e, vy = %e, vz = %e", ptclJ->PredVelocity[0], ptclJ->PredVelocity[1], ptclJ->PredVelocity[2]);
    //fflush(binout);

}








void Binary::IntegrateBinary(double next_time) {

    // variables for position, velocity prediction

    double dt, dt2, dt3;
    double rinv, rinv2, rinv3, rinv4, rinv5;

    double R[Dim], Rdot[Dim];
    double Rinv;
    double ratioM;

    double dtau, dtau2, dtau3, dtau4, dtau5, dtau6;
    double dtau_temp;


    // variables for perturbing acceleration calculation

    double P[Dim], Pdot[Dim];  // perturbing acceleration and jerk
    double dxi[Dim],dvi[Dim],dxj[Dim],dvj[Dim];
    double dr2i, _dr3i, dr2j, _dr3j, dxdvi, dxdvj;


    // variables for calculating correction factors

    double dh0, dh;
    double h3dot_hermite, h4dot_hermite; 

    double Q_pred[4], Qdot_pred[4], Q2dot_pred[4];
    double L[3][4], Ldot[3][4], L2dot[3][4]; 

    double r_pred, rdot_pred, hdot_pred, h2dot_pred;
    double u2dot_pred[4], u3dot_pred[4];
    double u4dot_hermite[4], u5dot_hermite[4];

    double z;

    double mdot = 0.0;

    Particle* ptclI;
    Particle* ptclJ;

    // variables for updating time step
    
    int TimeLevelTmp;
    double TimeStepTmp;

    //fprintf(binout, "\nIn KSRegularlizationIntegration.cpp, IntegrateBinary\n");


    ptclI = ptclCM->BinaryParticleI;
    ptclJ = ptclCM->BinaryParticleJ;


    dtau = dTau;
    dtau2 = dtau*dtau;
    dtau3 = dtau2*dtau;
    dtau4 = dtau3*dtau;
    dtau5 = dtau4*dtau;
    dtau6 = dtau5*dtau;


    // refer to kspred.f

    // first predict the positions of neighbors

    // initialize the perturbing acceleration and jerk

    for (int dim=0; dim<Dim; dim++) {
        P[dim] = 0.0;
        Pdot[dim] = 0.0;
    }


    if (ptclCM->NumberOfAC >1) {

        for (Particle* ptcl: ptclCM->ACList) {
            ptcl->predictParticleSecondOrder(next_time);
        }

        // predict the positions of binary pair particles to highest order possible
        // using stumpff coefficients


        for (int dimu=0; dimu<4; dimu++) {
            u_pred[dimu]    = u[dimu] + udot[dimu]*dtau + u2dot[dimu]*dtau2/2 + u3dot[dimu]*dtau3/6 \
                            + cn[4]*u4dot[dimu]*dtau4/24 + cn[5]*u5dot[dimu]*dtau5/120;
            udot_pred[dimu] = udot[dimu] + u2dot[dimu]*dtau + u3dot[dimu]*dtau2/2 \
                            + cn[4]*u4dot[dimu]*dtau3/6+ cn[5]*u5dot[dimu]*dtau4/24;
        }


        R[0]   = u_pred[0]*u_pred[0] - u_pred[1]*u_pred[1] - u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3];
        R[1]   = 2*(u_pred[0]*u_pred[1] - u_pred[2]*u_pred[3]);
        R[2]   = 2*(u_pred[0]*u_pred[2] + u_pred[1]*u_pred[3]);
        ratioM = ptclJ->Mass/ptclCM->Mass;

        // this was already done in predictBinary routine!
        // for (int dim=0; dim<Dim; dim++) {
        //     ptclI->PredPosition[dim] = PredPosition[dim] + ratioM*R[dim];
        //     ptclJ->PredPosition[dim] = ptclI->PredPosition[dim] - R[dim];
        // }

        generate_Matrix(u_pred,L);

        r_pred = (u_pred[0]*u_pred[0] + u_pred[1]*u_pred[1] + u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3]) ;


        for (int dim=0; dim<Dim; dim++) {

            Rdot[dim] = 0.0;

            for (int dimu=0; dimu<4; dimu++) {
                Rdot[dim] += 2*L[dim][dimu]*udot_pred[dim]/r_pred;
            }
        }

        // for (int dim=0; dim<Dim; dim++) {
        //     ptclI->PredVelocity[dim] = PredVelocity[dim] + ratioM*Rdot[dim];
        //     ptclJ->PredVelocity[dim] = ptclI->PredVelocity[dim] - Rdot[dim];
        // } // -> already done at predictBinary


        // calculate perturbation from single particles

        for (Particle* ptcl: ptclCM->ACList) {

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

        
        for (int dim=0; dim<Dim; dim++) {
            Pdot[dim] *= r_pred;
        }

    }

    // calculating perturbation from particle systems
    // add later on

    gamma = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2])*r_pred*r_pred/ptclCM->Mass;



    // refer to kscorr.f

    dh = hdot*dtau + h2dot*dtau2/2 + h3dot*dtau3/6 + h4dot*dtau4/24;
    dh0 = hdot*dtau + h2dot*dtau2/2;


    // caculation of lower order derivatives

    // initialize the variables

    for(int dimu=0; dimu<4; dimu++) {
        Q_pred[dimu] = 0.0;
        Qdot_pred[dimu] = 0.0;
        u2dot_pred[dimu] = 0.0;
        u3dot_pred[dimu] = 0.0;
        u4dot_hermite[dimu] = 0.0;
        u5dot_hermite[dimu] = 0.0;
    }

    rdot_pred = 0.0; 
    hdot_pred = 0.0;
    h2dot_pred = 0.0;

    generate_Matrix(u_pred,L);
    generate_Matrix(udot_pred,Ldot);

     
    
    for (int dimu=0; dimu<4; dimu++) {

        Q_pred[dimu]     = L[0][dimu]*P[0] + L[1][dimu]*P[1] + L[2][dimu]*P[2];
        Qdot_pred[dimu]  = Ldot[0][dimu]*P[0] + Ldot[1][dimu]*P[1] + Ldot[2][dimu]*P[2] \
                         + L[0][dimu]*Pdot[0] + L[1][dimu]*Pdot[1] + L[2][dimu]*Pdot[2];

        u2dot_pred[dimu] = 0.5*(h+dh)*u_pred[dimu] + 0.5*r_pred*Q_pred[dimu];

        rdot_pred       += 2*u_pred[dimu]*udot_pred[dimu];
        hdot_pred       += 2*udot_pred[dimu]*Q_pred[dimu];
        h2dot_pred      += 2*(u2dot_pred[dimu]*Q_pred[dimu] + udot_pred[dimu]*Qdot_pred[dimu]);
        
    }

    for (int dimu = 0; dimu<4; dimu++) {
        u3dot_pred[dimu] = 0.5*((h+dh)*udot_pred[dimu] + hdot_pred*u_pred[dimu]) \
                         + 0.5*(rdot_pred*Q_pred[dimu] + r_pred*Qdot_pred[dimu]);
    }


    // apply the hermite corrector method to calculate 4th and 5th derivatives of u
    // and also the 3rd and 4th derivatives of h

    for (int dimu = 0; dimu<4; dimu++) {
        u4dot_hermite[dimu] = - 6*(u2dot[dimu] - u2dot_pred[dimu])/dtau2 \
                              - 2*(2*u3dot[dimu] + u3dot_pred[dimu])/dtau;
        u5dot_hermite[dimu] =   12*(u2dot[dimu] - u2dot_pred[dimu])/dtau3 \
                              + 6*(u3dot[dimu] + u3dot_pred[dimu])/dtau2;
    }

    h3dot_hermite = -6*(hdot-hdot_pred)/dtau2 -2*(2*h2dot + h2dot_pred)/dtau;
    h4dot_hermite = 12*(hdot-hdot_pred)/dtau3 + 6*(h2dot + h2dot_pred)/dtau2;
        // for (int dim=0; dim<Dim; dim++) {
        //     ptclI->PredVelocity[dim] = PredVelocity[dim] + ratioM*Rdot[dim];
        //     ptclJ->PredVelocity[dim] = ptclI->PredVelocity[dim] - Rdot[dim];
        // }

    // based on the higher-order derivatives, correct the values of u, udot and r

    for (int dimu=0; dimu<4; dimu++) {
        u_pred[dimu]    = u[dimu] + udot[dimu]*dtau + u2dot[dimu]*dtau2/2 + u3dot[dimu]*dtau3/6 \
                        + cn[4]*u4dot_hermite[dimu]*dtau4/24 + cn[5]*u5dot_hermite[dimu]*dtau5/120;
        udot_pred[dimu] = udot[dimu] + u2dot[dimu]*dtau + u3dot[dimu]*dtau2/2 \
                        + cn[4]*u4dot_hermite[dimu]*dtau3/6+ cn[5]*u5dot_hermite[dimu]*dtau4/24;           
    }

    r_pred = (u_pred[0]*u_pred[0] + u_pred[1]*u_pred[1] + u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3]);


    // intialize and save the updated values of all relevant derivatives

    r = r_pred;
    rdot = 0.0;
    r2dot = 0.0;
    r3dot = 0.0;
    r4dot = 0.0;
    r5dot = 0.0;

    h = h + dh0 + h3dot_hermite*dtau3/6 + h4dot_hermite*dtau4/24;
    hdot = hdot_pred;
    h2dot = h2dot_pred;
    h3dot = h3dot_hermite + h4dot_hermite*dtau;
    h4dot = h4dot_hermite;

    for (int dimu=0; dimu<4 ;dimu++) {

        Q[dimu] = Q_pred[dimu];
        Qdot[dimu] = Qdot_pred[dimu];

        u[dimu] = u_pred[dimu];
        udot[dimu] = udot_pred[dimu];
        u2dot[dimu] = 0.5*h*u_pred[dimu] + 0.5*r_pred*Q_pred[dimu];
        u3dot[dimu] = 0.5*(h*udot_pred[dimu] + hdot_pred*u_pred[dimu]) \
                    + 0.5*(rdot_pred*Q_pred[dimu] + r_pred*Qdot_pred[dimu]);
        u4dot[dimu] = u4dot_hermite[dimu] + u5dot_hermite[dimu]*dTau;
        u5dot[dimu] = u5dot_hermite[dimu];
    }


    for (int dimu=0; dimu<4; dimu++) {
        rdot      += 2*u[dimu]*udot[dimu];
        r2dot     += 2*(udot[dimu]*udot[dimu] + u[dimu]*u2dot[dimu]);
        r3dot     += 2*(3*u2dot[dimu]*udot[dimu] + u3dot[dimu]*u[dimu]);
        r4dot     += 2*(u4dot[dimu]*u[dimu] + 4*u3dot[dimu]*udot[dimu] + 3*u2dot[dimu]*u2dot[dimu]);
        r5dot     += 2*(u5dot[dimu]*u[dimu] + 5*u4dot[dimu]*udot[dimu] + 10*u3dot[dimu]*u2dot[dimu]);
    }

    tdot = r;
    t2dot = rdot;
    t3dot = r2dot;
    t4dot = r3dot;
    t5dot = r4dot;
    t6dot = r5dot;


    // update the time


    // generate new stumpff coefficients
    // refer to stumpff.for

    z = -0.5*h*dtau;
    getStumpffCoefficients(z);

    // TimeStep = tdot*dtau + t2dot*dtau2/2 + t3dot*dtau3/6 + t4dot*dtau4/24 \
    //          + t5dot*cn_4z[5]*dtau5/120 + t6dot*cn_4z[6]*dtau6/720;

    // generate new list of perturbers

    CurrentTime = next_time;
    CurrentTau += dtau;

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

//    TimeLevel = TimeLevelTmp;
    dTau = dtau;


    //fprintf(binout, "\nPerturbing term Q - Q0:%e, Q1:%e, Q2:%e, Q3: %e \n", Q_pred[0], Q_pred[1], Q_pred[2], Q_pred[3]);
    //fprintf(binout, "Perturbing term Qdot - Qdot0:%e, Qdot1:%e, Qdot2:%e, Qdot3: %e \n", Qdot_pred[0], Qdot_pred[1], Qdot_pred[2], Qdot_pred[3]);
    //fflush(binout);

    //fprintf(binout, "derivatives of r: r0 = %e r1 = %e r2 = %e r3 = %e r4 = %e r5 = %e \n", r, rdot, r2dot, r3dot, r4dot, r5dot);
    //fprintf(binout, "derivatives of h: h0 = %e h1 = %e h2 = %e h3 = %e h4 = %e \n", h, hdot, h2dot, h3dot, h4dot);
    //fflush(binout);
    
    //fprintf(binout, "derivatives of u[0]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[0], udot[0], u2dot[0], u3dot[0], u4dot[0], u5dot[0]);
    //fprintf(binout, "derivatives of u[1]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[1], udot[1], u2dot[1], u3dot[1], u4dot[1], u5dot[1]);
    //fprintf(binout, "derivatives of u[2]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[2], udot[2], u2dot[2], u3dot[2], u4dot[2], u5dot[2]);
    //fprintf(binout, "derivatives of u[3]: u0 = %e u1 = %e u2 = %e u3 = %e u4 = %e u5 = %e \n", u[3], udot[3], u2dot[3], u3dot[3], u4dot[3], u5dot[3]);
    //fflush(binout);

    //fprintf(binout,"Stumpff Coefficinets : c1 = %e c2 = %e c3 = %e c4 = %e c5 = %e \n \n ", cn[1], cn[2], cn[3], cn[4], cn[5]);
    //fprintf(binout,"dTau = %e, TimeStep = %e, gamma = %e\n\n", dTau, TimeStep, gamma);
    //fflush(binout);
    

}
