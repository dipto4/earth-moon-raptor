/******************************************
 * RAPTOR-solar-system
 * Apr 2, 2023
 * Author: Dipto Mukherjee
 * NOTE:
 * (1) pos, vel refer to positions and velocities in cartesian coordinates
 * but xyz, vxyz refers to positions and velocites in jacobi coordinates
 * (2) N: number of bodies including the test object
 *
 *
 * *****************************************/

//#include "universal.h"
#include "wh_kepler.h"
#define VSCALE 6.2830666414874994
//#define ETA_SYM 0.00625
#define ETA_SYM 0.025
#define ETA 0.001
#define NMAX 20
#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


/*void evolve_kepler_2body(double mass[2], double pos[2][3], double vel[2][3], double dt) {
    double xb[3] = {pos[0][0], pos[0][1], pos[0][2]}; 
    double xc[3] = {pos[1][0], pos[1][1], pos[1][2]}; 
    
    double vxb[3] = {vel[0][0], vel[0][1], vel[0][2]}; 
    double vxc[3] = {vel[1][0], vel[1][1], vel[1][2]}; 
    
    double dpos[3],dpos0[3];
    double pos_cm[3];
    double dvel[3],dvel0[3];
    double vel_cm[3];
    double m1 = mass[0];
    double m2 = mass[1];
    double mtot = (m1+m2);
    double f1 = m2 / mtot;
    double f2 = m1 / mtot;

    for(int d=0; d<3; d++){
        dpos0[d] = xb[d]-xc[d];
        dvel0[d] = vxb[d]-vxc[d];
    }

    for(int d=0; d<3; d++){
        pos_cm[d] = (m1*xb[d] + m2*xc[d]) / mtot;
        vel_cm[d] = (m1*vxb[d] + m2*vxc[d]) / mtot;
    }
    // evolve center of mass for dt

    // call the Kepler solver
    Binary state0 = {dpos0[0],dpos0[1],dpos0[2], dvel0[0],dvel0[1],dvel0[2], 0, 0, 0};
    Binary state1;

    kepler_step(mtot, dt, &state0, &state1);

    // translate coordinates from 2-body frame to original frame
    xb[0] = pos_cm[0]+f1*state1.x;
    xb[1] = pos_cm[1]+f1*state1.y;
    xb[2] = pos_cm[2]+f1*state1.z;
    vxb[0] = vel_cm[0]+f1*state1.xd;
    vxb[1] = vel_cm[1]+f1*state1.yd;
    vxb[2] = vel_cm[2]+f1*state1.zd;
    xc[0] = pos_cm[0]-f2*state1.x;
    xc[1] = pos_cm[1]-f2*state1.y;
    xc[2] = pos_cm[2]-f2*state1.z;
    vxc[0] = vel_cm[0]-f2*state1.xd;
    vxc[1] = vel_cm[1]-f2*state1.yd;
    vxc[2] = vel_cm[2]-f2*state1.zd;
    
    for(int d=0;d<3;d++) {
        pos[0][d] = xb[d];
        pos[1][d] = xc[d];
        vel[0][d]= vxb[d];
        vel[1][d]= vxc[d];
    }


}*/

void debug_print(double xyz[NMAX][3],double vxyz[NMAX][3],int N) {
    for(size_t i=8;i<9;i++) {
        for(size_t d=0;d<3;d++)
            printf("%21.16f ",xyz[i][d]);
        for(size_t d=0;d<3;d++)
            printf("%21.16f ",vxyz[i][d]);
        printf("\n");
    }
}

void cart2jacobi(double xyz[NMAX][3], double vxyz[NMAX][3], double mass[NMAX],int N,
        double x_star[NMAX][3], double v_star[NMAX][3]) {
    double eta[NMAX];
    eta[0] = mass[0];
    double rcm[3];
    double vcm[3];
    for(size_t i=1; i<N;i++) {
        eta[i] = eta[i-1] + mass[i];
    }
    
    for(size_t d=0;d<3;d++) {
        rcm[d] = mass[0] * xyz[0][d];
        vcm[d] = mass[0] * vxyz[0][d];
    }
    for(size_t i=1;i<N;i++) {
        double xi_star[3];
        double vi_star[3];

        for(size_t d=0;d<3;d++) {
            xi_star[d] = xyz[i][d] - rcm[d]     / eta[i-1];
            vi_star[d] = vxyz[i][d] - vcm[d]    / eta[i-1];
        }

        for(size_t d=0;d<3;d++) {
            rcm[d] = rcm[d]*(1+mass[i]/eta[i-1]) + mass[i] * xi_star[d];
            vcm[d] = vcm[d]*(1+mass[i]/eta[i-1]) + mass[i] * vi_star[d];
        }

        for(size_t d=0;d<3;d++) {
            x_star[i][d] = xi_star[d];
            v_star[i][d] = vi_star[d];
        }
     }
    for(size_t d=0;d<3;d++) {
        x_star[0][d] = rcm[d] / eta[N-1];
        v_star[0][d] = vcm[d] / eta[N-1];
    }


}
/*taken from ABIE */
void helio2jacobi_nonABIE(double mass[NMAX], double xyz[NMAX][3], double vxyz[NMAX][3], int N,
        double eta[NMAX], double x_star[NMAX][3], double v_star[NMAX][3]) {
    eta[0] = mass[0];
    double rcm[3];
    double vcm[3];
    for(size_t i=1; i<N;i++) {
        eta[i] = eta[i-1] + mass[i];
    }
    
    for(size_t d=0;d<3;d++) {
        x_star[0][d] = 0.0;
        v_star[0][d] = 0.0;
        x_star[1][d] = xyz[1][d];
        v_star[1][d] = vxyz[1][d];
    }

    for(size_t d=0;d<3;d++) {
        rcm[d] = mass[1] * xyz[1][d];
        vcm[d] = mass[1] * vxyz[1][d];
    }
    for(size_t i=2;i<N;i++) {
        double xi_star[3];
        double vi_star[3];

        for(size_t d=0;d<3;d++) {
            xi_star[d] = xyz[i][d] - rcm[d]     / eta[i-1];
            vi_star[d] = vxyz[i][d] - vcm[d]    / eta[i-1];
        }

        for(size_t d=0;d<3;d++) {
            rcm[d] = rcm[d]*(1+mass[i]/eta[i-1]) + mass[i] * xi_star[d];
            vcm[d] = vcm[d]*(1+mass[i]/eta[i-1]) + mass[i] * vi_star[d];
        }

        for(size_t d=0;d<3;d++) {
            x_star[i][d] = xi_star[d];
            v_star[i][d] = vi_star[d];
        }
     }
    /*for(size_t d=0;d<3;d++) {
        x_star[0][d] = rcm[d] / eta[N-1];
        v_star[0][d] = vcm[d] / eta[N-1];
    }*/


}

void jacobi2cart(double x_star[NMAX][3], double v_star[NMAX][3], double mass[NMAX],int N,
        double xyz[NMAX][3], double vxyz[NMAX][3]) {

    double eta[NMAX];
    eta[0] = mass[0];
    double rcm[3];
    double vcm[3];
    for(size_t i=1; i<N;i++) {
        eta[i] = eta[i-1] + mass[i];
    }
    

    for(size_t d=0;d<3;d++) {
        rcm[d] = eta[N-1] * x_star[0][d];
        vcm[d] = eta[N-1] * v_star[0][d];
    }

    for(size_t i=N-1;i>0;i--) {
        for(size_t d=0;d<3;d++) {
            rcm[d] = (rcm[d] - mass[i]*x_star[i][d]) / eta[i];
            vcm[d] = (vcm[d] - mass[i]*v_star[i][d]) / eta[i];
        }

        for(size_t d=0;d<3;d++) {
            xyz[i][d] =     x_star[i][d] + rcm[d];
            vxyz[i][d] =    v_star[i][d] + vcm[d];
        }

        for(size_t d=0;d<3;d++) {
            rcm[d] = eta[i-1] * rcm[d];
            vcm[d] = eta[i-1] * vcm[d];
        }
    }
    for(size_t d=0;d<3;d++) {
        xyz[0][d] = rcm[d] / mass[0];
        vxyz[0][d] = vcm[d] / mass[0];
    }
}

void jacobi2helio_nonABIE(double eta[NMAX], double x_star[NMAX][3], double v_star[NMAX][3], int N,
        double mass[NMAX], double xyz[NMAX][3], double vxyz[NMAX][3]) {

    double rcm[3];
    double vcm[3];

    for(size_t d=0;d<3;d++) {
        rcm[d] = eta[N-1] * x_star[0][d];
        vcm[d] = eta[N-1] * v_star[0][d];
    }

    for(size_t i=N-1;i>0;i--) {
        for(size_t d=0;d<3;d++) {
            rcm[d] = (rcm[d] - mass[i]*x_star[i][d]) / eta[i];
            vcm[d] = (vcm[d] - mass[i]*v_star[i][d]) / eta[i];
        }

        for(size_t d=0;d<3;d++) {
            xyz[i][d] =     x_star[i][d] + rcm[d];
            vxyz[i][d] =    v_star[i][d] + vcm[d];
        }

        for(size_t d=0;d<3;d++) {
            rcm[d] = eta[i-1] * rcm[d];
            vcm[d] = eta[i-1] * vcm[d];
        }
    }
    for(size_t d=0;d<3;d++) {
        xyz[0][d] = rcm[d] / mass[0];
        vxyz[0][d] = vcm[d] / mass[0];
    }
}

void helio2jacobi(double xyz[NMAX][3], double vxyz[NMAX][3], double mass[NMAX], int N, double x_star[NMAX][3], double v_star[NMAX][3]) {
//void helio2jacobi(double *pos, double *vel, double *masses, size_t nbodies, double *jacobi_pos, double *jacobi_vel) {


    // Compute eta (interior masses):
    double eta[N];
    eta[0] = mass[0];
    for (size_t i = 1; i < N; i++) {
        eta[i] = mass[i] + eta[i - 1];
    }

    // Assume central body at rest:
    x_star[0][0] = 0.0; x_star[0][1] = 0.0; x_star[0][2] = 0.0;
    v_star[0][0] = 0.0; v_star[0][1] = 0.0; v_star[0][2] = 0.0;

    // the jacobi coordinates for the second body is the same as the heliocentric
    x_star[1][0] = xyz[1][0]; x_star[1][1] = xyz[1][1]; x_star[1][2] = xyz[1][2];
    v_star[1][0] = vxyz[1][0]; v_star[1][1] = vxyz[1][1]; v_star[1][2] = vxyz[1][2];

    // Jacobi coordinates of first body coincide with heliocentric, leave as they are.

    // Compute internal c.o.m. and momentum:
    double auxR[3], auxV[3], Ri[3], Vi[3];
    auxR[0] = mass[1] * xyz[1][0];
    auxR[1] = mass[1] * xyz[1][1];
    auxR[2] = mass[1] * xyz[1][2];
    auxV[0] = mass[1] * vxyz[1][0];
    auxV[1] = mass[1] * vxyz[1][1];
    auxV[2] = mass[1] * vxyz[1][2];
    Ri[0] = auxR[0] / eta[1];
    Ri[1] = auxR[1] / eta[1];
    Ri[2] = auxR[2] / eta[1];
    Vi[0] = auxV[0] / eta[1];
    Vi[1] = auxV[1] / eta[1];
    Vi[2] = auxV[2] / eta[1];
    for (size_t i = 2; i < N; i++) {
        for(size_t d=0;d<3;d++) {
            x_star[i][d] = xyz[i][d] - Ri[d];
            v_star[i][d] = vxyz[i][d] - Vi[d];
        }
        /*jacobi_pos[i][0] = pos[i][0] - Ri[0];
        jacobi_pos[3 * i + 1] = pos[3 * i + 1] - Ri[1];
        jacobi_pos[3 * i + 2] = pos[3 * i + 2] - Ri[2];
        jacobi_vel[3 * i] = vel[3 * i] - Vi[0];
        jacobi_vel[3 * i + 1] = vel[3 * i + 1] - Vi[1];
        jacobi_vel[3 * i + 2] = vel[3 * i + 2] - Vi[2];
        */
        // Compute the next internal c.o.m. and momentum of the sequence:
        if (i < (N - 1)) {
            for(size_t d=0;d<3;d++) {
                auxR[d] += (mass[i] * xyz[i][d]);
                auxV[d] += (mass[i] * vxyz[i][d]);
            }
            /*auxR[0] += (masses[i] * pos[3 * i]);
            auxR[1] += (masses[i] * pos[3 * i + 1]);
            auxR[2] += (masses[i] * pos[3 * i + 2]);
            auxV[0] += (masses[i] * vel[3 * i]);
            auxV[1] += (masses[i] * vel[3 * i + 1]);
            auxV[2] += (masses[i] * vel[3 * i + 2]);*/
            for(size_t d=0;d<3;d++) {
                Ri[d] = auxR[d] / eta[i];
                Vi[d] = auxV[d] / eta[i];
            }
            /*Ri[0] = auxR[0] / eta[i];
            Ri[1] = auxR[1] / eta[i];
            Ri[2] = auxR[2] / eta[i];
            Vi[0] = auxV[0] / eta[i];
            Vi[1] = auxV[1] / eta[i];
            Vi[2] = auxV[2] / eta[i];*/
        }
    }
    return;
}
void jacobi2helio(double x_star[NMAX][3], double v_star[NMAX][3], double mass[NMAX], int N, double xyz[NMAX][3], double vxyz[NMAX][3]) {
//void jacobi2helio(double *jacobi_pos, double *jacobi_vel, double *masses, size_t nbodies, double *pos, double *vel) {

    // the helio coordinates for the second body is the same as the jacobi
    xyz[1][0] = x_star[1][0]; xyz[1][1] = x_star[1][1]; xyz[1][2] = x_star[1][2];
    vxyz[1][0] = v_star[1][0]; vxyz[1][1] = v_star[1][1]; vxyz[1][2] = v_star[1][2];

    // Compute etas (interior masses):
    double eta[N];
    eta[0] = mass[0];
    for (size_t i = 1; i < N; i++) {
        eta[i] = mass[i] + eta[i - 1];
    }

    // Assume central body at rest:
    xyz[0][0] = 0.0;
    xyz[0][1] = 0.0;
    xyz[0][2] = 0.0;
    vxyz[0][0] = 0.0;
    vxyz[0][1] = 0.0;
    vxyz[0][2] = 0.0;

    // Heliocentric coordinates of first body coincide with Jacobi, leave as they are.

    // Compute internal c.o.m. and momentum:
    double Ri[3], Vi[3];
    Ri[0] = mass[1] * xyz[1][0] / eta[1];
    Ri[1] = mass[1] * xyz[1][1] / eta[1];
    Ri[2] = mass[1] * xyz[1][2] / eta[1];
    Vi[0] = mass[1] * vxyz[1][0] / eta[1];
    Vi[1] = mass[1] * vxyz[1][1] / eta[1];
    Vi[2] = mass[1] * vxyz[1][2] / eta[1];
    for (size_t i = 2; i < N; i++) {
        for(size_t d=0;d<3;d++) {
            xyz[i][d] = x_star[i][d] + Ri[d];
            vxyz[i][d] = v_star[i][d] + Vi[d];
        }
        /*pos[3 * i] = jacobi_pos[3 * i] + Ri[0];
        pos[3 * i + 1] = jacobi_pos[3 * i + 1] + Ri[1];
        pos[3 * i + 2] = jacobi_pos[3 * i + 2] + Ri[2];
        vel[3 * i] = jacobi_vel[3 * i] + Vi[0];
        vel[3 * i + 1] = jacobi_vel[3 * i + 1] + Vi[1];
        vel[3 * i + 2] = jacobi_vel[3 * i + 2] + Vi[2];
        */
        // Compute the next internal c.o.m. and momentum of the sequence:
        if (i < (N - 1)) {
            for(size_t d=0;d<3;d++) {
                Ri[d] += mass[i] * x_star[i][d] / eta[i];
                Vi[d] += mass[i] * v_star[i][d] / eta[i];
            }
            /*Ri[0] += masses[i] * jacobi_pos[3 * i] / eta[i];
            Ri[1] += masses[i] * jacobi_pos[3 * i + 1] / eta[i];
            Ri[2] += masses[i] * jacobi_pos[3 * i + 2] / eta[i];
            Vi[0] += masses[i] * jacobi_vel[3 * i] / eta[i];
            Vi[1] += masses[i] * jacobi_vel[3 * i + 1] / eta[i];
            Vi[2] += masses[i] * jacobi_vel[3 * i + 2] / eta[i];*/
        }
        //for (size_t i = 0; i < 3; i++) {
        //}
    }
    return;
}

void bary2helio(double pos[NMAX][3], double vel[NMAX][3], int N) {
    for(size_t i=0;i<N;i++) {
        for(size_t d=0;d<3;d++) {
            pos[i][d] -= pos[0][d];
            vel[i][d] -= vel[0][d];
        }
    }
}

void helio2bary(double mass[NMAX], double pos[NMAX][3], double vel[NMAX][3], int N) {
    double rcm[3] = {0.0,0.0,0.0};
    double vcm[3] = {0.0, 0.0, 0.0};
    double mtot = 0.0;
    for(size_t i=0;i<N;i++) {
        for(size_t d=0;d<3;d++) {
            rcm[d] += pos[i][d] * mass[i];
            vcm[d] += vel[i][d] * mass[i];
        }
        mtot += mass[i];
    }
    
    for(size_t d=0;d<3;d++) {
        rcm[d] /= mtot;
        vcm[d] /= mtot;
    }
    for(size_t i=0;i<N;i++) {
        for(size_t d=0;d<3;d++) {
            pos[i][d] -= rcm[d];
            vel[i][d] -= vcm[d];
        }
    }

}





//These functions are mainly used for bug fixing
void Py_cart2jacobi(double *mass, double *x, double *y, double *z,
        double *vx, double *vy, double *vz, int N,
        double *eta, double *x_star, double *y_star, double *z_star,
        double *vx_star, double *vy_star, double *vz_star) {
//void Py_cart2jacobi() {
    double xyz[NMAX][3];
    double vxyz[NMAX][3];
    double eeta[NMAX];
    double mmass[NMAX];
    for(size_t i=0;i<N;i++) {
        //printf("%f\n",x[i]);
        mmass[i] = mass[i];
        xyz[i][0] = x[i];
        xyz[i][1] = y[i];
        xyz[i][2] = z[i];

        vxyz[i][0] = vx[i];
        vxyz[i][1] = vy[i];
        vxyz[i][2] = vz[i];

    }
    double xx_star[NMAX][3];
    double vv_star[NMAX][3];

    //bary2helio(xyz,vxyz,N);
    //helio2jacobi(mmass,xyz,vxyz,N,eeta,xx_star,vv_star);
    helio2jacobi(xyz,vxyz,mmass,N,xx_star,vv_star);

    for(size_t i=0;i<N;i++) {
        eta[i] = eeta[i];
        x_star[i] = xx_star[i][0];
        y_star[i] = xx_star[i][1];
        z_star[i] = xx_star[i][2];

        vx_star[i] = vv_star[i][0];
        vy_star[i] = vv_star[i][1];
        vz_star[i] = vv_star[i][2];
    }
}

void Py_jacobi2cart(double *mass, double *x, double *y, double *z,
        double *vx, double *vy, double *vz, int N,
        double *eta, double *x_star, double *y_star, double *z_star,
        double *vx_star, double *vy_star, double *vz_star) {
    double xx_star[NMAX][3];
    double vv_star[NMAX][3];
    double xyz[NMAX][3];
    double vxyz[NMAX][3];
    double eeta[NMAX];
    double mmass[NMAX];
    for(size_t i=0;i<N;i++) {
        mmass[i] = mass[i];
        eeta[i] = eta[i];
        xx_star[i][0] = x_star[i];
        xx_star[i][1] = y_star[i];
        xx_star[i][2] = z_star[i];

        vv_star[i][0] = vx_star[i];
        vv_star[i][1] = vy_star[i];
        vv_star[i][2] = vz_star[i];
    }
    jacobi2helio(xx_star, vv_star, mmass,N, xyz, vxyz);
    //jacobi2helio(eta, xx_star, vv_star, N,mmass, xyz, vxyz);
    
    for(size_t i=0;i<N;i++) {
        x[i] = xyz[i][0];
        y[i] = xyz[i][1];
        z[i] = xyz[i][2];
        
        vx[i] = vxyz[i][0];
        vy[i] = vxyz[i][1];
        vz[i] = vxyz[i][2];
    }

}
void timestep(double mass[NMAX], double pos[NMAX][3], double vel[NMAX][3], int N,double *dt, int neigh[3]) {
    double stepi;
    stepi = 1.e38;
    for (int i=0; i<3; i++) {
        int j = neigh[i];
        double dx = pos[N-1][0] - pos[j][0];
        double dy = pos[N-1][1] - pos[j][1];
        double dz = pos[N-1][2] - pos[j][2];
        double r2 = dx*dx + dy*dy + dz*dz;
        double r = sqrt(r2);

        double dvx = vel[N-1][0] - vel[j][0];
        double dvy = vel[N-1][1] - vel[j][1];
        double dvz = vel[N-1][2] - vel[j][2];

        double v2 = dvx * dvx + dvy * dvy + dvz * dvz ; 


        double vdotdr2 = (dx*dvx+dy*dvy+dz*dvz) / r2;

        double tau = ETA_SYM / M_SQRT2 * r * sqrt(r/(mass[N-1]+mass[j]));
        double dtau = 3 * tau * vdotdr2 / 2;

        dtau = (dtau < 1) ? dtau : 1;
        tau /= (1 - dtau / 2);

        stepi = min(tau, stepi);
        tau = ETA_SYM * r / sqrt(v2);
        dtau = tau * vdotdr2 * (1 + (mass[N-1]+mass[j]) / (v2 * r));
        dtau = (dtau < 1)? dtau: 1;

        tau /= (1 - dtau / 2);
        stepi = min(tau, stepi);

    }
    *dt = stepi;
}


/* works ONLY in Jacobi coordinates!
 * Translation of wh_drift from ABIE */
void evolve_kepler(double mass[NMAX], double xyz[NMAX][3], double vxyz[NMAX][3], 
        int N, double dt) {
    //double mi_star[NMAX];
    //double Mi_star[NMAX];
    //double temp_xyz[NMAX];
    /*for(size_t i=1;i<N;i++) {
        mi_star[i] = (eta[i-1]/eta[i]) * mass[i];
        Mi_star[i] = (eta[i] / eta[i-1]) * mass[0];
    }
    mi_star[0] = eta[N-1];
    Mi_star[0] = 0.0;
    */
    double Np = N; //last object is test particle
    double x_star[NMAX][3];
    double v_star[NMAX][3];
    helio2jacobi(xyz,vxyz,mass,Np,x_star,v_star);
    //cart2jacobi(xyz,vxyz,mass,Np,x_star,v_star);


    double eta0 = mass[0];
    //IMPORTANT: starting index NEEDS to be 1
    for(size_t i=1;i<Np;i++) {
        double eta = eta0 + mass[i];
        double gm = mass[0] * eta / eta0;
        //printf("from evolve_kepler\n");
        //printf("%21.16f %21.16f %21.16f\n",x_star[i][0],x_star[i][1],x_star[i][2]);
        //printf("%21.16f %21.16f %21.16f\n",v_star[i][0],v_star[i][1],v_star[i][2]);
        propagate_kepler_wh(x_star,v_star,gm,dt, Np, i);
        eta0 = eta;
    }
        jacobi2helio(x_star,v_star,mass,Np,xyz,vxyz);
    //jacobi2cart(x_star,v_star,mass,Np,xyz,vxyz);
    /*for(size_t d=0;d<3;d++) {
        xyz[0][d] += vxyz[0][d] * dt;
    }*/

}
void scale2nbody(double vel[NMAX][3], int N) {
    for(size_t i=0; i <N; i++) {
        for(size_t d=0;d<3;d++) {
            vel[i][d] = vel[i][d] / VSCALE;
        }
    }

    /*for(size_t d=0; d<3;d++) {
        vel[0][d] = vel[0][d] / VSCALE;
        vel[1][d] = vel[1][d] / VSCALE;
        vel[2][d] = vel[2][d] / VSCALE;
        vel[3][d] = vel[3][d] / VSCALE;
    }*/
}
void scale2physical(double vel[NMAX][3], int N) {
    for(size_t i=0;i<N;i++) {
        for(size_t d=0;d<3;d++) {
            vel[i][d] = vel[i][d] * VSCALE;
        }
    }
    /*for(size_t d=0; d<3;d++) {
        vel[0][d] = vel[0][d] * VSCALE;
        vel[1][d] = vel[1][d] * VSCALE;
        vel[2][d] = vel[2][d] * VSCALE;
        vel[3][d] = vel[3][d] * VSCALE;
    }*/
}

void force_and_jerk(double pos[NMAX][3], double vel[NMAX][3], double mass[NMAX], int N,double *tpos, double *tvel,
        double *ac, double *jc, int calcN, int neigh[3]) {
    //for calculating the nearest neighbors
    neigh[0] = 0;
    double first = 1.e38;
    double second = 1.e38;

    double accx =0.0;
    double accy = 0.0;
    double accz = 0.0;

    double jerkx = 0.0;
    double jerky = 0.0;
    double jerkz = 0.0;
    
    for(int j=0; j<(N-1);j++) {
        double dx = tpos[0] - pos[j][0];
        double dy = tpos[1] - pos[j][1];
        double dz = tpos[2] - pos[j][2];
        double dist = sqrt(dx*dx + dy*dy + dz*dz );

        double dist2 = dist*dist;
        double dist3 = dist*dist*dist;

        double mass_over_dist3 = mass[j]/dist3;
        accx += -mass_over_dist3*dx;
        accy += -mass_over_dist3*dy;
        accz += -mass_over_dist3*dz;

        double dvx = tvel[0] - vel[j][0];
        double dvy = tvel[1] - vel[j][1];
        double dvz = tvel[2] - vel[j][2];
        double rrdot = (dx*dvx+dy*dvy+dz*dvz);


        jerkx += -mass[j]*(dvx/dist3 - 3*rrdot*dx/(dist3*dist2));
        jerky += -mass[j]*(dvy/dist3 - 3*rrdot*dy/(dist3*dist2));
        jerkz += -mass[j]*(dvz/dist3 - 3*rrdot*dz/(dist3*dist2));
        
        if(calcN==1 && j>0) {
            if(dist < first) {
                //int temp = neigh[1];
                neigh[2] = neigh[1];
                neigh[1] = j;
                second = first;
                first = dist;
                //second = dist;
            } else if (dist < second) {
                neigh[2] = j;
                second = dist;
            }
        }

    }
    ac[0] = accx;
    ac[1] = accy;
    ac[2] = accz;
    jc[0] = jerkx;
    jc[1] = jerky;
    jc[2] = jerkz;

}
void initial(double pos[NMAX][3],double vel[NMAX][3], double mass[NMAX],int N,double *ac, double *jc, double *dt)  {
    double tpos[3] = {pos[N-1][0],pos[N-1][1],pos[N-1][2]};
    double tvel[3] = {vel[N-1][0],vel[N-1][1],vel[N-1][2]};
    int neigh[3] = {0,-1,-1};
    force_and_jerk(pos,vel,mass,N,tpos,tvel,ac,jc,1,neigh);
    timestep(mass,pos,vel,N,dt,neigh);
}

void predictor(double pos[NMAX][3], double vel[NMAX][3], double mass[NMAX], int N,double *ac, double *jc, double dt,
        double *tpos, double *tvel) {
    
    evolve_kepler(mass,pos,vel,N-1,dt);

    double dt2_over_2 = dt*dt/2;
    double dt3_over_6 = dt*dt*dt/6;
    
    for(int d=0;d<3;d++) {
        tpos[d] = pos[N-1][d] + vel[N-1][d]*dt + ac[d]*dt2_over_2 + jc[d]*dt3_over_6;
        tvel[d] = vel[N-1][d] + ac[d]*dt + jc[d]*dt2_over_2;
    }
}

void corrector (double *ac, double *jc,
        double *tpos, double *tvel, double *tposc, double *tvelc,double *aci, double *jci, double *a2ci, double *a3ci,double dt) {
    double a2i[3], a3i[3];
    double step_i2 = dt * dt;
    double step_i3 =dt * dt * dt;

    for(int d=0;d<3;d++){
        a2i[d] = 2*(-3*(ac[d]-aci[d]) - dt*(2*jc[d] + 1*jci[d]))/(step_i2);
        a3i[d] = 6*(2*(ac[d]-aci[d]) + dt*(jc[d] + jci[d]))/(step_i3);
    }

    for(int d=0;d<3;d++) {
        tposc[d] = tpos[d] + (step_i2*step_i2)*(a2i[d])/24 + (step_i2*step_i3)*(a3i[d])/120;
        tvelc[d] = tvel[d] + (step_i3)*(a2i[d])/6 + (step_i2*step_i2)*(a3i[d])/24;
    }
    for(int d=0;d<3;d++) {
        a2ci[d] = a2i[d];
        a3ci[d] = a3i[d];
    }


}
void timestep_aarseth(double *aci, double *jci, double *a2ci, double *a3ci, double *dt) {
    double mod_f, mod_f1, mod_f2, mod_f3;
    double temp_t_i;
    double old = *dt;
    mod_f = sqrt(aci[0] * aci[0] + aci[1] * aci[1] + aci[2] * aci[2]);
    mod_f1 = sqrt(jci[0] * jci[0] + jci[1] * jci[1] + jci[2] * jci[2]);
    mod_f2 = sqrt(a2ci[0] * a2ci[0] + a2ci[1] * a2ci[1] + a2ci[2] * a2ci[2]);
    mod_f3 = sqrt(a3ci[0] * a3ci[0] + a3ci[1] * a3ci[1] + a3ci[2] * a3ci[2]);

//find an optimal timestep
    temp_t_i = sqrt(ETA* ((mod_f * mod_f2) + mod_f1*mod_f1 ) / (mod_f1*mod_f3 + mod_f2*mod_f2));

//check if timestep is not zero

    *dt = min((1.2*old), temp_t_i);
}


void evolve_hermite(double pos[NMAX][3], double vel[NMAX][3], double mass[NMAX],int N,double *ac, double *jc, double *dt, int neigh[3]) {
    double tpos[3], tvel[3];
    double tposc[3], tvelc[3];
    double aci[3], jci[3];
    double a2ci[3], a3ci[3];
    //int neigh[3] = {0,-1,-1};
    predictor(pos,vel,mass,N,ac,jc,*dt,tpos,tvel);

    force_and_jerk(pos,vel,mass,N,tpos,tvel,aci,jci,0,neigh);
    corrector(ac,jc,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);
    
    force_and_jerk(pos,vel,mass,N,tposc,tvelc,aci,jci,0,neigh);
    corrector(ac,jc,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);
    
    force_and_jerk(pos,vel,mass,N,tposc,tvelc,aci,jci,1,neigh);
    corrector(ac,jc,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);

    //force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    //corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,dt);

    /*for(int i=0;i<2;i++) {
        force_and_jerk(pos,vel,tpos,tvel,aci,jci);
        corrector(pos,vel,ac,jc,tpos,tvel,aci,jci,dt);
    }*/
    //printf("%i %i %i\n",neigh[0],neigh[1],neigh[2]);
    //exit(0);
    for(int d=0;d<3;d++) {
        pos[N-1][d] = tposc[d];
        vel[N-1][d] = tvelc[d];
        ac[d] = aci[d];
        jc[d] = jci[d];
        //a2ci[d] = a2ci[d] + *dt * a3ci[d];
    }
    //timestep_aarseth(aci,jci,a2ci,a3ci,dt);

}


int check_collision_and_bound(double pos[NMAX][3],double vel[NMAX][3], double mass[NMAX], double thresh[NMAX],int N) {
    double x = pos[N-1][0];
    double y = pos[N-1][1];
    double z = pos[N-1][2];

    double v2 = vel[N-1][0] * vel[N-1][0] + vel[N-1][1] * vel[N-1][1] + vel[N-1][2] * vel[N-1][2];
    
    double ene = 0.5 * v2;

    int flag = 0; //0 = no collision/remains bound, 1 = collision with sun/planet, 4 = becomes unbound
    //double mass[3] = {M_SUN,M_EARTH,M_MOON}; 
    //double thresh[3] = {R_SUN, R_EARTH, R_MOON};
    //int colcode[3] = {0,2,3};
    for(int j=0;j<(N-1);j++) {
        double dx = pos[j][0] - x;
        double dy = pos[j][1] - y;
        double dz = pos[j][2] - z;
        double d = sqrt(dx*dx + dy*dy + dz*dz);

        ene -= (mass[j] / d);
        if(ene <0) {
            break;
        }
        if(d <= thresh[j]) {
            flag = 1;//colcode[j];
            break;
        }
    }
    if(flag==0 && ene > 0) {
        flag = 4;
    }
    return flag;
}

int evolve_system(double pos[NMAX][3], double vel[NMAX][3], double mass[NMAX], double thresh[NMAX],int N,double *ac, double *jc, double *dt, double *t, double t_end) {
    double t_cur = *t;
    int flag = 0;
    
    //printf("%.7e %.7e\n",t_cur, *dt);
    while(t_cur < t_end) {
        if((t_cur+*dt) > t_end)
            *dt  = t_end-t_cur;
        t_cur += *dt;
        int neigh[3] = {0,-1,-1};
        evolve_hermite(pos,vel,mass,N,ac,jc,dt,neigh);
        //debug_print(pos,vel);
        flag = check_collision_and_bound(pos,vel,mass,thresh,N);
        /*if(*dp_min <= 0.01) {
            break;
        }*/
        if(flag > 0)  {
            break;
        }
        timestep(mass,pos,vel,N,dt,neigh);
    }
    return flag;
}




void integrate(double *x, double *y, double *z, double *vx, double *vy, double *vz, double *m, double *r,int N, double t_end, int* collb_flag) {
    double xyz[NMAX][3];
    double vxyz[NMAX][3];
    double thresh[NMAX];
    //double x_star[NMAX][3];
    //double v_star[NMAX][3];
    double mass[NMAX];
    double t = 0.0;
    double dt = 0.0;
    double ac[3], jc[3];
    for(size_t i=0;i<N;i++) {
        xyz[i][0] = x[i];
        xyz[i][1] = y[i];
        xyz[i][2] = z[i];
        vxyz[i][0] = vx[i];
        vxyz[i][1] = vy[i];
        vxyz[i][2] = vz[i];
        mass[i] = m[i];
        thresh[i] = r[i];
    }

    scale2nbody(vxyz,N);
    int flag = 0;
    initial(xyz,vxyz,mass,N,ac,jc,&dt);
    flag = evolve_system(xyz,vxyz,mass,thresh,N,ac,jc,&dt,&t,(t+t_end*VSCALE));
    scale2physical(vxyz,N);
    *collb_flag = flag;
    //printf("%i\n",flag);
    debug_print(xyz,vxyz,N);
    for(size_t i=0;i<N;i++) {
        x[i] = xyz[i][0];
        y[i] = xyz[i][1];
        z[i] = xyz[i][2];
        vx[i] = vxyz[i][0];
        vy[i] = vxyz[i][1];
        vz[i] = vxyz[i][2];
    }
    //scale2nbody(vxyz,N);
    //double dt = t_end / 100.0;
    //for(size_t i=0;i<100;i++) {
        
    //    scale2nbody(vxyz,N);
        //helio2jacobi(xyz,vxyz,mass,N,x_star,v_star);
        //debug_print(xyz,N);

    //    evolve_kepler(mass,xyz,vxyz,N,dt*VSCALE);
        //jacobi2helio(x_star,v_star,mass,N,xyz,vxyz);
    //    scale2physical(vxyz,N);
    //    debug_print(xyz,vxyz,N);
    //}
    //helio2jacobi(xyz,vxyz,mass,N,x_star,v_star);
    //evolve_kepler(mass,x_star,v_star,N,t_end*VSCALE);
    //jacobi2helio(x_star,v_star,mass,N,xyz,vxyz);
    //scale2physical(vxyz,N);

}
