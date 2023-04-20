/******************************************
 * RAPTOR-solar-system
 * Apr 2, 2023
 * Author: Dipto Mukherjee
 * NOTE:
 * (1) pos, vel refer to positions and velocities in cartesian coordinates
 * but xyz, vxyz refers to positions and velocites in jacobi coordinates
 * (2) N: number of bodies excluding the test object
 *
 *
 * *****************************************/

#include "universal.h"
#define VSCALE 6.2830666414874994
#define ETA_SYM 0.00625
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


void evolve_kepler_2body(double mass[2], double pos[2][3], double vel[2][3], double dt) {
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


}
void cart2jacobi(double mass[NMAX], double xyz[NMAX][3], double vxyz[NMAX][3], int N,
        double eta[NMAX], double x_star[NMAX][3], double v_star[NMAX][3]) {
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
void helio2jacobi(double mass[NMAX], double xyz[NMAX][3], double vxyz[NMAX][3], int N,
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

void jacobi2cart(double eta[NMAX], double x_star[NMAX][3], double v_star[NMAX][3], int N,
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

void jacobi2helio(double eta[NMAX], double x_star[NMAX][3], double v_star[NMAX][3], int N,
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

    helio2jacobi(mmass,xyz,vxyz,N,eeta,xx_star,vv_star);

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

    jacobi2helio(eta, xx_star, vv_star, N,mmass, xyz, vxyz);
    
    for(size_t i=0;i<N;i++) {
        x[i] = xyz[i][0];
        y[i] = xyz[i][1];
        z[i] = xyz[i][2];
        
        vx[i] = vxyz[i][0];
        vy[i] = vxyz[i][1];
        vz[i] = vxyz[i][2];
    }

}
void timestep(double mass[NMAX], double pos[NMAX][3], double vel[NMAX][3], int N,double *dt) {
    double stepi;
    stepi = 1.e38;
    for (int j=0; j<N; j++) {
        double dx = pos[N][0] - pos[j][0];
        double dy = pos[N][1] - pos[j][1];
        double dz = pos[N][2] - pos[j][2];
        double r2 = dx*dx + dy*dy + dz*dz;
        double r = sqrt(r2);

        double dvx = vel[N][0] - vel[j][0];
        double dvy = vel[N][1] - vel[j][1];
        double dvz = vel[N][2] - vel[j][2];

        double v2 = dvx * dvx + dvy * dvy + dvz * dvz ; 


        double vdotdr2 = (dx*dvx+dy*dvy+dz*dvz) / r2;

        double tau = ETA_SYM / M_SQRT2 * r * sqrt(r/(mass[N]+mass[j]));
        double dtau = 3 * tau * vdotdr2 / 2;

        dtau = (dtau < 1) ? dtau : 1;
        tau /= (1 - dtau / 2);

        stepi = min(tau, stepi);
        tau = ETA_SYM * r / sqrt(v2);
        dtau = tau * vdotdr2 * (1 + (mass[N]+mass[j]) / (v2 * r));
        dtau = (dtau < 1)? dtau: 1;

        tau /= (1 - dtau / 2);
        stepi = min(tau, stepi);

    }
    *dt = stepi;
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

/* works ONLY in Jacobi coordinates!*/
/*void evolve_kepler(double mass[NMAX], double xyz[NMAX][3], double vxyz[NMAX][3], double eta[NMAX],
        int N) {
    double mi_star[NMAX];
    double Mi_star[NMAX];
    double temp_xyz[NMAX];
    for(size_t i=1;i<N;i++) {
        mi_star[i] = (eta[i-1]/eta[i]) * mass[i];
        Mi_star[i] = (eta[i] / eta[i-1]) * mass[0];
    }
    mi_star[0] = eta[N-1];
    Mi_star[0] = 0.0;
    
    //IMPORTANT: starting index NEEDS to be 1
    for(size_t i=1;i<N;i++) {

    }
}*/
