#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Kb 8.6173333e-5     //Boltzmann constant (in eV/K)
#define c0 2997924.58       //Speed of light (in angstrom/ps)
#define utoev 9.3149e8      //Coefficient between u and eV

const double timestep=0.01;     //Delta t (in ps)
const int count=512;            //Number of particles
const double size=45;           //Side length of the box (in angstrom)
const double temp0=600;         //Initial temperature (in K)
const double armass=39.948;     //Particle (Ar) mass (in u)
const int stepnum=5000;         //Number of steps
const int framerate=10;         //Output rate
const double sigma=3.4;         //(in angstrom)
const double epsilon=1.03e-2;   //(in eV)
const double tau=100;           //(in ps^-1)
const double bathtemp=85;       //Heatbath temperature (in K)
const int qlength=1000;         //Queue length

void 
    initialp(double a[][3], double b[][3][2]),                                                                                                              //Initialize particle position
    initialv(double a[][3], double tem, double mass, int n),                                                                                                //Initialize particle velocity
    calspeedsq(double a[][3], double b[], int n),                                                                                                           //Calculate velocity squared
    output(double parpos[][3], double parv[][3], double force[][3], double virpos[][3][2], int n, double speedsq[], double temperature[], double mass),
    walk(double s[][3], double v[][3], double f[][3], double vp[][3][2], int n, double deltat, double bd, double mass),                                     //Move one step for all particle
    cmcorr(double a[][3], int n),                                                                                                                           //Center of mass correction
    calv(double v[][3], double f[][3], int n, double mass, double deltat),                                                                                  //Calculate new velocity
    calf(double s[][3], double f[][3], int n, double lamb, double eps),                                                                                     //Calculate force
    thermalbath(double parv[][3], double temperature, double t0, double n, double tau);                                                                     //Heat bath

double
    sumarray(double a[], int n);    //Sum an array
    
int main()
{
    double parpos[count][3], parv[count][3], speedsq[count], force[count][3], virpos[count][3][2], temperature[qlength], mass;
    int seed;

    printf("Seed: ");
    scanf("%d", &seed);
    srand(seed);
    mass=armass*utoev/c0/c0;                                                                                                    //Calculate the mass of a particle in eV
    initialp(parpos, virpos);
    initialv(parv, temp0, mass, count);
    cmcorr(parv, count);
    output(parpos, parv, force, virpos, count, speedsq, temperature, mass);
    return 0;
}

//Generate particles in a box with regular separation
void initialp(double a[][3], double b[][3][2])
{
    int i, j, k, m;

    m=0;
    for (i=5; i<=40; i=i+5)
    {
        for (j=5; j<=40; j=j+5)
        {
            for (k=5; k<=40; k=k+5)
            {
                a[m][0]=i;
                a[m][1]=j;
                a[m][2]=k;
                b[m][0][0]=i;
                b[m][1][0]=j;
                b[m][2][0]=k;
                b[m][0][1]=0;
                b[m][1][1]=0;
                b[m][2][1]=0;
                m=m+1;
            }
        }
    }
    return;
}

//Random velocity with respect to temperature
void initialv(double a[][3], double tem, double mass, int n)
{
    int i, j;
    double x1, x2, sf;

    sf=sqrt(Kb*tem/mass);
    for (i=0; i<n; ++i)
    {
        for (j=0; j<3; ++j)
        {
            x1=(rand() % 999+1)/1000.0;
            x2=(rand() % 999+1)/1000.0;
            a[i][j]=sqrt(-2*log(x1))*cos(2*M_PI*x2)*sf;         //Box-Muller transform
        }
    }
    return;
}

void calspeedsq(double a[][3], double b[], int n)
{
    int i;

    for (i=0; i<n; ++i)
    {
        b[i]=a[i][0]*a[i][0]+a[i][1]*a[i][1]+a[i][2]*a[i][2];
    }
    return;
}

void output(double parpos[][3], double parv[][3], double force[][3], double virpos[][3][2], int n, double speedsq[], double temperature[], double mass)
{
    int i, j, frame;
    double temp;
    FILE *particle, *virtuald, *velocity, *tempf;

    frame=0;
    calf(parpos, force, n, sigma, epsilon);
    particle=fopen("particle.xyz", "w");                                                                                                                 //Store particle position
    virtuald=fopen("virtualposition.txt", "w");
    velocity=fopen("velocity.txt", "w");
    tempf=fopen("temperature.txt", "w");
    fprintf(particle, "%d\n", n);
    fprintf(particle, "frame 0\n");
    for (j=0; j<n; ++j)
    {
        fprintf(particle, "Ar %lf\t%lf\t%lf\n", parpos[j][0], parpos[j][1], parpos[j][2]);                                                               //Output initial position
    }
    fprintf(virtuald, "%d\n", 0);
    for (j=0; j<n; ++j)
    {
        fprintf(virtuald, "%d\t%lf\t%lf\t%lf\n", j, virpos[j][0][0], virpos[j][1][0], virpos[j][2][0]);
    }
    fprintf(velocity, "0\n");
    for (j=0; j<n; ++j)
    {
        fprintf(velocity, "%d\t%lf\t%lf\t%lf\n", j, parv[j][0], parv[j][1], parv[j][2]);
    }
    for (i=1; i<=stepnum; ++i)
    {
        walk(parpos, parv, force, virpos, n, timestep, size, mass);
        calv(parv, force, n, mass, timestep);
        calf(parpos, force, n, sigma, epsilon);
        calv(parv, force, n, mass, timestep);
        calspeedsq(parv, speedsq, n);
        temperature[i % qlength]=mass*sumarray(speedsq, n)/n/3/Kb;
        thermalbath(parv, temperature[i % qlength], bathtemp, n, tau);
        fprintf(virtuald, "%d\n", i);
        for (j=0; j<n; ++j)
        {
            fprintf(virtuald, "%d\t%lf\t%lf\t%lf\n", j, virpos[j][0][0], virpos[j][1][0], virpos[j][2][0]);
        }
        fprintf(velocity, "%d\n", i);
        for (j=0; j<n; ++j)
        {
            fprintf(velocity, "%d\t%lf\t%lf\t%lf\n", j, parv[j][0], parv[j][1], parv[j][2]);
        }
        if (i>=qlength)
        {
            temp=sumarray(temperature, qlength)/qlength;
            fprintf(tempf, "%d\t%lf\n", i, temp);
        }
        if (i % framerate==0)
        {
            frame=frame+framerate;
            fprintf(particle, "%d\n", n);
            fprintf(particle, "frame %d\n", frame);
            for (j=0; j<n; ++j)
            {
                fprintf(particle, "Ar %lf\t%lf\t%lf\n", parpos[j][0], parpos[j][1], parpos[j][2]);              //Output particle position
            }
        }
    }
    fclose(particle);
    fclose(virtuald);
    fclose(velocity);
    fclose(tempf);
    return;
}

void walk(double s[][3], double v[][3], double f[][3], double vp[][3][2], int n, double deltat, double bd, double mass)
{
    int i, j;
    
    for (i=0; i<n; ++i)
    {
        for (j=0; j<3; ++j)
        {
            s[i][j]=s[i][j]+(v[i][j]+f[i][j]/2/mass*deltat)*deltat;
            if (vp[i][j][1]==0)
            {
                vp[i][j][0]=vp[i][j][0]+(v[i][j]+f[i][j]/2/mass*deltat)*deltat;
            }
            else
            {
                vp[i][j][0]=vp[i][j][0]-(v[i][j]+f[i][j]/2/mass*deltat)*deltat;
            }
            if (s[i][j]<0)
            {
                s[i][j]=-s[i][j];
                v[i][j]=-v[i][j];
                if (vp[i][j][1]==0)
                {
                    vp[i][j][1]=1;
                }
                else
                {
                    vp[i][j][1]=0;
                }
            }
            if (s[i][j]>bd)
            {
                s[i][j]=bd+bd-s[i][j];
                v[i][j]=-v[i][j];
                if (vp[i][j][1]==0)
                {
                    vp[i][j][1]=1;
                }
                else
                {
                    vp[i][j][1]=0;
                }
            }
        }        
    }
    return;
}

double sumarray(double a[], int n)
{
    int i;
    double temp;

    temp=0;
    for (i=0; i<n; ++i)
    {
        temp=temp+a[i];
    }
    return temp;
}

void cmcorr(double a[][3], int n)
{
    double sum, temp;
    int i, j;

    for (i=0; i<3; ++i)
    {
        sum=0;
        for (j=0; j<n; ++j)
        {
            sum=sum+a[j][i];
        }
        temp=sum/n;
        for (j=0; j<n; ++j)
        {
            a[j][i]=a[j][i]-temp;
        }
    }
    return;
}

void calf(double s[][3], double f[][3], int n, double sig, double eps)
{
    int i, j, k;
    double r[3], fmag, d, d2, sig2, temp2, temp6, ftemp;

    sig2=sig*sig;
    for (i=0; i<n; ++i)
    {
        for (k=0; k<3; ++k)
        {
            f[i][k]=0;
        }
    }
    for (i=0; i<n-1; ++i)
    {
        for (j=i+1; j<n; ++j)
        {
            for (k=0; k<3; ++k)
            {
                r[k]=s[i][k]-s[j][k];
            }
            d2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            d=sqrt(d2);
            temp2=sig2/d2;
            temp6=temp2*temp2*temp2;
            fmag=24*eps/d*(2*temp6-1)*temp6;                                //Calculate force
            for (k=0; k<3; ++k)
            {
                ftemp=fmag*r[k]/d;
                f[i][k]=f[i][k]+ftemp;
                f[j][k]=f[j][k]-ftemp;
            }
        }
    };
}

void calv(double v[][3], double f[][3], int n, double mass, double deltat)
{
    int i, j;

    for (i=0; i<n; ++i)
    {
        for (j=0; j<3; ++j)
        {
            v[i][j]=v[i][j]+f[i][j]/2/mass*deltat;
        }
    }
    return;
}

void thermalbath(double parv[][3], double temperature, double t0, double n, double tau)
{
    int i, j;

    for (i=0; i<n; ++i)
    {
        for (j=0; j<3; ++j)
        {
            parv[i][j]=parv[i][j]*sqrt(1+((t0/temperature-1)/tau));
        }
    }
    return;
}
