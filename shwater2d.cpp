/** shwater2d.cpp
 *  solves the two dimensional shallow water equations
 *  using the Lax-Friedrich's scheme
 */

#include <iostream>
#include <cmath> // isfinite()
//#include <omp.h>
#include <sys/time.h>

const int cell_size = 3;
const double xstart = 0.0;
const double ystart = 0.0;
const double xend = 4.0;
const double yend = 4.0;

extern void save_vtk( double* Q, double* x, double* y, int m, int n );

#define Q(i, j, k) Q[((k) + n * ((j) + m * (i)))]

/** Timing function with a microsecond resolution.
 */
inline double gettime( void )
{
    struct timeval tv;
    gettimeofday( &tv,NULL );
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

/** Check if the solution is finite.
 */
void validate( double* Q, int m, int n )
{
    for( int i = 0; i < n; ++i )
    {
        for( int j = 0; j < m; ++j )
        {
            for( int k = 0; k < cell_size; ++k )
            {
                if( ! std::isfinite( Q( k, j, i ) ) )
                {
                    std::cerr << "Invalid solution" << std::endl;
                    exit( -1 );
                }
            }
        }
    }
}

/** Flux function in the x-direction
 */
void fx( double* Q, double** fq, int m, int n, int j )
{
    const double g = 9.81;

    for( int i = 0; i < m; ++i )
    {
        fq[0][i] = Q(1, i, j);
        fq[1][i] = (pow(Q(1, i, j), 2) / Q(0, i, j))  +
                   (g * pow(Q(0, i, j), 2)) / 2.0;
        fq[2][i] = (Q(1, i, j) * Q(2, i, j)) / Q(0, i, j);
    }

}

/** Flux function in the y-direction
 */
void fy( double* Q, double** fq, int m, int n, int i)
{
    const double g = 9.81;

    for( int j = 0; j < n; ++j )
    {
        fq[0][j] = Q(2, i, j);
        fq[1][j] = (Q(1, i, j) * Q(2, i, j)) / Q(0, i, j);
        fq[2][j] = (pow(Q(2, i, j), 2) / Q(0, i, j))  +
                   (g * pow(Q(0, i, j), 2)) / 2.0;
    }
}

/** This is the Lax-Friedrich's scheme for updating volumes
 *  Try to parallelize it in an efficient way!
 */
void laxf_scheme_2d
(
    double* Q,
    double** ffx, double** ffy,
    double** nFx, double** nFy,
    int m, int n,
    double dx, double dy, double dt
)
{
    // Calculate and update fluxes in the x-direction
    //
    for( int i = 1; i < n; ++i )
    {
        fx(Q, ffx, m, n, i);
        for( int j = 1; j < m; ++j )
            for( int k = 0; k < cell_size; ++k )
                nFx[k][j] = 0.5 * ((ffx[k][j-1] + ffx[k][j]) -
                                   dx/dt * (Q(k, j, i) - Q(k, j-1, i)));
        for( int j = 1; j < m-1; ++j )
            for( int k = 0; k < cell_size; ++k )
                Q(k, j, i) = Q(k, j, i)  - dt/dx * ((nFx[k][j+1] - nFx[k][j]));

    }

    // Calculate and update fluxes in the y-direction
    //
    for( int i = 1; i < m; ++i )
    {
        fy( Q, ffy, m, n, i );
        for( int j = 1; j < n; ++j )
            for( int k = 0; k < cell_size; ++k )
                nFy[k][j] = 0.5 * ((ffy[k][j-1] + ffy[k][j]) -
                                   dy/dt * (Q(k, i, j) - Q(k, i, j -1)));
        for( int j = 1; j <  n-1; ++j )
            for( int k = 0; k < cell_size; ++k )
                Q(k,i,j) = Q(k,i,j) -  dt/dy * ((nFy[k][j+1]  -  nFy[k][j]));
    }
}

/** This is the main solver routine, parallelize this.
 *  But don't forget the subroutine laxf_scheme_2d
 */
void solver
(
    double* Q, double** ffx, double** ffy,
    double** nFx, double** nFy,
    int m, int n, double tend,
    double dx, double dy, double dt
)
{
    double bc_mask[3] = { 1.0, -1.0, -1.0 };

    int steps = std::ceil( tend / dt );
    double time = 0;
    for( int i = 0; i < steps; ++i, time += dt )
    {
        // Apply boundary condition
        //
        for( int j = 1; j < n - 1; ++j )
        {
            for( int k = 0; k < cell_size; ++k )
            {
                Q(k, 0, j) = bc_mask[k] *  Q(k, 1, j);
                Q(k, m-1, j) = bc_mask[k] *  Q(k, m-2, j);
            }
        }

        for( int j = 0; j < m; ++j )
        {
            for( int k = 0; k < cell_size; ++k )
            {
                Q(k, j, 0) = bc_mask[k] * Q(k, j, 1);
                Q(k, j, n-1) = bc_mask[k] * Q(k, j, n-2);
            }
        }

        // Update all volumes with the Lax-Friedrich's scheme
        //
        laxf_scheme_2d( Q, ffx, ffy, nFx, nFy, m, n, dx, dy, dt );
    }
}

/** This is the main routine of the program, which allocates memory
 *  and setup all parameters for the problem.
 *
 *  You don't need to parallelize anything here!
 *
 *  However, it might be useful to change the m and n parameters
 *  during debugging
*/
void shwater2d ()
{
    int    m     = 1000;                   // Use m volumes in the x-direction
    int    n     = 1000;                   // Use n volumes in the y-direction

    double epsi  = 2.0;                    // Parameter used for the initial condition
    double delta = 0.5;                    // Parameter used for the initial condition

    double dx   = ( xend - xstart ) / m;   // Distance between two volumes (x-direction)
    double dy   = ( yend - ystart ) / n;   // Distance between two volumes (y-direction)
    double dt   = dx / sqrt( 9.81 * 5.0 ); // Time step
    double tend = 0.1;                     // End time

    // Add two ghost volumes at each side of the domain
    //
    m = m + 2;
    n = n + 2;

    // Allocate memory for the domain
    //
    double* Q = new double[ m * n * cell_size ];
    double* x = new double[ m ];
    double* y = new double[ n ];

    // Allocate memory for fluxes
    //
    double** ffx = new double*[ cell_size ];
    double** ffy = new double*[ cell_size ];
    double** nFx = new double*[ cell_size ];
    double** nFy = new double*[ cell_size ];

    ffx[0] = new double[ cell_size * m ];
    ffy[0] = new double[ cell_size * m ];
    nFx[0] = new double[ cell_size * n ];
    nFy[0] = new double[ cell_size * n ];

    for( int i = 0; i < cell_size; ++i )
    {
        ffx[i] =  ffx[0] + i * m;
        nFx[i] =  nFx[0] + i * m;
        ffy[i] =  ffy[0] + i * n;
        nFy[i] =  nFy[0] + i * n;
    }

    double tmp = xstart - dx/2;
    for( int i = 0; i < m; ++i ) {
        x[i] = tmp + i * dx;
    }

    tmp = ystart - dy/2;
    for( int i = 0; i < n; ++i  ) {
        y[i] = tmp + i * dy;
    }

    // Set the initial Gauss hump
    //
    for( int i = 0; i < m; ++i )
    {
        for( int j = 0; j < n; ++j )
        {
            Q( 0, i, j ) = 4.0;
            Q( 1, i, j ) = 0.0;
            Q( 2, i, j ) = 0.0;
        }
    }

    for( int i = 1; i < m-1; ++i )
    {
        for( int j = 1; j < n-1; ++j )
        {
            Q(0, i, j) = 4.0 + epsi * exp(-(pow(x[i] - xend / 4.0, 2) + pow(y[j] - yend / 4.0, 2)) /
                                          (pow(delta, 2)));
        }
    }

    double stime = gettime ();

    solver(Q, ffx, ffy, nFx, nFy, m, n, tend, dx, dy, dt);

    double etime = gettime ();

    validate( Q, m, n );

    std::cout << "Solver took " << etime - stime << " seconds" << std::endl;

    // Uncomment this line if you want visualize the result in ParaView
    //
    save_vtk( Q, x, y, m, n );

    delete[] Q;
    delete[] x;
    delete[] y;

    delete[] ffx[0];
    delete[] ffy[0];
    delete[] nFx[0];
    delete[] nFy[0];

    delete[] ffx;
    delete[] ffy;
    delete[] nFx;
    delete[] nFy;
}

int main( int argc, char** argv )
{
    shwater2d ();
    return 0;
}
