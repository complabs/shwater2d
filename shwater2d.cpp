#include <iostream>
#include <cmath> // isfinite()
#include <string>
#include <sys/time.h>
#include <omp.h>

/** class ShallowWater2D
 *  enacapsulates a solver for the two dimensional shallow water equations
 *  using the Lax-Friedrich's scheme.
 */
class ShallowWater2D
{
    const double g      = 9.81;  //!< Standard acceleration of the Earth
    const double xstart = 0.0;   //!< The domain starts here in the x-direction
    const double xend   = 4.0;   //!< The domain ends here in the x-direction
    const double ystart = 0.0;   //!< THe domain starts here in the y-direction
    const double yend   = 4.0;   //!< The domain ends here in the y-direction

    int m;  //!< Use m volumes in the x-direction
    int n;  //!< Use n volumes in the y-direction

    double tend;  //!< The end time of the integration

    double dx;  //!< The distance between two volumes in the x-direction
    double dy;  //!< The distance between two volumes in the y-direction
    double dt;  //!< The time step

    double* grid;  //!< The data grid (with the water height)
    double* x;     //!< The x-coordinates of the domain
    double* y;     //!< The y coordinates of the domain

    const int cell_size = 3;  //!< A gridpoint has three cells

    /** Access the gridpoint data.
     *  The cells are: 0 = water depth, 1 = x-velocity, 2 = y-velocity.
     */
    inline double& Q( int cell, int i, int j )
    {
        return grid[ ( cell * m + i ) * n + j ];
    }

    /** The timing function with a microsecond resolution.
     */
    static inline double gettime( void )
    {
        struct timeval tv;
        gettimeofday( &tv,NULL );
        return tv.tv_sec + 1e-6 * tv.tv_usec;
    }

public:

    /** Constructs an instance of the solver (allocates the memory
     *  and sets up all parameters for the problem).
     */
    ShallowWater2D( int mSize, int nSize, double tEnd )
        : m( mSize ), n( nSize ), tend( tEnd )
    {
        dx = ( xend - xstart ) / m;  // Distance between two volumes (x-direction)
        dy = ( yend - ystart ) / n;  // Distance between two volumes (y-direction)
        dt = dx / sqrt( g * 5.0 );   // Time step

        // Add two ghost volumes at each side of the domain
        //
        m = m + 2;
        n = n + 2;

        // Allocate memory for the domain
        //
        grid = new double[ m * n * cell_size ];

        // x coordinates
        //
        x = new double[ m ];
        for( int i = 0; i < m; ++i ) {
            x[i] = xstart - dx / 2 + i * dx;
        }

        // y coordinates
        //
        y = new double[ n ];
        for( int j = 0; j < n; ++j  ) {
            y[j] = ystart - dy / 2 + j * dy;
        }

        std::cout << "Shallow water 2D solver"
            << ": m = " << m << ", n = " << n << ", tEnd = " << tend << std::endl;
    }

    /** Destructs an object (it only frees the allocated memory).
     */
    ~ShallowWater2D ()
    {
        delete[] grid;
        delete[] x;
        delete[] y;
    }

    /** Square of a number.
     */
    static inline double pow2( const double& x ) { return x * x; }

    /** Constructs a Gauss hump as the initial data.
     */
    void InitialData
    (
        double height = 4.0, //!< A nominal height
        double eps    = 2.0, //!< A height of the bump (added to the height)
        double delta  = 0.5  //!< The bump is in fact divided by delta^2
    )
    {
        for( int i = 0; i < m; ++i )
        {
            for( int j = 0; j < n; ++j )
            {
                Q( 0, i, j ) = height;
                Q( 1, i, j ) = 0.0;
                Q( 2, i, j ) = 0.0;
            }
        }

        for( int i = 1; i < m-1; ++i )
        {
            for( int j = 1; j < n-1; ++j )
            {
                Q( 0, i, j ) = height
                    + eps * exp( - (   pow2( x[i] - xend / 4.0 )
                                     + pow2( y[j] - yend / 4.0 )
                                   ) / pow2( delta ) );
            }
        }
    }

    /** Checks if the solution is finite.
     */
    void Validate ()
    {
        for( int i = 0; i < m; ++i )
        {
            for( int j = 0; j < n; ++j )
            {
                for( int k = 0; k < cell_size; ++k )
                {
                    if( ! std::isfinite( Q( k, i, j ) ) )
                    {
                        std::cerr << std::endl << "Invalid solution!" 
                            << std::endl << std::endl;
                        exit( -1 );
                    }
                }
            }
        }
    }

    /** Exports results in the VTK format
     */
    void SaveData( const std::string filename = "result.vtk" );

    /** Solves the two dimensional shallow water equations.
     */
    void Solver ();

private:

    /** Updates the boundary conditions.
     */
    void update_bc ()
    {
        // Mask for the reflective boundary conditions
        //
        const double bc_mask[3] = { 1.0, -1.0, -1.0 }; 

        #pragma omp for collapse(2)
        for( int j = 1; j < n - 1; ++j )
        {
            for( int k = 0; k < cell_size; ++k )
            {
                Q( k, 0,   j ) = bc_mask[k] * Q( k, 1, j );
                Q( k, m-1, j ) = bc_mask[k] * Q( k, m-2, j );
            }
        }

        #pragma omp for collapse(2)
        for( int j = 0; j < m; ++j )
        {
            for( int k = 0; k < cell_size; ++k )
            {
                Q( k, j, 0   ) = bc_mask[k] * Q( k, j, 1   );
                Q( k, j, n-1 ) = bc_mask[k] * Q( k, j, n-2 );
            }
        }
    }

    /** Calculates the flux function in the x-direction.
     */
    void calculate_ffx( double** ffx, int j )
    {
        for( int i = 0; i < m; ++i )
        {
            ffx[0][i] = Q(1, i, j);

            ffx[1][i] = pow2( Q(1, i, j) ) / Q(0, i, j)
                        + g * pow2( Q(0, i, j) ) / 2.0;

            ffx[2][i] = Q(1, i, j) * Q(2, i, j) / Q(0, i, j);
        }
    }

    /** Calculates the flux function in the y-direction.
     */
    void calculate_ffy( double** ffy, int i )
    {
        for( int j = 0; j < n; ++j )
        {
            ffy[0][j] = Q(2, i, j);

            ffy[1][j] = Q(1, i, j) * Q(2, i, j) / Q(0, i, j);

            ffy[2][j] = pow2( Q(2, i, j) ) / Q(0, i, j) 
                        + g * pow2( Q(0, i, j) ) / 2.0;
        }
    }

    /** The Lax-Friedrich's scheme for updating volumes.
     */
    void laxf_scheme_2d( double** ffx, double** ffy, double** nFx, double** nFy );
};

/** Saves the result in the vtk fileformat.
 */
void ShallowWater2D::SaveData( const std::string filename )
{
    FILE *fp = fopen( filename.c_str(), "w" );

    // Write the vtk-file header
    //
    fprintf( fp, "# vtk DataFile Version 2.0\n"
                 "VTK\n"
                 "ASCII\n"
                 "DATASET POLYDATA\n" );

    // Store the water height as polydata
    //
    fprintf( fp, "\nPOINTS %d double\n", m * n );

    for( int j = 0; j < n; ++j )
    {
        for( int i = 0; i < m; ++i ) {
            fprintf( fp, "%e %e %e\n", x[i], y[j], Q(0, i, j) );
        }
    }

    fprintf( fp, "\nVERTICES %d %d\n", n, n * ( m + 1 ) );

    for( int j = 0; j < n; ++j )
    {
        fprintf(fp, "%d ", m);

        for( int i = 0; i < m; ++i ) {
            fprintf( fp, "%d ", i + j * m );
        }

        fprintf( fp, "\n" );
    }

    // Store the lookup table
    //
    fprintf( fp,
        "POINT_DATA %d\n"
        "SCALARS height double 1\n"
        "LOOKUP_TABLE default\n", m * n);

    for( int j = 0; j < n; ++j )
    {
        for( int i = 0; i < m; ++i ) {
            fprintf( fp, "%e\n", Q( 0, i, j ) );
        }
    }
    fclose( fp );
}

/** Uses the Lax-Friedrich's scheme for updating the volumes.
 */
void ShallowWater2D::laxf_scheme_2d( double** ffx, double** ffy, double** nFx, double** nFy )
{
    // Calculate and update the fluxes in the x-direction
    //
    #pragma omp for
    for( int j = 1; j < n; ++j )
    {
        calculate_ffx( ffx, j ); // calculates ffx from Q

        for( int i = 1; i < m; ++i ) 
        {
            for( int k = 0; k < cell_size; ++k ) 
            {
                nFx[k][i] = 0.5 * ( ( ffx[k][i-1] + ffx[k][i] ) 
                                   - dx/dt * ( Q(k, i, j) - Q(k, i-1, j) ) );
            }
        }

        for( int i = 1; i < m-1; ++i ) 
        {
            for( int k = 0; k < cell_size; ++k ) 
            {
                Q(k, i, j) = Q(k, i, j) - dt/dx * ( nFx[k][i+1] - nFx[k][i] );
            }
        }
    }

    // Calculate and update the fluxes in the y-direction
    //
    #pragma omp for
    for( int i = 1; i < m; ++i )
    {
        calculate_ffy( ffy, i ); // calculates ffy from Q

        for( int j = 1; j < n; ++j ) 
        {
            for( int k = 0; k < cell_size; ++k ) 
            {
                nFy[k][j] = 0.5 * ( ( ffy[k][j-1] + ffy[k][j] ) 
                                   - dy/dt * ( Q(k, i, j) - Q(k, i, j-1) ) );
            }
        }

        for( int j = 1; j <  n-1; ++j ) 
        {
            for( int k = 0; k < cell_size; ++k ) 
            {
                Q(k,i,j) = Q(k,i,j) - dt/dy * ( nFy[k][j+1] - nFy[k][j] );
            }
        }
    }
}

/** This is the main solver routine, parallelize this.
 *  But don't forget the subroutine laxf_scheme_2d.
 */
void ShallowWater2D::Solver ()
{
    double stime = gettime ();

    int steps = std::ceil( tend / dt );

    #pragma omp parallel
    {
        #pragma omp master
        std::cout << "Running; num_threads = " << omp_get_num_threads () << std::endl;

        double time = 0;

        // Allocate memory for the fluxes in the x- and y-direction
        //
        double** ffx = new double*[ cell_size ];    ffx[0] = new double[ cell_size * m ];
        double** nFx = new double*[ cell_size ];    nFx[0] = new double[ cell_size * m ];

        double** ffy = new double*[ cell_size ];    ffy[0] = new double[ cell_size * n ];
        double** nFy = new double*[ cell_size ];    nFy[0] = new double[ cell_size * n ];

        for( int k = 1; k < cell_size; ++k )
        {
            ffx[k] =  ffx[0] + k * m;      ffy[k] =  ffy[0] + k * n;
            nFx[k] =  nFx[0] + k * m;      nFy[k] =  nFy[0] + k * n;
        }

        // Apply the boundary conditions then 
        // update all the volumes using the Lax-Friedrich's scheme.
        //
        for( int i = 0; i < steps; ++i, time += dt )
        {
            update_bc ();
            laxf_scheme_2d( ffx, ffy, nFx, nFy );
        }

        // Free the memory reserved for the fluxes
        //
        delete[] ffx[0];   delete[] ffx;
        delete[] nFx[0];   delete[] nFx;

        delete[] ffy[0];   delete[] ffy;
        delete[] nFy[0];   delete[] nFy;
    }

    double etime = gettime ();
    std::cout << "Solver took " << etime - stime << " seconds" << std::endl;
}

/** This is the main routine of the program, which allocates memory
 *  and setup all parameters for the problem.
 *
 *  You don't need to parallelize anything here!
 *
 *  However, it might be useful to change the problem volume size (m,n)
 *  during debugging.
*/
int main( int argc, char** argv )
{
    int    nThreads = argc >= 2 ? atoi   ( argv[1]       ) : -1   ;
    int    mSize    = argc >= 3 ? atoi   ( argv[2]       ) : 1024 ;
    int    nSize    = argc >= 4 ? atoi   ( argv[3]       ) : 1024 ;
    double tEnd     = argc >= 5 ? strtod ( argv[4], NULL ) : 0.1  ;
    int    useVTK   = argc >= 6 ? atoi   ( argv[5]       ) : 0    ;

    if( nThreads > 0 ) {
        omp_set_num_threads( nThreads );
    }

    ShallowWater2D shWater( mSize, nSize, tEnd );

    shWater.InitialData ();

    shWater.Solver ();

    shWater.Validate ();

    if( useVTK ) {
        shWater.SaveData ();
    }

    return 0;
}
