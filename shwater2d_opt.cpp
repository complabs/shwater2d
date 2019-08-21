#include <iostream>
#include <string>
#include <cmath>       // needs isfinite()
#include <cstdio>      // used by SaveData
#include <sys/time.h>  // used by gettime
#include <omp.h>

namespace GF
{
    /** The grid functions (they are naming the cells at a gridpoint).
     */
    enum GridFunctions
    { 
        h = 0,  //!< Water height (or depth)
        u,      //!< x-velocity
        v,      //!< y-velocity
        count   //!< Number of grid functions at a gridpoint
    };
}

/** class ShallowWater2D
 *  enacapsulates a solver for the two dimensional shallow water equations
 *  using the Lax-Friedrich's scheme.
 */
class ShallowWater2D
{
    ///////////////////////////////////////////////////////////////////////////
    // Constants
    ///////////////////////////////////////////////////////////////////////////

    const double g = 9.81;  //!< The acceleration due to gravity 

    ///////////////////////////////////////////////////////////////////////////
    // The grid
    ///////////////////////////////////////////////////////////////////////////

    int m;  //!< Number of volumes (grid points) in the x-direction
    int n;  //!< Number of volumes (grid points) in the y-direction

    double* grid;  //!< The data grid

    double* x;  //!< The x-coordinates of the domain, dim = `m`
    double* y;  //!< The y coordinates of the domain, dim = `n`

    // Grid functions
    //
    inline double& Qh( int i, int j ) { return Q( GF::h, i, j ); }
    inline double& Qu( int i, int j ) { return Q( GF::u, i, j ); }
    inline double& Qv( int i, int j ) { return Q( GF::v, i, j ); }

    /** Q(gf,i,j), a generic grid function at the gridpoint (i,j).
     *
     *  @warning The structure of the multidimensional grid array affects 
     *           the optimization and the for-loop performance.
     */
    inline double& Q( int gf, int i, int j )
    {
        return grid[ ( i * n + j ) * GF::count + gf ];     // Fast
        // return grid[ ( j * m + i ) * GF::count + gf ];  // Fast
        // return grid[ ( gf * m + i ) * n + j ];          // Slow
    }

    ///////////////////////////////////////////////////////////////////////////
    // Integration variabls
    ///////////////////////////////////////////////////////////////////////////

    double tend;  //!< End time of integration

    double dx;    //!< Distance between two volumes in the x-direction
    double dy;    //!< Distance between two volumes in the y-direction
    double dt;    //!< Time step

    ///////////////////////////////////////////////////////////////////////////
    // Public methods
    ///////////////////////////////////////////////////////////////////////////

public:

    /** Constructs the solver: allocates the memory and sets up all
     *  parameters for the problem.
     */
    ShallowWater2D
    ( 
        int    mSize  = 1024,  //!< Use m volumes (grid points) in the x-direction
        int    nSize  = 1024,  //!< Use n volumes (grid points) in the y-direction
        double tEnd   = 0.1,   //!< The integration end time
        double xstart = 0.0,   //!< The domain starts here along the x-direction
        double xend   = 4.0,   //!< The domain ends here along the x-direction
        double ystart = 0.0,   //!< The domain starts here along the y-direction
        double yend   = 4.0    //!< The domain ends here along the y-direction
    )   
        : m( mSize ), n( nSize ), tend( tEnd )
    {
        std::cout << m << ", " << n << ", " << tend << ", " << std::flush;

        dx = ( xend - xstart ) / m;  // Grid spacing in the x-direction
        dy = ( yend - ystart ) / n;  // Grid spacing in the y-direction
        dt = dx / sqrt( g * 5.0 );   // Time step

        // Add two ghost volumes at each side of the domain
        //
        m = m + 2;
        n = n + 2;

        // Allocate memory for the domain
        //
        grid = new double[ m * n * GF::count ];

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
    }

    /** Destructs an object (it only frees the allocated memory).
     */
    ~ShallowWater2D ()
    {
        delete[] grid;
        delete[] x;
        delete[] y;
    }

    /** Constructs a Gauss hump as the initial data.
     */
    void InitialData
    (
        double height = 4.0, //!< A nominal height
        double eps    = 2.0, //!< A height of the bump (added to the height)
        double delta  = 0.5, //!< The bump is in fact divided by delta^2
        double xpos   = 1.0, //!< The bump x-position
        double ypos   = 1.0  //!< The bump y-position
    )
    {
        #pragma omp parallel
        {
            #pragma omp for collapse(2)
            for( int i = 0; i < m; ++i )
            {
                for( int j = 0; j < n; ++j )
                {
                    Qh( i, j ) = height;
                    Qu( i, j ) = 0.0;
                    Qv( i, j ) = 0.0;
                }
            }

            #pragma omp for collapse(2)
            for( int i = 1; i < m - 1; ++i )
            {
                for( int j = 1; j < n - 1; ++j )
                {
                    Qh( i, j ) = height
                        + eps * exp( - (   SQR( x[i] - xpos )
                                         + SQR( y[j] - ypos )
                                       ) / SQR( delta ) );
                }
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
                for( int gf = 0; gf < GF::count; ++gf )
                {
                    if( ! std::isfinite( Q( gf, i, j ) ) )
                    {
                        std::cerr << "Data validation failed." << std::endl;
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

    ///////////////////////////////////////////////////////////////////////////
    // Private methods
    ///////////////////////////////////////////////////////////////////////////

private:

    /** Square of a number.
     */
    static inline double SQR( const double& x ) { return x * x; }

    /** Returns the current time in seconds with a microsecond resolution.
     */
    static inline double gettime( void )
    {
        struct timeval tv;
        gettimeofday( &tv, NULL );
        return tv.tv_sec + 1e-6 * tv.tv_usec;
    }

    /** Updates the boundary conditions.
     */
    void apply_boundary_conditions ()
    {
        // Mask for the reflective boundary conditions
        //
        const double bc_mask[ GF::count ] = { 1.0, -1.0, -1.0 }; 

        #pragma omp for collapse(2)
        for( int j = 1; j < n - 1; ++j )
        {
            for( int gf = 0; gf < GF::count; ++gf )
            {
                Q( gf, 0,   j ) = bc_mask[gf] * Q( gf, 1,   j );
                Q( gf, m-1, j ) = bc_mask[gf] * Q( gf, m-2, j );
            }
        }

        #pragma omp for collapse(2)
        for( int j = 0; j < m; ++j )
        {
            for( int gf = 0; gf < GF::count; ++gf )
            {
                Q( gf, j, 0   ) = bc_mask[gf] * Q( gf, j, 1   );
                Q( gf, j, n-1 ) = bc_mask[gf] * Q( gf, j, n-2 );
            }
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
    if( fp == NULL ) {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return;
    }

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
            fprintf( fp, "%e %e %e\n", x[i], y[j], Qh( i, j ) );
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
            fprintf( fp, "%e\n", Qh( i, j ) );
        }
    }

    fclose( fp );

    std::cout << "Data saved to " << filename << std::endl;
}

/** Uses the Lax-Friedrich's scheme for updating the volumes.
 */
void ShallowWater2D::laxf_scheme_2d
    ( 
        double** ffx, double** ffy, 
        double** nFx, double** nFy 
    )
{
    const double half_dt_over_dx = 0.5 * dt/dx; 
    const double half_dt_over_dy = 0.5 * dt/dy;

    ///////////////////////////////////////////////////////////////////////////
    // Loop along the y-direction
    //
    #pragma omp for
    for( int j = 1; j < n; ++j )
    {
        // Calculate the flux functions in the x-direction.
        //
        for( int i = 0; i < m; ++i )
        {
            ffx[GF::h][i] = Qu( i, j );

            ffx[GF::u][i] = Qu( i, j ) * Qu( i, j ) / Qh( i, j )
                            + g * SQR( Qh(i, j ) ) / 2.0;

            ffx[GF::v][i] = Qu( i, j ) * Qv( i, j ) / Qh( i, j );

            nFx[GF::h][i] = Qh( i, j );
            nFx[GF::u][i] = Qu( i, j );
            nFx[GF::v][i] = Qv( i, j );
        }

        // Update the volume Q with the fluxes
        //
        for( int i = 1; i < m - 1; ++i ) 
        {
            for( int gf = 0; gf < GF::count; ++gf ) 
            {
                Q( gf, i, j ) = 0.5 * ( nFx[gf][i+1] + nFx[gf][i-1] )
                                + half_dt_over_dx * ( ffx[gf][i-1] - ffx[gf][i+1] );
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Loop along the x-direction
    //
    #pragma omp for
    for( int i = 1; i < m; ++i )
    {
        // Calculate the flux functions in the y-direction.
        //
        for( int j = 0; j < n; ++j )
        {
            ffy[GF::h][j] = Qv( i, j );

            ffy[GF::u][j] = Qu( i, j ) * Qv( i, j ) / Qh( i, j );

            ffy[GF::v][j] = Qv( i, j ) * Qv( i, j ) / Qh( i, j ) 
                            + g * SQR( Qh( i, j ) ) / 2.0;

            nFy[GF::h][j] = Qh( i, j );
            nFy[GF::u][j] = Qu( i, j );
            nFy[GF::v][j] = Qv( i, j );
        }

        // Update the volume Q with the fluxes
        //
        for( int j = 1; j < n - 1; ++j ) 
        {
            for( int gf = 0; gf < GF::count; ++gf ) 
            {
                Q( gf, i, j ) = 0.5 * ( nFy[gf][j+1] + nFy[gf][j-1] )
                                + half_dt_over_dy * ( ffy[gf][j-1] - ffy[gf][j+1] );
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
        std::cout <<  omp_get_num_threads () << ", " << std::flush;

        // Allocate memory for the fluxes in the x- and y-direction
        //
        double** ffx = new double*[ GF::count ];
        double** nFx = new double*[ GF::count ];
        double** ffy = new double*[ GF::count ];
        double** nFy = new double*[ GF::count ];

        ffx[0] = new double[ GF::count * m ];
        nFx[0] = new double[ GF::count * m ];
        ffy[0] = new double[ GF::count * n ];
        nFy[0] = new double[ GF::count * n ];

        for( int gf = 1; gf < GF::count; ++gf )
        {
            ffx[gf] = ffx[0] + gf * m;
            nFx[gf] = nFx[0] + gf * m;
            ffy[gf] = ffy[0] + gf * n;
            nFy[gf] = nFy[0] + gf * n;
        }

        // Apply the boundary conditions then 
        // update all the volumes using the Lax-Friedrich's scheme.
        //
        double time = 0;
        for( int i = 0; i < steps; ++i, time += dt )
        {
            apply_boundary_conditions ();
            laxf_scheme_2d( ffx, ffy, nFx, nFy );
        }

        // Free the memory reserved for the fluxes
        //
        delete[] ffx[0];   delete[] ffx;
        delete[] nFx[0];   delete[] nFx;

        delete[] ffy[0];   delete[] ffy;
        delete[] nFy[0];   delete[] nFy;
    }

    std::cout << ( gettime() - stime ) << std::endl;
}

/** The main program entry which instantiates the ShallowWater2D object,
 *  loads the initial data, runs the solver using the Lax-Friedrich's scheme,
 *  verifies the result, and optionally saves the vtk file.
*/
int main( int argc, char** argv )
{
    int    nThreads = argc >= 2 ? atoi   ( argv[1]       ) : -1   ;
    int    mSize    = argc >= 3 ? atoi   ( argv[2]       ) : 1000 ;
    int    nSize    = argc >= 4 ? atoi   ( argv[3]       ) : 1000 ;
    double tEnd     = argc >= 5 ? strtod ( argv[4], NULL ) : 0.1  ;
    int    useVTK   = argc >= 6 ? atoi   ( argv[5]       ) : 0    ;

    if ( argc == 3 ) { // If given m, but not n
        nSize = mSize;
    }

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
