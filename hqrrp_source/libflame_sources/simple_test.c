#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FLAME.h"
#include "FLA_HQRRP_UT_blk_var2.h"

#define PRINT_DATA

enum Dtype { FLOAT, INT };

static void matrix_generate( FLA_Obj A );

static void init_pvt( FLA_Obj p );

static void write_matrix_to_file( FLA_Obj A, char *filename, enum Dtype dtype ) {
  FILE *fp;
  int   m_A, n_A, ldim_A;
  void * buff_A;
  int   i, j;
  
  buff_A = ( void * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );

  fp = fopen( filename, "w" );
  if ( fp == NULL ) {
    printf( "Error opening file!\n" );
    return;
  }

  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i < m_A; i++ ) {
      if ( dtype == INT ) {
        fprintf( fp, "%d ", ((int *)buff_A)[ i + j * ldim_A ] );
      } else if ( dtype == FLOAT ) {
        fprintf( fp, "%f ", ((double *)buff_A)[ i + j * ldim_A ] );
      }
    }
    fprintf( fp, "\n" );
  }

  fclose( fp );
}

// ============================================================================
int main( int argc, char *argv[] ) {
  int      m_A, n_A;
  FLA_Obj  A, p, tau;

  // Initialize FLAME.
  FLA_Init();

  // Some initializations.
  m_A = 40;
  n_A = 40;

  // Create FLAME objects, and attach buffers.
  FLA_Obj_create( FLA_DOUBLE, m_A, n_A, 0, 0, & A );
  FLA_Obj_create( FLA_INT,    n_A, 1,   0, 0, & p );
  FLA_Obj_create( FLA_DOUBLE, n_A, 1,   0, 0, & tau );

  // Generate matrix.
  //// FLA_Random_matrix( A );
  matrix_generate( A );

  write_matrix_to_file(A, "../../project_implementation/A.dat", FLOAT);

  // Initialize vector with pivots.
  init_pvt( p );

  // Print initial data.
#ifdef PRINT_DATA
//  FLA_Obj_show( " Ai = [ ", A, "%le", " ];" );
//  FLA_Obj_show( " pi = [ ", p, "%d", " ];" );
//  FLA_Obj_show( " taui = [ ", tau, "%le", " ];" );
#endif

  // Factorize matrix.
  printf( "%% Just before computing factorization.\n" );
  // New factorization.
  FLA_HQRRP_UT_blk_var2( A, p, tau, 64, 10, 1 );
  printf( "%% Just after computing factorization.\n" );

  // Print results.
#ifdef PRINT_DATA
//  FLA_Obj_show( " Af = [ ", A, "%le", " ];" );
write_matrix_to_file(A, "../../project_implementation/A_out.dat", FLOAT);
FLA_Obj_show( " pf = [ ", p, "%d", " ];" );
write_matrix_to_file(p, "../../project_implementation/p_out.dat", INT);
FLA_Obj_show( " tauf = [ ", tau, "%le", " ];" );
write_matrix_to_file(tau, "../../project_implementation/tau_out.dat", FLOAT);
#endif

  // Free objects.
  FLA_Obj_free( & A );
  FLA_Obj_free( & p );
  FLA_Obj_free( & tau );

  // Finalize FLAME.
  printf( "%% End of Program\n" );
  FLA_Finalize();


  return 0;
}

// ============================================================================
static void matrix_generate( FLA_Obj A ) {
  double  * buff_A;
  int     m_A, n_A, ldim_A;
  int     i, j, num;

  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );

  //
  // Matrix with integer values.
  // ---------------------------
  //
  if( ( m_A > 0 )&&( n_A > 0 ) ) {
    num = 1;
    for ( j = 0; j < n_A; j++ ) {
      for ( i = ( j % m_A ); i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
      for ( i = 0; i < ( j % m_A ); i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
    if( ( m_A > 0 )&&( n_A > 0 ) ) {
      buff_A[ 0 + 0 * ldim_A ] = 1.2;
    }
#if 0
    // Scale down matrix.
    if( num == 0.0 ) {
      rnum = 1.0;
    } else {
      rnum = 1.0 / num;
    }
    for ( j = 0; j < n_A; j++ ) {
      for ( i = 0; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] *= rnum;
      }
    }
#endif
  }
}

// ============================================================================
static void init_pvt( FLA_Obj p ) {
  int  * buff_p, n_p;
  int  i;

  buff_p = ( int * ) FLA_Obj_buffer_at_view( p );
  n_p    = FLA_Obj_length( p );
  for( i = 0; i < n_p; i++ ) {
    buff_p[ i ] = ( i + 1 );
  }
}


