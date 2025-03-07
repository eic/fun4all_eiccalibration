#include <iostream>
#include <fstream>
#include <TMath.h>

using namespace std;

void makeMapping_FHCAL_EIC( TString setup="asymmetric_XL" )
{

  /* Global detector position / transformation */
  double femc_x0 = 0.0; // cm,
  double femc_y0 = 0.0; // cm,
  double femc_z0 = 400.0; // cm,

  double femc_rot_x0 = 0.0;
  double femc_rot_y0 = 0.0;
  double femc_rot_z0 = 0.0;

  /* Detector envelope size (cone shape) */
  double femc_rmin1 = 0; // cm
  double femc_rmax1 = 262.; // cm
  double femc_rmin2 = 0; // cm
  double femc_rmax2 = 262.; // cm
  double femc_dz = 100; // cm

  /* Tower parameters */
  double tower_dx = 10.0; // cm
  double tower_dy = 10.0; // cm
  double tower_dz = 100.0; // cm

  double offset_rmin_x = 0.0; // cm
  double offset_rmin_y = 0.0; // cm

  cout << "Setup selected: " << setup << endl;

  if ( setup == "2x" ){
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 14; // cm
      femc_rmin2 = 14; // cm

      /* Tower parameters */
      tower_dx = tower_dx/2.; // cm
      tower_dy = tower_dy/2.; // cm
  } else if ( setup == "4x" ) {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 14; // cm
      femc_rmin2 = 14; // cm
      femc_dz = 100; // cm

      /* Tower parameters */
      tower_dx = tower_dx/4.; // cm
      tower_dy = tower_dy/4.; // cm
  } else if ( setup == "2xdefeta" ){
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 23; // cm
      femc_rmin2 = 23; // cm

      /* Tower parameters */
      tower_dx = tower_dx/2.; // cm
      tower_dy = tower_dy/2.; // cm
  } else if ( setup == "4xdefeta" ){
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 23; // cm
      femc_rmin2 = 23; // cm
      /* Tower parameters */
      tower_dx = tower_dx/4.; // cm
      tower_dy = tower_dy/4.; // cm
      tower_dz = 100.0; // cm
  } else if ( setup == "default" ) {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 14; // cm
      femc_rmin2 = 14; // cm
  } else if ( setup == "wDR" ){
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 50; // cm
      femc_rmin2 = 50; // cm
  } else if ( setup == "asymmetric" ){
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 17; // cm
      femc_rmin2 = 17; // cm

      offset_rmin_x = 10.0; // cm
  } else if ( setup == "IP6-asymmetric" ){
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 17; // cm
      femc_rmin2 = 17; // cm

      offset_rmin_x = -10.0; // cm
  } else if ( setup == "asymmetric_XL" ){
      femc_z0 = 420.0; // cm,

      /* Detector envelope size (cone shape) */
      femc_rmin1 = 17; // cm
      femc_rmin2 = 17; // cm

      /* Tower parameters */
      tower_dz = 140.0; // cm

      offset_rmin_x = 10.0; // cm
  }

  // NOTE: code below assumes tower_dx = tower_dy
  // Will need to be updated if that's not the case JGL 12/27/2015
  unsigned n_towers_j = 2 * ( (unsigned)( ( (femc_rmax1+0.5*tower_dx) /tower_dx) )) + 1;
  unsigned n_towers_k = n_towers_j;
  // unsigned n_towers_offset_x = ( (unsigned)( ( (offset_rmin_x+0.5*tower_dx) /tower_dx) )) + 1;

  // double xpos_j0_k0 = (-1 * ( (double)( n_towers_j - 1 ) / 2 ) * tower_dx) + 0.5*tower_dx;
  double xpos_j0_k0 = -1 * ( (float)( n_towers_j - 1 ) / 2 ) * tower_dx - tower_dx;
  double ypos_j0_k0 = xpos_j0_k0;

  cout << "n_towers_j = " << n_towers_j << endl;
  cout << "n_towers_k = " << n_towers_k << endl;
  cout << "xpos_j0_k0 = " << xpos_j0_k0 << endl;
  cout << "ypos_j0_k0 = " << ypos_j0_k0 << endl;

  // create map
  ofstream fout(Form("towerMap_FHCAL_%s.txt",setup.Data()));

  /* Global detector transformation */
  fout << "#Global detector geometry and transforamtion; lengths given in cm" << endl;
  fout << "Gtype " << 1 << endl;
  fout << "Gr1_inner " << femc_rmin1-0.05 << endl;
  fout << "Gr1_outer " << femc_rmax1+0.1 << endl;
  fout << "Gr2_inner " << femc_rmin2-0.05 << endl;
  fout << "Gr2_outer " << femc_rmax2+0.1 << endl;
  fout << "Gdz " << femc_dz << endl;
  fout << "Gx0 " << femc_x0 << endl;
  fout << "Gy0 " << femc_y0 << endl;
  fout << "Gz0 " << femc_z0 << endl;
  fout << "Grot_x " << femc_rot_x0 << endl;
  fout << "Grot_y " << femc_rot_y0 << endl;
  fout << "Grot_z " << femc_rot_z0 << endl;
  fout << "Gtower_dx " << tower_dx << endl;
  fout << "Gtower_dy " << tower_dy << endl;
  fout << "Gtower_dz " << tower_dz << endl;
  fout << "xoffset " << offset_rmin_x << endl;
  fout << "yoffset " << offset_rmin_y << endl;

  /* Tower mapping */
  fout << "#Tower type,idx_j,idx_k,idx_l,x[cm],y[cm],z[cm],dx[cm],dy[cm],dz[cm],rot_x,rot_y,rot_z" << endl;

  unsigned int twr_count = 0; 
  unsigned idx_l = 0;

  for (int idx_j = 0; idx_j < n_towers_j; idx_j++)
  {
    for (int idx_k = 0; idx_k < n_towers_k; idx_k++)
    {

      /* Calculate center position for tower */
      double xpos = xpos_j0_k0 + idx_j * tower_dx;
      double ypos = ypos_j0_k0 + idx_k * tower_dy;
      double zpos = 0;

      // check if all four corners are within envelope volume
      double r_corner_1 = sqrt( pow( xpos + tower_dx/2. , 2 ) + pow( ypos + tower_dy/2. , 2 ) );
      double r_corner_2 = sqrt( pow( xpos - tower_dx/2. , 2 ) + pow( ypos + tower_dy/2. , 2 ) );
      double r_corner_3 = sqrt( pow( xpos + tower_dx/2. , 2 ) + pow( ypos - tower_dy/2. , 2 ) );
      double r_corner_4 = sqrt( pow( xpos - tower_dx/2. , 2 ) + pow( ypos - tower_dy/2. , 2 ) );
      // check if all four corners are within envelope volume
      double r_corner_1_offset = sqrt( pow( xpos - offset_rmin_x + tower_dx/2. , 2 ) + pow( ypos + tower_dy/2. , 2 ) );
      double r_corner_2_offset = sqrt( pow( xpos - offset_rmin_x - tower_dx/2. , 2 ) + pow( ypos + tower_dy/2. , 2 ) );
      double r_corner_3_offset = sqrt( pow( xpos - offset_rmin_x + tower_dx/2. , 2 ) + pow( ypos - tower_dy/2. , 2 ) );
      double r_corner_4_offset = sqrt( pow( xpos - offset_rmin_x - tower_dx/2. , 2 ) + pow( ypos - tower_dy/2. , 2 ) );

      if ( r_corner_1 > femc_rmax1 ||
          r_corner_2 > femc_rmax1 ||
          r_corner_3 > femc_rmax1 ||
          r_corner_4 > femc_rmax1 )
        continue;

      if ( r_corner_1_offset < femc_rmin1 ||
          r_corner_2_offset < femc_rmin1 ||
          r_corner_3_offset < femc_rmin1 ||
          r_corner_4_offset < femc_rmin1 )
        continue;

      fout << "Tower " << 0 << "\t" << idx_j << "\t" << idx_k << "\t" << idx_l << "\t" << xpos << "\t" << ypos << "\t" << zpos << "\t" << tower_dx << "\t" << tower_dy << "\t" << tower_dz << "\t0\t0\t0" << endl;
      
      twr_count++; 
      
    }

  }

  fout.close();

  cout << "Placed " << twr_count << " towers in mapping file." << endl; 

}
