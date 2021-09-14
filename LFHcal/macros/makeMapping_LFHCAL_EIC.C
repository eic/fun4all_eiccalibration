// configuration macro to create tower maps for foward hcal with different setups
// provided by Nicolas Schmidt (nschmidt@cern.ch)

#include <iostream>
#include <fstream>
#include <TMath.h>

using namespace std;

void makeMapping_LFHCAL_EIC( TString setup="2x" )
{

  /* Global detector position / transformation */
  double lhcal_x0 = 0.0; // cm,
  double lhcal_y0 = 0.0; // cm,
  double lhcal_z0 = 400.0; // cm,

  double lhcal_rot_x0 = 0.0;
  double lhcal_rot_y0 = 0.0;
  double lhcal_rot_z0 = 0.0;

  /* Detector envelope size (cone shape) */
  double lhcal_rmin1 = 0; // cm
  double lhcal_rmax1 = 262; // cm
  double lhcal_rmin2 = 0; // cm
  double lhcal_rmax2 = 262; // cm
  double lhcal_dz = 100; // cm
      

  /* Tower parameters */
  double tower_dx = 5.0; // cm
  double tower_dy = 5.0; // cm
  double tower_dz = 100.0; // cm

  double scintthick = 0.4;
  double absthick   = 1.6;
  int nLayerPerSeg  = 10;
  
  double offset_rmin_x = 0.0; // cm
  double offset_rmin_y = 0.0; // cm

  cout << "Setup selected: " << setup << endl;

  if ( setup.Contains("2x") ){
      /* Detector envelope size (cone shape) */
      lhcal_rmin1 = 14; // cm
      lhcal_rmin2 = 14; // cm

      /* Tower parameters */
      tower_dx = 2.5; // cm
      tower_dy = 2.5; // cm
  } else if ( setup.Contains("default") ) {

      /* Detector envelope size (cone shape) */
      lhcal_rmin1 = 14; // cm
      lhcal_rmin2 = 14; // cm

  } else if ( setup.Contains("FwdConfig") ) {

      /* Detector envelope size (cone shape) */
      lhcal_rmin1 = 50; // cm
      lhcal_rmin2 = 50; // cm
      tower_dx = 5.535; // cm
      tower_dy = 5.535; // cm

  }  else if ( setup.Contains("wDR") ){

      /* Detector envelope size (cone shape) */
      lhcal_rmin1 = 50; // cm
      lhcal_rmin2 = 50; // cm
  } else if ( setup.Contains("IP6-asymmetric") ) {
      /* Detector envelope size (cone shape) */
      lhcal_rmin1 = 17; // cm
      lhcal_rmin2 = 17; // cm
      
      offset_rmin_x = -10.0; // cm
      offset_rmin_y = 0.0; // cm
  } else if ( setup.Contains("asymmetric") ) {
      /* Detector envelope size (cone shape) */
      lhcal_rmin1 = 17; // cm
      lhcal_rmin2 = 17; // cm
      
      offset_rmin_x = 10.0; // cm
      offset_rmin_y = 0.0; // cm
  }

  if (setup.Contains("long")){
      lhcal_z0    = 420.0; // cm,
      lhcal_dz    = 140.; // cm
      tower_dz    = 140.; //cm
  }

  // NOTE: code below assumes tower_dx = tower_dy
  // Will need to be updated if that's not the case JGL 12/27/2015
  unsigned n_towers_j = 2 * ( (unsigned)( ( (lhcal_rmax1+0.5*tower_dx) /tower_dx) )) + 1;
  unsigned n_towers_k = n_towers_j;

  // double xpos_j0_k0 = (-1 * ( (double)( n_towers_j - 1 ) / 2 ) * tower_dx) + 0.5*tower_dx;
  double xpos_j0_k0 = -1 * ( (float)( n_towers_j - 1 ) / 2 ) * tower_dx - tower_dx;
  double ypos_j0_k0 = xpos_j0_k0;
  
  cout << "n_towers_j = " << n_towers_j << endl;
  cout << "n_towers_k = " << n_towers_k << endl;
  cout << "xpos_j0_k0 = " << xpos_j0_k0 << endl;
  cout << "ypos_j0_k0 = " << ypos_j0_k0 << endl;

  // create map
  ofstream fout(Form("towerMap_LFHCAL_%s.txt",setup.Data()));

  /* Global detector transformation */
  fout << "#Global detector geometry and transforamtion; lengths given in cm" << endl;
  fout << "Gtype " << 1 << endl;
  fout << "Gr1_inner " << lhcal_rmin1 << endl;
  fout << "Gr1_outer " << lhcal_rmax1+0.1 << endl;
  fout << "Gr2_inner " << lhcal_rmin2 << endl;
  fout << "Gr2_outer " << lhcal_rmax2+0.1 << endl;
  fout << "Gdz " << lhcal_dz+0.1 << endl;
  fout << "Gx0 " << lhcal_x0 << endl;
  fout << "Gy0 " << lhcal_y0 << endl;
  fout << "Gz0 " << lhcal_z0 << endl;
  fout << "Grot_x " << lhcal_rot_x0 << endl;
  fout << "Grot_y " << lhcal_rot_y0 << endl;
  fout << "Grot_z " << lhcal_rot_z0 << endl;
  fout << "Gtower_dx " << tower_dx << endl;
  fout << "Gtower_dy " << tower_dy << endl;
  fout << "Gtower_dz " << tower_dz << endl;
  fout << "# Layer settings for tower and longitudinal segmentation" << endl;
  fout << "nlayerspertowerseg " << nLayerPerSeg << endl;
  fout << "thickness_absorber " << absthick << endl;
  fout << "thickness_scintillator " << scintthick << endl;
  fout << "xoffset " << offset_rmin_x << endl;
  fout << "yoffset " << offset_rmin_y << endl;
  
  /* Tower mapping */
  fout << "#Tower type,idx_j,idx_k,idx_l,x[cm],y[cm],z[cm],dx[cm],dy[cm],dz[cm],rot_x,rot_y,rot_z" << endl;

  unsigned int twr_count = 0; 
  unsigned idx_l = 0;

  for (int idx_j = 0; idx_j < n_towers_j; idx_j++){
    for (int idx_k = 0; idx_k < n_towers_k; idx_k++){

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

      if ( r_corner_1 > lhcal_rmax1 ||
          r_corner_2 > lhcal_rmax1 ||
          r_corner_3 > lhcal_rmax1 ||
          r_corner_4 > lhcal_rmax1 )
        continue;

      if ( r_corner_1_offset < lhcal_rmin1 ||
          r_corner_2_offset < lhcal_rmin1 ||
          r_corner_3_offset < lhcal_rmin1 ||
          r_corner_4_offset < lhcal_rmin1 )
        continue;
      
      fout << "Tower " << 0 << "\t" << idx_j << "\t" << idx_k << "\t" << idx_l << "\t" << xpos << "\t" << ypos << "\t" << zpos << "\t" << tower_dx << "\t" << tower_dy << "\t" << tower_dz << "\t0\t0\t0" << endl;
      
      twr_count++; 
      
    }
  }

  fout.close();

  cout << "Placed " << twr_count << " towers in mapping file." << endl; 

}
