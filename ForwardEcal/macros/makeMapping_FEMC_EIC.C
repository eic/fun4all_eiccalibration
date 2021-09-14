#include <iostream>
#include <fstream>
#include <TMath.h>

using namespace std;

void makeMapping_FEMC_EIC( TString setup="asymmetric_ROS" )
{

  /* Global detector position / transformation */
  double femc_x0 = 0.0; // cm,
  double femc_y0 = 0.0; // cm,
  double femc_z0 = 310.0; // cm,

  double femc_rot_x0 = 0.0;
  double femc_rot_y0 = 0.0;
  double femc_rot_z0 = 0.0;

  /* Detector envelope size (cone shape) */
  double femc_rmin1 = 0; // cm
  double femc_rmax1 = 182.655; // cm
  double femc_rmin2 = 0.; // cm
  double femc_rmax2 = 182.655; // cm
  double femc_dz    = 36.5; // cm

  /* Tower parameters */
  double tower_dx = 5.535; // cm // From PHENIX EMCal JGL 12/27/2015
  double tower_dy = 5.535; // cm // From PHENIX EMCal JGL 12/27/2015
  double tower_dz = 36.3; // cm
  double tower_ROS = 1.0; // cm

  double offset_rmin_x = 0.0; // cm
  double offset_rmin_y = 0.0; // cm

  cout << "Setup selected: " << setup << endl;

  if ( setup == "PHENIXEMCal" )
    {

      /* Detector envelope size (cone shape) */
      femc_rmin1 = 11.0; // cm - accomodate Mar 2020 EIC beam pipe // was 20cm
      femc_rmin2 = 12.0; // cm - accomodate Mar 2020 EIC beam pipe // was 20cm
    }
  else if ( setup == "PHENIXEMCal_2x" )
    {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 16.0; // cm - accomodate Mar 2020 EIC beam pipe // was 20cm
      femc_rmin2 = 16.0; // cm - accomodate Mar 2020 EIC beam pipe // was 20cm

      /* Tower parameters */
      // From PHENIX EMCal JGL 12/27/2015
      tower_dx = tower_dx/2; // cm
      tower_dy = tower_dy/2; // cm
    }
  else if ( setup == "PHENIXEMCal" )
    {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 11.0; // cm - accomodate Mar 2020 EIC beam pipe
      femc_rmin2 = 12.0; // cm - accomodate Mar 2020 EIC beam pipe

      /* Tower parameters */
      // From PHENIX EMCal JGL 12/27/2015
      tower_dz = 36.3; // cm (FULL SIZE)
    }
  else if ( setup == "IP6-asymmetric" )
    {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 11.5; // cm - accomodate Mar 2020 EIC beam pipe
      femc_rmin2 = 11.5; // cm - accomodate Mar 2020 EIC beam pipe

      offset_rmin_x = -7.1; // cm
      offset_rmin_y = 0.0; // cm
    }
  else if ( setup == "asymmetric" )
    {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 11.5; // cm - accomodate Mar 2020 EIC beam pipe
      femc_rmin2 = 11.5; // cm - accomodate Mar 2020 EIC beam pipe

      offset_rmin_x = 7.1; // cm
      offset_rmin_y = 0.0; // cm
    }
  else if ( setup == "IP6-asymmetric_ROS" )
    {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 11.; // cm - accomodate Mar 2020 EIC beam pipe
      femc_rmin2 = 11.; // cm - accomodate Mar 2020 EIC beam pipe

      /* Tower parameters */
      // From PHENIX EMCal JGL 12/27/2015
      tower_ROS = 3.; // readout-separated

      offset_rmin_x = -6.3; // cm
      offset_rmin_y = 0.0; // cm
    }
  else if ( setup == "asymmetric_ROS" )
    {
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 11.5; // cm - accomodate Mar 2020 EIC beam pipe
      femc_rmin2 = 11.5; // cm - accomodate Mar 2020 EIC beam pipe

      /* Tower parameters */
      // From PHENIX EMCal JGL 12/27/2015
      tower_ROS = 3.; // readout-separated

      offset_rmin_x = 7.1; // cm
      offset_rmin_y = 0.0; // cm
    }
  else if ( setup == "wDR" )
    {
     
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 50.0; // cm - accomodate Mar 2020 EIC beam pipe
      femc_rmin2 = 50.0; // cm - accomodate Mar 2020 EIC beam pipe

    }
  else if ( setup == "EIC_v1" )
    {
     
      /* Detector envelope size (cone shape) */
      femc_rmin1 = 11; // cm
      femc_rmax1 = 225; // cm
      femc_rmin2 = 12; // cm
      femc_rmax2 = 246; // cm
      femc_dz = 17; // cm

      /* Tower parameters */
      tower_dx = 3.0; // cm
      tower_dy = 3.0; // cm
      tower_dz = 17.0; // cm
    }

  // create map
  ofstream fout(Form("towerMap_FEMC_%s.txt",setup.Data()));

  /* Global detector transformation */
  fout << "#Global detector geometry and transforamtion; lengths given in cm" << endl;
  fout << "Gtype " << 2 << endl;
  fout << "Gr1_inner " << femc_rmin1 << endl;
  fout << "Gr1_outer " << Form("%5.5f",femc_rmax1+0.1) << endl;
  fout << "Gr2_inner " << femc_rmin2 << endl;
  fout << "Gr2_outer " << Form("%5.5f",femc_rmax2+0.1) << endl;
  fout << "Gdz " << femc_dz << endl;
  fout << "Gx0 " << femc_x0 << endl;
  fout << "Gy0 " << femc_y0 << endl;
  fout << "Gz0 " << femc_z0 << endl;
  fout << "Grot_x " << femc_rot_x0 << endl;
  fout << "Grot_y " << femc_rot_y0 << endl;
  fout << "Grot_z " << femc_rot_z0 << endl;
  if (tower_ROS > 1){
    fout << "Gtower2_dx " << tower_dx/tower_ROS << endl;
    fout << "Gtower2_dy " << tower_dy/tower_ROS << endl;
  } else {
    fout << "Gtower2_dx " << tower_dx << endl;
    fout << "Gtower2_dy " << tower_dy << endl;    
  }
  fout << "Gtower2_dz " << tower_dz << endl;
  fout << "xoffset " << offset_rmin_x << endl;
  fout << "yoffset " << offset_rmin_y << endl;

  /* Tower mapping */
  fout << "#Tower type,idx_j,idx_k,idx_l,x[cm],y[cm],z[cm],dx[cm],dy[cm],dz[cm],rot_x,rot_y,rot_z" << endl;

  // NOTE: code below assumes tower_dx = tower_dy
  // Will need to be updated if that's not the case JGL 12/27/2015
  unsigned n_towers_j = 2 * ( (unsigned)( (femc_rmax1/tower_dx) ));
  unsigned n_towers_k = 2 * ( (unsigned)( (femc_rmax1/tower_dx) ));

  double xpos_j0_k0 = (-1 * ( (double)( n_towers_j - 1 ) / 2 ) * tower_dx) + 0.5*tower_dx;
  double ypos_j0_k0 = (-1 * ( (double)( n_towers_j - 1 ) / 2 ) * tower_dx) + 0.5*tower_dx;

  cout << "n_towers_j = " << n_towers_j << endl;
  cout << "n_towers_k = " << n_towers_k << endl;
  cout << "xpos_j0_k0 = " << xpos_j0_k0 << endl;
  cout << "ypos_j0_k0 = " << ypos_j0_k0 << endl;

  unsigned int twr_count = 0;
  unsigned idx_l = 0;
  unsigned idtx_j = 0;
  unsigned idtx_k = 0;
  double halfNtowersize = tower_dy/2;
  double towersize    = tower_dy/tower_ROS;
  
  for (int idx_j = 0; idx_j < n_towers_j; idx_j++){
    for (int idx_k = 0; idx_k < n_towers_k; idx_k++){
      /* Calculate center position for tower */
      double xpos = xpos_j0_k0 + idx_j * tower_dx;
      double ypos = ypos_j0_k0 + idx_k * tower_dy;
      double zpos = 0;
      
      
      // check if all four corners are within envelope volume
      double r_corner_1 = sqrt( pow( xpos + halfNtowersize , 2 ) + pow( ypos + halfNtowersize , 2 ) );
      double r_corner_2 = sqrt( pow( xpos - halfNtowersize , 2 ) + pow( ypos + halfNtowersize , 2 ) );
      double r_corner_3 = sqrt( pow( xpos + halfNtowersize , 2 ) + pow( ypos - halfNtowersize , 2 ) );
      double r_corner_4 = sqrt( pow( xpos - halfNtowersize , 2 ) + pow( ypos - halfNtowersize , 2 ) );
      double r_corner_1_offset = sqrt( pow( xpos - offset_rmin_x + tower_dx/2. , 2 ) + pow( ypos + halfNtowersize , 2 ) );
      double r_corner_2_offset = sqrt( pow( xpos - offset_rmin_x - tower_dx/2. , 2 ) + pow( ypos + halfNtowersize , 2 ) );
      double r_corner_3_offset = sqrt( pow( xpos - offset_rmin_x + tower_dx/2. , 2 ) + pow( ypos - halfNtowersize , 2 ) );
      double r_corner_4_offset = sqrt( pow( xpos - offset_rmin_x - tower_dx/2. , 2 ) + pow( ypos - halfNtowersize , 2 ) );
      
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
      if(tower_ROS==1){
        fout << "Tower " << 2 << " " << idx_j << " " << idx_k << " " << idx_l << " " << Form("%4.5f",xpos) << " " << Form("%4.5f",ypos) << " " << zpos << " " << tower_dx << " " << tower_dy << " " << tower_dz << " 0 0 0" << endl;
        twr_count++;
      } else {
        for(int isubj=0; isubj<tower_ROS; isubj++){
          idtx_j = idx_j*tower_ROS + isubj;
          double xpos_sub = (double)(xpos-halfNtowersize) + (double)((isubj+0.5) * towersize);
          
          for(int isubk=0; isubk<tower_ROS; isubk++){
            double ypos_sub = (double)(ypos-halfNtowersize) + (double)((isubk+0.5) * towersize);
            idtx_k = idx_k*tower_ROS + isubk;
//             if (idtx_j == 167 && (idtx_k == 43 || idtx_k == 44)) 
//             cout << xpos_sub << "\t" << idx_k << "\t" << isubk << "\t"  << ypos << "\t"<< Form("%4.5f",(double)(ypos-(tower_dy/2.)))  << "\t"<< Form("%4.5f",ypos_sub) << "\t" << (double)(isubk * towersize)<<endl;
            double r_corner_5 = sqrt( pow( xpos_sub + towersize/2 , 2 ) + pow( ypos_sub + towersize/2 , 2 ) );
            double r_corner_6 = sqrt( pow( xpos_sub - towersize/2 , 2 ) + pow( ypos_sub + towersize/2 , 2 ) );
            double r_corner_7 = sqrt( pow( xpos_sub + towersize/2 , 2 ) + pow( ypos_sub - towersize/2 , 2 ) );
            double r_corner_8 = sqrt( pow( xpos_sub - towersize/2 , 2 ) + pow( ypos_sub - towersize/2 , 2 ) );
            if ( r_corner_5 > femc_rmax1 ||
                r_corner_6 > femc_rmax1 ||
                r_corner_7 > femc_rmax1 ||
                r_corner_8 > femc_rmax1 )
              cout << "bananas ... "<< idtx_j << "\t" << idtx_k << "\t" << r_corner_5 << "\t" << r_corner_6 << "\t" << r_corner_7 << "\t" << r_corner_8 << "\t" << femc_rmax1 << endl;

            
            if (idtx_j == 112 && (idtx_k == 0 )) { // || idtx_k == 1)
              cout << "****************************" << endl;
            cout << xpos_sub << "\t" << idx_k << "\t" << isubk << "\t"  << ypos << "\t"<< Form("%4.5f",(double)(ypos-(halfNtowersize)))  << "\t"<< Form("%4.5f",ypos_sub) << "\t" << (double)(isubk * towersize)<<endl;
              cout << xpos << "\t" << ypos << "\t" << tower_dx/2. << "\t" << r_corner_1 << "\t" << r_corner_2 << "\t" << r_corner_3 << "\t" << r_corner_4 << "\t" << femc_rmax1 << endl;
              cout << "Tower " << 2 << " " << idtx_j << " " << idtx_k << " " << idx_l << " " << Form("%4.5f",xpos_sub) << " " << Form("%4.5f",ypos_sub) << " " << zpos << " " << Form("%4.5f",towersize) << " " << Form("%4.5f",towersize) << " " << tower_dz << " 0 0 0" << endl;
              
            }
            fout << "Tower " << 2 << " " << idtx_j << " " << idtx_k << " " << idx_l << " " << Form("%4.5f",xpos_sub) << " " << Form("%4.5f",ypos_sub) << " " << zpos << " " << Form("%4.5f",towersize) << " " << Form("%4.5f",towersize) << " " << tower_dz << " 0 0 0" << endl;
            twr_count++;
          }
        }
      }
    }
	}

  fout.close();

  cout << "Placed " << twr_count << " towers in mapping file." << endl; 

}
