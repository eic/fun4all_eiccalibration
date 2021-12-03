#include <fstream>
#include <iostream>
using namespace std;

void tower_pos_gen_crystal()
{

  // The Default unit if the fun4all is cm
  // Here are all the parameters you would like to change
  // One thing need to be noticed here, when you want to change the tower width,
  // no matter the inner or outer tower, please keep their width ratio always at 2 (outer/inner)
  // Or there will be an overlap problem.
  // material para: 0.5[G4_PbWO4],
  //                1.5[G4_GLASS_LEAD],
  //                2.5[G4_BARIUM_SULFATE],
  //                3.5[G4_CESIUM_IODIDE],
  //                4.5[Sci-glass-selfdefined]
  //                5.5[heavier Sci-glass-selfdefined]
  //                6.5[Sci-glass(lead)-selfdefined]
  //                7.5[Sci-glass Nathaly]
  //                8.5[Sci-glass g4e]
  
  double PWO_carbon_Gap = 0.025,
         PWO_air_Gap = 0.,
         PWO_reffoil_Gap = 0.0066,
         PWO_tedlar_Gap = 0.0051,
         PWO_Gap = PWO_carbon_Gap + PWO_air_Gap + PWO_reffoil_Gap + PWO_tedlar_Gap,
         PWO_frame_depth = 2.0,
         PWO_face_lip = 0.1,
         PWO_Width = 2.0,
         PWO_tower_width = (PWO_Width + 2. * PWO_Gap),
         PWO_Thickness = 20,
         PWO_OuterR = 61., // 200cm version, pure PbW04 SciGlass
         PWO_InnerR = 7.,
         PWO_Frame_OuterR = 64.,
         PWO_Frame_InnerR = 55.,
         Inner_mat = 0.5; 
    
  int PWO_frame_style = 1; // 0 = carbon wrap around each tower, 1 = spacers on front and back
double inner_glo_zpos = -185; // 250.cm alternative
// double inner_glo_zpos = -200., // 250.cm alternative
 

  
  // CRYSTAL
  double diameter = 2 * PWO_OuterR;
  double leftTowerPos, topTowerPos;
  // How many towers do we have per row/column?
  // Add a gap + diameter as if we have N towers, we have N-1 gaps;
   int towersInRow = std::floor((diameter + 2. * PWO_Gap) /  PWO_tower_width);
  //  cout << "towersInRow: " << towersInRow << endl;
  
  // Is it odd or even number of towersInRow

  if(towersInRow%2)
    {
      //             |
      //      [ ][ ][ ][ ][ ]
      //       ^     |
      int towersInHalfRow = std::ceil(towersInRow/2.0);
           leftTowerPos = topTowerPos = -towersInHalfRow * PWO_tower_width;
    }
  else
    {
      //               |
      //      [ ][ ][ ][ ][ ][ ]
      //       ^      |
      int towersInHalfRow = towersInRow/2;
      leftTowerPos = topTowerPos = -(towersInHalfRow - 0.5) * PWO_tower_width;
    }

  std::printf("\nCE EMCAL PWO SQUARE START\n");
  std::printf("PWO_Thickness  = %.5f cm;\n", PWO_Thickness );
  std::printf("PWO_tower_width= %.5f cm;\n", PWO_tower_width );
  std::printf("PWO_Width      = %.5f cm;\n", PWO_Width );
  std::printf("PWO_total_Gap  = %.5f cm;\n", PWO_Gap );
  std::printf("PWO_carbon_Gap = %.5f cm;\n", PWO_carbon_Gap );
  std::printf("PWO_reflec_Gap = %.5f cm;\n", PWO_reffoil_Gap );
  std::printf("PWO_tedlar_Gap = %.5f cm;\n", PWO_tedlar_Gap );
  std::printf("PWO_air_Gap    = %.5f cm;\n", PWO_air_Gap );
  std::printf("PWO_InnerR     = %.5f cm;\n", PWO_InnerR );
  std::printf("PWO_OuterR     = %.5f cm;\n", PWO_OuterR );
  //  std::printf("PWO_PosZ      = %.1f cm;\n", PWO_PosZ );
  std::printf("Towers in Row/Col   = %d;\n", towersInRow);
  std::printf("Top left tower pos  = %.5f, %.5f cm;\n", leftTowerPos , topTowerPos );
  std::printf("Top left tower pos  = %.5f, %.5f cm;\n", leftTowerPos , topTowerPos );


  int towerIndex = 0, count = 0;
  ofstream ofs;

  
  // =========================================
  // Generating the crystal position mapping
  // =========================================
  
  ofs.open("tower_map_purecrystal_185cm_EEEMCcarbon.txt");
  if (!ofs.is_open()) 
    cout << "Failed to open file." << endl;

  //The prefix parameter of the mapping
  ofs << "Gtype 1" << endl;
  ofs << "Gr1_inner " << PWO_InnerR << endl;
  ofs << "Gr1_outer " << PWO_OuterR << endl;
  ofs << "Gr2_inner " << PWO_Frame_InnerR << endl;
  ofs << "Gr2_outer " << PWO_Frame_OuterR << endl;
  ofs << "Gdx " << PWO_tower_width << endl;
  ofs << "Gdy " << PWO_tower_width << endl;
  ofs << "Gdz " << PWO_Thickness << endl;
  ofs << "Gx0 " << 0 << endl;
  ofs << "Gy0 " << 0 << endl;
  ofs << "Gz0 " << inner_glo_zpos << endl;
  ofs << "Grot_x " << 0 << endl;
  ofs << "Grot_y " << 0 << endl;
  ofs << "Grot_z " << 0 << endl;
  ofs << "Gcrystal_dx " << PWO_Width << endl;
  ofs << "Gcrystal_dy " << PWO_Width << endl;
  ofs << "Gcrystal_dz " << PWO_Thickness << endl;
  ofs << "Gcarbon_gap " << PWO_carbon_Gap << endl;
  ofs << "Greflective_foil_thickness " << PWO_reffoil_Gap << endl;
  ofs << "Gtedlar_thickness " << PWO_tedlar_Gap << endl;
  ofs << "carbon_frame_style " << PWO_frame_style << endl;
  ofs << "carbon_frame_depth " << PWO_frame_depth << endl;
  ofs << "carbon_face_lip " << PWO_face_lip << endl;
  ofs << "Gair_gap " << PWO_air_Gap << endl;
  ofs << "Gmaterial " << Inner_mat << endl;
  ofs << "Gcolor_R " << 0.0 << endl;
  ofs << "Gcolor_G " << 0.0 << endl;
  ofs << "Gcolor_B " << 1.0 << endl;
  ofs << "#Tower: type, idx_j, idx_k, idx_l, x[cm], y[cm], z[cm], dx[cm], dy[cm], dz[cm], rot_x, rot_y, rot_z" << endl;
  
  for(int colIndex = 0 ; colIndex < towersInRow ; colIndex++){
    for(int rowIndex = 0 ; rowIndex < towersInRow ; rowIndex++)	{
      double x = leftTowerPos + colIndex * PWO_tower_width;
      double y = topTowerPos + rowIndex * PWO_tower_width;
      double z = inner_glo_zpos;
      //	  double z = -250.;
      double Cry_R = std::sqrt(x * x + y * y);
      double inner_r = 0.;
      
      //	  if ( (std::abs(y) < PWO_OuterR && std::abs(x) < PWO_OuterR) && (std::abs(y) > PWO_InnerR || std::abs(x) > PWO_InnerR) )
      //	  if ( (Cry_R < PWO_OuterR) && (std::abs(y) > (PWO_InnerR + 2.) || std::abs(x) > (PWO_InnerR + 2.) ) )
      if ( (Cry_R < PWO_OuterR) && (Cry_R > PWO_InnerR) ){
        int index_dis = 0;
        bool accepted = true;
        for(int i = 0 ; i < 2 ; i++) {
          for(int j = 0 ; j < 2 ; j++) {
            double x_cry, y_cry;

            if( (i==0) && (j==0) ) {
              x_cry = x - PWO_tower_width / 2.;
              y_cry = y - PWO_tower_width / 2.;
              inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
              if(inner_r < PWO_InnerR){
                accepted = false;
                continue;
              }
            } else if( (i==0) && (j==1) ) {
              x_cry = x - PWO_tower_width / 2.;
              y_cry = y + PWO_tower_width / 2.;
              inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
              if(inner_r < PWO_InnerR){
                accepted = false;
                continue;
              }
            } else if( (i==1) && (j==0) ){
              x_cry = x + PWO_tower_width / 2.;
              y_cry = y - PWO_tower_width / 2.;
              inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
              if(inner_r < PWO_InnerR){
                accepted = false;
                continue;
              }
            } else if( (i==1) && (j==1) ) {
              x_cry = x + PWO_tower_width / 2.;
              y_cry = y + PWO_tower_width / 2.;
              inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
              if(inner_r < PWO_InnerR){
                accepted = false;
                continue;
              }
            }
          }
        }
        if(!accepted) continue;
        // int code = 1000 * cry_col_index + cry_row_index;

        ofs << "Tower " << 0 << " ";
        ofs << colIndex << " " << rowIndex << " " << 0 << " ";
        ofs << x << " " << y << " " << z << " ";
        ofs << PWO_tower_width << " " << PWO_tower_width << " " << PWO_Thickness << " ";
        ofs << 0 << " " << 0 << " " << 0 << endl;
      
        towerIndex++;
        count++;
            
      } else if (Cry_R <= PWO_InnerR) { // inner range 
        cout << "rejected inner: "<<  colIndex << "\t"<< rowIndex << "\t" << Cry_R << endl;
      } else if (Cry_R >= PWO_OuterR) { // inner range 
        cout << "rejected outer: "<<  colIndex << "\t"<< rowIndex << "\t" << Cry_R << endl;
        
      }
    }
  }

  ofs.close();
  towerIndex = 0;
  cout << "Inner use : " << count << " towers to construct!!!" << endl << endl;
  count = 0;



}
