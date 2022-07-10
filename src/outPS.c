/*
 
 <outPS.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 functions for writing PostScript file.

 <PostScript Basics>

 *COORDINATE

  Y
 /| 
  |
  |
  |
  o--------> X

  Unit of length(point) : 1/72 inch
  1 inch = 25.4 mm = 2.54 cm

  1 point  = 1/72 inch = 0.35 mm

 *PAPER

 The default paper is "Letter" size.

 Letter size:  
    X: 215.9 mm  :  8.5 inch -> 612 
    Y: 279.4 mm  : 11.0 inch -> 792 
 
 A4 paper:
    X: 210mm : 8.27 inch  -> 595.44  
    Y: 297mm : 11.69 inch -> 841.68

 
 *LINE
 
 newpath
 [start_x] [start_y] moveto
 [end_x]   [ned_y]   lineto
 stroke

 *FILLED CIRCLE
 
 newpath
 [center_x] [center_y] [radius] 0 360 arc
 closepath
 fill
 stroke

 *TEXT

 /Helvetica findfont
 [font_size (point)] scalefont
 setfont
 [pos_x] [pos_y] moveto
 ([Text String]) show


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "ioSDF.h"
#include "match.h"
#include "qRMS.h"
#include "2dRMS.h"
#include "molprop.h"
#include "vector3.h"
#include "options.h"
#include "PCAfit.h"


/*** FUNCTIONS (GLOBAL) ***/
void Write_Molecule_in_PostScript();
void Write_MATCH_in_PostScript();

/*** FUNCTIONS (LOCAL) ***/
static void write_one_molecule_in_postscript();
static void get_xyps_from_atom_pos();
static void get_min_max_XY_molecule();
static void get_most_empty_2D_point();
static void RGB_from_element();
static float get_scale();

void Write_Molecule_in_PostScript(ofname,mol)
  char *ofname;
  struct MOLECULE *mol;
{
  FILE *fpo;
  float paperXi,paperYi;     /* size of paper (inch) */
  float scale;

  printf("#Write_Molecule_in_PostScript()-->'%s'\n",ofname);
  fpo = fopen(ofname,"w");
  paperXi =  8.5;
  paperYi = 11.0;
  scale = get_scale(mol,paperXi/2,paperYi);
  write_one_molecule_in_postscript(fpo,mol,0.0,paperYi/2,paperXi/2,paperYi,scale);
  fprintf(fpo,"/Helvetica findfont\n");
  fprintf(fpo,"8 scalefont\n");
  fprintf(fpo,"setfont\n");
  fprintf(fpo,"0 0 0 setrgbcolor\n");
  fprintf(fpo,"%f %f moveto\n",72*0.1,      72*paperYi/2);
  fprintf(fpo,"(%s) show\n",mol->name);
  fprintf(fpo,"showpage\n"); 
  fclose(fpo);

} /* end of Write_Molecule_in_PostScript() */



void Write_MATCH_in_PostScript(ofname,m,molA,molB)
  char *ofname;
  struct MATCH *m;
  struct MOLECULE *molA,*molB;
{
  FILE *fpo;
  float paperXi,paperYi;     /* size of paper (inch) */
  int a,b,i; 
  double g1[3],g2[3],Rmat[3][3];
  float scaleA,scaleB,scale; 
  printf("#Write_MATCH_in_PostScript()-->'%s'\n",ofname);
  paperXi =  8.5;
  paperYi = 11.0;

  /** [1] Rotate molB in XY-flat, and superimpose molA onto molB  **/ 
  if (PAR.ps_planeB=='T'){ Rotate_Molecule_XY_flat_plane_by_PCA(molB);}
  
  if (PAR.ps_supAonB=='T'){
    Calculate_CRMS_MATCH_rotation_Z(m,molA,molB,g1,g2,Rmat,"AonB");
    Rotate_Molecule(molA,g1,g2,Rmat);
  }

  /** [2] Mark matched atoms **/ 
  for (a=0;a<molA->Natom;++a){molA->atoms[a].mark = molA->atoms[a].subnum = 0;}
  for (b=0;b<molB->Natom;++b){molB->atoms[b].mark = molB->atoms[b].subnum = 0;}

  /* mark   :  matched atom number */
  /* subnum :  label nuber */

  for (i=0;i<m->Npair;++i){
    molA->atoms[m->anumA[i]].mark = i+1;
    molB->atoms[m->anumB[i]].mark = i+1;
  }

 
  if (PAR.ps_label == 'M'){ 
    for (a=0;a<molA->Natom;++a){ molA->atoms[a].subnum = molA->atoms[a].mark;}
    for (b=0;b<molB->Natom;++b){ molB->atoms[b].subnum = molB->atoms[b].mark;}
  }

  if (PAR.ps_label == 'A'){ 
    /*
    for (a=0;a<molA->Natom;++a) molA->atoms[a].subnum =  molA->atoms[a].num_in_file;
    */
    for (i=0;i<m->Npair;++i){
      molA->atoms[m->anumA[i]].subnum =  molA->atoms[m->anumA[i]].num_in_file;
      molB->atoms[m->anumB[i]].subnum =  molA->atoms[m->anumA[i]].num_in_file;
    }
  }

  if (PAR.ps_label == 'B'){ 
    for (b=0;b<molB->Natom;++b){ molB->atoms[b].subnum =  molB->atoms[b].num_in_file;}
    for (i=0;i<m->Npair;++i){
      molA->atoms[m->anumA[i]].subnum =  molB->atoms[m->anumB[i]].num_in_file;
    }
  }


  /** [3] Write PostScript file **/ 
  scaleA = get_scale(molA,paperXi/2,paperYi);
  scaleB = get_scale(molB,paperXi/2,paperYi);
  if (scaleA<scaleB) {scale = scaleA;} else {scale = scaleB;}
  printf("#scale molA %f molB %f --> %f\n" ,scaleA,scaleB,scale);
  fpo = fopen(ofname,"w");
  write_one_molecule_in_postscript(fpo,molA,0.0,      paperYi/2,paperXi/2,paperYi,scale);
  write_one_molecule_in_postscript(fpo,molB,paperXi/2,paperYi/2,paperXi/2,paperYi,scale);

  fprintf(fpo,"/Helvetica findfont\n");
  fprintf(fpo,"8 scalefont\n");
  fprintf(fpo,"setfont\n");
  fprintf(fpo,"0 0 0 setrgbcolor\n");
  fprintf(fpo,"%f %f moveto\n",72*0.1,      72*paperYi/2);
  fprintf(fpo,"(%s) show\n",molA->name);
  fprintf(fpo,"%f %f moveto\n",72*(paperXi/2+0.1),72*paperYi/2);
  fprintf(fpo,"(%s) show\n",molB->name);

  fprintf(fpo,"showpage\n"); 
  fclose(fpo);

} /* end of Write_MATCH_in_PostScript() */







void write_one_molecule_in_postscript(fpo,mol,oriXi,oriYi,widthXi,widthYi,sc)
  FILE *fpo;
  struct MOLECULE *mol;
  float oriXi,  oriYi;   /* origin (left bottom ) point in paper (inch) */
  float widthXi,widthYi; /* width of papers     (inch) */
  float sc;              /* scale (inch/angstrom) */
{
  int a,b;
  float minXa,maxXa,minYa,maxYa; /* min, max of molecule (angstrom) */
  float oriXa,oriYa;             /* origin of molecule (angstrom)   */ 
  float margin_a;       /* margin of molecule (angstrom) */
  float xps,yps;                 /* 2D position in post script (72 inch) */
  float radius,pos[3],posA[3],posB[3],bvec[2],nvec[2],len;
  float double_bondWa; /* (angstrom) */
  float R,G,B;
  int text_pointsize;
  char string[128];
  float sc_circleL,sc_circleS;

  margin_a = 2.0;
  double_bondWa = 0.08; 
  text_pointsize = 8;
  pos[2] = posA[2] = posB[2] = 0.0;

  get_min_max_XY_molecule(mol,&minXa,&maxXa,&minYa,&maxYa);
  oriXa = minXa - margin_a;
  oriYa = minYa - margin_a;
 
  sc_circleL = 0.35;
  sc_circleS = 0.25;


/*
  printf("#minXa %f maxXa %f minYa %f maxYa %f scX %f scY %f sc %f\n",minXa,maxXa,minYa,maxYa,scX,scY,sc);
*/

  /*** [1] Draw Bond Lines ***/
  for (a=0;a<mol->Natom;++a){
    for (b=a+1;b<mol->Natom;++b){
      if ((mol->conmap.map[a][b]!='0') &&
          ((PAR.ps_hydro=='T')||(mol->atoms[a].one_char_ele!='H')) &&
          ((PAR.ps_hydro=='T')||(mol->atoms[b].one_char_ele!='H')) ){
        bvec[0] = mol->atoms[b].Pos[0] - mol->atoms[a].Pos[0];
        bvec[1] = mol->atoms[b].Pos[1] - mol->atoms[a].Pos[1];
        len = sqrt(bvec[0]*bvec[0] + bvec[1]*bvec[1]);
        nvec[0] =  bvec[1]/len;
        nvec[1] = -bvec[0]/len;

        if ((PAR.ps_care_bondtype=='T')&&(mol->conmap.map[a][b]=='2')){
          posA[0] = mol->atoms[a].Pos[0] + nvec[0] * double_bondWa;
          posA[1] = mol->atoms[a].Pos[1] + nvec[1] * double_bondWa;
          posB[0] = mol->atoms[b].Pos[0] + nvec[0] * double_bondWa;
          posB[1] = mol->atoms[b].Pos[1] + nvec[1] * double_bondWa;
          fprintf(fpo,"newpath\n");
          fprintf(fpo,"0 0 0 setrgbcolor\n");
          get_xyps_from_atom_pos(&xps,&yps,posA,sc,oriXa,oriYa,oriXi,oriYi);
          fprintf(fpo,"%f %f moveto\n",xps,yps);
          get_xyps_from_atom_pos(&xps,&yps,posB,sc,oriXa,oriYa,oriXi,oriYi);
          fprintf(fpo,"%f %f lineto\n",xps,yps);
          fprintf(fpo,"stroke\n");

          posA[0] = mol->atoms[a].Pos[0] - nvec[0] * double_bondWa;
          posA[1] = mol->atoms[a].Pos[1] - nvec[1] * double_bondWa;
          posB[0] = mol->atoms[b].Pos[0] - nvec[0] * double_bondWa;
          posB[1] = mol->atoms[b].Pos[1] - nvec[1] * double_bondWa;
          fprintf(fpo,"newpath\n");
          fprintf(fpo,"0 0 0 setrgbcolor\n");
          get_xyps_from_atom_pos(&xps,&yps,posA,sc,oriXa,oriYa,oriXi,oriYi);
          fprintf(fpo,"%f %f moveto\n",xps,yps);
          get_xyps_from_atom_pos(&xps,&yps,posB,sc,oriXa,oriYa,oriXi,oriYi);
          fprintf(fpo,"%f %f lineto\n",xps,yps);
          fprintf(fpo,"stroke\n");
        }
        else{
          fprintf(fpo,"newpath\n");
          fprintf(fpo,"0 0 0 setrgbcolor\n");
          get_xyps_from_atom_pos(&xps,&yps,mol->atoms[a].Pos,sc,oriXa,oriYa,oriXi,oriYi);
          fprintf(fpo,"%f %f moveto\n",xps,yps);
          get_xyps_from_atom_pos(&xps,&yps,mol->atoms[b].Pos,sc,oriXa,oriYa,oriXi,oriYi);
          fprintf(fpo,"%f %f lineto\n",xps,yps);
          fprintf(fpo,"stroke\n");
        }
      }
    } 
  }

 
  /*** [2] Draw Atom Circles ***/
  for (a=0;a<mol->Natom;++a){
    if ((PAR.ps_hydro=='T')||(mol->atoms[a].one_char_ele!='H')){
      get_xyps_from_atom_pos(&xps,&yps,mol->atoms[a].Pos,sc,oriXa,oriYa,oriXi,oriYi);
      if (PAR.ps_circle=='M'){
        if (mol->atoms[a].mark>0){ radius = 72 * sc * sc_circleL; }
                             else{ radius = 72 * sc * sc_circleS; }
      }
      else if (PAR.ps_circle=='H'){
        if (mol->atoms[a].one_char_ele == 'H'){ radius = 36 * sc * sc_circleL; }
                                        else  { radius = 72 * sc * sc_circleS; }
      }
      else if (PAR.ps_circle=='A'){ radius = 72 * sc * sc_circleL; } 
      else                        { radius = 72 * sc * sc_circleS; } 

      fprintf(fpo,"newpath\n");
      /* fprintf(fpo,"0 0 0 setrgbcolor\n"); */
      fprintf(fpo,"%f %f %f 0.0 360.0 arc\n",xps,yps,radius);
      fprintf(fpo,"closepath\n");
      RGB_from_element(mol->atoms[a].element,&R,&G,&B);
      fprintf(fpo,"%f %f %f setrgbcolor\n",R,G,B);
      fprintf(fpo,"fill\n");
      fprintf(fpo,"stroke\n");
      if (PAR.ps_circ_outline=='T'){
        fprintf(fpo,"newpath\n");
        fprintf(fpo,"%f %f %f 0.0 360.0 arc\n",xps,yps,radius);
        fprintf(fpo,"closepath\n");
        RGB_from_element(mol->atoms[a].element,&R,&G,&B);
        fprintf(fpo,"%f %f %f setrgbcolor\n",0.0,0.0,0.0);
        fprintf(fpo,"stroke\n");
      }
    }
  } 

  /*** [3] Draw Text on the atoms ***/
  fprintf(fpo,"/Helvetica findfont\n");
  fprintf(fpo,"%d scalefont\n",text_pointsize);
  fprintf(fpo,"setfont\n");
  fprintf(fpo,"0 0 0 setrgbcolor\n");
  for (a=0;a<mol->Natom;++a){
    if ((PAR.ps_hydro=='T')||(mol->atoms[a].one_char_ele!='H')){
       string[0] = '\0';
      if (((PAR.ps_label=='M')||(PAR.ps_label=='A')||(PAR.ps_label=='B')) && (mol->atoms[a].subnum>0)) sprintf(string,"%d",mol->atoms[a].subnum);
      if (PAR.ps_label=='a')  sprintf(string,"%s",mol->atoms[a].atomname);
      if (PAR.ps_label=='E')  sprintf(string,"%s",mol->atoms[a].element);
      if (PAR.ps_label=='T')  sprintf(string,"%s",mol->atoms[a].atomtype);
      if (PAR.ps_label=='N')  sprintf(string,"%d",mol->atoms[a].num_in_file);

      if (strlen(string)>0){
        get_most_empty_2D_point(&(bvec[0]),&(bvec[1]),mol,a);
        pos[0] = mol->atoms[a].Pos[0] + 0.6 * bvec[0];
        pos[1] = mol->atoms[a].Pos[1] + 0.6 * bvec[1];
        pos[2] = 0.0;
        get_xyps_from_atom_pos(&xps,&yps,pos,sc,oriXa,oriYa,oriXi,oriYi);
        xps = xps - strlen(string)/2.0 * text_pointsize/2.0;
        yps = yps - text_pointsize/2.0;
        fprintf(fpo,"%f %f moveto\n",xps,yps);
        fprintf(fpo,"(%s) show\n",string);
      }
    }
  }


} /* end of write_one_molecule_in_postscript() */


void get_most_empty_2D_point(empXa,empYa,mol,anum)
  float *empXa,*empYa; /* X,Y coordinate of empty points (angstrom; to be calculated) */
  struct MOLECULE *mol;
  int anum;  /* atom number for the focused atom */
{
  int i,j,Nnei,Ndiv;
  float *neiX,*neiY,theta,posX,posY,DD,minDD,max_minDD;
  float len;

  Ndiv = 12;
  Nnei = mol->atoms[anum].Nneighbor;

  neiX = (float *)malloc(sizeof(float)*mol->atoms[anum].Nneighbor);
  neiY = (float *)malloc(sizeof(float)*mol->atoms[anum].Nneighbor);

  Nnei = 0;
  for (i=0;i<mol->atoms[anum].Nneighbor;++i){
     j = mol->atoms[anum].neighbors[i];
     if ((PAR.ps_hydro=='T')||(mol->atoms[j].one_char_ele!='H')){
       neiX[Nnei] = mol->atoms[j].Pos[0] - mol->atoms[anum].Pos[0];
       neiY[Nnei] = mol->atoms[j].Pos[1] - mol->atoms[anum].Pos[1];
       len = neiX[Nnei]*neiX[Nnei] + neiY[Nnei]*neiY[Nnei];
       if (len>0.0) len = sqrt(len);
       neiX[Nnei] /= len;
       neiY[Nnei] /= len;
       Nnei += 1;
     }
  }

  max_minDD = -1.0;
  for (i=0;i<Ndiv;++i){
    theta = 2.0*M_PI*i/Ndiv; 
    posX = cos(theta);
    posY = sin(theta);

    minDD = -1.0;
    for (j=0;j<Nnei;++j){
      DD = (neiX[j]-posX)*(neiX[j]-posX) +  (neiY[j]-posY)*(neiY[j]-posY);
      if ((minDD<0.0)||(DD<minDD)) minDD = DD;
    }

    if ((max_minDD<0.0)||(minDD>max_minDD)){
      max_minDD = minDD;
      *empXa = posX;
      *empYa = posY;
    } 
  }

  free(neiX);
  free(neiY);

} /* end of get_most_empty_2D_point() */



void RGB_from_element(ele,r,g,b)
  char *ele;
  float *r, *g, *b;
{
   if (PAR.ps_color=='M'){
         if (strcmp(ele,"C")==0) { *r = 1.00; *g = 1.00; *b = 1.00;}
    else if (strcmp(ele,"N")==0) { *r = 0.70; *g = 0.70; *b = 0.70;}
    else if (strcmp(ele,"O")==0) { *r = 0.00; *g = 0.00; *b = 0.00;}
    else { *r = 0.80; *g = 0.80; *b = 0.80;}
   }
  else{ 
 /*
  These colors are taken from RasMol manual.
 */
         if (strcmp(ele,"C")==0) { *r = 0.40; *g = 0.40; *b = 0.40;}
    /* else if (strcmp(ele,"N")==0) { *r = 0.40; *g = 0.40; *b = 1.00;} */
    else if (strcmp(ele,"N")==0) { *r = 0.00; *g = 0.00; *b = 1.00;} 
    else if (strcmp(ele,"O")==0) { *r = 1.00; *g = 0.00; *b = 0.00;}
    else if (strcmp(ele,"H")==0) { *r = 0.80; *g = 0.80; *b = 0.80;}
    else if (strcmp(ele,"S")==0) { *r = 0.90; *g = 0.80; *b = 0.00;}
    else if (strcmp(ele,"P")==0) { *r = 1.00; *g = 0.66; *b = 0.00;}
    else if (strcmp(ele,"F")==0) { *r = 0.85; *g = 0.64; *b = 0.13;}
    else if (strcmp(ele,"I")==0) { *r = 0.62; *g = 0.12; *b = 0.94;}
    else if (strcmp(ele,"B")==0) { *r = 0.00; *g = 1.00; *b = 0.0;}
    else if ((strcmp(ele,"HE")==0)||(strcmp(ele,"He")==0)) { *r = 1.00; *g = 0.75; *b = 0.79;}
    else if ((strcmp(ele,"NA")==0)||(strcmp(ele,"Na")==0)) { *r = 0.00; *g = 0.00; *b = 1.00;}
    else if ((strcmp(ele,"LI")==0)||(strcmp(ele,"Li")==0)) { *r = 0.69; *g = 0.13; *b = 0.13;}
    else if ((strcmp(ele,"CL")==0)||(strcmp(ele,"Cl")==0)) { *r = 0.00; *g = 1.00; *b = 0.00;}
    else if ((strcmp(ele,"BR")==0)||(strcmp(ele,"Br")==0)) { *r = 0.50; *g = 0.15; *b = 0.15;}
    else if ((strcmp(ele,"MG")==0)||(strcmp(ele,"Mg")==0)) { *r = 0.13; *g = 0.54; *b = 0.13;}
    else if ((strcmp(ele,"FE")==0)||(strcmp(ele,"Fe")==0)) { *r = 1.00; *g = 0.66; *b = 0.00;}
    else if ((strcmp(ele,"BA")==0)||(strcmp(ele,"Ba")==0)) { *r = 1.00; *g = 0.66; *b = 0.00;}
    else if ((strcmp(ele,"CA")==0)||(strcmp(ele,"Ca")==0)) { *r = 0.41; *g = 0.41; *b = 0.41;}
    else if ((strcmp(ele,"MN")==0)||(strcmp(ele,"Mn")==0)) { *r = 0.41; *g = 0.41; *b = 0.41;}
    else if ((strcmp(ele,"AL")==0)||(strcmp(ele,"Al")==0)) { *r = 0.41; *g = 0.41; *b = 0.41;}
    else if ((strcmp(ele,"TI")==0)||(strcmp(ele,"Ti")==0)) { *r = 0.41; *g = 0.41; *b = 0.41;}
    else if ((strcmp(ele,"CR")==0)||(strcmp(ele,"Cr")==0)) { *r = 0.41; *g = 0.41; *b = 0.41;}
    else if ((strcmp(ele,"AG")==0)||(strcmp(ele,"Ag")==0)) { *r = 0.41; *g = 0.41; *b = 0.41;}
    else if ((strcmp(ele,"NI")==0)||(strcmp(ele,"Ni")==0)) { *r = 0.50; *g = 0.15; *b = 0.15;}
    else if ((strcmp(ele,"CU")==0)||(strcmp(ele,"Cu")==0)) { *r = 0.50; *g = 0.15; *b = 0.15;}
    else if ((strcmp(ele,"ZN")==0)||(strcmp(ele,"Zn")==0)) { *r = 0.50; *g = 0.15; *b = 0.15;}
    else { *r = 0.98; *g = 0.08; *b = 0.56;}
  }
}




void get_xyps_from_atom_pos(Xps,Yps,AtomPos,scale,oriXa,oriYa,oriXi,oriYi)
  float *Xps,*Yps;   /* 2D coordiantes in postscript (72 inch ) */
  float AtomPos[3];     /* atom position (angstrom) */
  float scale;       /* scale (inch/angstrom)  */
  float oriXa,oriYa; /* origin position (angstrom) */
  float oriXi,oriYi; /* origin position (inch) */
{
  *Xps = 72*(scale*(AtomPos[0]-oriXa)+oriXi);
  *Yps = 72*(scale*(AtomPos[1]-oriYa)+oriYi);
}


void get_min_max_XY_molecule(mol,minXa,maxXa,minYa,maxYa)
  struct MOLECULE *mol;
  float *minXa,*maxXa;
  float *minYa,*maxYa;
{
  int a;
  
  *minXa = *maxXa = mol->atoms[0].Pos[0]; 
  *minYa = *maxYa = mol->atoms[0].Pos[1]; 

  for (a=1;a<mol->Natom;++a){
    if (mol->atoms[a].Pos[0] < (*minXa)){*minXa = mol->atoms[a].Pos[0];}
    if (mol->atoms[a].Pos[0] > (*maxXa)){*maxXa = mol->atoms[a].Pos[0];}
    if (mol->atoms[a].Pos[1] < (*minYa)){*minYa = mol->atoms[a].Pos[1];}
    if (mol->atoms[a].Pos[1] > (*maxYa)){*maxYa = mol->atoms[a].Pos[1];}
  }

} /* end of get_min_max_XY_molecule() */



float get_scale(mol,widthXi,widthYi)
  struct MOLECULE *mol;
  float widthXi,widthYi; /* width of papers     (inch) */
{
  float scX,scY,sc;
  float minXa,maxXa,minYa,maxYa; /* min, max of molecule (angstrom) */
  float margin_a;       /* margin of molecule (inch) */
  float margin_i;       /* margin of paper (inch) */


  margin_i = 0.5;
  margin_a = 2.0;
  get_min_max_XY_molecule(mol,&minXa,&maxXa,&minYa,&maxYa);
  scX = (widthXi - 2.0*margin_i)/(maxXa-minXa+2*margin_a);
  scY = (widthYi - 2.0*margin_i)/(maxYa-minYa+2*margin_a);
  if (scX<scY) sc = scX; else sc = scY;
  return(sc);
}
