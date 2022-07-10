/*

< steep_descen.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


functions for steepest decscent method
by energy minimizatioin procedure.


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string.h>
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "globalvar.h"
#include "molprop.h"
#include "transform.h"
#include "vector3.h"
#include "energies.h"
#include "stamp_transf.h"
#include "PCAfit.h"
#include "qRMS.h"

/*** FUNCTIONS (GLOBAL) ***/
int  steepest_descent_golden_linear_search();
void steepest_descent_search_using_only_Eselfclash();


/*** FUNCTIONS (LOCAL) ***/
static void Transform_MOLECULE_by_Trans_Rot_RotBond_with_step();





void Transform_MOLECULE_by_Trans_Rot_RotBond_with_step(tra,mol,step,FixStamp,FixPlane)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *mol;
  float  step;  /* step > 0.0 */
  char   FixStamp;         /* 'F':fixing stamped angles 'R'otating stamped angle. */
  char   FixPlane;         /* 'F':fixing plane angles   'R'otating plane angle. */
{
  float dT[3], rotangle,rbangle;
  int r,k;

  for (k=0;k<3;++k){
    dT[k] = step * PAR.maxTRA * tra->Force[k];
  }

  rotangle = step * M_PI*PAR.maxROT/180.0; 
  Transform_MOLECULE_by_Translation_Rotation(dT,tra->Torque,rotangle,tra,mol);

  for (r=0;r<tra->Nrbond;++r){
     /* printf("[RBOND %d] rbStamp %c rbPlane %c\n",r,tra->rbStamp[r],tra->rbPlane[r]); */
    if ( ((FixStamp=='R')||(tra->rbStamp[r]!='S'))  &&
         ((FixPlane=='R')||(tra->rbPlane[r]!='P'))  ){
      /*
      rbangle = -step * tra->maxRB * tra->Nrbond * tra->dEd_rbAngle[r];
      */
      rbangle = -step * M_PI * PAR.maxRB/180.0  *  tra->dEd_rbAngle[r];
      /* printf("[RBOND %d] dEd_rbAngle %f rbangle %f\n",r,tra->dEd_rbAngle[r],rbangle); */
      Transform_MOLECULE_by_RotBond(r,rbangle,tra,mol);
    }
  }

} /* end of Transform_MOLECULE_by_Trans_Rot_RotBond_with_step() */






int steepest_descent_golden_linear_search(traT, molT, molR, molP,M,FixStamp,FixPlane)
  struct TRANSFORM_VARIABLE *traT;
  struct MOLECULE *molT, *molR, *molP;
  struct MATCH *M;
  char   FixStamp;         /* 'F':fixing stamped angles 'R'otating stamped angle. */
  char   FixPlane;         /* 'F':fixing plane angles   'R'otating plane angle. */
{
  float gold;
  struct MOLECULE ma,mb,mc,mx;
  float a, b, c, x;
  float Ea, Eb, Ec, Ex;
  char init;
  char shorter_edge;  /* 'a':ab is shorter than bc, 'c':bc is shorter than ab */
  int  nrep, Nrepeat_max;

 /*
  printf("#golden_linear_search(traT, molT, molR,M)\n");
 */
  gold = 0.38196601125010515180;  /* (3-sqrt(5))/2) */

  Initialize_MOLECULE(&ma);
  Initialize_MOLECULE(&mb);
  Initialize_MOLECULE(&mc);
  Initialize_MOLECULE(&mx);

  Malloc_MOLECULE(&ma,molT->Natom,molT->Nbond); 
  Malloc_MOLECULE(&mb,molT->Natom,molT->Nbond); 
  Malloc_MOLECULE(&mc,molT->Natom,molT->Nbond); 
  Malloc_MOLECULE(&mx,molT->Natom,molT->Nbond); 

  /** [1] Find Initial (a,b,c) **/
  /*
  [case 'G']: --> Do golden division
   Ea         |        Ec
          Ec  | Ea
      Eb      |     Eb

  [case 'F']: --> finish. Take c and escape from this function. 
   Ea         |    Eb
      Eb      | Ea
          Ec  |       Ec

  [case 'I']: try again by narrowing the region (a',b',c') = (a,x,b)
         Ec   |    Eb  
      Eb      |       Ec 
   Ea         | Ea 

  */ 
  Nrepeat_max = 10; 
  a = 0.0; b = gold; c = 1.0; 
  shorter_edge = 'a';
  init = 'N';
  Copy_MOLECULE(&ma,molT);
  Transform_MOLECULE_by_Trans_Rot_RotBond_with_step(traT,&ma,a,FixStamp,FixPlane);
  Ea = energy_total(&ma,molR,molP,M);

  Copy_MOLECULE(&mb,molT);
  Transform_MOLECULE_by_Trans_Rot_RotBond_with_step(traT,&mb,b,FixStamp,FixPlane);
  Eb = energy_total(&mb,molR,molP,M);
  
  Copy_MOLECULE(&mc,molT);
  Transform_MOLECULE_by_Trans_Rot_RotBond_with_step(traT,&mc,c,FixStamp,FixPlane);
  Ec = energy_total(&mc,molR,molP,M);
  nrep = 0;
  do{
    /* printf("I[%d] a %f b %f c %f Ea %f Eb %f Ec %f\n",nrep,a,b,c,Ea,Eb,Ec); */
 
    if ((Eb<Ea) && (Eb<Ec)){
      init = 'G';      /* try golden division */
    }

    if ((Ec<Ea) && (Ec<Eb)){
      init = 'F';      /* finish everything. (escape immediately) */
      Copy_MOLECULE(molT,&mc);
    }

    if ((Ea<Eb) && (Ea<Ec)){
      c = b;
      Copy_MOLECULE(&mc,&mb);
      Ec = Eb;
      b = c * gold;
      Copy_MOLECULE(&mb,molT);
      Transform_MOLECULE_by_Trans_Rot_RotBond_with_step(traT,&mb,b,FixStamp,FixPlane);
      Eb = energy_total(&mb,molR,molP,M);
    }

    nrep += 1;

  } while ((init=='N') && (nrep <Nrepeat_max));


  if (init=='F'){
    /*
    printf("#init F (take c ) (monotonous decreasing)\n"); fflush(stdout);
    */
    Copy_MOLECULE(molT,&mc);
    Free_MOLECULE(&mx);
    Free_MOLECULE(&mc);
    Free_MOLECULE(&mb);
    Free_MOLECULE(&ma);
    return(1);
  }
  else if (init!='G'){
    /*
    printf("#WOOPS (init '%c') (converged) !!\n",init);
     */
    return(0);
  }

  /** [2] Golden division (a,b,c) **/
  nrep = 0;
  x = 0.0; 
  do{
    if (shorter_edge=='a'){
      x = a+(c-a)*(1.0-gold);
    }
    else if (shorter_edge=='c'){
      x = a+(c-a)*gold;
    }
    Copy_MOLECULE(&mx,molT);
    Transform_MOLECULE_by_Trans_Rot_RotBond_with_step(traT,&mx,x,FixStamp,FixPlane);
    Ex = energy_total(&mx,molR,molP,M);

    /*
    printf("G[%d](%c) a %f b %f c %f x %f Ea %f Eb %f Ec %f Ex %f\n",nrep,shorter_edge,a,b,c,x,Ea,Eb,Ec,Ex);
     */

    if (shorter_edge=='a'){
      if (Eb<Ex){
        c = x;
        Ec = Ex;
        Copy_MOLECULE(&mc,&mx);
        shorter_edge = 'c';
      }
      else if (Ex<=Eb){
        a = b;
        Ea = Eb;
        Copy_MOLECULE(&ma,&mb);
        b = x;
        Eb = Ex;
        Copy_MOLECULE(&mb,&mx);
        shorter_edge = 'a';
      }
    }

    else if (shorter_edge=='c'){
      if (Eb<Ex){
        a = x;
        Ea = Ex;
        Copy_MOLECULE(&ma,&mx);
        shorter_edge = 'a';
      }
      else if (Ex<=Eb){
        c = b;
        Ec = Eb;
        Copy_MOLECULE(&mc,&mb);
        b = x;
        Eb = Ex;
        Copy_MOLECULE(&mb,&mx);
        shorter_edge = 'c';
      }
    }

    nrep += 1;
  } while (nrep<Nrepeat_max);

  /*** (3) Ending procedure ***/
 /*
 printf("G[%d](%c) a %f b %f c %f x %f Ea %f Eb %f Ec %f Ex %f\n",nrep,shorter_edge,a,b,c,x,Ea,Eb,Ec,Ex);
  */
       if ((Ea<=Eb)&&(Ea<=Ec)){ Copy_MOLECULE(molT,&ma); }
  else if ((Eb<=Ea)&&(Eb<=Ec)){ Copy_MOLECULE(molT,&mb); }
  else if ((Ec<=Ea)&&(Ec<=Eb)){ Copy_MOLECULE(molT,&mc); }
  
  Free_MOLECULE(&mx);
  Free_MOLECULE(&mc);
  Free_MOLECULE(&mb);
  Free_MOLECULE(&ma);
  return(2);

} /* end of steepest_descent_golden_linear_search() */





void steepest_descent_search_using_only_Eselfclash(Niterate,traT,molT,molR,molP,M,FixStamp,FixPlane)
  int    Niterate;
  struct TRANSFORM *traT;
  struct MOLECULE *molT; /* target molecule */
  struct MOLECULE *molR; /* reference molecule (actually, not used in this function)*/
  struct MOLECULE *molP; /* protein receptor molecule (actually, not used in this function)*/
  struct MATCH    *M;    /* MATCH (actually, not used in this function) */
  char   FixStamp;         /* 'F':fixing stamped angles 'R'otating stamped angle. */
  char   FixPlane;         /* 'F':fixing plane angles   'R'otating plane angle. */
{
  int k,ok;
  float WEatommatch_orig;
  float WEselfclash_orig;
  float WEprotclash_orig;
  float WEprotatrct_orig;
  float WEvolmovlap_orig;
  float WEtpdisrest_orig;

  printf("#steepest_descent_search_using_only_Eselfclash(Niterate:%d traT,molT,molR,molP,M,FixStamp,FixPlane)\n",Niterate);
  WEatommatch_orig = PAR.WEatommatch;
  WEselfclash_orig = PAR.WEselfclash;
  WEprotclash_orig = PAR.WEprotclash;
  WEprotatrct_orig = PAR.WEprotatrct;
  WEvolmovlap_orig = PAR.WEvolmovlap;
  WEtpdisrest_orig = PAR.WEtpdisrest;

  PAR.WEatommatch = 0.0;
  PAR.WEselfclash = 1.0;
  PAR.WEprotclash = 0.0;
  PAR.WEprotatrct = 0.0;
  PAR.WEvolmovlap = 0.0;
  PAR.WEtpdisrest = 0.0;

  for (k=0;k<Niterate;++k){
    cal_Force_and_Torque(traT, molT, molR, molP,M);
    cal_dEd_rbAngle(traT,molT,molR, molP, M);
    ok = steepest_descent_golden_linear_search(traT, molT, molR, molP, M, FixStamp,FixPlane);
    if (ok==0) k = 1000000;
    /*
    RMSD = Calculate_CRMS_MATCH_Quaternion(m_model,&(molTtrial[nn]),&molR,g1,g2,Rmat,"AonB");
    Rotate_Molecule(molT,g1,g2,Rmat);
    */
  } 

  PAR.WEatommatch = WEatommatch_orig;
  PAR.WEselfclash = WEselfclash_orig;
  PAR.WEprotclash = WEprotclash_orig;
  PAR.WEprotatrct = WEprotatrct_orig;
  PAR.WEvolmovlap = WEvolmovlap_orig;
  PAR.WEtpdisrest = WEtpdisrest_orig;

} /* end of steepest_descent_search_using_only_Eselfclash() */

