#!/usr/bin/env python
import sys
import molecule
import math

def cal_distance_bwn_heavy_atoms(mol):
  Dmin = 10000.0
  Rvdw = {}
  Rvdw['C'] = 1.87 
  Rvdw['O'] = 1.40 
  Rvdw['N'] = 1.65 
  str = ''
  for a in range (mol.Natom):
    if (mol.atoms[a].element !='H'):
      #print "%s '%s'"%(mol.atoms[a].atomname,mol.atoms[a].element)
      for b in range (a+1,mol.Natom):
        if (mol.atoms[b].element !='H') and (mol.contable[a][b]=='0'):
          hit = 0
          for n in (mol.atoms[a].neighbors):
            for m in (mol.atoms[n].neighbors):
              if (m==b):
                hit = 1
          if (hit==0):
            dx = mol.atoms[a].x - mol.atoms[b].x
            dy = mol.atoms[a].y - mol.atoms[b].y
            dz = mol.atoms[a].z - mol.atoms[b].z
            D = math.sqrt(dx*dx + dy*dy + dz*dz)
            Dvdw = D - Rvdw[mol.atoms[a].element] - Rvdw[mol.atoms[b].element]  
          #print "D %f Dvdw %f r %f %f"%(D,Dvdw,Rvdw[mol.atoms[a].element],Rvdw[mol.atoms[b].element])  
            if (Dvdw<Dmin):
              Dmin = Dvdw
              str = "%d '%s' vs %s '%s'"%(mol.atoms[a].num_in_file,mol.atoms[a].atomname,mol.atoms[b].num_in_file,mol.atoms[b].atomname) 
  return([Dmin,str])


##############
#### MAIN ####
##############


if (len(sys.argv)<2):
  print "phipsi.py [sdf or pdb or kcf or mol2 file]"
  sys.exit(1)
m = molecule.Molecule()
if (sys.argv[1].endswith(".sdf")):
  m.read_in_sdf(sys.argv[1])
elif (sys.argv[1].endswith(".pdb")):
  m.read_in_pdb(sys.argv[1])
  m.guess_bonds_from_atom_xyz()
else:
  print "#no input file"
  sys.exit(1)

print "Natom %d"%(m.Natom)
m.count_neighbor_atoms()

m.write_in_pdb("out0.pdb")
m.set_internal_coodinate_topology()
m.set_internal_coodinate_angle()
m.write_in_mopin("out0.mopin")

m.set_xyz_coodinates_from_internal()
m.write_in_pdb("out1.pdb")

m.set_internal_coodinate_angle()
m.write_in_mopin("out2.mopin")
m.set_xyz_coodinates_from_internal()
m.write_in_pdb("out2.pdb")

phi0 = m.atoms[12].torsion_ang
psi0 = m.atoms[20].torsion_ang 

Ndiv = 36
of = open(sys.argv[1]+".phipsi","w") 
of.write("[phi] [psi] [vdwOverlap] [crashatom]\n")
for i in range(Ndiv):
  phi = -180 + 360*float(i)/float(Ndiv)
  for j in range(Ndiv):
    psi = -180 + 360*float(j)/float(Ndiv)
    ofname = "tri%.0f_%.0f.pdb"%(phi,psi) 
    if (sys.argv[1]=='AAA.pdb'):
      m.atoms[12].torsion_ang = phi
      m.atoms[14].torsion_ang = phi - 120
      m.atoms[16].torsion_ang = phi + 120
      m.atoms[20].torsion_ang = psi
      m.atoms[13].torsion_ang = psi + 180
    if (sys.argv[1]=='AVA.pdb'):
      m.atoms[12].torsion_ang = phi
      m.atoms[14].torsion_ang = phi - 120
      m.atoms[18].torsion_ang = phi + 120
      m.atoms[26].torsion_ang = psi
      m.atoms[13].torsion_ang = psi + 180
    if (sys.argv[1]=='AGA.pdb'):
      m.atoms[12].torsion_ang = phi
      m.atoms[14].torsion_ang = phi - 120
      m.atoms[16].torsion_ang = phi + 120
      m.atoms[17].torsion_ang = psi
      m.atoms[13].torsion_ang = psi + 180




    m.set_xyz_coodinates_from_internal()
    m.set_internal_coodinate_angle()
    m.write_in_pdb("tmpout/"+ofname)
    ofname = "tri%.0f_%.0f.sdf"%(phi,psi) 
    m.write_in_sdf("tmpout/"+ofname)
    [dmin,crashstr] = cal_distance_bwn_heavy_atoms(m)
    print "phi %f psi %f dmin %f '%s'"%(phi,psi,dmin,ofname)
    of.write("%f %f %f %s\n"%(phi,psi,-dmin,crashstr))
#m.write_in_mopin("out.mopin")
  of.write("\n")
print "phi0 %f psi0 %f"%(phi0,psi0)
of.close()
