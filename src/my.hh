#ifndef MY_HH
#define MY_HH

#include <G4GenericMessenger.hh>
#include <G4SystemOfUnits.hh>
#include <G4ParticleGun.hh>
#include <G4Types.hh>
#include <Randomize.hh>

struct my {
  G4double lab_size                  =    4      *  m;
  G4double d2o_z                     =    1      *  m;
  G4double d2o_r                     =    0.62   *  m;
  G4double h2o_z_delta               =    1.15    *  m;
  G4double h2o_r_delta               =    0.67   *  m;
  G4double vessel_out_delta          =    1      * cm;
  G4double teflon_delta              =    0.5    * mm;
  G4double vessel_in_r_delta         =    1.    *cm;
  G4double vessel_in_z_delta         =    1.    *cm;
  G4double pmt_thickness             =   100      * mm;
  G4double pmt_radius                =   202.    * mm;
  G4double particle_energy           =   30    * MeV;
  G4double scint_yield               = 3200    / MeV;
  G4int    physics_verbose           =    0;
  G4int         em_verbose           =    0;
  G4int    optical_verbose           =    0;
  std::unique_ptr<G4ParticleGun> gun{};
  G4String particle = "e-";
  void set_seed(G4long seed){G4Random::setTheSeed(seed);}
  G4double full_r() const  {return d2o_r +   vessel_in_r_delta +   h2o_r_delta + teflon_delta +   vessel_out_delta;}
  G4double full_z() const  {return d2o_z + 2*vessel_in_z_delta + 2*h2o_z_delta + teflon_delta + 2*vessel_out_delta;}

  my() : msngr{new G4GenericMessenger{nullptr, "/my/", "docs: bla bla bla"}} {
    // The trailing slash after '/my' is CRUCIAL: without it, the msngr
    // violates the principle of least surprise.
    msngr -> DeclarePropertyWithUnit("lab_size"           ,   "m", lab_size);
    msngr -> DeclarePropertyWithUnit("d2o_r"              ,   "m", d2o_r);
    msngr -> DeclarePropertyWithUnit("d2o_z"              ,   "m", d2o_z);
    msngr -> DeclarePropertyWithUnit("h2o_r_delta"        ,  "cm", h2o_r_delta);
    msngr -> DeclarePropertyWithUnit("h2o_z_delta"        ,  "cm", h2o_z_delta);
    msngr -> DeclarePropertyWithUnit("vessel_in_r_delta"  ,  "cm", vessel_in_r_delta);
    msngr -> DeclarePropertyWithUnit("vessel_in_z_delta"  ,  "cm", vessel_in_z_delta);
    msngr -> DeclarePropertyWithUnit("vessel_out_delta"   ,  "cm", vessel_out_delta);
    msngr -> DeclarePropertyWithUnit("teflon_delta"       ,  "mm", teflon_delta);
    msngr -> DeclarePropertyWithUnit("pmt_thickness"      ,  "mm", pmt_thickness);
    msngr -> DeclarePropertyWithUnit("pmt_diameter"       ,  "mm", pmt_radius);
    msngr -> DeclarePropertyWithUnit("particle_energy"    , "MeV", particle_energy);
    msngr -> DeclareProperty        ("physics_verbose"    ,        physics_verbose);
    msngr -> DeclareProperty        ("em_verbose"         ,        em_verbose);
    msngr -> DeclareProperty        ("optical_verbose"    ,        optical_verbose);
    msngr -> DeclareProperty        ("scint_yield"        ,        scint_yield);
    msngr -> DeclareProperty        ("particle"           ,        particle);
    msngr -> DeclareMethod          ("seed"               ,        &my::set_seed);
  }
  std::unique_ptr<G4GenericMessenger> msngr;
};

#endif // MY_HH
