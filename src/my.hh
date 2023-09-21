#ifndef MY_HH
#define MY_HH

#include <G4GenericMessenger.hh>
#include <G4SystemOfUnits.hh>
#include <G4ParticleGun.hh>

struct my {
  G4double lab_size         =    3    *  m;
  G4double detector_length  =    1    *  m;
  G4double detector_radius  =    0.56 *  m;
  G4double vessel_thickness =    1    * cm;
  G4double teflon_thickness =    0.5  * mm;
  G4double particle_energy  =   30    * MeV;
  G4double scint_yield      = 3200    / MeV;
  std::unique_ptr<G4ParticleGun> gun{};
  G4String particle = "e-";

  my() : msngr{new G4GenericMessenger{nullptr, "/my/", "docs: bla bla bla"}} {
    // The trailing slash after '/my' is CRUCIAL: without it, the msngr
    // violates the principle of least surprise.
    msngr -> DeclarePropertyWithUnit("lab_size"        ,   "m", lab_size);
    msngr -> DeclarePropertyWithUnit("detector_radius" ,   "m", detector_radius);
    msngr -> DeclarePropertyWithUnit("detector_length" ,   "m", detector_length);
    msngr -> DeclarePropertyWithUnit("vessel_thickness",   "m", vessel_thickness);
    msngr -> DeclarePropertyWithUnit("teflon_thickness",   "m", teflon_thickness);
    msngr -> DeclarePropertyWithUnit("particle_energy" , "MeV", particle_energy);
    msngr -> DeclareProperty        ("scint_yield"     ,        scint_yield);
    msngr -> DeclareProperty        ("particle"        ,        particle);
  }
  std::unique_ptr<G4GenericMessenger> msngr;
};

#endif // MY_HH
