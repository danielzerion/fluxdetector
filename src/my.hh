#ifndef MY_HH
#define MY_HH

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
  std::unique_ptr<G4ParticleGun> gun;
  G4String particle = "e-";
};

#endif // MY_HH
