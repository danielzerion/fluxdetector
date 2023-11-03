#ifndef MATERIALS_HH
#define MATERIALS_HH

#include <G4Material.hh>
#include <G4Types.hh>
#include <G4PVPlacement.hh>

G4Material* d2o_csi_hybrid_FIXME_with_properties(G4double scint_yield);
G4Material * acrylic_with_properties();
G4Material * h2o_with_properties();
G4Material* air_with_properties();
G4Material* teflon_with_properties();
void place_D2O_teflon_border_surface_between(G4PVPlacement* one, G4PVPlacement* two); 

#endif // MATERIALS_HH
