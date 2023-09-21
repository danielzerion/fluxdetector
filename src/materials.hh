#ifndef MATERIALS_HH
#define MATERIALS_HH


#include <G4Material.hh>
#include <G4Types.hh>

G4Material* d2o_csi_hybrid_FIXME_with_properties(G4double scint_yield);
G4Material* air_with_properties();
G4Material* teflon_with_properties();

#endif // MATERIALS_HH
