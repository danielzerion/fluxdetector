#include "materials.hh"

#include <n4-material.hh>
#include <n4-utils.hh>
#include <n4-place.hh>

#include <G4SystemOfUnits.hh>

using vec_double = std::vector<G4double>;
using vec_int    = std::vector<G4int>;

const vec_double OPTPHOT_ENERGY_RANGE{1*eV, 8.21*eV};
const G4double hc = CLHEP::h_Planck * CLHEP::c_light;

auto D2O_without_properties() {
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D -> AddIsotope(H2, 100*perCent);

  return nain4::material_from_elements_N(
    "D2O", 1.107*g/cm3,
    {.state=kStateLiquid, .temp=293.15*kelvin, .pressure = 1*atmosphere},
    {{D, 2}, {"O", 1}});
}

// TODO: this gives us an idea of what is likely to be needed for D2O
G4Material* d2o_csi_hybrid_FIXME_with_properties(G4double scint_yield) {
  std::cout << "-------------------- SCINT YIELD: " << scint_yield / (1/MeV) << std::endl;

  auto csi = D2O_without_properties();
  // csi_rindex: values taken from "Optimization of Parameters for a CsI(Tl) Scintillator Detector Using GEANT4-Based Monte Carlo..." by Mitra et al (mainly page 3)
  //  csi_scint: Fig. 2 in the paper
  // must be in increasing ENERGY order (decreasing wavelength) for scintillation to work properly
  auto     csi_energies = n4::scale_by(hc*eV, {1/0.9, 1/0.7, 1/0.54, 1/0.35}); // denominator is wavelength in micrometres
  vec_double csi_rindex =                     {1.79 , 1.79 , 1.79  , 1.79  };  //vec_double csi_rindex = {2.2094, 1.7611};
  vec_double  csi_scint =                     {0.0  , 0.1  , 1.0   , 0.0   };
  auto    csi_abslength = n4::scale_by(m    , {5    , 5    , 5     , 5     });
  // Values from "Temperature dependence of pure CsI: scintillation light yield and decay time" by Amsler et al
  // "cold" refers to ~77K, i.e. liquid nitrogen temperature
  G4double csi_scint_yield      =  scint_yield;
  G4double csi_scint_yield_cold = 50000 / MeV;
  G4double csi_time_fast        =     6 * ns;
  G4double csi_time_slow        =    28 * ns;
  G4double csi_time_fast_cold   =  1015 * ns; // only one component at cold temps!
  G4double csi_time_slow_cold   =  1015 * ns;
  G4MaterialPropertiesTable *csi_mpt = n4::material_properties()
    .add("RINDEX"                 , csi_energies, csi_rindex)
    .add("SCINTILLATIONCOMPONENT1", csi_energies,  csi_scint)
    .add("SCINTILLATIONCOMPONENT2", csi_energies,  csi_scint)
    .add("ABSLENGTH"              , csi_energies, csi_abslength)
    .add("SCINTILLATIONTIMECONSTANT1", csi_time_fast)
    .add("SCINTILLATIONTIMECONSTANT2", csi_time_slow)
    .add("SCINTILLATIONYIELD"        , csi_scint_yield)
    .add("SCINTILLATIONYIELD1"       ,     0.57   )
    .add("SCINTILLATIONYIELD2"       ,     0.43   )
    .add("RESOLUTIONSCALE"           ,     1.0    )
    .done();
  csi -> GetIonisation() -> SetBirksConstant(0.00152 * mm/MeV);
  csi -> SetMaterialPropertiesTable(csi_mpt);
  return csi;
}

G4Material* air_with_properties() {
    auto air = n4::material("G4_AIR");
    G4MaterialPropertiesTable *mpt_air = n4::material_properties()
        .add("RINDEX", OPTPHOT_ENERGY_RANGE, {1, 1})
        .done();
    air -> SetMaterialPropertiesTable(mpt_air);
    return air;
}

G4Material* teflon_with_properties() {
    auto teflon = n4::material("G4_TEFLON");
    // Values could be taken from "Optical properties of Teflon AF amorphous fluoropolymers" by Yang, French & Tokarsky (using AF2400, Fig.6)
    // but are also stated in the same paper as above
    G4MaterialPropertiesTable *mpt_teflon = n4::material_properties()
        .add("RINDEX", OPTPHOT_ENERGY_RANGE, {1.35, 1.35})
        .done();
    teflon -> SetMaterialPropertiesTable(mpt_teflon);
    return teflon;
}

// --------------------------------------------------------------------------------

// TODO: we haven't used this at all yet, but the reflection already works!
#include <G4LogicalBorderSurface.hh>
#include <G4OpticalSurface.hh>

void place_D2O_teflon_border_surface_between(G4PVPlacement* one, G4PVPlacement* two) {
  static G4OpticalSurface* D2O_teflon_surface = nullptr;
  auto name = "D2O-Teflon-Surface";
  if (! D2O_teflon_surface) {
    D2O_teflon_surface = new G4OpticalSurface(name);
    // TODO: Values copied blindly from double-sipm, need verification
    D2O_teflon_surface -> SetType(dielectric_dielectric);
    D2O_teflon_surface -> SetModel(unified);
    D2O_teflon_surface -> SetFinish(groundfrontpainted);
    D2O_teflon_surface -> SetSigmaAlpha(0.0);

    vec_double pp = {2.038*eV, 4.144*eV};
    // According to the docs, for UNIFIED, dielectric_dielectric surfaces only the Lambertian reflection is turned on
    D2O_teflon_surface -> SetMaterialPropertiesTable(
      n4::material_properties{}
      .add("REFLECTIVITY", pp, {0.98 , 0.98})
      .done());
  }
  new G4LogicalBorderSurface(name, one, two, D2O_teflon_surface);
}
