#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4_ui.hh"
#include "n4-utils.hh"
#include "n4-volumes.hh"


#include <G4EmStandardPhysics_option4.hh>
#include <G4GenericMessenger.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4OpticalPhysics.hh>
#include <G4OpticalSurface.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4PrimaryParticle.hh>
#include <G4String.hh>
#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include <G4Event.hh>           // needed to inject primary particles into an event
#include <G4Box.hh>             // for creating shapes in the geometry
#include <G4Sphere.hh>          // for creating shapes in the geometry
#include <FTFP_BERT.hh>         // our choice of physics list
#include <G4RandomDirection.hh> // for launching particles in random directions

#include <G4ThreeVector.hh>
#include <Randomize.hh>
#include <cstdlib>

using vec_double = std::vector<G4double>;
using vec_int    = std::vector<G4int>;

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

auto my_generator(my& my) {
  my.gun.reset(new G4ParticleGun{n4::find_particle(my.particle)});

  return [&my](G4Event *event) {
    auto random_position_in_detector = [&my] {
      auto [x, z] = n4::random::random_on_disc(my.detector_radius);
      auto y      = n4::random::uniform(
        -my.detector_length / 2,
         my.detector_length / 2
      );
      return G4ThreeVector{x, y, z};
    };
    my.gun -> SetParticlePosition(random_position_in_detector());
    my.gun -> SetParticleEnergy(my.particle_energy);
    my.gun -> SetParticleMomentumDirection(G4RandomDirection());
    my.gun -> GeneratePrimaryVertex(event);
  };
}

n4::actions* create_actions(my& my, unsigned& n_event) {
  auto my_stepping_action = [&] (const G4Step* step) {
    auto pt = step -> GetPreStepPoint();
    auto volume_name = pt -> GetTouchable() -> GetVolume() -> GetName();
    if (volume_name == "Detector") {
      auto pos = pt -> GetPosition();
      auto p = step -> GetTrack() -> GetParticleDefinition() -> GetParticleName();
      std::cout << "pos: " << pos << "   " << p << std::endl;
    }
  };

  auto my_event_action = [&] (const G4Event*) {
     n_event++;
     std::cout << "end of event " << n_event << std::endl;
  };

  return (new n4::        actions{my_generator(my)  })
 -> set( (new n4::   event_action{                  }) -> end(my_event_action) )
 -> set(  new n4::stepping_action{my_stepping_action});
}

// TODO: we haven't used this at all yet, but the reflection already works!
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

auto D2O_without_properties() {
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D -> AddIsotope(H2, 100*perCent);

  return nain4::material_from_elements_N(
    "D2O", 1.107*g/cm3,
    {.state=kStateLiquid, .temp=293.15*kelvin, .pressure = 1*atmosphere},
    {{D, 2}, {"O", 1}});
}

const vec_double OPTPHOT_ENERGY_RANGE{1*eV, 8.21*eV};
const G4double hc = CLHEP::h_Planck * CLHEP::c_light;

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
    auto air = air_with_properties();
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

auto my_geometry(const my& my) {
  // Heavy water

  auto D2O = d2o_csi_hybrid_FIXME_with_properties(my.scint_yield);

  auto air    = n4::material("G4_AIR");
  auto Al     = n4::material("G4_Al");
  auto teflon = teflon_with_properties();
  auto world  = n4::box("World").cube(my.lab_size).volume(air);

  auto vessel = n4::tubs("Vessel")
    .r(my.detector_radius + my.vessel_thickness + my.teflon_thickness)
    .z(my.detector_length + my.vessel_thickness + my.teflon_thickness)
    .place(Al)
    .in(world).at_z((my.vessel_thickness+my.teflon_thickness)/2).rotate_x(90*deg).now();

  auto reflector = n4::tubs("Reflector")
    .r(my.detector_radius + my.teflon_thickness)
    .z(my.detector_length + my.teflon_thickness)
    .place(teflon)
    .in(vessel).at_z(-my.vessel_thickness/2).now();

  n4::tubs("Detector")
    .r(my.detector_radius)
    .z(my.detector_length)
    .place(D2O)
    .in(reflector).at_z(-my.teflon_thickness/2).now();

  return n4::place(world).now();
}

auto my_physics_list() {
    G4int verbosity;
    auto physics_list = new FTFP_BERT{verbosity = 0};
    physics_list ->  ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    return physics_list;
}

int main(int argc, char* argv[]) {
  unsigned n_event = 0;

  my my;

  // The trailing slash after '/my' is CRUCIAL: without it, the messenger
  // violates the principle of least surprise.
  auto messenger = new G4GenericMessenger{nullptr, "/my/", "docs: bla bla bla"};
  messenger -> DeclarePropertyWithUnit("lab_size"        , "m"  , my.lab_size  );
  messenger -> DeclarePropertyWithUnit("detector_radius" , "m"  , my.detector_radius);
  messenger -> DeclarePropertyWithUnit("detector_length" , "m"  , my.detector_length);
  messenger -> DeclarePropertyWithUnit("vessel_thickness", "m"  , my.vessel_thickness);
  messenger -> DeclarePropertyWithUnit("teflon_thickness", "m"  , my.teflon_thickness);
  messenger -> DeclarePropertyWithUnit("particle_energy" , "MeV", my.particle_energy);
  messenger -> DeclareProperty("scint_yield", my.scint_yield);
  messenger -> DeclareProperty("particle"   , my.particle);

    n4::run_manager::create()
    .ui("fluxdetector", argc, argv)
    .macro_path("macs")
    .apply_command("/my/straw_radius 0.5 m")
    .apply_early_macro("early-hard-wired.mac")
    .apply_cli_early_macro() // CLI --early-macro executed at this point
    // .apply_command(...) // also possible after apply_early_macro

    .physics(my_physics_list)
    .geometry([&] { return my_geometry(my); })
    .actions(create_actions(my, n_event))

    .apply_command("/my/particle e-")
    .apply_late_macro("late-hard-wired.mac")
    .apply_cli_late_macro() // CLI --late-macro executed at this point
    // .apply_command(...) // also possible after apply_late_macro

    .run();
}
