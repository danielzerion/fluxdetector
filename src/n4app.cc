#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4_ui.hh"
#include "n4-utils.hh"
#include "n4-volumes.hh"

#include <G4GenericMessenger.hh>

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

void verify_number_of_args(int argc){
  if (argc != 2) {
    std::cerr << "Wrong number of arguments: " << argc
              << "\nUsage:\n./n4app <number of events>" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

struct my {
  G4double lab_size         =  3    * m;
  G4double detector_length  =  1    * m;
  G4double detector_radius  =  0.56 * m;
  G4double vessel_thickness =  1    *cm;
  G4double particle_energy  = 30    * MeV;
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
    my.gun->SetParticlePosition(random_position_in_detector());
    my.gun->SetParticleEnergy(my.particle_energy);
    my.gun->SetParticleMomentumDirection(G4RandomDirection());
    my.gun->GeneratePrimaryVertex(event);
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

auto my_geometry(const my& my) {

  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);

  // Heavy water
  auto D2O = nain4::material_from_elements_N(
    "D2O", 1.107*g/cm3,
    {.state=kStateLiquid, .temp=293.15*kelvin, .pressure = 1*atmosphere},
    {{D, 2}, {"O", 1}});

  auto air    = n4::material("G4_AIR");
  auto Al     = n4::material("G4_Al");
  auto world  = n4::box("World").cube(my.lab_size).volume(air);
  auto vessel = n4::tubs("Vessel")
    .r(my.detector_radius + my.vessel_thickness)
    .z(my.detector_length + my.vessel_thickness)
    .place(Al)
    .in(world).at_z(my.vessel_thickness/2).rotate_x(90*deg).now();
  auto det    = n4::tubs("Detector")
    .r(my.detector_radius)
    .z(my.detector_length)
    .place(D2O)
    .in(vessel).at_z(-my.vessel_thickness/2).now();
  return n4::place(world).now();
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
  messenger -> DeclarePropertyWithUnit("particle_energy" , "MeV", my.particle_energy);
  messenger -> DeclareProperty("particle", my.particle);

    n4::run_manager::create()
    .ui("fluxdetector", argc, argv)
    .macro_path("macs")
    .apply_command("/my/straw_radius 0.5 m")
    .apply_early_macro("early-hard-wired.mac")
    .apply_cli_early_macro() // CLI --early-macro executed at this point
    // .apply_command(...) // also possible after apply_early_macro

    .physics<FTFP_BERT>(0) // verbosity 0
    .geometry([&] { return my_geometry(my); })
    .actions(create_actions(my, n_event))

    .apply_command("/my/particle e-")
    .apply_late_macro("late-hard-wired.mac")
    .apply_cli_late_macro() // CLI --late-macro executed at this point
    // .apply_command(...) // also possible after apply_late_macro

    .run();
}
