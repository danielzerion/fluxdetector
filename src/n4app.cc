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
  G4double lab_size = 3*m;
  G4double detector_length = 1*m;
  G4double detector_radius = 56*cm;
  std::unique_ptr<G4ParticleGun> gun;
  G4int particles_per_event = 1000;

  G4double      bubble_radius{0.2 * m};
  G4double      socket_rot   {-90 * deg};
  G4double      particle_energy{511 * keV};
  G4ThreeVector particle_dir {};
};

auto my_generator(my& my) {
  my.gun.reset(new G4ParticleGun{});
  std::vector<G4ParticleDefinition*> nu {
      n4::find_particle(     "nu_e"),
      n4::find_particle(     "nu_mu"),
      n4::find_particle("anti_nu_mu")
  };

  std::vector<G4double> weights{1,2,3};
  auto pick = n4::random::biased_choice{weights};

  return [nu = std::move(nu), pick = std::move(pick), &my](G4Event* event) {
    for (size_t i=0; i<my.particles_per_event; i++) {
      my.gun -> SetParticleDefinition(nu[pick()]);
      G4double x = my.lab_size * (G4UniformRand() - 0.5);
      G4double y = my.lab_size * (G4UniformRand() - 0.5);
      my.gun -> GeneratePrimaryVertex(event);
      my.gun -> SetParticlePosition({x, y, -my.lab_size/2});
      my.gun -> SetParticleEnergy(30 * MeV);
      my.gun -> SetParticleMomentumDirection({0,0,1});
      my.gun -> GeneratePrimaryVertex(event);
    }
  };
}

n4::actions* create_actions(my& my, unsigned& n_event) {
  auto my_stepping_action = [&] (const G4Step* step) {
    auto pt = step -> GetPreStepPoint();
    auto volume_name = pt -> GetTouchable() -> GetVolume() -> GetName();
    if (volume_name == "Detector") {
      auto pos = pt -> GetPosition();
      std::cout << volume_name << " " << pos << std::endl;
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
  auto world  = n4::box("World").cube(my.lab_size).volume(air);
  auto det    = n4::tubs("Detector").r(my.detector_radius).z(my.detector_length).place(D2O).in(world).rotate_x(90*deg).now();
  return n4::place(world).now();
}

int main(int argc, char* argv[]) {
  unsigned n_event = 0;

  my my;

  // The trailing slash after '/my_geometry' is CRUCIAL: without it, the
  // messenger violates the principle of least surprise.
  auto messenger = new G4GenericMessenger{nullptr, "/my/", "docs: bla bla bla"};
  messenger -> DeclarePropertyWithUnit("bubble_radius"     , "m"  , my.bubble_radius  );
  messenger -> DeclarePropertyWithUnit("socket_rot"        , "deg", my.socket_rot     );
  messenger -> DeclarePropertyWithUnit("particle_energy"   , "keV", my.particle_energy);
  messenger -> DeclareProperty("n_particles_per_event", my.particles_per_event);
  messenger -> DeclareProperty("particle_direction"   , my.particle_dir       );

    n4::run_manager::create()
    .ui("my-program-name", argc, argv)
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
