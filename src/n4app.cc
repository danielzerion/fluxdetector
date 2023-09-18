#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4_ui.hh"
#include "n4-volumes.hh"

#include <G4GenericMessenger.hh>

#include <G4PrimaryParticle.hh>
#include <G4String.hh>
#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include <G4Event.hh>           // needed to inject primary particles into an event
#include <G4Box.hh>             // for creating shapes in the geometry
#include <G4Sphere.hh>          // for creating shapes in the geometry
#include <FTFP_BERT.hh>         // our choice of physics list
#include <G4RandomDirection.hh> // for launching particles in random directions


#include <G4ThreeVector.hh>
#include <cstdlib>

void verify_number_of_args(int argc){
  if (argc != 2) {
    std::cerr << "Wrong number of arguments: " << argc
              << "\nUsage:\n./n4app <number of events>" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

struct my {
  G4double       straw_radius{0.1 * m};
  G4double      bubble_radius{0.2 * m};
  G4double      socket_rot   {-90 * deg};
  G4String      particle_name{"geantino"};
  G4double      particle_energy{511 * keV};
  G4ThreeVector particle_dir {};
};

auto my_generator(const my& my) {
  return [&](G4Event* event) {
    auto particle_type = n4::find_particle(my.particle_name);
    auto vertex = new G4PrimaryVertex();
    auto r = my.particle_dir.mag2() > 0 ? my.particle_dir : G4RandomDirection();
    vertex -> SetPrimary(new G4PrimaryParticle(
                           particle_type,
                           r.x(), r.y(), r.z(),
                           my.particle_energy
                         ));
    event  -> AddPrimaryVertex(vertex);
  };
}

n4::actions* create_actions(my& my, unsigned& n_event) {
  auto my_stepping_action = [&] (const G4Step* step) {
    auto pt = step -> GetPreStepPoint();
    auto volume_name = pt -> GetTouchable() -> GetVolume() -> GetName();
    if (volume_name == "straw" || volume_name == "bubble") {
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

  // Pure Geant4 definition of heavy water
  // G4int ncomponents, natoms;
  // G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  // G4Element* O  = new G4Element("Oxygen"        ,"O" , 8., 16.00*g/mole);
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);
  // G4Material* D2O = new G4Material("D2O", 1.107*g/cm3, ncomponents=2,
  //                                  kStateLiquid, 293.15*kelvin, 1*atmosphere);
  // D2O->AddElement(D, natoms=2);
  // D2O->AddElement(O, natoms=1);

  // Working towards a nain4 definition of heavy water
  auto D2O_new = nain4::material_from_elements_N(
    "D2O", 1.107*g/cm3,
    {.state=kStateLiquid, .temp=293.15*kelvin, .pressure = 1*atmosphere},
    {{D, 2}, {"O", 1}});

  auto air    = n4::material("G4_AIR");
  auto world  = n4::box("World").cube(2*m).x(3*m).volume(air);
  auto det    = n4::tubs("Detector").r(0.56*m).z(1*m).place(D2O_new).in(world).rotate_x(90*deg).now();
  return n4::place(world).now();
}

int main(int argc, char* argv[]) {
    unsigned n_event = 0;

  my my;

  // The trailing slash after '/my_geometry' is CRUCIAL: without it, the
  // messenger violates the principle of least surprise.
  auto messenger = new G4GenericMessenger{nullptr, "/my/", "docs: bla bla bla"};
  messenger -> DeclarePropertyWithUnit("straw_radius"      , "m"  , my. straw_radius  );
  messenger -> DeclarePropertyWithUnit("bubble_radius"     , "m"  , my.bubble_radius  );
  messenger -> DeclarePropertyWithUnit("socket_rot"        , "deg", my.socket_rot     );
  messenger -> DeclarePropertyWithUnit("particle_energy"   , "keV", my.particle_energy);
  messenger -> DeclareProperty        ("particle"          ,        my.particle_name  );
  messenger -> DeclareProperty        ("particle_direction",        my.particle_dir   );

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
