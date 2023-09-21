#include "my.hh"
#include "materials.hh"

#include <nain4.hh>
#include <g4-mandatory.hh>
#include <n4_ui.hh>
#include <n4-utils.hh>
#include <n4-volumes.hh>


#include <G4EmStandardPhysics_option4.hh>
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

auto PMT(const my& my) {
  static auto dummy_material = n4::material("G4_Cu");
  static auto unnnumbered_pmt = n4::tubs("PMT")
    .z(my.pmt_thickness).r(my.pmt_radius)
    .place(dummy_material)
    .at_z((my.pmt_thickness-my.detector_length)/2);
  return unnnumbered_pmt;
}

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

  auto water = n4::tubs("Water")
    .r(my.detector_radius)
    .z(my.detector_length)
    .place(D2O)
    .in(reflector).at_z(-my.teflon_thickness/2).now();

  auto one_pmt = PMT(my).in(water);


  // TODO: steal ideas for PMT distribution from this code:

  // void jaszczak_phantom::rod_sector(unsigned long n, G4double r,
  //                                   G4LogicalVolume* body, G4Material* material) const {
  //   auto d = 2 * r;
  //   auto z = (height_rods - height_body) / 2;
  //   G4RotationMatrix around_z_axis{{0,0,1}, n*pi/3};

  //   // Sector displacement from centre, to accommodate gap between sectors
  //   auto dx = gap * cos(pi/6);
  //   auto dy = gap * sin(pi/6);
  //   // Displacement of first rod WRT sector corner
  //   dx += r * sqrt(3);
  //   dy += r;
  //   // Basis vectors of rod lattice
  //   const auto Ax = 2.0, Ay = 0.0;
  //   const auto Bx = 1.0, By = sqrt(3);
  //   auto a = 0;
  //   for (bool did_b=true ; did_b; a+=1) {
  //     did_b = false;
  //     for (auto b = 0; /*break in body*/; b+=1, did_b = true) {
  //       auto x = (a*Ax + b*Bx) * d + dx;
  //       auto y = (a*Ay + b*By) * d + dy;
  //       if (sqrt(x*x + y*y) + r + margin >= radius_body) { break; }
  //       auto label = std::string("Rod-") + std::to_string(n);
  //       auto rod = volume<G4Tubs>(label, material, 0.0, r, height_rods/2, 0.0, twopi);
  //       place(rod).in(body).at(x,y,z).rotate(around_z_axis).now();
  //     }
  //   }
  // }


  n4::place(one_pmt).copy_no(0)                                            .now();
  n4::place(one_pmt).copy_no(1).at_x(  my.detector_radius - my.pmt_radius ).now();
  n4::place(one_pmt).copy_no(2).at_x(-(my.detector_radius - my.pmt_radius)).now();
  n4::place(one_pmt).copy_no(3).at_y(  my.detector_radius - my.pmt_radius ).now();
  n4::place(one_pmt).copy_no(4).at_y(-(my.detector_radius - my.pmt_radius)).now();


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
