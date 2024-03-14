#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // 10% loss of the scintillation internal surface reflection
  if ( volume->GetName() == "Surface" )
  {
  auto currentTrack = step->GetTrack();
  auto name = currentTrack->GetParticleDefinition()->GetParticleName();
    if ( name == "opticalphoton" )
    {
      G4double random = G4UniformRand();
      if(random < 0.1)
      {
        currentTrack->SetTrackStatus(fStopAndKill);
      }
    }
  }

  // make the photons that enter SiPM region disappear
  if ( volume->GetName() == "SiPM" )
  {
    auto currentTrack = step->GetTrack();
    auto name = currentTrack->GetParticleDefinition()->GetParticleName();
    if ( name == "opticalphoton" )
    {
      currentTrack->SetTrackStatus(fStopAndKill); 
    } 
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

