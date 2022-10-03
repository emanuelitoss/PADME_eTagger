//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "../include/RunAction.hh"
#include "../include/PrimaryGeneratorAction.hh"
#include "../include/DetectorConstruction.hh"
#include "../include/RunData.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "g4root.hh"


RunAction::RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{
  // set printing event number per each event
  // G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleMerging(true);
  // Note: merging ntuples is available only with Root output

  // Creating histograms
  analysisManager->CreateH1("EnergyBGO","Energy deposited in BGO crystal", 800, 0., 100, "MeV");
  //analysisManager->CreateH1("EnergyPlastic1","Energy deposited in first plastic scintillator", 800, 0., 100*MeV);
  //analysisManager->CreateH1("EnergyPlastic2","Energy deposited in second plastic scintillator", 800, 0., 100*MeV);
  analysisManager->CreateH1("Cherenkov","Cherenkov energy production in BGO", 200, 0., 100*keV);
  analysisManager->CreateH1("Scintillation","Scintillation energy production in BGO", 200, 0., 100*keV);
  
  // Creating ntuple
  analysisManager->CreateNtuple("eTagDataTuples", "EnengyDeposit");
  analysisManager->CreateNtupleDColumn("EnergyInBGO");
  //analysisManager->CreateNtupleDColumn("EnergyInScintillator1");
  //analysisManager->CreateNtupleDColumn("EnergyInScintillator2");
  analysisManager->CreateNtupleDColumn("EnergyInBGO_Cherenkov");
  analysisManager->CreateNtupleDColumn("EnergyInBGO_Scintillation");
  analysisManager->FinishNtuple();


  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
  
}

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

G4Run* RunAction::GenerateRun()
{
  return (new RunData);
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  G4String fileName = "../eTag_data/eTag_data.root";
  analysisManager->OpenFile(fileName);
  
  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  // reset number of detected particles
  detectedParticles = 0;

}

void RunAction::EndOfRunAction(const G4Run* run){

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents; // non dovrebbe essere (edep2 - edep*edep)/nofEvents ???
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const DetectorConstruction* detectorConstruction
   = static_cast<const DetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const PrimaryGeneratorAction* generatorAction
    = static_cast<const PrimaryGeneratorAction*>
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  // print histogram statistics
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl;
    }
  
    G4cout << " EnergyBGO : mean = "
      << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

    /*
    G4cout << " EnergyScintillator1 : mean = "
      << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

    G4cout << " EnergyScintillator2 : mean = "
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
    */
    
    G4cout << " EnergyCherenkovBGO : mean = "
      << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
    G4cout << " EnergyScintillationBGO : mean = "
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
    /*
    G4cout << " NumberCherenkovBGO : mean = "
      << analysisManager->GetH1(5)->mean() << " Number of photons"
      << " rms = "
      << analysisManager->GetH1(5)->rms() <<  " Number of photons" << G4endl;

    G4cout << " NumberScintillationBGO : mean = "
      << analysisManager->GetH1(6)->mean() << " number of photons"
      << " rms = "
      << analysisManager->GetH1(6)->rms() <<  " number of photons" << G4endl;
    */
  }

  // number od detected events
  if(!IsMaster())
    std::cout << OBOLDWHITE << "Number of detected events: " << detectedParticles << ORESET << std::endl;

  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

}

void RunAction::AddEdep(G4double edep){

  fEdep  += edep;
  fEdep2 += edep*edep;

}
