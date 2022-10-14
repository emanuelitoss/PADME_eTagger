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

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunData.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "g4root.hh"

#include <string>


RunAction::RunAction(double X, double Y)
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{
  // set printing event number per each event
  // G4RunManager::GetRunManager()->SetPrintProgress(1);

  fPercentX = X;
  fPercentY = Y;

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleMerging(true);

  // Creating histograms
  // virtual G4int 	CreateH1 (const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName="none", const G4String &fcnName="none")=0
  analysisManager->CreateH1("EnergyPlasticScintillator","Energy deposited in plastic scintillator", 100, 0., 100, "MeV");
  analysisManager->CreateH1("PhotonsTime[1]","Arival time of photons[1]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[2]","Arival time of photons[2]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[3]","Arival time of photons[3]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[4]","Arival time of photons[4]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[5]","Arival time of photons[5]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[6]","Arival time of photons[6]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[7]","Arival time of photons[7]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[8]","Arival time of photons[8]", 100, 0, 30, "ns");

  // Creating ntuples
  analysisManager->CreateNtuple("energyLoss", "Energy loss in each event in the tagger");
  analysisManager->CreateNtupleDColumn("EnergyPlasticScintillator");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("eTagDataTuples", "Tuples of arrival time of photons");
  analysisManager->CreateNtupleDColumn("PhotonsTime[1]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[2]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[3]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[4]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[5]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[6]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[7]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[8]");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("firstTimes","Tuples of first arrival time for each event and SiPM");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[1]");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[2]");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[3]");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[4]");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[5]");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[6]");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[7]");
  analysisManager->CreateNtupleDColumn("1stTimeSiPM[8]");
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
  // create the right string
  G4String fileName = "../data_eTag/data_eTag.root";
  analysisManager->OpenFile(fileName);
  
  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

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
  
    G4cout << " EnergyPlastic : mean = "
      << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

  }

  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

}

void RunAction::AddEdep(G4double edep){

  fEdep  += edep;
  fEdep2 += edep*edep;

}
