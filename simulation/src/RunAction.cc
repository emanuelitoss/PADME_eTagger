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
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "g4root.hh"

#include <string>

RunAction::RunAction()
: G4UserRunAction()
{

  DefineCommands();

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleMerging(true);

  // Create histograms of photons detection time
  analysisManager->CreateH1("PhotonsTime[1]","Arival time of photons[1]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[2]","Arival time of photons[2]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[3]","Arival time of photons[3]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[4]","Arival time of photons[4]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[5]","Arival time of photons[5]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[6]","Arival time of photons[6]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[7]","Arival time of photons[7]", 100, 0, 30, "ns");
  analysisManager->CreateH1("PhotonsTime[8]","Arival time of photons[8]", 100, 0, 30, "ns");

  // Create Ntuple of photons detection time
  analysisManager->CreateNtuple("generic_times", "Tuples of arrival time of photons");
  analysisManager->CreateNtupleDColumn("PhotonsTime[1]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[2]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[3]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[4]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[5]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[6]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[7]");
  analysisManager->CreateNtupleDColumn("PhotonsTime[8]");
  analysisManager->FinishNtuple();

  // Create Ntuple of minimum photons detection time
  analysisManager->CreateNtuple("arrival_times","Tuples of first arrival time for each event and SiPM");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[1]");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[2]");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[3]");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[4]");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[5]");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[6]");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[7]");
  analysisManager->CreateNtupleDColumn("arrival_times_SiPM[8]");
  analysisManager->FinishNtuple();

  // Create Ntuple of signals charge
  analysisManager->CreateNtuple("charges","Total charge for each SiPM for event");
  analysisManager->CreateNtupleDColumn("Charges[1]");
  analysisManager->CreateNtupleDColumn("Charges[2]");
  analysisManager->CreateNtupleDColumn("Charges[3]");
  analysisManager->CreateNtupleDColumn("Charges[4]");
  analysisManager->CreateNtupleDColumn("Charges[5]");
  analysisManager->CreateNtupleDColumn("Charges[6]");
  analysisManager->CreateNtupleDColumn("Charges[7]");
  analysisManager->CreateNtupleDColumn("Charges[8]");
  analysisManager->FinishNtuple();

  // Create Ntuple of hitting position
  analysisManager->CreateNtuple("positions","Positions (x,y) of the incident electrons (or positrons)");
  analysisManager->CreateNtupleDColumn("X_position");
  analysisManager->CreateNtupleDColumn("Y_position");
  analysisManager->FinishNtuple();
  
}

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
  delete fMessenger;
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

  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  G4String OutputCompleteName = G4String("../../data_analysis/data_eTag") + OutputFileName + G4String(".root");
  analysisManager->OpenFile(OutputCompleteName);
}

void RunAction::EndOfRunAction(const G4Run* run){

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

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
     << G4endl;

  // print histogram statistics
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

}

void RunAction::DefineCommands() {
  
  fMessenger = new G4GenericMessenger(this, "/PadmETag/OutputFile/", "Output file");

  // position X
  auto& setFileNameCommand
  = fMessenger->DeclareProperty("SetName", OutputFileName);
  setFileNameCommand.SetParameterName("Output File Name", true);
  setFileNameCommand.SetDefaultValue("");

}
