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

#include "include/DetectorConstruction.hh"
#include "include/ActionInitialization.hh"
#include "include/RunData.hh"

#include "G4RunManagerFactory.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "FTFP_BERT.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "Randomize.hh"
#include "G4StepLimiterPhysics.hh"

#include <time.h>

int main(int argc, char** argv)
{

  // storing runtime of the code:
  clock_t tStart = clock();

  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  // G4Random::setTheEngine(new CLHEP::RanecuEngine);
  double myseed = round(tStart);
  G4Random::setTheSeed(myseed);

  // Construct the default run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Detector construction
  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);
  runManager->SetNumberOfThreads(1);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  auto actionInitialization = new ActionInitialization(detConstruction);
  runManager->SetUserInitialization(actionInitialization);
  
  // Initialize visualization
  // G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else{
    // interactive mode if no macro file is involved
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    /* vedi esempio B5
    if (ui->IsGUI()) {
         UImanager->ApplyCommand("/control/execute gui.mac");
    }
    */
    ui->SessionStart();
    delete ui;
  }

  // evaluate running time of the program
  int execution_time_s = round((clock() - tStart)/CLOCKS_PER_SEC);
  int remainder_s = execution_time_s%60;
  int execution_time_min = (int)((execution_time_s-remainder_s)/60);

  std::cout << OBOLDWHITE << "Execution time of the code: "
    << execution_time_min << " min " << remainder_s << " s" << ORESET << std::endl;

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  delete visManager;
  delete runManager;
}