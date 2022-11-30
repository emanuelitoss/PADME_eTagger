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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "OutputColors.hh"

#include "G4RunManager.hh" 
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

G4double conversionFactor = 1239.8*nm*eV;

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fTagger(nullptr),
  fStepLimit(nullptr),
  fScoringVolume(0)
{}

DetectorConstruction::~DetectorConstruction(){
  delete fStepLimit;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){ 

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // Dimensions of my solids
  const G4double shape_taggerX = 600.*mm, shape_taggerY = 44.*mm, shape_taggerZ = 5.*mm;
  const G4double shape_SiPMX = 1.*mm, shape_SiPMY = 3.*mm, shape_SiPMZ = shape_SiPMY;
  G4double deltaShape = 5.*cm;
  const G4double shapeEnvelopeX = 0.5*shape_taggerX + deltaShape, shapeEnvelopeY = 0.5*shape_taggerY + deltaShape,
    shapeEnvelopeZ = 0.5*shape_taggerZ + deltaShape;

  /*************** WORLD ***************/

  G4Material* VacuumMaterial = CreateLowVacuumAir();

  G4Box* solidWorld = new G4Box("World", 1.1*shapeEnvelopeX, 1.1*shapeEnvelopeY, 1.1*shapeEnvelopeZ);
  
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        VacuumMaterial,      //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      false);                //overlaps checking
  
  /*************** ENVELOPE ***************/

  G4Box* solidEnv = new G4Box("Envelope", shapeEnvelopeX, shapeEnvelopeY, shapeEnvelopeZ);
  
  G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv, VacuumMaterial, "Envelope");

  G4VPhysicalVolume* physEnvelope
    = new G4PVPlacement(0, G4ThreeVector(), logicEnv, "Envelope", logicWorld, false, 0, checkOverlaps);

  /*************** PLASTIC SCINTILLATOR ***************/

  G4Material* plastic_material = this->CreatePlasticMaterial();

  G4Box* TaggerShape = new G4Box("Tagger Box", 0.5*shape_taggerX, 0.5*shape_taggerY, 0.5*shape_taggerZ);

  G4LogicalVolume* Tagger_LogicalVolume 
    = new G4LogicalVolume(TaggerShape, plastic_material, "Tagger LV");

  fTagger = 
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Tagger_LogicalVolume, "Tagger PV", logicEnv, false, 0, checkOverlaps);

  OpticalSurfaceTagger_Vacuum(fTagger, physEnvelope);

  /*************** SILICON PHOTOMULTIPLIERS ***************/

  G4Material* borosilicate = this->CreatePyrex();

  G4Box* SiPM_Shape = new G4Box("SiPM Box", 0.5*shape_SiPMX, 0.5*shape_SiPMY, 0.5*shape_SiPMZ);
  
  G4LogicalVolume* SiPM_LogicalVolume 
    = new G4LogicalVolume(SiPM_Shape, borosilicate, "SiPM Logical Volume");
  
  SiPM_LogicalVolume->SetUserLimits(new G4UserLimits(shape_SiPMX/2));

  G4int counter = 0;
  G4ThreeVector positions_sipm;
  G4double starting_position = 0.5*(shape_SiPMY - shape_taggerY);
  G4double deltaPos = (shape_taggerY - 4.*shape_SiPMY)/3 + shape_SiPMY;
  for( auto& sipm : fSiPMs )
  {
    if (counter < 4)
    {
      positions_sipm = G4ThreeVector(0.5*(shape_taggerX+shape_SiPMX), starting_position + (counter%4)*deltaPos, 0);
    }
    else
    {
      positions_sipm = G4ThreeVector(-0.5*(shape_taggerX+shape_SiPMX), starting_position + (counter-4)*deltaPos, 0);
    }

    sipm = new G4PVPlacement(0, positions_sipm, SiPM_LogicalVolume, "SiPM detector", logicEnv, false, 0, checkOverlaps); 
    OpticalSurfaceTagger_SiPM(fTagger, sipm);
    
    counter++;
  }

  /*************** END ***************/

  fScoringVolume = logicEnv;

  return physWorld;

}

G4Material* DetectorConstruction::CreatePlasticMaterial() const {
  
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* plastic_basic_material = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material* plastic_material = new G4Material("PlasticScintillator", 1.023*g/cm3, plastic_basic_material);

  // photon energy = hPlanck * (speed of light) / wavelength
  G4double energies_photons[10] = {conversionFactor/(380.*nm),
                                   conversionFactor/(390.*nm),
                                   conversionFactor/(400.*nm),
                                   conversionFactor/(410.*nm),
                                   conversionFactor/(420.*nm),
                                   conversionFactor/(430.*nm),
                                   conversionFactor/(440.*nm),
                                   conversionFactor/(460.*nm),
                                   conversionFactor/(480.*nm),
                                   conversionFactor/(500.*nm)};

  G4double rindex[10] = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58};
  G4double absorption[10] = {140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm};
  G4double scintillation_spectrum[10] = {0.04, 0.3, 0.87, 0.92, 0.55, 0.45, 0.34, 0.13, 0.05, 0.01};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
 
  // properties independent of energy
  G4double light_yield_Anthracene = 17400./MeV;
  MPT->AddConstProperty("SCINTILLATIONYIELD", 0.68*light_yield_Anthracene);
  MPT->AddConstProperty("FASTTIMECONSTANT", 1.8*ns);
  MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);

  // properties that depend on energy
  MPT->AddProperty("RINDEX", energies_photons, rindex, 10)->SetSpline(true);
  MPT->AddProperty("ABSLENGTH", energies_photons, absorption, 10)->SetSpline(true);
  MPT->AddProperty("FASTCOMPONENT", energies_photons, scintillation_spectrum, 10)->SetSpline(true);

  plastic_material->SetMaterialPropertiesTable(MPT);
  plastic_material->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  return plastic_material;

}

G4Material* DetectorConstruction::CreateLowVacuumAir() const {

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* basic_air = nist->FindOrBuildMaterial("G4_AIR");

  double p_atm = 1.013, p_vacuum = 0.001; // pressure [bar]
  G4double density_atm = 1.204*kg/m3;
  G4double density_vacuum = density_atm*p_vacuum/p_atm; // scaling the desity according to ideal gas law

  G4Material* vacuum_air = new G4Material("Low Vacuum Air", density_vacuum, basic_air);

  // photon energy = hPlanck * (speed of light) / wavelength
  G4double photonenergy[3] = {conversionFactor/(380.*nm),
                              conversionFactor/(440.*nm), 
                              conversionFactor/(500.*nm)};
  G4double rindex[3] = {1.000293, 1.000293, 1.000293};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
  MPT->AddProperty("RINDEX", photonenergy, rindex, 3);

  vacuum_air->SetMaterialPropertiesTable(MPT);

  return vacuum_air;

}

G4Material* DetectorConstruction::CreatePyrex() const {
  
  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* pyrex_basic = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* pyrex = new G4Material("Borosilicate_Glass", 2.23*g/cm3, pyrex_basic);
  
  // photon energy = hPlanck * (speed of light) / wavelength
  G4double photonenergy[3] = {conversionFactor/(320.*nm),
                              conversionFactor/(400.*nm),  
                              conversionFactor/(480.*nm)};
                            
  G4double rindex[3] = {1.471, 	1.471, 	1.471};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  MPT->AddProperty("RINDEX", photonenergy, rindex, 3);

  pyrex->SetMaterialPropertiesTable(MPT);

  return pyrex;
}

void DetectorConstruction::OpticalSurfaceTagger_SiPM(G4VPhysicalVolume* Tagger_PV, G4VPhysicalVolume* TheOtherPV) const {

  G4OpticalSurface* opPlasticSurface = new G4OpticalSurface("Plastic scintillator surface");
  opPlasticSurface->SetModel(unified);
  opPlasticSurface->SetType(dielectric_dielectric);
  opPlasticSurface->SetFinish(polished);

  G4LogicalBorderSurface* LogicalPlasticSurface = new G4LogicalBorderSurface(
    "Plastic Surface", Tagger_PV, TheOtherPV, opPlasticSurface);
  
  // // To re-get this surface use the following:
  // G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
  //  LogicalPlasticSurface->GetSurface(Tagger_PV, TheOtherPV)->GetSurfaceProperty());

  G4MaterialPropertiesTable* SurfaceTable = new G4MaterialPropertiesTable();

  // SiPM model Hamamatsu 13360-3050PE: see datasheet
  G4double photonenergy[10] = {conversionFactor/(320.*nm),
                               conversionFactor/(350.*nm),
                               conversionFactor/(400.*nm),
                               conversionFactor/(450.*nm),
                               conversionFactor/(500.*nm),
                               conversionFactor/(550.*nm),
                               conversionFactor/(600.*nm),
                               conversionFactor/(650.*nm),
                               conversionFactor/(700.*nm),
                               conversionFactor/(800.*nm)};

  G4double transmissiom[10] = {0.025, 0.18, 0.35, 0.4, 0.38, 0.33, 0.27, 0.20, 0.15, 0.08};
  SurfaceTable->AddProperty("TRANSMISSION", photonenergy, transmissiom, 10)->SetSpline(true);
  
  // if(opticalSurface) opticalSurface->DumpInfo();
  // SurfaceTable->DumpTable();
  opPlasticSurface->SetMaterialPropertiesTable(SurfaceTable);

}

void DetectorConstruction::OpticalSurfaceTagger_Vacuum(G4VPhysicalVolume* Tagger_PV, G4VPhysicalVolume* Vacuum_PV) const {

  G4OpticalSurface* opPlasticSurface = new G4OpticalSurface("Scintillator-Vacuum(Air) surface");

  G4LogicalBorderSurface* LogicalPlasticSurface = new G4LogicalBorderSurface(
    "Plastic Surface", Tagger_PV, Vacuum_PV, opPlasticSurface);
  
  G4double sigma_alpha = 0.001;
  opPlasticSurface->SetModel(unified);
  opPlasticSurface->SetType(dielectric_dielectric);
  opPlasticSurface->SetFinish(polishedbackpainted); //polishedfrontpainted - polishedbackpainted
  opPlasticSurface->SetSigmaAlpha(sigma_alpha);

  // // To re-get this surface use the following:
  // G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
  //   LogicalPlasticSurface->GetSurface(Tagger_PV, Vacuum_PV)->GetSurfaceProperty());

  G4MaterialPropertiesTable* SurfaceTable = new G4MaterialPropertiesTable();

  // photon energy = hPlanck * (speed of light) / wavelength
  G4double energies_photons[10] = {conversionFactor/(380.*nm),
                                   conversionFactor/(390.*nm),
                                   conversionFactor/(400.*nm),
                                   conversionFactor/(410.*nm),
                                   conversionFactor/(420.*nm),
                                   conversionFactor/(430.*nm),
                                   conversionFactor/(440.*nm),
                                   conversionFactor/(460.*nm),
                                   conversionFactor/(480.*nm),
                                   conversionFactor/(500.*nm)};

  G4double reflectivity[10] = {0.68, 0.84, 0.9, 0.93, 0.95, 0.955, 0.96, 0.96, 0.96, 0.96};
  SurfaceTable->AddProperty("RINDEX", energies_photons, reflectivity, 10)->SetSpline(true);

  // if(opticalSurface) opticalSurface->DumpInfo();
  // SurfaceTable->DumpTable();
  opPlasticSurface->SetMaterialPropertiesTable(SurfaceTable);
  
}
