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

#include "../include/DetectorConstruction.hh"
#include "../include/OutputColors.hh"

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
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

const double hPlanck = 4.135655e-15; // [eV*s]
const double c_light = 3e+8; // [m/s]
const double meters_to_nanometers = 1e9; // 1 m = 1e9 nm

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScintillator(nullptr),
  fCerenkovPMT(nullptr),
  fScintillatorPMT(nullptr),
  fStepLimit(nullptr),
  fScoringVolume(0)
{}

DetectorConstruction::~DetectorConstruction(){
  delete fStepLimit;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){ 
  
  // Get nist material manager:
  // (per prendere materiali già preparati da geant. Ma tanto noi li costruiamo)
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // Dimensions of BGO
  const G4double shape_bgoX = 600.*mm, shape_bgoY = 44.*mm, shape_bgoZ = 5.*mm;
  const G4double shape_PMT_X = 5.*mm, shape_PMT_Y = shape_bgoY, shape_PMT_Z = shape_bgoZ;

  // envelope and world radius
  if (shape_bgoZ > shape_bgoY)
  {
    if (shape_bgoZ > shape_bgoX) minimal_radius = shape_bgoZ;
    else minimal_radius = shape_bgoX;
  }
  else if (shape_bgoY > shape_bgoX) minimal_radius = shape_bgoY;
  else minimal_radius = shape_bgoX;

  G4double radius_sphere = minimal_radius;

  // create vacuum material
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  G4Material* VacuumMaterial = new G4Material("LowDensityVacuum", atomicNumber, massOfMole, density, kStateGas,temperature, pressure);

  //
  // World
  //
  
  G4Sphere* solidWorld = new G4Sphere("World", 0, sqrt(3)*radius_sphere, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
      
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
                      checkOverlaps);        //overlaps checking
                   
  //     
  // Envelope
  //  
  G4Sphere* solidEnv = new G4Sphere("World", 0, radius_sphere, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
      
  G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv, VacuumMaterial, "Envelope");

  
  G4VPhysicalVolume* physEnvelope
    = new G4PVPlacement(0, G4ThreeVector(), logicEnv, "Envelope", logicWorld, false, 0, checkOverlaps);
  

  // limit stepLength to 4.8cm in order to optimize rejecting logic of non-detected events 
  // (see SteppingAction.cc SteppingAction::EstinguishParticleIfNotTrigger )
  G4double maxStep = 4.8*cm;
  fStepLimit = new G4UserLimits(maxStep);
  logicEnv->SetUserLimits(fStepLimit);

  //     
  // BGO Crystal + Photomultipliers
  //  
  G4Material* plastic_material = this->CreatePlasticMaterial();
  G4Material* borosilicate = this->CreatePyrex();

  // position: choosen in order to minimize the envelope sphere
  G4ThreeVector bgo_position = G4ThreeVector(0, 0, 0);
  G4ThreeVector pmt1_position = G4ThreeVector(0.5*(shape_bgoX+shape_PMT_X), 0, 0);
  G4ThreeVector pmt2_position = G4ThreeVector(-0.5*(shape_bgoX+shape_PMT_X), 0, 0);
        
  // BGO & PMTs shape
  G4Box* BGOShape = new G4Box("BGO Box", 0.5*shape_bgoX, 0.5*shape_bgoY, 0.5*shape_bgoZ);
  G4Box* PMTsShape = new G4Box("PMT Box", 0.5*shape_PMT_X, 0.5*shape_PMT_Y, 0.5*shape_PMT_Z);

  // logical volumes
  G4LogicalVolume* BGO_LogicalVolume 
    = new G4LogicalVolume(BGOShape, plastic_material, "BGO Crystal");
  
  G4LogicalVolume* PMT_LogicalVolume 
    = new G4LogicalVolume(PMTsShape, borosilicate, "PMT made up of borosilicate glass");

  // in order to ensure at least one step in the PMT
  PMT_LogicalVolume->SetUserLimits(new G4UserLimits(shape_PMT_X/2));

  // delta vector useful to translate PMTs
  G4ThreeVector original_vector = G4ThreeVector(0, 0.5*(shape_bgoY+shape_PMT_Y), 0);

  fScintillator = 
    new G4PVPlacement(0, bgo_position, BGO_LogicalVolume, "BGO Crystal", logicEnv, false, 0, checkOverlaps);
  
  fCerenkovPMT 
   = new G4PVPlacement(0, pmt1_position, PMT_LogicalVolume, "Cherenkov PMT", logicEnv, false, 0, checkOverlaps);

  fScintillatorPMT = new G4PVPlacement(0, pmt2_position, PMT_LogicalVolume, "Scintill. PMT", logicEnv, false, 0, checkOverlaps);

  // optical properties of BGO surface
  OpticalSurfacePlastic_SiPM(fScintillator, fCerenkovPMT);
  OpticalSurfacePlastic_SiPM(fScintillator, fScintillatorPMT);

  fScoringVolume = logicEnv;

  //always return the physical World
  return physWorld;
}

// optical between of BGO and PMTs
void DetectorConstruction::OpticalSurfacePlastic_SiPM(G4VPhysicalVolume* BGO_PV, G4VPhysicalVolume* TheOtherPV) const {

  G4OpticalSurface* opBGOSurface = new G4OpticalSurface("BGO Surface");
  // see here: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/BackupVersions/V10.7/html/TrackingAndPhysics/physicsProcess.html
  opBGOSurface->SetType(dielectric_LUTDAVIS);
  opBGOSurface->SetModel(DAVIS);
  opBGOSurface->SetFinish(PolishedESRGrease_LUT); // Using optical grease between PMTs and BGO

  G4LogicalBorderSurface* LogicalBGOSurface = new G4LogicalBorderSurface(
    "BGO Surface", BGO_PV, TheOtherPV, opBGOSurface);
  
  G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
    LogicalBGOSurface->GetSurface(BGO_PV, TheOtherPV)->GetSurfaceProperty());

  if(opticalSurface) opticalSurface->DumpInfo();

  // Generate & Add Material Properties Table attached to the optical surface
  // energy = hPlanck * (light speed) / wavelength
  // reference for reflectivity https://aip.scitation.org/doi/10.1063/1.3272909
  // reference for transmittance https://www.researchgate.net/publication/264828576_The_Radiation_Hard_BGO_Crystals_for_Astrophysics_Applications
  G4int n10 = 10;
  // G4int n6 = 6;
  G4double energies_photons[] = { 2.*eV, 2.5*eV, 3.*eV, 3.5*eV, 4.*eV,
                                4.5*eV, 4.75*eV, 4.9*eV, 5.*eV, 5.4*eV };
  G4double reflectivity[] = { 0.125, 0.13, 0.14, 0.155, 0.175,
                              0.25, 0.29, 0.243, 0.21, 0.22 };

  G4MaterialPropertiesTable* SurfaceTable = new G4MaterialPropertiesTable();

  SurfaceTable->AddProperty("REFLECTIVITY", energies_photons, reflectivity, n10);
  
  SurfaceTable->DumpTable();
  opBGOSurface->SetMaterialPropertiesTable(SurfaceTable);

}

// BGO - material
G4Material* DetectorConstruction::CreatePlasticMaterial() const {
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  /*
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  G4Material* VacuumMaterial = new G4Material("LowDensityVacuum", atomicNumber, massOfMole, density, kStateGas,temperature, pressure);
  */

  // Composition
  G4Material* plastic_basic_material = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material* plastic_material = new G4Material("BismuthGermaniumOxygen Crystal", 1.023*g/cm3, plastic_basic_material);

  // energy = hPlanck * (light speed) / wavelength
  G4double energies_photons[10] = {hPlanck*c_light*meters_to_nanometers/380.*eV,
                                    hPlanck*c_light*meters_to_nanometers/390.*eV,
                                    hPlanck*c_light*meters_to_nanometers/400.*eV,
                                    hPlanck*c_light*meters_to_nanometers/410.*eV,
                                    hPlanck*c_light*meters_to_nanometers/420.*eV,
                                    hPlanck*c_light*meters_to_nanometers/430.*eV,
                                    hPlanck*c_light*meters_to_nanometers/440.*eV,
                                    hPlanck*c_light*meters_to_nanometers/460.*eV,
                                    hPlanck*c_light*meters_to_nanometers/480.*eV,
                                    hPlanck*c_light*meters_to_nanometers/500.*eV};

  // from our documents
  G4double rindex[10] = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58};
  G4double absorption[10] = {140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm, 140.*cm};
  G4double scintillation_spectrum[10] = {0.04, 0.3, 0.87, 0.92, 0.55, 0.45, 0.34, 0.13, 0.05, 0.01};

  // new instance of Material Properties
  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
 
  // properties independent of energy
  MPT->AddConstProperty("FASTTIMECONSTANT", 1.8*ns);
  MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  G4double light_yield_Anthracene = 17400./MeV;
  MPT->AddConstProperty("SCINTILLATIONYIELD", 0.68*light_yield_Anthracene);

  // properties that depend on energy
  MPT->AddProperty("RINDEX", energies_photons, rindex, 10)->SetSpline(true);
  MPT->AddProperty("ABSLENGTH", energies_photons, absorption, 10);
  MPT->AddProperty("FASTCOMPONENT", energies_photons, scintillation_spectrum, 10)->SetSpline(true);

  // bgo material
  plastic_material->SetMaterialPropertiesTable(MPT);
  plastic_material->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  return plastic_material;
}

// Birosilicate glass - material
G4Material* DetectorConstruction::CreatePyrex() const {
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // BGO material definition
  G4Material* pyrex_basic = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* pyrex = new G4Material("Borosilicate_Glass", 2.23*g/cm3, pyrex_basic);
  
  const G4int n3 = 3;
  // energy = hPlanck * (light speed) / wavelength
  G4double photonenergy[n3] = {hPlanck*c_light*meters_to_nanometers/320.*eV,  // lower wavelength cutoff 320 nm
                            hPlanck*c_light*meters_to_nanometers/400.*eV,     // intermediate energy
                            hPlanck*c_light*meters_to_nanometers/480.*eV};    // maximum emission at 480 nm
  G4double rindex[n3]     = {1.471, 	1.471, 	1.471};                   // google
  // new instance of Material Properties
  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  // properties that depend on energy
  MPT->AddProperty("RINDEX", photonenergy, rindex, n3);

  // material
  pyrex->SetMaterialPropertiesTable(MPT);

  return pyrex;
}