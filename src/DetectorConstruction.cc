#include "DetectorConstruction.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"

#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"

#include "G4SubtractionSolid.hh"

#include "G4SystemOfUnits.hh"

#include "math.h"
#include "G4VisAttributes.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{ 
  fBC420 = fAir = fSiPM = fsurface = fPMMA = fPethylene1 = fFe = fAl = nullptr;
  fN = fO = fC = fH = nullptr;
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

// define the materials
void DetectorConstruction::DefineMaterials()
{
    G4double a;  // atomic mass
    G4double z;  // atomic number
    G4double density;

    G4NistManager* nist_manager = G4NistManager::Instance();

    //define elements and element composition
    fH = new G4Element("H", "H", z = 1., a = 1.01 * g / mole);
    fC = new G4Element("C", "C", z = 6., a = 12.01 * g / mole);
    fN = new G4Element("N", "N", z = 7., a = 14.01 * g / mole);
    fO = new G4Element("O", "O", z = 8., a = 16.00 * g / mole);

    fFe = nist_manager->FindOrBuildMaterial("G4_Fe");
    fAl = nist_manager->FindOrBuildMaterial("G4_Al");
    fAir = new G4Material("Air", density = 1.29 * mg / cm3, 2);
    fAir->AddElement(fN, 70 * perCent);
    fAir->AddElement(fO, 30 * perCent);

    fBC420 = new G4Material("BC420", density = 1.032 * g / cm3, 2);
    fBC420->AddElement(fC, 10);
    fBC420->AddElement(fH, 11);

    fsurface = new G4Material("surface", density = 1.032 * g / cm3, 2);
    fsurface->AddElement(fC, 10);
    fsurface->AddElement(fH, 11);

    fSiPM = new G4Material("SiPM", density = 1.29 * mg / cm3, 2);
    fSiPM->AddElement(fN, 70 * perCent);
    fSiPM->AddElement(fO, 30 * perCent);

    G4int polyPMMA = 1;
    G4int nC_PMMA = 3+2*polyPMMA;
    G4int nH_PMMA = 6+2*polyPMMA;

    G4int polyeth = 1;
    G4int nC_eth = 2*polyeth;
    G4int nH_eth = 4*polyeth;

    fPMMA = new G4Material("PMMA", density=1190*kg/m3,3);
    fPMMA->AddElement(fH,nH_PMMA);
    fPMMA->AddElement(fC,nC_PMMA);
    fPMMA->AddElement(fO,2);
    
    fPethylene1 = new G4Material("Pethylene1", density=1200*kg/m3,2);
    fPethylene1->AddElement(fH,nH_eth);
    fPethylene1->AddElement(fC,nC_eth);


    //define the parameters of the material
    G4double BC420_Energy[] = { 2.38 * eV, 2.88 * eV, 3.45 * eV };
    const G4int nEntries = sizeof(BC420_Energy) / sizeof(G4double);
    G4double BC420_RefractionIndex[] = { 1.58, 1.58, 1.58 };
    assert(sizeof(BC420_RefractionIndex) == sizeof(BC420_Energy));
    G4double BC420_AbsorptionLength[] = { 140. * cm, 140. * cm, 140. * cm };
    assert(sizeof(BC420_AbsorptionLength) == sizeof(BC420_Energy));
    G4double BC420_ScintilFast[] = {1.0,1.0,1.0 };
    assert(sizeof(BC420_ScintilFast) == sizeof(BC420_Energy));
    G4double wavelength[] = { 0., 1, 0. };
    assert(sizeof(wavelength) == sizeof(BC420_Energy));

    BC420MPT = new G4MaterialPropertiesTable();
    BC420MPT->AddProperty("RINDEX", BC420_Energy, BC420_RefractionIndex, nEntries);
    BC420MPT->AddProperty("ABSLENGTH", BC420_Energy, BC420_AbsorptionLength, nEntries);
    BC420MPT->AddProperty("WLSCOMPONENT", BC420_Energy, wavelength, nEntries);
    BC420MPT->AddConstProperty("SCINTILLATIONYIELD", 10240./ MeV);
    BC420MPT->AddConstProperty("RESOLUTIONSCALE", 2.5);
    BC420MPT->AddConstProperty("FASTTIMECONSTANT", 0.9 * ns);
    BC420MPT->AddConstProperty("SLOWTIMECONSTANT", 2.1 * ns);
    fBC420->SetMaterialPropertiesTable(BC420MPT);
    fBC420->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

    G4double vacuum_Energy[] = { 2.0 * eV,7.0 * eV,7.14 * eV };
    G4double AirRefractiveIndex[] = { 1.00, 1.00, 1.00 };
    assert(sizeof(AirRefractiveIndex) == sizeof(vacuum_Energy));
    G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
    vacuum_mt->AddProperty("RINDEX", vacuum_Energy, AirRefractiveIndex, nEntries);
    fAir->SetMaterialPropertiesTable(vacuum_mt);

    G4double SiPM_RefractionIndex[] = { 1.57, 1.57, 1.57 };
    assert(sizeof(SiPM_RefractionIndex) == sizeof(BC420_Energy));
    G4MaterialPropertiesTable* SiPM_mt = new G4MaterialPropertiesTable();
    SiPM_mt->AddProperty("RINDEX", BC420_Energy, SiPM_RefractionIndex, nEntries);
    fSiPM->SetMaterialPropertiesTable(SiPM_mt);

    G4double wls_Energy[] = {2.00*eV,2.87*eV,2.90*eV,3.47*eV};
    const G4int wlsnum = sizeof(wls_Energy)/sizeof(G4double);

    G4double RefractiveIndexFiber[]={ 1.60, 1.60, 1.60, 1.60};
    assert(sizeof(RefractiveIndexFiber) == sizeof(wls_Energy));
    G4double AbsFiber[]={9.00*m,9.00*m,0.1*mm,0.1*mm};
    assert(sizeof(AbsFiber) == sizeof(wls_Energy));
    G4double EmissionFib[]={1.0, 1.0, 0.0, 0.0};
    assert(sizeof(EmissionFib) == sizeof(wls_Energy));
    G4MaterialPropertiesTable* fiberProperty = new G4MaterialPropertiesTable();
    fiberProperty->AddProperty("RINDEX",wls_Energy,RefractiveIndexFiber,wlsnum);
    fiberProperty->AddProperty("WLSABSLENGTH",wls_Energy,AbsFiber,wlsnum);
    fiberProperty->AddProperty("WLSCOMPONENT",wls_Energy,EmissionFib,wlsnum);
    fiberProperty->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
    fPMMA->SetMaterialPropertiesTable(fiberProperty);

    G4double RefractiveIndexClad1[]={ 1.49, 1.49, 1.49, 1.49};
    assert(sizeof(RefractiveIndexClad1) == sizeof(wls_Energy));
    G4MaterialPropertiesTable* clad1Property = new G4MaterialPropertiesTable();
    clad1Property->AddProperty("RINDEX",wls_Energy,RefractiveIndexClad1,wlsnum);
    clad1Property->AddProperty("ABSLENGTH",wls_Energy,AbsFiber,wlsnum);
    fPethylene1->SetMaterialPropertiesTable(clad1Property);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  G4bool checkOverlaps = true;    //check for overlaps
  G4int strip_num[12] = { 30, 100, 40, 100, 55, 100, 70, 100, 80, 100, 95, 100};    //the number of stripes in each layer
  G4RotationMatrix* rm_Fe = new G4RotationMatrix;
  rm_Fe->rotateX(90 * deg);
  G4double x_a = 105 * cm;
  G4double x_b = x_a + 210 * sqrt(3) * cm;
  G4VisAttributes* blank = new G4VisAttributes(false);

  auto solidworld = new G4Box( "World", 20 * m , 20 * m , 20 * m );
  auto logicworld = new G4LogicalVolume( solidworld, fAir, "World" );
  auto physworld = new G4PVPlacement( nullptr, G4ThreeVector(), logicworld, "World", 0, false, 0, checkOverlaps);
  logicworld->SetVisAttributes(blank);
  //Fe frame
  for ( G4int i5 = 0; i5 < 2; i5 ++)
  {
    for ( G4int i4 = 0; i4 < 12; i4 ++ )
    {
      G4RotationMatrix* rm_env = new G4RotationMatrix;
      rm_env->rotateZ( i4 * 30 * deg);
      G4double env_sizeX = 210 * ( 1 + sqrt(3) ) * cm;
      G4double env_sizeY = 210 * ( 3 + 1.5 * sqrt(3) ) * cm;
      G4double env_posZ = 202.5 * ( 2 * i5 - 1 ) * cm;
      auto solidenv = new G4Box( "Envelope", env_sizeX , env_sizeY , 202.5 * cm );
      auto logicenv = new G4LogicalVolume( solidenv, fAir, "Envelope" );
      auto physenv = new G4PVPlacement( rm_env, G4ThreeVector(0,0,env_posZ), logicenv, "Envelope", logicworld, false, i4, checkOverlaps);
      logicenv->SetVisAttributes(blank);

      G4double Fe_posX = -1 * 105 * cm;
      G4double Fe_posY = -1 * 105 * ( 2.5 + 1.5 * sqrt(3) ) * cm;
      auto solidFe = new G4Trd( "Fe",  0.5 * x_a, 0.5 * x_b, 202.5 * cm, 202.5 * cm, 52.5 * cm);
      auto logicFe = new G4LogicalVolume( solidFe, fFe, "Fe" );
      auto physFe = new G4PVPlacement( rm_Fe, G4ThreeVector( Fe_posX, Fe_posY, 0 ), logicFe, "Fe", logicenv, false, i4, checkOverlaps);

      G4Box* solidSiPM = new G4Box("SiPM", 3 * mm, 0.005 * cm, 3 * mm);

      //place the scintillator
      for ( G4int i1 = 0; i1 < 6; i1 ++ )
      {
        G4int i7 = 2 * i1;
        G4double layer_sizeX =  ( 4 * strip_num[i7] + 0.2 ) * cm;
        G4double layer_posy = ( 1 - 2 * i5 ) * 2.4 * cm;
        G4double layer_posZ = ( i1 * 15 - 35 ) * cm;
        auto solidlayer = new G4Box("Layer", 0.5 * layer_sizeX, 200.1 * cm, 2.1 * cm);
        auto logiclayer =
              new G4LogicalVolume(solidlayer,
                                 fAir,
                                  "Layer");
              new G4PVPlacement(nullptr,
                                G4ThreeVector( 0, layer_posy, layer_posZ ),
                                logiclayer,
                                "Layer",
                                logicFe,
                                false,
                                i1,
                                checkOverlaps);
        
        G4double Al_sizeX = ( 4 * strip_num[i7] + 0.1 ) * cm;
        auto solidAl = new G4Box("Al", 0.5 * Al_sizeX, 200.05 * cm, 2.05 * cm);
        auto logicAl =
              new G4LogicalVolume(solidAl,
                                 fAl,
                                  "Al");
              new G4PVPlacement(nullptr,
                                G4ThreeVector(),
                                logicAl,
                                "Al",
                                logiclayer,
                                false,
                                i1,
                                checkOverlaps);
        for ( G4int i6 = 0; i6 < 2; i6 ++ )
        {
          G4int i3 = 2 * i1 + i6;
          for( G4int i2 = 0; i2 < strip_num[i3]; i2 ++ )
          {
            G4RotationMatrix* rm = new G4RotationMatrix;
            G4RotationMatrix* rm_fiber = new G4RotationMatrix;
            G4RotationMatrix* rm_cut3 = new G4RotationMatrix;
            rm_fiber->rotateX( 90 * deg);
            rm_cut3->rotateX( 90 * deg);
            G4double strip_posX;
            G4double strip_posY;
            G4double strip_posZ;
            G4double strip_sizeY;
            //horizontal and vertical alternate arrangement
            if( i6 == 1 )
            {
              rm->rotateZ(90 * deg);
              strip_sizeY = 2 * strip_num[i7] * cm;
              strip_posX = 0;
              strip_posY = ( 4 * i2 - 198 ) * cm;
              strip_posZ = 0.5 * cm;
            }
            else
            {
              strip_sizeY = 200 * cm;
              strip_posX = ( 2 + 4 * i2 - 2 * strip_num[i3] ) * cm;
              strip_posY = 0;
              strip_posZ = -0.5 * cm;
            } 
            G4double surface_sizeY = strip_sizeY;
            G4double BC420_sizeY = strip_sizeY - 0.01 * cm;
            G4double cut1_sizeY = BC420_sizeY;
            G4double cut2_sizeY = BC420_sizeY;
            G4double cut3_sizeZ = BC420_sizeY;
            G4double Cladding_sizeZ = BC420_sizeY;
            G4double Core_sizeZ = BC420_sizeY;
            G4double SiPM_posY = strip_sizeY - 0.005 * cm;
            G4Box* solidstrip = new G4Box("Strip", 2 * cm , strip_sizeY, 0.5 * cm);
            G4Box* solidsurface= new G4Box("Surface", 2 * cm, surface_sizeY, 0.5 * cm);
            G4Box* solidBC420 = new G4Box("BC420", 1.99 * cm, BC420_sizeY, 0.49 * cm);
            G4Box* solidcut1 = new G4Box("Cut1", 1.1 * mm, cut1_sizeY, 0.05 * mm);
            G4Box* solidcut2 = new G4Box("Cut2", 1.1 * mm, cut2_sizeY, 2.4 * mm);
            G4Tubs* solidcut3 = new G4Tubs("Cut3", 0, 1.1 * mm, cut3_sizeZ, 0, 180 * deg);
            G4Tubs* solidCladding = new G4Tubs("Cladding", 0.95 * mm , 1 * mm,  Cladding_sizeZ, 0, 360 * deg);
            G4Tubs* solidCore = new G4Tubs("Core", 0, 0.95 * mm, Core_sizeZ, 0, 360 * deg);
            auto logicstrip =
              new G4LogicalVolume(solidstrip,
                                 fAir,
                                  "Strip");
              new G4PVPlacement(rm,
                                G4ThreeVector( strip_posX, strip_posY, strip_posZ ),
                                logicstrip,
                                "Strip",
                                logicAl,
                                false,
                                i3,
                                checkOverlaps);

            auto logicsurface =
              new G4LogicalVolume(solidsurface,
                                  fsurface,
                                  "Surface");
            G4PVPlacement* physsurface =
              new G4PVPlacement(nullptr,
                                G4ThreeVector(),
                                logicsurface,
                                "Surface",
                                logicstrip,
                                false,
                                i2,
                                checkOverlaps);

            auto logiccut1 = 
              new G4LogicalVolume(solidcut1,
	      	                fAir,
	                        "Cut1");
              new G4PVPlacement(nullptr,
	  	              G4ThreeVector( 0, 0, 4.95 * mm),
		                logiccut1,
		                "Cut1",
		                logicsurface,
		                false,
		                i2,
	                  checkOverlaps);
  
            auto logicSiPM = new G4LogicalVolume(solidSiPM, fSiPM, "SiPM");
            for ( G4int j = 0; j < 2; j++ )
            {
              G4PVPlacement* physSiPM = new G4PVPlacement(nullptr,
                                                          G4ThreeVector( 0, ( 2 * j - 1 ) * SiPM_posY, 0),
                                                          logicSiPM,
                                                          "SiPM",
                                                          logicsurface,
                                                          false,
                                                          j,
                                                          checkOverlaps);
            }

            auto logicBC420 =
              new G4LogicalVolume(solidBC420,
                                  fBC420,
                                  "BC420");
            G4PVPlacement* physBC420 =
              new G4PVPlacement(nullptr,
                                G4ThreeVector(),
                                logicBC420,
                                "BC420",
                                logicsurface,
                                false,
                                i2,
                                checkOverlaps);  
  
          auto logiccut2 = 
            new G4LogicalVolume(solidcut2,
	  	                fAir,
 	  	         	      "Cut2");
            new G4PVPlacement(nullptr,
	                            G4ThreeVector( 0, 0, 2.5 * mm),
		                          logiccut2,
		                          "Cut2",
		                          logicBC420,
		                          false,
		                          i2,
	                            checkOverlaps);

	        auto logiccut3 =
              new G4LogicalVolume(solidcut3,
                                  fAir,
                                  "Cut3");
              new G4PVPlacement(rm_cut3,
                                G4ThreeVector( 0, 0, 0.1 * mm),
                                logiccut3,
                                "Cut3",
                                logicBC420,
                                false,
                                i2,
                                checkOverlaps);


            G4LogicalVolume* logicCladding =
              new G4LogicalVolume(solidCladding,
                                  fPethylene1,
                                  "Cladding");
            G4PVPlacement* physCladding =
              new G4PVPlacement(rm_fiber,
                                G4ThreeVector(),
                                logicCladding,
                                "Cladding",
				logicBC420,
				false,
                                i3,
                                checkOverlaps);

            G4LogicalVolume* logicCore =
              new G4LogicalVolume(solidCore,
                                  fPMMA,
                                  "Core");
            G4PVPlacement* physCore =
              new G4PVPlacement(rm_fiber,
                                G4ThreeVector(),
                                logicCore,
                                "Core",
                                logicBC420,
                                false,
                                i2,
                                checkOverlaps);
  
            //define surface
            G4OpticalSurface* Surface = new G4OpticalSurface("Surface");
            new G4LogicalBorderSurface("Surface", physBC420, physsurface, Surface);
            Surface->SetType(dielectric_metal);
            Surface->SetFinish(polished);
            Surface->SetModel(glisur);
            G4double sur_Energy[] = { 2.38 * eV, 2.88 * eV, 3.45 * eV };
            const G4int num = sizeof(sur_Energy) / sizeof(G4double);
            G4double sur_RefractionIndex[] = { 1.58, 1.58, 1.58 };
            assert(sizeof(sur_RefractionIndex) == sizeof(sur_Energy));
            G4MaterialPropertiesTable* SURMPT = new G4MaterialPropertiesTable();
            SURMPT->AddProperty("RINDEX", sur_Energy, sur_RefractionIndex,num);
            Surface->SetMaterialPropertiesTable(SURMPT);

            G4OpticalSurface* Cladding = new G4OpticalSurface("Cladding");
            new G4LogicalBorderSurface("Surface", physCore, physCladding, Cladding);
            Cladding->SetType(dielectric_metal);
            Cladding->SetFinish(polished);
            Cladding->SetModel(glisur);
            G4double cladding_RefractionIndex[] = { 1.49, 1.49, 1.49 };
            assert(sizeof(sur_RefractionIndex) == sizeof(sur_Energy));
            G4MaterialPropertiesTable* CLAMPT = new G4MaterialPropertiesTable();
            CLAMPT->AddProperty("RINDEX", sur_Energy, cladding_RefractionIndex,num);
            Cladding->SetMaterialPropertiesTable(CLAMPT);
          }
        }
      }
    }
  }
  //
  //always return the physical World
  //
  return physworld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
