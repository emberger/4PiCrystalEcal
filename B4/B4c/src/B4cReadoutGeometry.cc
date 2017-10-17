#include "B4cReadoutGeometry.hh"
#include "B4cCalorimeterSD.hh"


#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include <iostream>

#include "B4cDetParams.hh"

MyRO::MyRO() : G4VReadOutGeometry(){
}

MyRO::MyRO(G4String aString) : G4VReadOutGeometry(aString){



}

MyRO::~MyRO(){
}





G4VPhysicalVolume* MyRO::Build(){

								// Geometry parameters also for Detector Construction
								GetInst().SetfNofLayers(50);
								GetInst().SetcalorSizeXY(2000); // in mm
								GetInst().SettileLenX(10); // in mm
								GetInst().SettileLenY(10); // in mm

								GetInst().SetabsoThickness(1.8); // in mm
								GetInst().SetgapThickness(10); // in mm

								GetInst().SetcrystaltileLen(10);  //in mm
								GetInst().SetcrystalThickness(10); //in mm
								GetInst().SetpcbThickness(1.5);//in mm

								GetInst().SetWorldMult(10.);//in mm
								GetInst().InitDet(); // dont forget!!!




//	   auto ROworldSizeXY = 1.2 * GetInst().GetcalorSizeXY();
//	    auto ROworldSizeZ  = 1.2 * GetInst().GetcalorThickness();

								auto fCheckOverlaps=true;

								// A dummy material is used to fill the volumes of the readout geometry.
								// ( It will be allowed to set a NULL pointer in volumes of such virtual
								// division in future, since this material is irrelevant for tracking.)
								G4Material * dummyMat  = new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);



								//Build Readout World

								auto ROWorldS = new G4Box("ROWorld",
																																		GetInst().GetWorldSizeXY()/2,
																																		GetInst().GetWorldSizeXY()/2,
																																		GetInst().GetWorldSizeZ()/2);

								auto ROWorldLog = new G4LogicalVolume(ROWorldS,
																																														dummyMat,
																																														"ROWorldLogical",
																																														0,
																																														0,
																																														0);


								G4PVPlacement *ROWorld = new G4PVPlacement(0,
																																																			G4ThreeVector(),
																																																			"ROWorldPhysical",
																																																			ROWorldLog,
																																																			0,
																																																			false,
																																																			0);


								//Place Calorimeter

								auto calorimeterS
																= new G4Box("Calorimeter", // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetcalorThickness()/2); // its size

								auto calorLV
																= new G4LogicalVolume(
																calorimeterS, // its solid
																dummyMat, // its material
																"Calorimeter"); // its name

								G4ThreeVector seg1(0*mm,0*mm,((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcalorThickness()/2)+GetInst().GetcrystallayerThickness())*mm);
								new G4PVPlacement(
																0,                             // no rotation
																seg1,               // at (0,0,0)
																calorLV,                       // its logical volume
																"Calorimeter",                 // its name
																ROWorldLog,                       // its mother  volume
																false,                         // no boolean operation
																1,                             // copy number
																fCheckOverlaps);            // checking overlaps

								G4RotationMatrix * myrotation2= new G4RotationMatrix();
								myrotation2->rotateX(180.*deg);
								myrotation2->rotateY(0.*deg);
								myrotation2->rotateZ(0.*deg);
								G4ThreeVector seg2(0*mm,0*mm,-((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcalorThickness()/2)+GetInst().GetcrystallayerThickness())*mm);
								new G4PVPlacement(
																myrotation2,                                             // no rotation
																seg2,                               // at (0,0,0)
																calorLV,                                       // its logical volume
																"Calorimeter",                                 // its name
																ROWorldLog,                                       // its mother  volume
																false,                                         // no boolean operation
																2,                                             // copy number
																fCheckOverlaps);                               // checking overlaps
								G4RotationMatrix * myrotation3= new G4RotationMatrix();
								myrotation3->rotateX(90.*deg);
								myrotation3->rotateY(0.*deg);
								myrotation3->rotateZ(0.*deg);
								G4ThreeVector seg3(0*mm,((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcalorThickness()/2)+GetInst().GetcrystallayerThickness())*mm,0*mm);
								new G4PVPlacement(
																myrotation3,                                                     // no rotation
																seg3,                                       // at (0,0,0)
																calorLV,                                               // its logical volume
																"Calorimeter",                                         // its name
																ROWorldLog,                                               // its mother  volume
																false,                                                 // no boolean operation
																3,                                                     // copy number
																fCheckOverlaps);                                       // checking overlaps
								G4RotationMatrix * myrotation4= new G4RotationMatrix();
								myrotation4->rotateX(270.*deg);
								myrotation4->rotateY(0.*deg);
								myrotation4->rotateZ(0.*deg);
								G4ThreeVector seg4(0*mm,-((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcalorThickness()/2)+GetInst().GetcrystallayerThickness())*mm,0*mm);
								new G4PVPlacement(
																myrotation4,                                                             // no rotation
																seg4,                                               // at (0,0,0)
																calorLV,                                                       // its logical volume
																"Calorimeter",                                                 // its name
																ROWorldLog,                                                       // its mother  volume
																false,                                                         // no boolean operation
																4,                                                             // copy number
																fCheckOverlaps);                                               // checking overlaps

								G4RotationMatrix * myrotation5= new G4RotationMatrix();
								myrotation5->rotateX(0.*deg);
								myrotation5->rotateY(270.*deg);
								myrotation5->rotateZ(0.*deg);
								G4ThreeVector seg5(((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcalorThickness()/2)+GetInst().GetcrystallayerThickness())*mm,0*mm,0*mm);
								new G4PVPlacement(
																myrotation5,                                                                     // no rotation
																seg5,                                                       // at (0,0,0)
																calorLV,                                                               // its logical volume
																"Calorimeter",                                                         // its name
																ROWorldLog,                                                               // its mother  volume
																false,                                                                 // no boolean operation
																5,                                                                     // copy number
																fCheckOverlaps);                                                       // checking overlaps

								G4RotationMatrix * myrotation6= new G4RotationMatrix();
								myrotation6->rotateX(0.*deg);
								myrotation6->rotateY(90.*deg);
								myrotation6->rotateZ(0.*deg);
								G4ThreeVector seg6(-((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcalorThickness()/2)+GetInst().GetcrystallayerThickness())*mm,0*mm,0*mm);
								new G4PVPlacement(
																myrotation6,                                                                             // no rotation
																seg6,                                                               // at (0,0,0)
																calorLV,                                                                       // its logical volume
																"Calorimeter",                                                                 // its name
																ROWorldLog,                                                                       // its mother  volume
																false,                                                                         // no boolean operation
																6,                                                                             // copy number
																fCheckOverlaps);                                                               // checking overlaps



								// Build Layers

								auto ROlayerS
																= new G4Box("Layer",  // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetlayerThickness()/2); //its size

								auto ROlayerLV
																= new G4LogicalVolume(
																ROlayerS,               // its solid
																dummyMat,      // its material
																"Layer");             // its name

								new G4PVReplica(
																"Layer",              // its name
																ROlayerLV,              // its logical volume
																calorLV,              // its mother
																kZAxis,               // axis of replication
																GetInst().GetfNofLayers(),            // number of replica
																GetInst().GetlayerThickness());      // witdth of replica

								//-------------------------------
								//build calorimeter readout cells
								//-------------------------------

								//build Gap
								auto GapS
																=new G4Box("Gap", GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetgapThickness()/2);

								auto GapLV
																=new G4LogicalVolume(GapS,
																																					dummyMat,
																																					"Gap");
								new G4PVPlacement(0,
																										G4ThreeVector(0., 0., -GetInst().GetabsoThickness()/2 ),
																										GapLV,
																										"Gap",
																										ROlayerLV,
																										false,
																										0,
																										true);
								//build strips parallel to X axis

								auto StripS
																=new G4Box("Strip", GetInst().GetcalorSizeXY()/2,GetInst().GettileLenY()/2,GetInst().GetgapThickness()/2);

								auto StripLV
																=new G4LogicalVolume(StripS,
																																					dummyMat,
																																					"Strip");

								new G4PVReplica("Strip",
																								StripLV, //logical volume
																								GapLV, //mother volume
																								kYAxis,
																								GetInst().GetnofTilesY(),
																								GetInst().GettileLenY());


								// build cells in Xaxis strip

								auto CellS
																=new G4Box("Cell", GetInst().GettileLenX()/2, GetInst().GettileLenY()/2, GetInst().GetgapThickness()/2);

								auto CellLV
																=new G4LogicalVolume(CellS,
																																					dummyMat,
																																					"Cell");

								new G4PVReplica("Cell",
																								CellLV, //logical volume
																								StripLV, //mother volume
																								kXAxis,
																								GetInst().GetnofTilesX(),
																								GetInst().GettileLenX());




								auto RoCrystalLayerS
																= new G4Box("RoCrystalLayer",GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetcrystallayerThickness()/2);

								auto RoCrystalLayerLV
																=new G4LogicalVolume(
																RoCrystalLayerS,
																dummyMat,
																"RoCrystalLayer"
																);
								G4ThreeVector crystal1(0*mm,0*mm,((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcrystallayerThickness()/2))*mm);
								new G4PVPlacement(
																0,                                             // no rotation
																crystal1,                               // at (0,0,0)
																RoCrystalLayerLV,                                       // its logical volume
																"RoCrystalP",                                 // its name
																ROWorldLog,                                       // its mother  volume
																false,                                         // no boolean operation
																1,                                             // copy number
																fCheckOverlaps);                            // checking overlaps


								G4ThreeVector crystal2(0*mm,0*mm,-((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcrystallayerThickness()/2))*mm);
								new G4PVPlacement(
																myrotation2,                                                     // no rotation
																crystal2,                                       // at (0,0,0)
																RoCrystalLayerLV,                                               // its logical volume
																"RoCrystalP",                                         // its name
																ROWorldLog,                                               // its mother  volume
																false,                                                 // no boolean operation
																2,                                                     // copy number
																fCheckOverlaps);                                    // checking overlaps

								G4ThreeVector crystal3(0*mm,((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcrystallayerThickness()/2))*mm,0*mm);
								new G4PVPlacement(
																myrotation3,                                                             // no rotation
																crystal3,                                               // at (0,0,0)
																RoCrystalLayerLV,                                                       // its logical volume
																"RoCrystalP",                                                 // its name
																ROWorldLog,                                                       // its mother  volume
																false,                                                         // no boolean operation
																3,                                                             // copy number
																fCheckOverlaps);                                            // checking overlaps

								G4ThreeVector crystal4(0*mm,-((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcrystallayerThickness()/2))*mm,0*mm);
								new G4PVPlacement(
																myrotation4,                                                                     // no rotation
																crystal4,                                                       // at (0,0,0)
																RoCrystalLayerLV,                                                               // its logical volume
																"RoCrystalP",                                                         // its name
																ROWorldLog,                                                               // its mother  volume
																false,                                                                 // no boolean operation
																4,                                                                     // copy number
																fCheckOverlaps);                                                    // checking overlaps

								G4ThreeVector crystal5(((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcrystallayerThickness()/2))*mm,0*mm,0*mm);
								new G4PVPlacement(
																myrotation5,                                                     // no rotation
																crystal5,                                       // at (0,0,0)
																RoCrystalLayerLV,                                               // its logical volume
																"RoCrystalP",                                         // its name
																ROWorldLog,                                               // its mother  volume
																false,                                                 // no boolean operation
																5,                                                     // copy number
																fCheckOverlaps);                                    // checking overlaps

								G4ThreeVector crystal6(-((GetInst().GetcalorSizeXY()/2)+(GetInst().GetcrystallayerThickness()/2))*mm,0*mm,0*mm);
								new G4PVPlacement(
																myrotation6,                                                             // no rotation
																crystal6,                                               // at (0,0,0)
																RoCrystalLayerLV,                                                       // its logical volume
																"RoCrystalP",                                                 // its name
																ROWorldLog,                                                       // its mother  volume
																false,                                                         // no boolean operation
																6,                                                             // copy number
																fCheckOverlaps);                                            // checking overlaps


								auto RoCrystalS
																= new G4Box("RoCrystal", GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetcrystalThickness()/2);
								auto RoCrystalLV
																= new G4LogicalVolume(
																RoCrystalS,
																dummyMat,
																"RoCrystal"
																);

								new G4PVPlacement(
																0,                                             // no rotation
																G4ThreeVector(0., 0., -GetInst().GetpcbThickness()/2),                              // its position
																RoCrystalLV,                                    // its logical volume
																"RoCrystal",                                        // its name
																RoCrystalLayerLV,                                       // its mother  volume
																false,                                         // no boolean operation
																0,                                             // copy number
																fCheckOverlaps);                               // checking overlaps


								auto CrystalStripS
																=new G4Box("CrystalStrip", GetInst().GetcalorSizeXY()/2,GetInst().GetcrystaltileLen()/2,GetInst().GetcrystalThickness()/2);

								auto CrystalStripLV
																=new G4LogicalVolume(CrystalStripS,
																																					dummyMat,
																																					"CrystalStrip");

								new G4PVReplica("CrystalStrip",
																								CrystalStripLV,         //logical volume
																								RoCrystalLV,         //mother volume
																								kYAxis,
																								GetInst().GetnofcrystalTilesY(),
																								GetInst().GetcrystaltileLen());


								// build cells in Xaxis strip

								auto CrystalCellS
																=new G4Box("CrystalCell", GetInst().GetcrystaltileLen()/2,GetInst().GetcrystaltileLen()/2, GetInst().GetgapThickness()/2);

								auto CrystalCellLV
																=new G4LogicalVolume(CrystalCellS,
																																					dummyMat,
																																					"CrystalCell");

								new G4PVReplica("Cell",
																								CrystalCellLV,         //logical volume
																								CrystalStripLV,         //mother volume
																								kXAxis,
																								GetInst().GetnofcrystalTilesX(),
																								GetInst().GetcrystaltileLen());

								auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
								simpleBoxVisAtt->SetVisibility(true);
								CellLV->SetVisAttributes(simpleBoxVisAtt);



								auto dummy
																= new B4cCalorimeterSD("dummy", "dummyCollection",0,0,0,0);

								CellLV->SetSensitiveDetector(dummy);
								CrystalCellLV->SetSensitiveDetector(dummy);


								return ROWorld;




}
