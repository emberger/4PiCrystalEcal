#include "B4cDetParams.hh"



DetParams::~DetParams(){
}

void DetParams::InitDet(){

								nofTilesX=calorSizeXY/tileLenX;
								nofTilesY=calorSizeXY/tileLenY;

								nofcrystalTilesX=calorSizeXY/crystaltileLen;
								nofcrystalTilesY=calorSizeXY/crystaltileLen;

								layerThickness=absoThickness+gapThickness;
								crystallayerThickness=pcbThickness+crystalThickness;

								calorThickness=fNofLayers*layerThickness;

								tilesPerLayer=nofTilesX*nofTilesY;

								WorldSizeXY=4.0 * calorSizeXY;
								WorldSizeZ=4.0 * calorSizeXY;
}

void DetParams::SetfNofLayers(G4double nla){
								fNofLayers=nla;
}

G4double DetParams::GetfNofLayers(){
								return fNofLayers;
}

void DetParams::SettileLenX(G4double tx){
								tileLenX=tx;
}

G4double DetParams::GettileLenX(){
								return tileLenX;
}

void DetParams::SettileLenY(G4double ty){
								tileLenY=ty;
}

G4double DetParams::GettileLenY(){
								return tileLenY;
}

void DetParams::SetabsoThickness(G4double abso){
								absoThickness=abso;
}

G4double DetParams::GetabsoThickness(){
								return absoThickness;
}

void DetParams::SetgapThickness(G4double gap){
								gapThickness=gap;
}

G4double DetParams::GetgapThickness(){
								return gapThickness;
}


G4double DetParams::GetcrystallayerThickness(){
								return crystallayerThickness;
}

void DetParams::SetcrystaltileLen(G4double ctile){
								crystaltileLen=ctile;
}

G4double DetParams::GetcrystaltileLen(){
								return crystaltileLen;
}

void DetParams::SetcrystalThickness(G4double crystal){
								crystalThickness=crystal;
}

G4double DetParams::GetcrystalThickness(){
								return crystalThickness;
}

void DetParams::SetpcbThickness(G4double pcb){
								pcbThickness=pcb;
}

G4double DetParams::GetpcbThickness(){
								return pcbThickness;
}

void DetParams::SetcalorSizeXY(G4double cs){
								calorSizeXY=cs;
}

G4double DetParams::GetcalorSizeXY(){
								return calorSizeXY;
}

void DetParams::SetWorldMult(G4double wm){
								WorldMult=wm;
}



G4double DetParams::GetnofTilesX(){
								return nofTilesX;
}
G4double DetParams::GetnofTilesY(){
								return nofTilesY;
}
G4double DetParams::GetnofcrystalTilesX(){
								return nofcrystalTilesX;
}
G4double DetParams::GetnofcrystalTilesY(){
								return nofcrystalTilesY;
}
G4double DetParams::GetlayerThickness(){
								return layerThickness;
}
G4double DetParams::GetcalorThickness(){
								return calorThickness;
}
G4double DetParams::GettilesPerLayer(){
								return tilesPerLayer;
}
G4double DetParams::GetWorldSizeXY(){
								return WorldSizeXY;
}
G4double DetParams::GetWorldSizeZ(){
								return WorldSizeZ;
}
