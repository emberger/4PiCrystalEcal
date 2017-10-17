// Implementation file for Pi0Event and components
// create dictionary with rootcint RootClasses_dict.cpp -c include/B4ROOTEvent.hh include/LinkDef.h
//

#include "B4ROOTEvent.hh"
//#include "Stiostream.h"




TClonesArray* B4ROOTEvent::aHits = 0;

// TPi0Event class

B4ROOTEvent::B4ROOTEvent() {
								if (!aHits) aHits = new TClonesArray("B4ROOTHit", 1000);
								m_Hits = aHits;
								m_EventNo = 0;
								m_GapEnergy = 0;
								m_NHits = 0;
								m_GapThickness=0;
								m_AbsoThickness=0;
								m_crystalsizeXY=0;
								m_crystalThickness=0;
								m_pcbThickness=0;
								m_Layerno=0;
								m_tilesizeX=0;
								m_tilesizeY=0;
								m_calsizeXY=0;
								m_EnergyPhoton1=0;
								m_EnergyPhoton2=0;
								m_MomentumPh1.SetXYZ(0.,0.,0.);
								m_MomentumPh2.SetXYZ(0.,0.,0.);
								m_GunPos.SetXYZ(0.,0.,0.);

}

B4ROOTEvent::B4ROOTEvent(const B4ROOTEvent& origEvent) : TObject(origEvent) {

								//std::cout<<"CC called"<<std::endl;
								m_Hits=origEvent.m_Hits;
								m_EventNo=origEvent.m_EventNo;
								m_GapEnergy=origEvent.m_GapEnergy;
								m_NHits=origEvent.m_NHits;
								m_GapThickness=origEvent.m_GapThickness;
								m_AbsoThickness=origEvent.m_AbsoThickness;
								m_Layerno=origEvent.m_Layerno;
								m_tilesizeX=origEvent.m_tilesizeX;
								m_tilesizeY=origEvent.m_tilesizeY;
								m_calsizeXY=origEvent.m_calsizeXY;

								m_EnergyPhoton1=origEvent.m_EnergyPhoton1;
								m_EnergyPhoton2=origEvent.m_EnergyPhoton2;
								m_MomentumPh1=origEvent.m_MomentumPh1;
								m_MomentumPh2=origEvent.m_MomentumPh2;
								m_GunPos=origEvent.m_GunPos;

								m_segmentPh1=origEvent.m_segmentPh1;
								m_segmentPh2=origEvent.m_segmentPh2;


}

B4ROOTEvent::~B4ROOTEvent() {
								Clear();
}

void B4ROOTEvent::Clear(const Option_t*) {
								//std::cout << "Clearing Event..." << std::endl;
								m_EventNo = 0;
								m_GapEnergy = 0;
								m_NHits=0;
								m_GapThickness=0;
								m_AbsoThickness=0;
								m_Layerno=0;
								m_tilesizeX=0;
								m_tilesizeY=0;
								m_calsizeXY=0;
								m_EnergyPhoton1=0;
								m_EnergyPhoton2=0;
								m_MomentumPh1.SetXYZ(0.,0.,0.);
								m_MomentumPh2.SetXYZ(0.,0.,0.);
								m_GunPos.SetXYZ(0.,0.,0.);
								m_Hits->Clear();
								//std::cout << "Clear done!" << std::endl;
}
void B4ROOTEvent::Reset(Option_t * /*option*/)
{
								// Static function to reset all static objects for this event
								delete aHits; aHits = 0;
}




B4ROOTHit* B4ROOTEvent::AddHit(B4ROOTHit& cand) {
								//cout << "Adding Pi0 candidate after " << m_NPi0Candidates << endl;
								TClonesArray &aCand = *m_Hits;
								B4ROOTHit *piCand = new B4ROOTHit(cand);
								aCand[m_NHits++]=piCand;
								//std::cout<<"AddHitCand: "<<cand.PhotNr()<<" AddHitPiCand: "<<piCand->PhotNr()<<std::endl;
								return piCand;
}


B4ROOTHit::B4ROOTHit() : TObject() {
								m_X = 0;
								m_Y = 0;
								m_Z = 0;
								m_EnergyDeposit = 0;
								m_PhotonNr=0;
								m_calorimeterSegment=0;

								m_CrystalHit=false;

}

B4ROOTHit::B4ROOTHit(const B4ROOTHit& orig) : TObject(orig) {
								m_X = orig.m_X;
								m_Y = orig.m_Y;
								m_Z = orig.m_Z;
								//std::cout<<"orig.Phot: "<<orig.m_PhotonNr<<std::endl;
								m_PhotonNr = orig.m_PhotonNr;
								m_EnergyDeposit = orig.m_EnergyDeposit;
								m_calorimeterSegment=orig.m_calorimeterSegment;
								m_CrystalHit=orig.m_CrystalHit;
}

B4ROOTHit::~B4ROOTHit(){
}
