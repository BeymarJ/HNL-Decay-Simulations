
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include <cmath> 

void HNL_GenLF() {

   Double_t m_N = 150.0;    //HNL mass
   
   // Simulate HNL from KDAR
   Double_t m_K = 493.677;
   Double_t m_mu = 105.658;
  
  
   Double_t m_z = 91187.6;
   Double_t m_e = 0.511;
   
   Double_t EN = (m_K*m_K - m_mu*m_mu + m_N*m_N)/(2*m_K);
   Double_t pzN = sqrt(EN*EN - m_N*m_N);
   
   
   TLorentzVector pN(0.0, 0.0, pzN, EN); 
   TLorentzVector *p1 = new TLorentzVector(0.0, 0.0, pzN, EN);
   Double_t masses[3] = {0.0, m_e, m_e};
   TVector3 z(0.0, 0.0, 1.0);
   
   TGenPhaseSpace event;
   event.SetDecay(pN, 3, masses);
   
   TH1D *hEv = new TH1D("hEv","Energy nu",250,0.0,250);
   TH1D *hEe3 = new TH1D("hEe3","Energy e+",250,0.0,250); //e3 positron
   TH1D *hEe4 = new TH1D("hEe4","Energy e-",250,0.0,250); //e4 electron
   
   TH2D *htL_tll = new TH2D("htL_tll","tll vs tLead", 180, 0.0, 180.0, 180, 0.0, 180.0);
   // leadind particle angle vs angle between e+ e-
   
   
   Double_t sw2 = 0.22290;
   Double_t a = 4.0*sw2 - 1;
   
   Double_t Ev, E3, E4, Num, Den, weight;
   
   for (Int_t n=0; n<100000;n++) {
      
      Double_t W0 = event.Generate();
      TLorentzVector *pNu = event.GetDecay(0);
      TLorentzVector *pE3 = event.GetDecay(1);
      TLorentzVector *pE4 = event.GetDecay(2);
      
      Ev = pNu->E();
      E3 = pE3->E();
      E4 = pE4->E();
      
      Double_t tll = pE3->Angle(pE4->Vect()) * 180./3.141592;
      Double_t tLead0;
      if ( E3 > E4) {
         tLead0 = pE3->Angle(z);
      } else {
         tLead0 = pE4->Angle(z);
      }
      Double_t tLead = tLead0 * 180./3.141592;
      
      Den = (m_N*m_N - 2*p1->Dot(*pNu) - m_z*m_z)*(m_N*m_N - 2*p1->Dot(*pNu) - m_z*m_z);
      Num = m_e*m_e*(a*a-1)*(p1->Dot(*pNu)) + (a*a+1)*((pNu->Dot(*pE4))*(p1->Dot(*pE3))+(pNu->Dot(*pE3))*(p1->Dot(*pE4)));
      
      
      
      weight = W0 * Num / Den; 
      
      
      
      hEv->Fill(Ev, weight);
      hEe3->Fill(E3, weight);
      hEe4->Fill(E4, weight);
      htL_tll->Fill(tll,tLead,weight);
      
   }
   
   
   TFile *g = new TFile("HNL150LF-Sim.root","recreate");
   
   hEv->Write();
   hEe3->Write();
   hEe4->Write();
   htL_tll->Write();
   
   g->Close();

}
