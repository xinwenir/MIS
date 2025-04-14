using namespace std;
//#include "Hists.h"
//#include "TVirtualPad.h"

void CanvasPartition(TCanvas *C,Int_t Nx,Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
   for (Int_t i=0;i<Nx;i++) {
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
      for (Int_t j=0;j<Ny;j++) {
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
         C->cd(0);
         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
         pad->Draw();
      }
   }
}

void canvas_multiple_pad(TCanvas* C, int nx, int ny, TList *lis_h)
{
   C->SetFillStyle(4000);
   // Number of PADS
   const Int_t Nx = nx;
   const Int_t Ny = ny;


   TPad *pad[Nx][Ny];
   for (Int_t i=0;i<Nx;i++) {
      for (Int_t j=0;j<Ny;j++) {
         C->cd(0);
         TH1F *h = (TH1F*) lis_h->At(i*Ny+j);
         // Get the pads previously created.
         char pname[16];
         sprintf(pname,"pad_%i_%i",i,j);
         pad[i][j] = (TPad*) gROOT->FindObject(pname);
         if ((i*Ny+j) == 0) pad[i][j]->SetPad(0., 0.3, 0.8, 1.);
         if ((i*Ny+j) == 1) pad[i][j]->SetPad(0., 0., 1., 0.3);

         pad[i][j]->Draw();
         pad[i][j]->SetFillStyle(4000);
         pad[i][j]->SetFrameFillStyle(4000);
         pad[i][j]->cd();

         // Size factors
//         Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
//         Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();
//         char hname[16];
//         sprintf(hname,"h_%i_%i",i,j);
//         TH1F *hFrame = (TH1F*) h->Clone(hname);
//         h->Reset();
         h->Draw();
//         // y axis range
//         h->GetYaxis()->SetRangeUser(0.0001,1.2*h->GetMaximum());
//         // Format for y axis
//         h->GetYaxis()->SetLabelFont(43);
//         h->GetYaxis()->SetLabelSize(16);
//         h->GetYaxis()->SetLabelOffset(0.02);
//         h->GetYaxis()->SetTitleFont(43);
//         h->GetYaxis()->SetTitleSize(16);
//         h->GetYaxis()->SetTitleOffset(5);
//         h->GetYaxis()->CenterTitle();
//         h->GetYaxis()->SetNdivisions(505);
//         // TICKS Y Axis
//         h->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
//         // Format for x axis
//         h->GetXaxis()->SetLabelFont(43);
//         h->GetXaxis()->SetLabelSize(16);
//         h->GetXaxis()->SetLabelOffset(0.02);
//         h->GetXaxis()->SetTitleFont(43);
//         h->GetXaxis()->SetTitleSize(16);
//         h->GetXaxis()->SetTitleOffset(5);
//         h->GetXaxis()->CenterTitle();
//         h->GetXaxis()->SetNdivisions(505);
//         // TICKS X Axis
//         h->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
//         h->Draw("");
      }
   }
   C->cd();
}