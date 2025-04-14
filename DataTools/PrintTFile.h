void print_tfile(string fin_name, string pdf_name) {
  TString pdf_name1 = Form("%s(", pdf_name.c_str());
  TString pdf_name2 = Form("%s]", pdf_name.c_str());

  TFile * fin = new TFile(fin_name.c_str(), "read");
  TKey *key;
  TIter next(fin->GetListOfKeys());
  TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1500);
   while ((key = (TKey*)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (cl->InheritsFrom("TCanvas"))
      {
        TCanvas *c = (TCanvas*)key->ReadObj();
        c->Draw();
        c->Print(pdf_name1.Data());
      }
      else if (cl->InheritsFrom("TH1F"))
      {
        c1->cd();
        cl->Draw("hist e");
        c1->Print(pdf_name1.Data());
      }
      else if (cl->InheritsFrom("TH2F"))
      {
        c1->cd();
        cl->Draw("colz");
        c1->Print(pdf_name1.Data());
      }
      else
        continue;
   }
   c1->Print(pdf_name2.Data());

}
