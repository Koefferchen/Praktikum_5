{
	UInt_t messungen = 4;
	TString histogrammName = "Driftzeiten";
	TString dateien[4] = { "nachmittags/2600V1FDDC10THRhist.root", "nachmittags/2700V1FDDC10THRhist.root", "nachmittags/2750V1FDDC10THRhist.root", "nachmittags/2850V1FDDC10THRhist.root"};
 	TString titel[4] = { "2600V10THR", "2700V10THR", "2750V10THR", "2850V10THR"};

	TLegend* leg = new TLegend(0.6, 0.5, 0.9, 0.7);
	leg->SetHeader("Spannungen");

	Bool_t first=true;
	UInt_t num = messungen;
	do {
		--num;
		TFile::Open(dateien[num]);
		TH1* plot = static_cast<TH1*>(gDirectory->FindObjectAny(histogrammName));
		plot->SetLineColor(num+1);
		if (first) {
			plot->Draw();
			first=false;
		} else {
			plot->Draw("same");
		}
		leg->AddEntry(plot, titel[num], "lep");
	} while (num != 0);

	leg->Draw("SAME");
}
