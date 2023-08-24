{
	/*
	  This code is for plotting the kinematics of 
	  the output file NOT in the sampling mode!!!
	*/

	gStyle->SetOptStat(0);
	gStyle->SetNdivisions(505,"XY");


	TFile file("DEMP-pion-pole-at-EicC.root");

	TCanvas c1("c1","Q2 vs xB, unweighted",700,550);
	tree->Draw("Q2:xB>>hQ2xB","","colz");
	gPad->SetLogz();
	hQ2xB->SetTitle("");
	hQ2xB->GetXaxis()->SetTitle("x_{B}");
	hQ2xB->GetXaxis()->SetTitleSize(0.06);
	hQ2xB->GetXaxis()->CenterTitle();
	hQ2xB->GetXaxis()->SetTitleOffset(1.05);
	hQ2xB->GetXaxis()->SetLabelSize(0.06);
	hQ2xB->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
	hQ2xB->GetYaxis()->SetTitleSize(0.06);
	hQ2xB->GetYaxis()->CenterTitle();
	hQ2xB->GetYaxis()->SetTitleOffset(1.05);
	hQ2xB->GetYaxis()->SetLabelSize(0.06);
	hQ2xB->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c2("c2","-t vs xB, unweighted",700,550);
	tree->Draw("-t:xB>>htxB","","colz");
	gPad->SetLogz();
	htxB->SetTitle("");
	htxB->GetXaxis()->SetTitle("x_{B}");
	htxB->GetXaxis()->SetTitleSize(0.06);
	htxB->GetXaxis()->CenterTitle();
	htxB->GetXaxis()->SetTitleOffset(1.05);
	htxB->GetXaxis()->SetLabelSize(0.06);
	htxB->GetYaxis()->SetTitle("-t (GeV^{2})");
	htxB->GetYaxis()->SetTitleSize(0.06);
	htxB->GetYaxis()->CenterTitle();
	htxB->GetYaxis()->SetTitleOffset(1.05);
	htxB->GetYaxis()->SetLabelSize(0.06);
	htxB->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c3("c3","Q2 vs xB, xsec. weighted",700,550);
	tree->Draw("Q2:xB>>hQ2xB_2","d4sigma*1000","colz");
	gPad->SetLogz();
	hQ2xB_2->SetTitle("");
	hQ2xB_2->GetXaxis()->SetTitle("x_{B}");
	hQ2xB_2->GetXaxis()->SetTitleSize(0.06);
	hQ2xB_2->GetXaxis()->CenterTitle();
	hQ2xB_2->GetXaxis()->SetTitleOffset(1.05);
	hQ2xB_2->GetXaxis()->SetLabelSize(0.06);
	hQ2xB_2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
	hQ2xB_2->GetYaxis()->SetTitleSize(0.06);
	hQ2xB_2->GetYaxis()->CenterTitle();
	hQ2xB_2->GetYaxis()->SetTitleOffset(1.05);
	hQ2xB_2->GetYaxis()->SetLabelSize(0.06);
	hQ2xB_2->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c4("c4","-t vs xB, xsec. weighted",700,550);
	tree->Draw("-t:xB>>htxB_2","d4sigma*1000","colz");
	gPad->SetLogz();
	htxB_2->SetTitle("");
	htxB_2->GetXaxis()->SetTitle("x_{B}");
	htxB_2->GetXaxis()->SetTitleSize(0.06);
	htxB_2->GetXaxis()->CenterTitle();
	htxB_2->GetXaxis()->SetTitleOffset(1.05);
	htxB_2->GetXaxis()->SetLabelSize(0.06);
	htxB_2->GetYaxis()->SetTitle("-t (GeV^{2})");
	htxB_2->GetYaxis()->SetTitleSize(0.06);
	htxB_2->GetYaxis()->CenterTitle();
	htxB_2->GetYaxis()->SetTitleOffset(1.05);
	htxB_2->GetYaxis()->SetLabelSize(0.06);
	htxB_2->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c5("c5","electron P vs theta, xsec. weighted",700,550);
	tree->Draw("elec_out.P():elec_out.Theta()*TMath::RadToDeg()>>helec","d4sigma*1000","colz");
	gPad->SetLogz();
	helec->SetTitle("");
	helec->GetXaxis()->SetTitle("#theta_e (#circ)");
	helec->GetXaxis()->SetTitleSize(0.06);
	helec->GetXaxis()->CenterTitle();
	helec->GetXaxis()->SetTitleOffset(1.05);
	helec->GetXaxis()->SetLabelSize(0.06);
	helec->GetYaxis()->SetTitle("P_e (GeV/c)");
	helec->GetYaxis()->SetTitleSize(0.06);
	helec->GetYaxis()->CenterTitle();
	helec->GetYaxis()->SetTitleOffset(1.05);
	helec->GetYaxis()->SetLabelSize(0.06);
	helec->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c6("c6","pion plus P vs theta, xsec. weighted",700,550);
	tree->Draw("pip_out.P():pip_out.Theta()*TMath::RadToDeg()>>hpip","d4sigma*1000","colz");
	gPad->SetLogz();
	hpip->SetTitle("");
	hpip->GetXaxis()->SetTitle("#theta of #pi+ (#circ)");
	hpip->GetXaxis()->SetTitleSize(0.06);
	hpip->GetXaxis()->CenterTitle();
	hpip->GetXaxis()->SetTitleOffset(1.05);
	hpip->GetXaxis()->SetLabelSize(0.06);
	hpip->GetYaxis()->SetTitle("P_#pi+ (GeV/c)");
	hpip->GetYaxis()->SetTitleSize(0.06);
	hpip->GetYaxis()->CenterTitle();
	hpip->GetYaxis()->SetTitleOffset(1.05);
	hpip->GetYaxis()->SetLabelSize(0.06);
	hpip->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c7("c7","neutron, P vs theta, xsec. weighted",700,550);
	tree->Draw("neut_out.P():neut_out.Theta()*TMath::RadToDeg()>>hneut","d4sigma*1000","colz");
	gPad->SetLogz();
	hneut->SetTitle("");
	hneut->GetXaxis()->SetTitle("#theta_n (#circ)");
	hneut->GetXaxis()->SetTitleSize(0.06);
	hneut->GetXaxis()->CenterTitle();
	hneut->GetXaxis()->SetTitleOffset(1.05);
	hneut->GetXaxis()->SetLabelSize(0.06);
	hneut->GetYaxis()->SetTitle("P_n (GeV/c)");
	hneut->GetYaxis()->SetTitleSize(0.06);
	hneut->GetYaxis()->CenterTitle();
	hneut->GetYaxis()->SetTitleOffset(1.05);
	hneut->GetYaxis()->SetLabelSize(0.06);
	hneut->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);





}
