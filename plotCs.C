
{
    TFile f("~/data/hists2012.root");
    ProtonCalibratedSource0_0->SetLineColor(1);
    ProtonCalibratedSource1_0->SetLineColor(2);
    ProtonCalibratedSource3_0->SetLineColor(3);

    ProtonCalibratedSource0_0->Draw();
    ProtonCalibratedSource1_0->Draw("same");
    ProtonCalibratedSource3_0->Draw("same");












}
