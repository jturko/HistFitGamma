
{
    TFile f("~/data/hists2012.root");
    ProtonCalibratedSource0_1->SetLineColor(1);
    ProtonCalibratedSource1_1->SetLineColor(2);
    ProtonCalibratedSource1_2->SetLineColor(3);
    ProtonCalibratedSource2_0->SetLineColor(4);
    ProtonCalibratedSource2_1->SetLineColor(5);

    ProtonCalibratedSource0_1->Draw();
    ProtonCalibratedSource1_1->Draw("same");
    ProtonCalibratedSource1_2->Draw("same");
    ProtonCalibratedSource2_0->Draw("same");
    ProtonCalibratedSource2_1->Draw("same");












}
