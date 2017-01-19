
{
    TFile f("~/data/hists2012.root");
    ProtonCalibratedSource0_2->SetLineColor(1);
    ProtonCalibratedSource0_3->SetLineColor(2);
    ProtonCalibratedSource3_1->SetLineColor(3);
    ProtonCalibratedSource3_2->SetLineColor(4);

    ProtonCalibratedSource0_2->Draw();
    ProtonCalibratedSource0_3->Draw("same");
    ProtonCalibratedSource3_1->Draw("same");
    ProtonCalibratedSource3_2->Draw("same");












}
