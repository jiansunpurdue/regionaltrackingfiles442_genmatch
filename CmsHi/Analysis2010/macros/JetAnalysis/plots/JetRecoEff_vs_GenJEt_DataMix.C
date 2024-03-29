{
//=========Macro generated from canvas: c1/
//=========  (Mon Jan  3 08:16:49 2011) by ROOT version5.22/00d
   TCanvas *c1 = new TCanvas("c1", "",10,73,1250,500);
   gStyle->SetOptStat(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(0);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.17);
   c1->SetRightMargin(0.04);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.14);
   c1->SetFrameLineColor(0);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: p_0_0
   TPad *p_0_0 = new TPad("p_0_0", "p_0_0",0,0,0.388211,0.9831461);
   p_0_0->Draw();
   p_0_0->cd();
   p_0_0->Range(-38.97436,-0.2731707,220,1.434146);
   p_0_0->SetFillColor(0);
   p_0_0->SetBorderMode(0);
   p_0_0->SetBorderSize(0);
   p_0_0->SetTickx(1);
   p_0_0->SetTicky(1);
   p_0_0->SetLeftMargin(0.22);
   p_0_0->SetRightMargin(0);
   p_0_0->SetTopMargin(0.02);
   p_0_0->SetBottomMargin(0.16);
   p_0_0->SetFrameLineColor(0);
   p_0_0->SetFrameBorderMode(0);
   p_0_0->SetFrameLineColor(0);
   p_0_0->SetFrameBorderMode(0);
   Double_t xAxis1[18] = {18, 21, 24, 28, 32, 37, 43, 50, 56, 64, 74, 84, 97, 114, 133, 174, 220, 300}; 
   
   TH1 *hTmp = new TH1F("hTmp","",17, xAxis1);
   hTmp->SetMinimum(0);
   hTmp->SetMaximum(1.4);
   hTmp->SetDirectory(0);
   hTmp->SetStats(0);
   hTmp->SetLineStyle(0);
   hTmp->SetLineWidth(2);
   hTmp->SetMarkerStyle(20);
   hTmp->SetMarkerSize(1.3);
   hTmp->GetXaxis()->SetTitle("p_{T}^{GenJet} (GeV/c)");
   hTmp->GetXaxis()->SetRange(1,16);
   hTmp->GetXaxis()->CenterTitle(true);
   hTmp->GetXaxis()->SetNdivisions(505);
   hTmp->GetXaxis()->SetLabelFont(63);
   hTmp->GetXaxis()->SetLabelOffset(0.01);
   hTmp->GetXaxis()->SetLabelSize(24);
   hTmp->GetXaxis()->SetTitleSize(28);
   hTmp->GetXaxis()->SetTitleOffset(1.2);
   hTmp->GetXaxis()->SetTitleFont(63);
   hTmp->GetYaxis()->SetTitle("Eff. (Matched Reco/Gen)");
   hTmp->GetYaxis()->CenterTitle(true);
   hTmp->GetYaxis()->SetLabelFont(63);
   hTmp->GetYaxis()->SetLabelOffset(0.01);
   hTmp->GetYaxis()->SetLabelSize(24);
   hTmp->GetYaxis()->SetTitleSize(24);
   hTmp->GetYaxis()->SetTitleOffset(1.8);
   hTmp->GetYaxis()->SetTitleFont(63);
   hTmp->GetZaxis()->SetLabelFont(42);
   hTmp->GetZaxis()->SetTitleFont(42);
   hTmp->Draw("");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10);
   grae->SetName("");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.3);
   grae->SetPoint(0,53,0.9847215);
   grae->SetPointError(0,0,0,0.001403303,0.00132349);
   grae->SetPoint(1,60,0.993087);
   grae->SetPointError(1,0,0,0.0009460755,0.0008674256);
   grae->SetPoint(2,69,0.9972072);
   grae->SetPointError(2,0,0,0.0007361501,0.0006261397);
   grae->SetPoint(3,79,0.9983649);
   grae->SetPointError(3,0,0,0.0008593427,0.0006351309);
   grae->SetPoint(4,90.5,0.9989744);
   grae->SetPointError(4,0,0,0.0009502488,0.000582359);
   grae->SetPoint(5,105.5,1);
   grae->SetPointError(5,0,0,0.001119114,0);
   grae->SetPoint(6,123.5,1);
   grae->SetPointError(6,0,0,0.002441384,0);
   grae->SetPoint(7,153.5,1);
   grae->SetPointError(7,0,0,0.003994973,0);
   grae->SetPoint(8,197,1);
   grae->SetPointError(8,0,0,0.01865742,0);
   grae->SetPoint(9,260,1);
   grae->SetPointError(9,0,0,0.09129775,0);
   
   TH1 *Graph1 = new TH1F("Graph1","",100,32.3,280.7);
   Graph1->SetMinimum(0.8995725);
   Graph1->SetMaximum(1.00913);
   Graph1->SetDirectory(0);
   Graph1->SetStats(0);
   Graph1->SetLineStyle(0);
   Graph1->SetLineWidth(2);
   Graph1->SetMarkerStyle(20);
   Graph1->SetMarkerSize(1.3);
   Graph1->GetXaxis()->SetLabelFont(42);
   Graph1->GetXaxis()->SetLabelOffset(0.01);
   Graph1->GetXaxis()->SetTitleSize(0.045);
   Graph1->GetXaxis()->SetTitleOffset(1.2);
   Graph1->GetXaxis()->SetTitleFont(42);
   Graph1->GetYaxis()->SetLabelFont(42);
   Graph1->GetYaxis()->SetLabelOffset(0.01);
   Graph1->GetYaxis()->SetTitleSize(0.045);
   Graph1->GetYaxis()->SetTitleOffset(1.8);
   Graph1->GetYaxis()->SetTitleFont(42);
   Graph1->GetZaxis()->SetLabelFont(42);
   Graph1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph1);
   
   grae->Draw("p");
   
   grae = new TGraphAsymmErrors(17);
   grae->SetName("");
   grae->SetTitle("");
   grae->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(25);
   grae->SetMarkerSize(1.3);
   grae->SetPoint(0,19.5,0.1666667);
   grae->SetPointError(0,0,0,0.04400851,0.05126909);
   grae->SetPoint(1,22.5,0.2125984);
   grae->SetPointError(1,0,0,0.03461512,0.03760129);
   grae->SetPoint(2,26,0.3743719);
   grae->SetPointError(2,0,0,0.02400555,0.0244251);
   grae->SetPoint(3,30,0.5791411);
   grae->SetPointError(3,0,0,0.01734826,0.01721891);
   grae->SetPoint(4,34.5,0.7448609);
   grae->SetPointError(4,0,0,0.01081913,0.01062168);
   grae->SetPoint(5,40,0.8644918);
   grae->SetPointError(5,0,0,0.006141073,0.005988488);
   grae->SetPoint(6,46.5,0.9447906);
   grae->SetPointError(6,0,0,0.002794341,0.002708286);
   grae->SetPoint(7,53,0.9724518);
   grae->SetPointError(7,0,0,0.00213693,0.002034663);
   grae->SetPoint(8,60,0.9893255);
   grae->SetPointError(8,0,0,0.001627296,0.001478511);
   grae->SetPoint(9,69,0.9935345);
   grae->SetPointError(9,0,0,0.001647931,0.001409405);
   grae->SetPoint(10,79,0.9970037);
   grae->SetPointError(10,0,0,0.001795377,0.001279786);
   grae->SetPoint(11,90.5,0.9987421);
   grae->SetPointError(11,0,0,0.001883317,0.0009199984);
   grae->SetPoint(12,105.5,1);
   grae->SetPointError(12,0,0,0.002494391,0);
   grae->SetPoint(13,123.5,1);
   grae->SetPointError(13,0,0,0.005643393,0);
   grae->SetPoint(14,153.5,1);
   grae->SetPointError(14,0,0,0.01029665,0);
   grae->SetPoint(15,197,1);
   grae->SetPointError(15,0,0,0.03638148,0);
   grae->SetPoint(16,260,1);
   grae->SetPointError(16,0,0,0.2052842,0);
   
   TH1 *Graph2 = new TH1F("Graph2","",100,0,284.05);
   Graph2->SetMinimum(0.03492397);
   Graph2->SetMaximum(1.087734);
   Graph2->SetDirectory(0);
   Graph2->SetStats(0);
   Graph2->SetLineStyle(0);
   Graph2->SetLineWidth(2);
   Graph2->SetMarkerStyle(20);
   Graph2->SetMarkerSize(1.3);
   Graph2->GetXaxis()->SetLabelFont(42);
   Graph2->GetXaxis()->SetLabelOffset(0.01);
   Graph2->GetXaxis()->SetTitleSize(0.045);
   Graph2->GetXaxis()->SetTitleOffset(1.2);
   Graph2->GetXaxis()->SetTitleFont(42);
   Graph2->GetYaxis()->SetLabelFont(42);
   Graph2->GetYaxis()->SetLabelOffset(0.01);
   Graph2->GetYaxis()->SetTitleSize(0.045);
   Graph2->GetYaxis()->SetTitleOffset(1.8);
   Graph2->GetYaxis()->SetTitleFont(42);
   Graph2->GetZaxis()->SetLabelFont(42);
   Graph2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph2);
   
   grae->Draw("p");
   TLine *line = new TLine(18,1,300,1);
   line->SetLineStyle(2);
   line->Draw();
   TLatex *   tex = new TLatex(0.28,0.9,"CMS");
tex->SetNDC();
   tex->SetTextFont(63);
   tex->SetTextSize(18);
   tex->Draw();
      tex = new TLatex(0.54,0.51,"30-100%");
tex->SetNDC();
   tex->SetTextFont(63);
   tex->SetTextSize(18);
   tex->Draw();
   p_0_0->Modified();
   c1->cd();
  
// ------------>Primitives in pad: p_1_0
   p_1_0 = new TPad("p_1_0", "p_1_0",0.388211,0,0.6910157,0.9831461);
   p_1_0->Draw();
   p_1_0->cd();
   p_1_0->Range(18,-0.2731707,220,1.434146);
   p_1_0->SetFillColor(0);
   p_1_0->SetBorderMode(0);
   p_1_0->SetBorderSize(0);
   p_1_0->SetTickx(1);
   p_1_0->SetTicky(1);
   p_1_0->SetLeftMargin(0);
   p_1_0->SetRightMargin(0);
   p_1_0->SetTopMargin(0.02);
   p_1_0->SetBottomMargin(0.16);
   p_1_0->SetFrameLineColor(0);
   p_1_0->SetFrameBorderMode(0);
   p_1_0->SetFrameLineColor(0);
   p_1_0->SetFrameBorderMode(0);
   Double_t xAxis2[18] = {18, 21, 24, 28, 32, 37, 43, 50, 56, 64, 74, 84, 97, 114, 133, 174, 220, 300}; 
   
   TH1 *hTmp = new TH1F("hTmp","",17, xAxis2);
   hTmp->SetMinimum(0);
   hTmp->SetMaximum(1.4);
   hTmp->SetDirectory(0);
   hTmp->SetStats(0);
   hTmp->SetLineStyle(0);
   hTmp->SetLineWidth(2);
   hTmp->SetMarkerStyle(20);
   hTmp->SetMarkerSize(1.3);
   hTmp->GetXaxis()->SetTitle("p_{T}^{GenJet} (GeV/c)");
   hTmp->GetXaxis()->SetRange(1,16);
   hTmp->GetXaxis()->CenterTitle(true);
   hTmp->GetXaxis()->SetNdivisions(505);
   hTmp->GetXaxis()->SetLabelFont(63);
   hTmp->GetXaxis()->SetLabelOffset(0.01);
   hTmp->GetXaxis()->SetLabelSize(24);
   hTmp->GetXaxis()->SetTitleSize(28);
   hTmp->GetXaxis()->SetTitleOffset(1.2);
   hTmp->GetXaxis()->SetTitleFont(63);
   hTmp->GetYaxis()->SetTitle("Eff. (Matched Reco/Gen)");
   hTmp->GetYaxis()->CenterTitle(true);
   hTmp->GetYaxis()->SetLabelFont(63);
   hTmp->GetYaxis()->SetLabelOffset(0.01);
   hTmp->GetYaxis()->SetLabelSize(24);
   hTmp->GetYaxis()->SetTitleSize(24);
   hTmp->GetYaxis()->SetTitleOffset(1.8);
   hTmp->GetYaxis()->SetTitleFont(63);
   hTmp->GetZaxis()->SetLabelFont(42);
   hTmp->GetZaxis()->SetTitleFont(42);
   hTmp->Draw("");
   
   grae = new TGraphAsymmErrors(10);
   grae->SetName("");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.3);
   grae->SetPoint(0,53,0.9763214);
   grae->SetPointError(0,0,0,0.003267653,0.002998279);
   grae->SetPoint(1,60,0.9847799);
   grae->SetPointError(1,0,0,0.002624307,0.002357206);
   grae->SetPoint(2,69,0.9955752);
   grae->SetPointError(2,0,0,0.00176896,0.001396347);
   grae->SetPoint(3,79,0.9930314);
   grae->SetPointError(3,0,0,0.003276664,0.00249432);
   grae->SetPoint(4,90.5,0.9966102);
   grae->SetPointError(4,0,0,0.003129135,0.001922157);
   grae->SetPoint(5,105.5,1);
   grae->SetPointError(5,0,0,0.003735201,0);
   grae->SetPoint(6,123.5,1);
   grae->SetPointError(6,0,0,0.008350729,0);
   grae->SetPoint(7,153.5,1);
   grae->SetPointError(7,0,0,0.01480942,0);
   grae->SetPoint(8,197,1);
   grae->SetPointError(8,0,0,0.04491424,0);
   grae->SetPoint(9,260,1);
   grae->SetPointError(9,0,0,0.2496484,0);
   
   TH1 *Graph3 = new TH1F("Graph3","",100,32.3,280.7);
   Graph3->SetMinimum(0.7253868);
   Graph3->SetMaximum(1.024965);
   Graph3->SetDirectory(0);
   Graph3->SetStats(0);
   Graph3->SetLineStyle(0);
   Graph3->SetLineWidth(2);
   Graph3->SetMarkerStyle(20);
   Graph3->SetMarkerSize(1.3);
   Graph3->GetXaxis()->SetLabelFont(42);
   Graph3->GetXaxis()->SetLabelOffset(0.01);
   Graph3->GetXaxis()->SetTitleSize(0.045);
   Graph3->GetXaxis()->SetTitleOffset(1.2);
   Graph3->GetXaxis()->SetTitleFont(42);
   Graph3->GetYaxis()->SetLabelFont(42);
   Graph3->GetYaxis()->SetLabelOffset(0.01);
   Graph3->GetYaxis()->SetTitleSize(0.045);
   Graph3->GetYaxis()->SetTitleOffset(1.8);
   Graph3->GetYaxis()->SetTitleFont(42);
   Graph3->GetZaxis()->SetLabelFont(42);
   Graph3->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph3);
   
   grae->Draw("p");
   
   grae = new TGraphAsymmErrors(16);
   grae->SetName("");
   grae->SetTitle("");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(25);
   grae->SetMarkerSize(1.3);
   grae->SetPoint(0,19.5,0.3333333);
   grae->SetPointError(0,0,0,0.118443,0.1346547);
   grae->SetPoint(1,22.5,0.4166667);
   grae->SetPointError(1,0,0,0.07857571,0.08151994);
   grae->SetPoint(2,26,0.622807);
   grae->SetPointError(2,0,0,0.04574949,0.0443332);
   grae->SetPoint(3,30,0.7102804);
   grae->SetPointError(3,0,0,0.03155007,0.03024836);
   grae->SetPoint(4,34.5,0.7713098);
   grae->SetPointError(4,0,0,0.01950511,0.01875447);
   grae->SetPoint(5,40,0.8709016);
   grae->SetPointError(5,0,0,0.01098825,0.01048144);
   grae->SetPoint(6,46.5,0.9289617);
   grae->SetPointError(6,0,0,0.005872228,0.005587743);
   grae->SetPoint(7,53,0.96);
   grae->SetPointError(7,0,0,0.004866839,0.004515586);
   grae->SetPoint(8,60,0.9819234);
   grae->SetPointError(8,0,0,0.003830992,0.003363644);
   grae->SetPoint(9,69,0.9884615);
   grae->SetPointError(9,0,0,0.004291396,0.003443304);
   grae->SetPoint(10,79,0.9823678);
   grae->SetPointError(10,0,0,0.007517494,0.005868544);
   grae->SetPoint(11,90.5,1);
   grae->SetPointError(11,0,0,0.00520845,0);
   grae->SetPoint(12,105.5,0.9910714);
   grae->SetPointError(12,0,0,0.01311826,0.006496046);
   grae->SetPoint(13,123.5,1);
   grae->SetPointError(13,0,0,0.01835927,0);
   grae->SetPoint(14,153.5,1);
   grae->SetPointError(14,0,0,0.03229151,0);
   grae->SetPoint(15,197,1);
   grae->SetPointError(15,0,0,0.09917226,0);
   
   TH1 *Graph4 = new TH1F("Graph4","",100,1.75,214.75);
   Graph4->SetMinimum(0.1363794);
   Graph4->SetMaximum(1.078511);
   Graph4->SetDirectory(0);
   Graph4->SetStats(0);
   Graph4->SetLineStyle(0);
   Graph4->SetLineWidth(2);
   Graph4->SetMarkerStyle(20);
   Graph4->SetMarkerSize(1.3);
   Graph4->GetXaxis()->SetLabelFont(42);
   Graph4->GetXaxis()->SetLabelOffset(0.01);
   Graph4->GetXaxis()->SetTitleSize(0.045);
   Graph4->GetXaxis()->SetTitleOffset(1.2);
   Graph4->GetXaxis()->SetTitleFont(42);
   Graph4->GetYaxis()->SetLabelFont(42);
   Graph4->GetYaxis()->SetLabelOffset(0.01);
   Graph4->GetYaxis()->SetTitleSize(0.045);
   Graph4->GetYaxis()->SetTitleOffset(1.8);
   Graph4->GetYaxis()->SetTitleFont(42);
   Graph4->GetZaxis()->SetLabelFont(42);
   Graph4->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph4);
   
   grae->Draw("p");
   line = new TLine(18,1,300,1);
   line->SetLineStyle(2);
   line->Draw();
   
   TLegend *leg = new TLegend(0.31,0.8,0.72,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","PYTHIA+DATA","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("","Leading Jet","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.3);
   entry=leg->AddEntry("","SubLeading Jet","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.3);
   leg->Draw();
      tex = new TLatex(0.43,0.51,"10-30%");
tex->SetNDC();
   tex->SetTextFont(63);
   tex->SetTextSize(18);
   tex->Draw();
   p_1_0->Modified();
   c1->cd();
  
// ------------>Primitives in pad: p_2_0
   p_2_0 = new TPad("p_2_0", "p_2_0",0.6910157,0,1,0.9831461);
   p_2_0->Draw();
   p_2_0->cd();
   p_2_0->Range(18,-0.2731707,224.1224,1.434146);
   p_2_0->SetFillColor(0);
   p_2_0->SetBorderMode(0);
   p_2_0->SetBorderSize(0);
   p_2_0->SetTickx(1);
   p_2_0->SetTicky(1);
   p_2_0->SetLeftMargin(0);
   p_2_0->SetRightMargin(0.02);
   p_2_0->SetTopMargin(0.02);
   p_2_0->SetBottomMargin(0.16);
   p_2_0->SetFrameLineColor(0);
   p_2_0->SetFrameBorderMode(0);
   p_2_0->SetFrameLineColor(0);
   p_2_0->SetFrameBorderMode(0);
   Double_t xAxis3[18] = {18, 21, 24, 28, 32, 37, 43, 50, 56, 64, 74, 84, 97, 114, 133, 174, 220, 300}; 
   
   TH1 *hTmp = new TH1F("hTmp","",17, xAxis3);
   hTmp->SetMinimum(0);
   hTmp->SetMaximum(1.4);
   hTmp->SetStats(0);
   hTmp->SetLineStyle(0);
   hTmp->SetLineWidth(2);
   hTmp->SetMarkerStyle(20);
   hTmp->SetMarkerSize(1.3);
   hTmp->GetXaxis()->SetTitle("p_{T}^{GenJet} (GeV/c)");
   hTmp->GetXaxis()->SetRange(1,16);
   hTmp->GetXaxis()->CenterTitle(true);
   hTmp->GetXaxis()->SetNdivisions(505);
   hTmp->GetXaxis()->SetLabelFont(63);
   hTmp->GetXaxis()->SetLabelOffset(0.01);
   hTmp->GetXaxis()->SetLabelSize(24);
   hTmp->GetXaxis()->SetTitleSize(28);
   hTmp->GetXaxis()->SetTitleOffset(1.2);
   hTmp->GetXaxis()->SetTitleFont(63);
   hTmp->GetYaxis()->SetTitle("Eff. (Matched Reco/Gen)");
   hTmp->GetYaxis()->CenterTitle(true);
   hTmp->GetYaxis()->SetLabelFont(63);
   hTmp->GetYaxis()->SetLabelOffset(0.01);
   hTmp->GetYaxis()->SetLabelSize(24);
   hTmp->GetYaxis()->SetTitleSize(24);
   hTmp->GetYaxis()->SetTitleOffset(1.8);
   hTmp->GetYaxis()->SetTitleFont(63);
   hTmp->GetZaxis()->SetLabelFont(42);
   hTmp->GetZaxis()->SetTitleFont(42);
   hTmp->Draw("");
   
   grae = new TGraphAsymmErrors(10);
   grae->SetName("");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.3);
   grae->SetPoint(0,53,0.956245);
   grae->SetPointError(0,0,0,0.006021209,0.005536137);
   grae->SetPoint(1,60,0.9742619);
   grae->SetPointError(1,0,0,0.00460988,0.004129203);
   grae->SetPoint(2,69,0.9875566);
   grae->SetPointError(2,0,0,0.004135499,0.003390637);
   grae->SetPoint(3,79,0.9896907);
   grae->SetPointError(3,0,0,0.005371879,0.003988193);
   grae->SetPoint(4,90.5,0.9967532);
   grae->SetPointError(4,0,0,0.004837377,0.002371439);
   grae->SetPoint(5,105.5,0.9934211);
   grae->SetPointError(5,0,0,0.00972201,0.004794304);
   grae->SetPoint(6,123.5,1);
   grae->SetPointError(6,0,0,0.01462095,0);
   grae->SetPoint(7,153.5,1);
   grae->SetPointError(7,0,0,0.0283128,0);
   grae->SetPoint(8,197,1);
   grae->SetPointError(8,0,0,0.09129775,0);
   grae->SetPoint(9,260,1);
   grae->SetPointError(9,0,0,0.3181538,0);
   
   TH1 *Graph5 = new TH1F("Graph5","",100,32.3,280.7);
   Graph5->SetMinimum(0.6500308);
   Graph5->SetMaximum(1.031815);
   Graph5->SetDirectory(0);
   Graph5->SetStats(0);
   Graph5->SetLineStyle(0);
   Graph5->SetLineWidth(2);
   Graph5->SetMarkerStyle(20);
   Graph5->SetMarkerSize(1.3);
   Graph5->GetXaxis()->SetLabelFont(42);
   Graph5->GetXaxis()->SetLabelOffset(0.01);
   Graph5->GetXaxis()->SetTitleSize(0.045);
   Graph5->GetXaxis()->SetTitleOffset(1.2);
   Graph5->GetXaxis()->SetTitleFont(42);
   Graph5->GetYaxis()->SetLabelFont(42);
   Graph5->GetYaxis()->SetLabelOffset(0.01);
   Graph5->GetYaxis()->SetTitleSize(0.045);
   Graph5->GetYaxis()->SetTitleOffset(1.8);
   Graph5->GetYaxis()->SetTitleFont(42);
   Graph5->GetZaxis()->SetLabelFont(42);
   Graph5->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph5);
   
   grae->Draw("p");
   
   grae = new TGraphAsymmErrors(16);
   grae->SetName("");
   grae->SetTitle("");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(25);
   grae->SetMarkerSize(1.3);
   grae->SetPoint(0,19.5,0.1666667);
   grae->SetPointError(0,0,0,0.1081155,0.1696671);
   grae->SetPoint(1,22.5,0.3333333);
   grae->SetPointError(1,0,0,0.1081144,0.1214185);
   grae->SetPoint(2,26,0.6);
   grae->SetPointError(2,0,0,0.06335432,0.06119324);
   grae->SetPoint(3,30,0.6615385);
   grae->SetPointError(3,0,0,0.04204659,0.04040942);
   grae->SetPoint(4,34.5,0.7949791);
   grae->SetPointError(4,0,0,0.02686696,0.02522927);
   grae->SetPoint(5,40,0.8675373);
   grae->SetPointError(5,0,0,0.01509547,0.01418184);
   grae->SetPoint(6,46.5,0.9187063);
   grae->SetPointError(6,0,0,0.008329598,0.007841113);
   grae->SetPoint(7,53,0.9383562);
   grae->SetPointError(7,0,0,0.008469778,0.007801511);
   grae->SetPoint(8,60,0.9668175);
   grae->SetPointError(8,0,0,0.007453611,0.006510166);
   grae->SetPoint(9,69,0.983945);
   grae->SetPointError(9,0,0,0.006854615,0.005347492);
   grae->SetPoint(10,79,0.984456);
   grae->SetPointError(10,0,0,0.01092583,0.007440347);
   grae->SetPoint(11,90.5,1);
   grae->SetPointError(11,0,0,0.01020519,0);
   grae->SetPoint(12,105.5,0.984375);
   grae->SetPointError(12,0,0,0.02258393,0.01131564);
   grae->SetPoint(13,123.5,1);
   grae->SetPointError(13,0,0,0.03638148,0);
   grae->SetPoint(14,153.5,0.96);
   grae->SetPointError(14,0,0,0.05450613,0.02847887);
   grae->SetPoint(15,197,1);
   grae->SetPointError(15,0,0,0.2496484,0);
   
   TH1 *Graph6 = new TH1F("Graph6","",100,1.75,214.75);
   Graph6->SetMinimum(0);
   Graph6->SetMaximum(1.094145);
   Graph6->SetDirectory(0);
   Graph6->SetStats(0);
   Graph6->SetLineStyle(0);
   Graph6->SetLineWidth(2);
   Graph6->SetMarkerStyle(20);
   Graph6->SetMarkerSize(1.3);
   Graph6->GetXaxis()->SetLabelFont(42);
   Graph6->GetXaxis()->SetLabelOffset(0.01);
   Graph6->GetXaxis()->SetTitleSize(0.045);
   Graph6->GetXaxis()->SetTitleOffset(1.2);
   Graph6->GetXaxis()->SetTitleFont(42);
   Graph6->GetYaxis()->SetLabelFont(42);
   Graph6->GetYaxis()->SetLabelOffset(0.01);
   Graph6->GetYaxis()->SetTitleSize(0.045);
   Graph6->GetYaxis()->SetTitleOffset(1.8);
   Graph6->GetYaxis()->SetTitleFont(42);
   Graph6->GetZaxis()->SetLabelFont(42);
   Graph6->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph6);
   
   grae->Draw("p");
   line = new TLine(18,1,300,1);
   line->SetLineStyle(2);
   line->Draw();
      tex = new TLatex(0.43,0.51,"0-10%");
tex->SetNDC();
   tex->SetTextFont(63);
   tex->SetTextSize(18);
   tex->Draw();
   p_2_0->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
