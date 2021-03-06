#include <string.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <stdbool.h>
#include <iostream>
#include <TColor.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TArrow.h>
#include <TMath.h>
#include <TFile.h> 
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <fstream>
#include "Header.h"

#if defined(__GLIBCXX__)
#pragma link C++ class MyOtherClass;
#endif


using namespace std;

void contours()
{
// we import the file (ppbarContours.root) with the contours for p+pbar produced in Pb-Pb collisions
	TFile *ppbarFile = new TFile("ppbarContours.root");
// we import the file (phiContours.root) with the contours for phi mesons produced in Pb-Pb collisions
	TFile *phiFile = new TFile("phiContours.root");

// we define a canvas
	TCanvas *T1 = new TCanvas("T1","Contours",200,200);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(0);
// il manque des options gstyle ou gpad

// we define a TGraph for each contour 
	TGraph *ppbar_1sigma[11];
	TGraph *ppbar_2sigma[11];
	TGraph *phi_1sigma[9];
	TGraph *phi_2sigma[9];

// we define the multigraph that will contain all the contours
	TMultiGraph *mg = new TMultiGraph();
	mg->GetXaxis()->SetTitle(" #beta_{S} ");
	mg->GetYaxis()->SetTitle("T_{kin} (GeV) ");

//------------------------------------------------------
//----------------- p+pbar -----------------------------	
//------------------------------------------------------
// we explorate the ppbar file in order to get the 9 first pairs of contours it contains
	char ppbar_1sigma_name[] = "ppbar_1sigma_contour_X";
	char ppbar_2sigma_name[] = "ppbar_2sigma_contour_X";
	for(int digit=1;digit<=9;digit++)
	{
	// we do it for the 1sigma contours
		ppbar_1sigma_name[21] = digit+'0';
		ppbar_1sigma[digit] = (TGraph*)ppbarFile->Get(Form("ppbar_1sigma_contour_%d",digit)) ; 
		ppbar_1sigma[digit]->Draw();   
		mg->Add(ppbar_1sigma[digit]);
	// we do it for the 2sigma contours
		ppbar_2sigma_name[21] = digit+'0';
		ppbar_2sigma[digit] = (TGraph*)ppbarFile->Get(Form("ppbar_2sigma_contour_%d",digit)) ; 
		ppbar_2sigma[digit]->Draw();   
		mg->Add(ppbar_2sigma[digit]);
	}
// we add the 10th pair of contours	
	ppbar_1sigma[10] = (TGraph*)ppbarFile->Get("ppbar_1sigma_contour_10") ; 
	ppbar_1sigma[10]->Draw();  
	mg->Add(ppbar_1sigma[10]);	
	ppbar_2sigma[10] = (TGraph*)ppbarFile->Get("ppbar_2sigma_contour_10") ; 
	ppbar_1sigma[10]->Draw();   
	mg->Add(ppbar_2sigma[10]);
	
	

//------------------------------------------------------
//------------------------------------------------------	
// we draw the multigraph
	
	mg->Draw("A");	  
	TLine *line = new TLine(0.76,0.175,0.92,0.175);
	line->SetLineColor(30);
	line->SetLineWidth(3);
	line->SetLineStyle(2);
	line->Draw("SAME");
// change the axis limits
	gPad->Modified();   
	mg->GetXaxis()->SetLimits(0.76,0.92); // we limit the x axis range
	// then we limit the y axis range
	mg->SetMinimum(0.08);
	mg->SetMaximum(0.18);
	
// we add a legend
	TLatex Tl;
	Tl.SetNDC();
	Tl.SetTextSize(0.040);
	Tl.SetTextAlign(13);  //align at top
	Tl.SetTextFont(42);
// we indicate the centrality
	Tl.DrawLatex(0.145,0.84, "80-90%");  // top left of the plot
	Tl.DrawLatex(0.85,0.15, "0-5%");    // bottom right of the plot
// we indicate the color code
	Tl.DrawLatex(0.5,0.85,"#color[30]{QGP Transition Temperature}");
	Tl.DrawLatex(0.15,0.2, "Blast-wave fit to p+#bar{p} (0.3-3 GeV/c)");
	Tl.DrawLatex(0.15,0.3, "#color[4]{1#sigma contour}");   // 4 is blue
	Tl.DrawLatex(0.15,0.25, "#color[2]{2#sigma contour}");   // 2 is blue
	
// we close the files
	ppbarFile->Close();
	phiFile->Close();
	
}
