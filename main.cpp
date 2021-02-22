
#include "Header.h"
#include "function.cpp"

#if defined(__GLIBCXX__)
#pragma link C++ class MyOtherClass;
#endif


using namespace std;


int main(){
	TString name ;
	char particule[]="  ", collision[]="  " ;
	double masse ;
	char *file ;
	char *table ;
	int number ;
	double SCALE[11] ;
	ofstream myfile ;
	ofstream myfile2 ;
	cout << " quelle particule étudier vous ? proton, mesonphi " << endl;
	cin >> particule ;
	cout << " quelle type de collision étudier vous ? pp, PbPb " << endl;
	cin >> collision ;
	if(particule[1]=='r')
	{
		masse=0.938 ;
		file="HEPData-1569102768-v1-root.root" ;
		if(collision[1]=='p'){
			table="Table 6";
			number = 1 ;
			myfile.open("resultat_proton_pp.txt");
	myfile << "Centrality class" << "	" << " DN/DY " << "		" << "Uncertainty" <<"	"<<"Fraction of yield at low Pt"<< "    " <<"Uncertainty" <<"	"<<"chi_2"<<"		"<<"ndf"<<"		"<<"chi2/ndf"<< endl;
	myfile << "  " << endl;
		}
		if(collision[1]=='b'){
			myfile.open("resultat_proton_PbPb.txt");
	myfile << "Centrality class" << "	" << " DN/DY " << "		" << "Uncertainty" <<"	"<<"Fraction of yield at low Pt"<< "    " <<"Uncertainty" <<"	"<<"chi_2"<<"		"<<"ndf"<<"		"<<"chi2/ndf"<< endl;
	myfile << "  " << endl;
	//myfile.close();
	myfile2.open("resultat_proton_PbPb_BL.txt");
	myfile2 << "Centrality class" << "	" << " DN/D_(eta) " << "	" << "<Beta_T>" <<"	"<<"<T_kin GeV"<< "    " <<"n" << endl;
	myfile2 << "  " << endl;
			table="Table 5";
			number = 10 ;
		}
		SCALE[1]=101276 ;
		SCALE[2]=73945 ;
		SCALE[3]=39170 ;	
		SCALE[4]=19521 ;
		SCALE[5]=8718 ;
		SCALE[6]=2840 ;
		SCALE[7]=945 ;
		SCALE[8]=182 ;
		SCALE[9]=30 ;
		SCALE[10]=6 ;
		name = "ppbar" ;
	}
	if(particule[0]=='m')
	{
		masse=1.019461 ;
		file="HEPData-ins1762368-v1-root.root" ;
		if(collision[1]=='p'){
			table="Table 4";
			number = 1 ;
			myfile.open("resultat_meson_pp.txt");
			myfile << "Centrality class" << "	" << " DN/DY " << "		" << "Uncertainty" <<"	"<<"Fraction of yield at low Pt"<< "    " <<"Uncertainty" <<"	"<<"chi_2"<<"		"<<"ndf"<<"		"<<"chi2/ndf"<< endl;
			myfile << "  " << endl;
		}
		if(collision[1]=='b'){
			myfile.open("resultat_meson_PbPb.txt");
	myfile << "Centrality class" << "	" << " DN/DY " << "		" << "Uncertainty" <<"	"<<"Fraction of yield at low Pt"<< "    " <<"Uncertainty" <<"	"<<"chi_2"<<"		"<<"ndf"<<"		"<<"chi2/ndf"<< endl;
			myfile << "  " << endl;
			myfile2.open("resultat_meson_PbPb_BL.txt");
			myfile2 << "Centrality class" << "	" << " DN/D_(eta) " << "	" << "<Beta_T>" <<"	"<<"<T_kin GeV"<< "    " <<"n" << endl;
			myfile2 << "  " << endl;
			table="Table 3";
			number = 8 ;
		}
		SCALE[1]=39827 ;
		SCALE[2]=24080 ;
		SCALE[3]=11821 ;	
		SCALE[4]=5520 ;
		SCALE[5]=2224 ;
		SCALE[6]=607 ;
		SCALE[7]=164 ;
		SCALE[8]=23.5 ;
		name = "phi" ;
	}


	char data[] = "Hist1D_yX" ;
	char error_sys[] = "Hist1D_yX_e2" ;
	char error_stat[] = "Hist1D_yX_e1" ;
	int digit, Nbinx ;
	double a , b ;
	double Scale_Histo[11], SCALE2[11], SCALE3[11], NSCALE[11], SCALE_BOL[11] , T_value[11], Beta_T_value[11],A_value[11],chi2_expo, ndf_expo, chindf_expo, ndf_levy, chindf_levy, Beta_S_value[11],Integral_value,Integral_error;
	
	

	
	SCALE2[1]=1943 ;
	SCALE2[2]=1587 ;
	SCALE2[3]=1180 ;
	SCALE2[4]=786 ;
	SCALE2[5]=512 ;
	SCALE2[6]=318 ;
	SCALE2[7]=183 ;
	SCALE2[8]=96 ;
	SCALE2[9]=44 ;
	SCALE2[10]=17 ;
	
	SCALE3[1]=74.56 ;
	SCALE3[2]=61.51 ;
	SCALE3[3]=47.40 ;
	SCALE3[4]=33.17 ;
	SCALE3[5]=22.51 ;
	SCALE3[6]=14.46 ;
	SCALE3[7]=8.71 ;
	SCALE3[8]=4.74 ;
	SCALE3[9]=2.30 ;
	SCALE3[10]=0.92 ;
	
	NSCALE[1]=0.735 ;
	NSCALE[2]=0.736;
	NSCALE[3]=0.739;
	NSCALE[4]=0.771;
	NSCALE[5]=0.828;
	NSCALE[6]=0.908;
	NSCALE[7]=1.052;
	NSCALE[8]=1.262;
	NSCALE[9]=1.678;
	NSCALE[10]=2.423 ;
	
	T_value[1]=0.090;
	T_value[2]=0.091;
	T_value[3]=0.094;
	T_value[4]=0.097;
	T_value[5]=0.101;
	T_value[6]=0.108;
	T_value[7]=0.115;
	T_value[8]=0.129;
	T_value[9]=0.147;
	T_value[10]=0.161;
	
	Beta_T_value[1]=0.663;
	Beta_T_value[2]=0.660;
	Beta_T_value[3]=0.655;
	Beta_T_value[4]=0.643;
	Beta_T_value[5]=0.622;
	Beta_T_value[6]=0.595;
	Beta_T_value[7]=0.557;
	Beta_T_value[8]=0.506;
	Beta_T_value[9]=0.435;
	Beta_T_value[10]=0.355;
	
	Beta_S_value[1]=0.906;
	Beta_S_value[2]=0.902;
	Beta_S_value[3]=0.897;
	Beta_S_value[4]=0.890;
	Beta_S_value[5]=0.879;
	Beta_S_value[6]=0.865;
	Beta_S_value[7]=0.849;
	Beta_S_value[8]=0.825;
	Beta_S_value[9]=0.799;
	Beta_S_value[10]=0.740;
	
	
	

	
	TFile *myFile = new TFile(file);
	TDirectoryFile* dirFile = (TDirectoryFile*)myFile->Get(table);
	hist=(TH1F*)dirFile->Get("Hist1D_y1");
	Nbinx = hist->GetNbinsX();
	
	double X[Nbinx];  // x axis of the plot
	double Y[Nbinx];  // y axis of the plot
	double Xerr[Nbinx];  // we take into account the width of a bin (this is not an error strictly speaking!)
	double Yerr[Nbinx];  // error on y axis
	double par[5],err[5], parlevy[4],errlevy[4];

	TCanvas *test1 = new TCanvas("DATA","DATA",0,0,700,500);
		test1->SetGrid();
		test1->SetLogy();
	TCanvas *ROMAIN = new TCanvas("FIT & DATA" , " FIT & DATA " , 200 , 200 );
		ROMAIN->SetGrid();
		ROMAIN->SetLogy();
	TCanvas *ROMAIN2 = new TCanvas("Mutli DATA" , " Multi DATA " , 200 , 200 );
		ROMAIN->SetGrid();
	TCanvas *DATA_MODEL = new TCanvas("DATA_MODEL","DATA_MODEL",200,10,600,400);
		DATA_MODEL->SetGrid();
		gStyle->SetPadBorderMode(0);
		gStyle->SetFrameBorderMode(0);
		gStyle->SetOptStat(0);
		float small = 1e-5;
		DATA_MODEL->Divide(1,2,small,small);
	TCanvas *CONTOUR = new TCanvas("CONTOUR" , " Contour " , 200 , 200 );
	TMultiGraph  *mg  = new TMultiGraph();
		mg->SetTitle(" ;<Beta_S>; T");  // we set the title of the graph, and the one of the axis
	TMultiGraph  *DM = new TMultiGraph();
		DM->SetTitle(" ;Pt (Gev/c) ; Data / BW-FIT ");

	  
	

//--------------------------------------------------Coliisions Pb-Pb--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	if(collision[1]=='b'){
	
	TF1 *func2 = new TF1("boltzmann",boltzmann,0.3,upper_limit,3);
		func2->SetParameter(0,Scale_Histo[1]);
		func2->SetParameter(1,masse);
		func2->SetParLimits(1,masse,masse);
		func2->SetParameter(2,0.2);
		func2->SetParNames("C", "m", "T");	
		func2->SetLineColor(kGreen);
	
	TF1 *func1 = new TF1("expo_law",expo_law,0.3,upper_limit,3);
		func1->SetParameter(0,Scale_Histo[1]);
		func1->SetParameter(1,masse);
		func1->SetParLimits(1,masse,masse);
		func1->SetParameter(2,0.1);
		func1->SetParNames("dyN","m","T");
		func1->SetLineColor(kViolet);
	
	TString contours_sigma_output=name+"Contours.root";
	TFile *Contours_sigma_output = new TFile(contours_sigma_output, "RECREATE");
	
	TString outputfilename=name+"_PbPb_Data_Model.root" ;
	TFile* OutputHisto = new TFile(outputfilename, "RECREATE");
	
	
	for(int digit=1 ; digit<= number ; digit++){
		data[8]=digit+'0' ;
		error_sys[8]=digit+'0' ;
		error_stat[8]=digit+'0' ;	
		hist=(TH1F*)dirFile->Get(data);
		TH1F* H1_pp_e1=(TH1F*)dirFile->Get(error_sys);
		TH1F* H1_pp_E1=(TH1F*)dirFile->Get(error_stat);
		if(digit==10){
			hist=(TH1F*)dirFile->Get("Hist1D_y10");
			H1_pp_e1=(TH1F*)dirFile->Get("Hist1D_y10_e2");
			H1_pp_E1=(TH1F*)dirFile->Get("Hist1D_y10_e1");
			//
			}
		hist->SetNameTitle(" Table 5" , "pT-distributions " );
		hist->SetXTitle("p_{T} [GeV/c]");
		hist->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
		Nbinx = hist->GetNbinsX();
		a= hist->GetBinLowEdge(1);
		b= hist->GetBinLowEdge(Nbinx);
		b= b+ hist->GetBinWidth(Nbinx);
		cout << " a et b = " << a << "    " << b  <<  endl;
		for(int i = 0; i <= Nbinx  ; i++){
			hist->SetBinError(i, sqrt( pow(H1_pp_e1->GetBinContent(i),2) + pow(H1_pp_E1->GetBinContent(i),2) ) ) ;
		}	
		Scale_Histo[1]=hist->Integral("width") ;
		auto fitResult_2 = hist->Fit(func2,"S R I 0");
		auto covMatrix_2 = fitResult_2->GetCovarianceMatrix();
		double chi2_bol = func2->GetChisquare();
		double ndf_bol = func2->GetNDF();
		double chindf_bol = chi2_bol/ndf_bol ;
		cout <<" Chi square boltz = "<< chi2_bol << endl;
		cout <<" NDF boltz = "<< ndf_bol << endl; 
		cout <<" Chi2 / ndf = " << chindf_bol << endl;
		cout <<" integral value " << func2->Integral(0.3,3);
	
		//hist->Fit("expo_law","+ I  ","SAMES HIST",0.3,upper_limit);  

		auto fitResult = hist->Fit(func1,"S R I 0");
		auto covMatrix = fitResult->GetCovarianceMatrix();
		double chi2_expo = func1->GetChisquare();
		double ndf_expo = func1->GetNDF();
		double chindf_expo = chi2_expo/ndf_expo ;
		cout <<" Chi square expo law = "<< chi2_expo << endl;
		cout <<" NDF expo law = "<< ndf_expo << endl; 
		cout <<" Chi2 / ndf = " << chindf_expo << endl;
		cout <<" integral value " << func1->Integral(0,20);
		cout <<" integral error value " << func1->IntegralError(0,20, fitResult->GetParams() , covMatrix.GetMatrixArray());
	
// FIT BLAST WAVE	
		TMinuit *gMinuit_BL = new TMinuit(5);
			gMinuit_BL->SetFCN(fcn_blast_wave);
			
			gMinuit_BL->DefineParameter(0, "A", SCALE[digit], 1, 0, 200000);
			gMinuit_BL->FixParameter(0);
			
			gMinuit_BL->DefineParameter(1, "m_0",masse, 0.01, 0.1,2);
			gMinuit_BL->FixParameter(1);
			
			gMinuit_BL->DefineParameter(2, "T", 0.1, 0.0001, 0, 0.3);
			//gMinuit_BL->FixParameter(2);
			
			gMinuit_BL->DefineParameter(3, "n", NSCALE[digit], 0.01, 0, 20);
			//gMinuit_BL->FixParameter(3);
				
			gMinuit_BL->DefineParameter(4, "beta_s", 0.6 , 0.0001, 0, 1);
			//gMinuit_BL->FixParameter(4);
		gMinuit_BL->Command("MIGRAD");
		gMinuit_BL->Command("MINOS");
		
		ROMAIN->cd();
		for(int i=0;i<5;i++) gMinuit_BL->GetParameter(i,par[i],err[i]);
		TH1F* curve = new TH1F("curve","curve",UpperBin_Fit(upper_limit),0,upper_limit);
			for(int i=1;i<=UpperBin_Fit(upper_limit);i++) 
			{		
				double *x = new double ;
				x[0] = curve->GetBinCenter(i);
				double f = blast_wave(x,par);
				curve->SetBinContent(i,f);
			} 
		
			
		A_value[digit]=par[0];
		curve->SetLineWidth(3);
		hist->Draw("same");
		curve->Draw("csame");	
		
// we get the covariant matrix from the Tminuit instance
		int npars = gMinuit_BL->GetNumPars();
		double *covar = new double[npars*npars];
		gMinuit_BL->mnemat(covar,npars);	
		
	   	gMinuit_BL->SetErrorDef(4); //note 4 and not 2!
	   	TGraph *gr2 = (TGraph*)gMinuit_BL->Contour(25,4,2);
	   	gr2->SetFillColor(42);
	   	gr2->SetLineColor(kRed);
	   	gr2->SetLineWidth(2);
//Get contour for parameter 2 versus parameter 4 for ERRDEF=1
	   	gMinuit_BL->SetErrorDef(1);
	   	TGraph *gr1 = (TGraph*)gMinuit_BL->Contour(25,4,2);
	   	gr1->SetFillColor(38);
	   	gr1->SetLineColor(kBlue);
	   	gr1->SetLineWidth(2);


		mg->Add(gr2);
		mg->Add(gr1);
	    
		gr2->SetName( Form(name+"_1sigma_contour_%d", digit) );
		gr1->SetName( Form(name+"_2sigma_contour_%d", digit) );
		Contours_sigma_output->cd();  // we save the 2-sigma contours 
		gr2->Write();
		gr1->Write();
		Data_Model(X,Y,Xerr,Yerr,blast_wave,par,par_size_blast_wave,covar);
		//Integral_value = Integral(levy,parlevy,par_size_levy,covar_levy,0,20) ;
		//Integral_error = Integral_E(levy,parlevy,par_size_levy,covar_levy,0,20) ;
		
		Integral_value = Integral(blast_wave,par,par_size_blast_wave,covar,0,20) ;
		Integral_error = Integral_E(blast_wave,par,par_size_blast_wave,covar,0,20) ;
		
		cout << " test integral " << Integral_value << endl;
// we create the graph to plot Data/Model fit with rectangle errors for phi mesons
		OutputHisto->cd();
		DATA_MODEL->cd(1);
		gPad->SetTickx(2);
		gPad->SetBottomMargin(small); 
		gPad->SetTickx();
//gPad->SetTicky(1);   	
// we create the graph to plot Data/Model fit with rectangle errors for p+pbar
		DATA_MODEL->cd(2); 
		gPad->SetTopMargin(small);
// we create the graph to plot Data/Model fit with rectangle errors
	   	TGraphErrors *dataModel = new TGraphErrors(UpperBin_Fit(upper_limit),X,Y,Xerr,Yerr);
		dataModel->SetName( Form(name+"_data_model_%d",digit) );
		dataModel->SetLineColor(digit);
		dataModel->SetLineWidth(2);
		dataModel->SetMarkerColor(digit);
		dataModel->SetMarkerSize(1);
		dataModel->SetMarkerStyle(21);
		dataModel->SetTitle(" ; p_T (GeV/c); Data / model fit");
		dataModel->SetLineWidth(2);
		dataModel->GetXaxis()->SetLimits(0.2,upper_limit+hist->GetBinWidth(UpperBin_Fit(upper_limit)));
		dataModel->Draw("A5 SAME");   // A has to be there (if not, nothing appears on the canvas) & 5 is to plot error rectangles
		dataModel->Draw("PX SAME");   // P is to put the chosen marker & X to remove the error bars (leave only the rectangles)
		dataModel->Write();
	// we plot a line at y=1 to show the deviation to Data = Model fit

		//We Store the graph in a Multigraph
		DM->Add(dataModel);
		//dataModel->Draw("A5");   // A has to be there (if not, nothing appears on the canvas) & 5 is to plot error rectangles
		//dataModel->Draw("PX");   // P is to put the chosen marker & X to remove the error bars (leave only the rectangles)
		// Export of the fit result into result.txt
		myfile <<"Centrality fichier "<< digit<<"	" << Integral_value << "		" << Integral_error << "	" << 100*Integral(blast_wave,par,par_size_blast_wave,covar,0,a)/Integral(blast_wave,par,par_size_blast_wave,covar,0,b) << "			" <<(100*Integral(blast_wave,par,par_size_blast_wave,covar,0,a)/Integral(blast_wave,par,par_size_blast_wave,covar,0,b)*TMath::Sqrt( TMath::Power( Integral_E(blast_wave,par,par_size_blast_wave,covar,0,a) / Integral(blast_wave,par,par_size_blast_wave,covar,0,a),2)  +  TMath::Power( Integral_E(blast_wave,par,par_size_blast_wave,covar,0,b) / Integral(blast_wave,par,par_size_blast_wave,covar,0,b) ,2) ) ) <<"	"<< chi2_blast_wave(par)<< "		"<<(Nbinx-gMinuit_BL->GetNumFreePars())<<"		"<<(chi2_blast_wave(par)/(Nbinx-gMinuit_BL->GetNumFreePars()))<<endl;
	// Export of the result of the fit into file result_BL.txt
		myfile2 << " Centrality fichier :"<<digit<<"	 "<<par[0]<<"		"<<(par[4]*(2/(2+par[3])))<<"	"<<par[2]<<"	"<<par[3]<<"	"<<endl;
		Contours_sigma_output -> cd();  // we save the 2-sigma contours 
		CONTOUR->cd();
		mg->Draw("AL");
		CONTOUR->Write();
		OutputHisto->cd();
		ROMAIN2->cd();
		DM->Draw("A5");
		DM->Draw("PX");
		TLine *line = new TLine(a,1,upper_limit+hist->GetBinWidth(UpperBin_Fit(upper_limit)),1);
		line->SetLineColor(30);
		line->SetLineWidth(3);
		line->SetLineStyle(2);
		line->Draw("SAME");
		ROMAIN2->Write();
		ROMAIN2->Close();
		//DATA_MODEL->Write();
 		test1->Close();
 		DATA_MODEL->Close();
		ROMAIN->Write();
		}
	}
//------------------------------------- Collision pp----------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
	if(collision[1]=='p'){
		TString outputfilename=name+"_pp_Data_Model.root" ;
		TFile* OutputHisto = new TFile(outputfilename, "RECREATE");
		for(int digit=1 ; digit<= number ; digit++)
		{
			data[8]=digit+'0' ;
			error_sys[8]=digit+'0' ;
			error_stat[8]=digit+'0' ;	
			hist=(TH1F*)dirFile->Get(data);
			TH1F* H1_pp_e1=(TH1F*)dirFile->Get(error_sys);
			TH1F* H1_pp_E1=(TH1F*)dirFile->Get(error_stat);
			if(digit==10)
			{
				hist=(TH1F*)dirFile->Get("Hist1D_y10");
				H1_pp_e1=(TH1F*)dirFile->Get("Hist1D_y10_e2");
				H1_pp_E1=(TH1F*)dirFile->Get("Hist1D_y10_e1");
				//
			}
			hist->SetNameTitle(" Table 5" , "" );
			hist->SetXTitle("p_{T} (GeV/c)");
			hist->SetYTitle("(1/N_{ev})*d^{2}(N)/dp_{T}dy  (Gev/c) ");
			Nbinx = hist->GetNbinsX();
			a= hist->GetBinLowEdge(1);
			b= hist->GetBinLowEdge(Nbinx);
			b= b+ hist->GetBinWidth(Nbinx);
			cout << " a et b = " << a << "    " << b  <<  endl;
			for(int i = 0; i <= Nbinx  ; i++){
				hist->SetBinError(i, sqrt( pow(H1_pp_e1->GetBinContent(i),2) + pow(H1_pp_E1->GetBinContent(i),2) ) ) ;
			}	
			Scale_Histo[1]=hist->Integral("width") ;	
		TF1 *func2 = new TF1("boltzmann",boltzmann,0.3,upper_limit,3);
			func2->SetParameter(0,Scale_Histo[1]);
			func2->SetParameter(1,masse);
			func2->SetParLimits(1,masse,masse);
			func2->SetParameter(2,0.2);
			func2->SetParNames("C", "m", "T");	
			func2->SetLineColor(kGreen);
	
		TF1 *func1 = new TF1("expo_law",expo_law,0.3,upper_limit,3);
			func1->SetParameter(0,Scale_Histo[1]);
			func1->SetParameter(1,masse);
			func1->SetParLimits(1,masse,masse);
			func1->SetParameter(2,0.1);
			func1->SetParNames("dyN","m","T");
			func1->SetLineColor(kViolet);
			
			
		TF1 *func3 = new TF1("levy",levy,0,upper_limit,4);
			func3->SetParameter(0,Scale_Histo[1]);
			//func3->SetParLimits(0,74.56,74.56);
			func3->SetParameter(1,masse);
			func3->SetParLimits(1,masse,masse);
			func3->SetParameter(2,0.6);
			func3->SetParameter(3,7.79);
			//func3->SetParLimits(3,7.79,7.79);
			//func3->SetParLimits(3,2,8);
			func3->SetParNames("C","m","T","n");	
			//func3->SetParLimits(3, 1.01 , 10);
			func3->SetLineColor(kCyan);
			
	
	
	
		auto fitResult_2 = hist->Fit(func2,"S R I 0");
			auto covMatrix_2 = fitResult_2->GetCovarianceMatrix();
			double chi2_bol = func2->GetChisquare();
			double ndf_bol = func2->GetNDF();
			double chindf_bol = chi2_bol/ndf_bol ;
			cout <<" Chi square boltz = "<< chi2_bol << endl;
			cout <<" NDF boltz = "<< ndf_bol << endl; 
			cout <<" Chi2 / ndf = " << chindf_bol << endl;
			cout <<" integral value " << func2->Integral(0.3,3);
	
		//hist->Fit("expo_law","+ I  ","SAMES HIST",0.3,upper_limit);  

		auto fitResult = hist->Fit(func1,"S R I 0");
			auto covMatrix = fitResult->GetCovarianceMatrix();
			double chi2_expo = func1->GetChisquare();
			double ndf_expo = func1->GetNDF();
			double chindf_expo = chi2_expo/ndf_expo ;
			cout <<" Chi square expo law = "<< chi2_expo << endl;
			cout <<" NDF expo law = "<< ndf_expo << endl; 
			cout <<" Chi2 / ndf = " << chindf_expo << endl;
			cout <<" integral value " << func1->Integral(0,20);
			cout <<" integral error value " << func1->IntegralError(0,20, fitResult->GetParams() , covMatrix.GetMatrixArray());
			
			
			
		auto fitResult_3 = hist->Fit(func3,"S R I 0");
			auto covMatrix_3 = fitResult_3->GetCovarianceMatrix();
			double chi2levy = func3->GetChisquare();
			double ndf_levy = func3->GetNDF();
			double chindf_levy = chi2levy/ndf_levy ;
			cout <<" Chi square levy = "<< chi2levy << endl;
			cout <<" NDF levy = "<< ndf_levy << endl; 
			cout <<" Chi2 / ndf = " << chindf_levy << endl;
			cout <<" integral value " << func3->Integral(0,20);
			cout <<" integral error value " << func3->IntegralError(0,20, fitResult_3->GetParams() , covMatrix_3.GetMatrixArray());
			cout <<" extrapolated yield = " <<  (100* (func3->Integral(0,0.3)/func3->Integral(0,20))) <<"		"<< (100*(func3->Integral(0,0.3)/func3->Integral(0,20))*TMath::Sqrt( TMath::Power( func3->IntegralError(0,a, fitResult_3->GetParams() , covMatrix_3.GetMatrixArray()) / func3->Integral(0,20),2)  +  TMath::Power( func3->IntegralError(0,b, fitResult_3->GetParams() , covMatrix_3.GetMatrixArray())/ func3->Integral(0,b) ,2) ) ) << endl;
			
			
			
// FIT LEVY
		TMinuit *gMinuit_LEVY = new TMinuit(4);
			gMinuit_LEVY->SetFCN(fcn_levy);
			gMinuit_LEVY->DefineParameter(0, "C", Scale_Histo[1], 0.1, 0, 2*Scale_Histo[1]);
			//gMinuit_LEVY->FixParameter(0);
			gMinuit_LEVY->DefineParameter(1, "m_0",masse, 0.01, 0.1,2);
			gMinuit_LEVY->FixParameter(1);
			gMinuit_LEVY->DefineParameter(2, "T", 0.6, 0.1, 0, 2);
			gMinuit_LEVY->DefineParameter(3, "n", 7.79, 0.1, 0, 20);
			//gMinuit_LEVY->FixParameter(3);	
			gMinuit_LEVY->Command("MIGRAD");
			gMinuit_LEVY->Command("MINOS");
			ROMAIN->cd();
			for(int i=0;i<5;i++) gMinuit_LEVY->GetParameter(i,parlevy[i],errlevy[i]);
			TH1F* curve1 = new TH1F("curve1","curve1",hist->GetNbinsX()*5,0,b);
				for(int i=1;i<=curve1->GetNbinsX();i++) 
				{		
					double *x = new double ;
					x[0] = curve1->GetBinCenter(i);
					double f = levy(x,parlevy);
					//double f = levy(x,par);
					curve1->SetBinContent(i,f);
				} 
		TLatex latex;
			latex.SetNDC();
			latex.SetTextSize(0.030);
			latex.SetTextAlign(13);  //align at top
			latex.SetTextFont(42);

		curve1->SetLineWidth(3);
		hist->Draw("same");
		curve1->Draw("csame");
		latex.DrawLatex(0.6,0.80, TString::Format("Parameter 1 : dN/dy = %s #pm 0.002",Form("%.3f",parlevy[0])) );
		latex.DrawLatex(0.6,0.75, TString::Format("Parameter 2 : m = %s GeV/c^{2}",Form("%.3f",parlevy[1])));
		latex.DrawLatex(0.6,0.70, TString::Format("Parameter 3 : T = (%s #pm 0.002) GeV", Form("%.3f",parlevy[2])));
		latex.DrawLatex(0.6,0.65, TString::Format("Parameter 4 : n = %s #pm 0.08", Form("%.2f",parlevy[3])));
		latex.DrawLatex(0.6,0.60, TString::Format("#chi^{2} = %s ", Form("%.1f",chi2_levy(parlevy))));
		latex.DrawLatex(0.6,0.55, TString::Format("ndf = %d ", (Nbinx-gMinuit_LEVY->GetNumFreePars())));
		latex.DrawLatex(0.6,0.50, TString::Format("#chi^{2}/ndf = %s ",	Form("%.2f",(chi2_levy(parlevy)/(Nbinx-gMinuit_LEVY->GetNumFreePars())))))	;

// we get the covariant matrix from the Tminuit instance

		int npars_levy = gMinuit_LEVY->GetNumPars();
		double *covar_levy = new double[npars_levy*npars_levy];
		gMinuit_LEVY->mnemat(covar_levy,npars_levy);	
		
		Data_Model(X,Y,Xerr,Yerr,levy,parlevy,par_size_levy,covar_levy);
		Integral_value = Integral(levy,parlevy,par_size_levy,covar_levy,0,20) ;
		Integral_error = Integral_E(levy,parlevy,par_size_levy,covar_levy,0,20) ;
		
		cout << " test integral " << Integral_value << endl;
// we create the graph to plot Data/Model fit with rectangle errors for phi mesons
		DATA_MODEL->cd(1);
		gPad->SetTickx(2);
		gPad->SetBottomMargin(small); 
		gPad->SetTickx();
// we create the graph to plot Data/Model fit with rectangle errors for p+pbar
		DATA_MODEL->cd(2); 
		gPad->SetTopMargin(small);
// we create the graph to plot Data/Model fit with rectangle errors
	   	TGraphErrors *dataModel = new TGraphErrors(Nbinx,X,Y,Xerr,Yerr);
		dataModel->SetName( Form(name+"_data_model_%d",digit) );
		dataModel->SetLineColor(digit);
		dataModel->SetLineWidth(2);
		dataModel->SetMarkerColor(digit);
		dataModel->SetMarkerSize(1);
		dataModel->SetMarkerStyle(21);
		dataModel->SetTitle(" ; p_T (GeV/c); Data / model fit");
		dataModel->SetLineWidth(2);
		dataModel->GetXaxis()->SetLimits(0.2,upper_limit+hist->GetBinWidth(UpperBin_Fit(upper_limit)));
		dataModel->Draw("A5 SAME");   // A has to be there (if not, nothing appears on the canvas) & 5 is to plot error rectangles
		dataModel->Draw("PX SAME");   // P is to put the chosen marker & X to remove the error bars (leave only the rectangles)
		dataModel->Write();	
//We Store the graph in a Multigraph
		DM->Add(dataModel);
// Export of the fit result into result.txt
		myfile <<"Centrality fichier "<< digit<<"	" << Integral_value << "		" << Integral_error << "	" << 100*Integral(levy,parlevy,par_size_levy,covar_levy,0,a)/Integral(levy,parlevy,par_size_levy,covar_levy,0,b) << "			" <<(100*Integral(levy,parlevy,par_size_levy,covar_levy,0,a)/Integral(levy,parlevy,par_size_levy,covar_levy,0,b)*TMath::Sqrt( TMath::Power( Integral_E(levy,parlevy,par_size_levy,covar_levy,0,a) / Integral(levy,parlevy,par_size_levy,covar_levy,0,a),2)  +  TMath::Power( Integral_E(levy,parlevy,par_size_levy,covar_levy,0,b) / Integral(levy,parlevy,par_size_levy,covar_levy,0,b) ,2) ) ) <<"		"<<chi2_levy(parlevy)<< "		"<<(Nbinx-gMinuit_LEVY->GetNumFreePars())<<"		"<<(chi2_levy(parlevy)/(Nbinx-gMinuit_LEVY->GetNumFreePars()))<<endl;
		OutputHisto->cd();
		ROMAIN2->cd();
		ROMAIN->Write();
		DM->Draw("A5");
		DM->Draw("PX");
		ROMAIN2->Write();
		TLine *line = new TLine(a,1,upper_limit+hist->GetBinWidth(UpperBin_Fit(upper_limit)),1);
		line->SetLineColor(30);
		line->SetLineWidth(3);
		line->SetLineStyle(2);
		line->Draw("SAME");
		ROMAIN2->Close();
		DATA_MODEL->Close();
		CONTOUR->Close();
		test1->Close();
		OutputHisto->Close();
		}
	}
// MAIN END
}
/*
// Création du fichier root de sortie et récupération du fichier root données
// Récupération des différents histogrammes 
// Changement des noms d'axes et titres
//Récupération du scaling en prépartion du fit
	double Scale_Histo[11] ;
	Scale_Histo[1]=H1_pp->Integral("width") ;
// Compilations des différentes erreurs et ajout sur l'histogramme PT
	int Nbinx = H1_pp->GetNbinsX();
	for(int i = 0; i <= Nbinx  ; i++){
		H1_pp->SetBinError(i, sqrt( pow(H1_pp_e1->GetBinContent(i),2) + pow(H1_pp_E1->GetBinContent(i),2) ) ) ;
	}		
	//H1_pp->SetAxisRange(0.3,3,"X");
	//H1_pp->SetAxisRange(0.3,3,"X");
	double a= H1_pp->GetBinLowEdge(1);
	double b= H1_pp->GetBinLowEdge(Nbinx);
	b += H1_pp->GetBinWidth(Nbinx);
// Création des fonctions modèle pour le fit
			
	TF1 *func1 = new TF1("expo_law",expo_law,0,b,3);
		func1->SetParameter(0,Scale_Histo[1]);
		func1->SetParameter(1,0.938);
		func1->SetParLimits(1,0.938,0.938);
		func1->SetParameter(2,0.1);
		func1->SetParNames("C","m","T");
		func1->SetLineColor(kViolet);
	TF1 *func2 = new TF1("boltzmann",boltzmann,0,b,3);
		func2->SetParameter(0,Scale_Histo[1]);
		func2->SetParameter(1,0.938);
		func2->SetParLimits(1,1.0194,1.0194);
		func2->SetParameter(2,0.2);
		func2->SetParNames("C", "m", "T");	
		func2->SetLineColor(kGreen);
	TF1 *func3 = new TF1("levy",levy,0,b,4);
		func3->SetParameter(0,0.2331);
		//func3->SetParLimits(0,74.56,74.56);
		func3->SetParameter(1,0.938);
		func3->SetParLimits(1,0.938,0.938);
		func3->SetParameter(2,0.6);
		func3->SetParameter(3,7.79);
		//func3->SetParLimits(3,7.79,7.79);
		//func3->SetParLimits(3,2,8);
		func3->SetParNames("C","m","T","n");	
		//func3->SetParLimits(3, 1.01 , 10);
		func3->SetLineColor(kCyan);
	TF1 *func4 = new TF1("power_law_Five", power_law_Five,0,b,3);
		func4->SetParameter(0,Scale_Histo[1]);
		func4->SetParameter(1,0.938);
		func4->SetParameter(2,4);
		func4->SetParNames("C","p0","n");	
		func4->SetLineColor(kMagenta);
	TF1 *func5 = new TF1("power_law_Seven", power_law_Seven ,0,b,3);
		func5->SetParameter(0,Scale_Histo[1]);
		func5->SetParameter(1,1.0194);
		func5->SetParameter(2,1);
		func5->SetParNames("C","p0","n");	
		func5->SetLineColor(kRed);
	TF1 *func6 = new TF1("blast_wave",blast_wave,0,b,5);
		func6->SetParameter(0,72.74);
		func6->SetParLimits(0,0.2331,0.2331);
		func6->SetParameter(1,0.938);
		func6->SetParLimits(1, 0.938,0.938);
		func6->SetParameter(2,0.09);
		func6->SetParameter(3,0.73);
		func6->SetParameter(4,0.66);
		func6->SetParNames("A","m","T","n","beta_s");
		func6->SetLineColor(kBlue);
TCanvas *test1 = new TCanvas("test1","The FillRandom example",0,0,700,500);
test1->SetFillColor(18);
test1->cd();
TH1F* HF = new TH1F("HF","TEST RANDOM BW", 250 , 0 , 20 );
HF->SetFillColor(45);
HF->FillRandom("blast_wave",10000);
HF->Draw();
 hist = (TH1F *)dirFile->Get("Hist1D_y1");
 HF->Scale(hist->Integral()/HF->Integral());
		TH1F* hist_e1=(TH1F*)dirFile->Get("Hist1D_y1_e1");
		TH1F* hist_E1=(TH1F*)dirFile->Get("Hist1D_y1_e2");
		hist->SetXTitle("p_{T} [GeV/c]");
		hist->SetYTitle("(1/Nev)*d^2(N)/dPtdYrap  [Gev/c] ");
	for(int i = 0; i <= Nbinx  ; i++)
	{
		HF->SetBinError(i, sqrt( pow(hist_e1->GetBinContent(i),2) + pow(hist_E1->GetBinContent(i),2) ) ) ;
		hist->SetBinError(i, sqrt( pow(hist_e1->GetBinContent(i),2) + pow(hist_E1->GetBinContent(i),2) ) ) ;
	}		
		//hist->GetXaxis()->SetRangeUser(0.3,3);
		
		
	TMinuit *gMinuit = new TMinuit(5);
		gMinuit->SetFCN(fcn);
// Paramètre blast wave
		gMinuit->DefineParameter(0, "A",72.74, 0.01, 0, 100);
		gMinuit->FixParameter(0);
		gMinuit->DefineParameter(1, "m_0",0.938, 0.01, 0.1,1);
		gMinuit->FixParameter(1);
		gMinuit->DefineParameter(2, "T", 0.09, 0.0001, 0, 1);
		gMinuit->DefineParameter(3, "n", 1, 0.01, 0, 20);
		//gMinuit->FixParameter(3);	
		gMinuit->DefineParameter(4, "beta_s", 0.66, 0.0001, 0, 1);
		gMinuit->DefineParameter(0, "C",Scale_Histo[1], 0.01, 0, 1);
		gMinuit->DefineParameter(1, "m_0",0.938, 0.01, 0.1,1);
		gMinuit->FixParameter(1);
		gMinuit->DefineParameter(2, "T", 0.3, 0.01, 0, 2);
		gMinuit->DefineParameter(3, "n", 7.79, 0.01, 0, 10);
		gMinuit->FixParameter(3);
		
		gMinuit->Command("MIGRAD");
		gMinuit->Command("MINOS");
		
		double par[5],err[5];
		for(int i=0;i<5;i++) gMinuit->GetParameter(i,par[i],err[i]);
		TH1F* curve = new TH1F("curve","curve",hist->GetNbinsX()*5,0,b);
			for(int i=1;i<=curve->GetNbinsX();i++) 
			{		
				double *x = new double ;
				x[0] = curve->GetBinCenter(i);
				double f = blast_wave(x,par);
				//double f = levy(x,par);
				curve->SetBinContent(i,f);
			} 
OutputHisto->cd();	
TCanvas *ROMAIN = new TCanvas("FIT BLAST-WAVE" , " PT DISTRIBUTION " , 200 , 200 );
ROMAIN->cd();
		curve->SetLineWidth(3);
		hist->Draw();
		curve->Draw("csame");
ROMAIN->Write();
// Contour
TCanvas *ROMAIN2 = new TCanvas("Contour" , " " , 200 , 200 );
ROMAIN2->cd();
		TGraph* graph1 =(TGraph*) gMinuit->Contour(25,4,2);
		gMinuit->SetErrorDef(4);
		TGraph* graph2 =(TGraph*) gMinuit->Contour(25,4,2);
		graph2->Draw("alp");
		graph1->Draw("same alf");
ROMAIN2->Write();
// Création des canvas
	TCanvas *T1 = new TCanvas("Canvas 1" , " CANVAS_1_PP1 " , 200 , 200 );
	T1->SetGrid();
	T1->SetLogy();
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(1111);
// Sauvegarde des canvas dans l'output-histo-file càd notre result.root
// We fit the data with the different fitting functions
	T1->cd();
	H1_pp->Draw("");
	//func6->Draw();
	H1_pp->Fit("expo_law","+ I  ","SAMES HIST",0,b);  
		double chi2_expo = func1->GetChisquare();
		double ndf_expo = func1->GetNDF();
		double chindf_expo = chi2_expo/ndf_expo ;
		cout <<" Chi square expo law = "<< chi2_expo << endl;
		cout <<" NDF expo law = "<< ndf_expo << endl; 
		cout <<" Chi2 / ndf = " << chindf_expo << endl;
	
	H1_pp->Fit("boltzmann","+ I ","SAMES HIST",0,b);
		double chi2_bol = func2->GetChisquare();
		double ndf_bol = func2->GetNDF();
		double chindf_bol = chi2_bol/ndf_bol ;
		cout <<" Chi square boltz = "<< chi2_bol << endl;
		cout <<" NDF boltz = "<< ndf_bol << endl; 
		cout <<" Chi2 / ndf = " << chindf_bol << endl;
	H1_pp->Fit("levy","I +","SAMES HIST",0,b);
		double chi2_levy = func3->GetChisquare();
		double ndf_levy = func3->GetNDF();
		double chindf_levy = chi2_levy/ndf_levy ;
		cout <<" Chi square levy = "<< chi2_levy << endl;
		cout <<" NDF levy = "<< ndf_levy << endl; 
		cout <<" Chi2 / ndf = " << chindf_levy << endl;
		
	H1_pp->Fit("power_law_Five"," I +","SAMES HIST",0,b);
		double chi2_PL_5 = func4->GetChisquare();
		double ndf_PL_5 = func4->GetNDF();
		double chindf_PL_5 = chi2_PL_5/ndf_PL_5 ;
		cout <<" Chi square power law = "<< chi2_PL_5 << endl;
		cout <<" NDF power law = "<< ndf_PL_5 << endl; 
		cout <<" Chi2 / ndf = " << chindf_PL_5 << endl;
	H1_pp->Fit("power_law_Seven","+ I ","SAMES HIST",0,b);
		double chi2_PL_7 = func5->GetChisquare();
		double ndf_PL_7 = func5->GetNDF();
		double chindf_PL_7 = chi2_PL_7/ndf_PL_7 ;
		cout <<" Chi square levy = "<< chi2_PL_7 << endl;
		cout <<" NDF levy = "<< ndf_PL_7 << endl; 
		cout <<" Chi2 / ndf = " << chindf_PL_7 << endl;
	//H1_pp -> Fit("blast_wave", " I V","SAMES HIST",a,b);
	H1_pp->GetFunction("expo_law")->SetLineColor(kViolet);
	H1_pp->GetFunction("boltzmann")->SetLineColor(kGreen);
	H1_pp->GetFunction("levy")->SetLineColor(kCyan);
	H1_pp->GetFunction("power_law_Five")->SetLineColor(kMagenta);
	//H1_pp->GetFunction("power_law_Seven")->SetLineColor(kRed);
	//H1_pp->GetFunction("blast_wave")->SetLineColor(kBlue);
// Ajout des résultats des fits sur le canvas
	TLatex latex;
		latex.SetNDC();
		latex.SetTextSize(0.020);
		latex.SetTextAlign(13);  //align at top
		latex.SetTextFont(60);
	   
	 
		latex.DrawLatex(0.5,0.65, TString::Format("Parametre Levy C= %g", func3->GetParameter(0)));
		latex.DrawLatex(0.5,0.625, TString::Format("Parametre Levy m= %g GeV/c", func3->GetParameter(1)));
		latex.DrawLatex(0.5,0.60, TString::Format("Parametre Levy T= %g GeV", func3->GetParameter(2)));
		latex.DrawLatex(0.5,0.575, TString::Format("Parametre Levy n= %g ", func3->GetParameter(3)));
		latex.DrawLatex(0.5,0.550, TString::Format("CHI2/NDF Levy = %g ", func3->GetChisquare()/func3->GetNDF()));
	  
		latex.DrawLatex(0.5,0.50, TString::Format("Parametre expo C= %g", func1->GetParameter(0)));
		latex.DrawLatex(0.5,0.475, TString::Format("Parametre expo m_0= %g GeV/c", func1->GetParameter(1)));
		latex.DrawLatex(0.5,0.45, TString::Format("Parametre expo T= %g GeV/c", func1->GetParameter(2)));
		latex.DrawLatex(0.5,0.425, TString::Format("CHI2/NDF expo = %g ", func1->GetChisquare()/func1->GetNDF()));
	 
		latex.DrawLatex(0.5,0.375, TString::Format("Parametre boltzmann C= %g", func2->GetParameter(0)));
		latex.DrawLatex(0.5,0.350, TString::Format("Parametre boltzmann m_0= %g GeV/c", func2->GetParameter(1)));
		latex.DrawLatex(0.5,0.325, TString::Format("Parametre boltzmann T= %g GeV/c", func2->GetParameter(2)));
		latex.DrawLatex(0.5,0.30, TString::Format("CHI2/NDF boltzmann = %g ", func2->GetChisquare()/func2->GetNDF()));
		
		latex.DrawLatex(0.5,0.25, TString::Format("Parametre power_law C= %g", func4->GetParameter(0)));
		latex.DrawLatex(0.5,0.225, TString::Format("Parametre power_law m_0= %g GeV/c", func4->GetParameter(1)));
		latex.DrawLatex(0.5,0.20, TString::Format("Parametre power_law T= %g GeV/c", func4->GetParameter(2)));
		latex.DrawLatex(0.5,0.175, TString::Format("CHI2/NDF power_law = %g ", func4->GetChisquare()/func4->GetNDF()));
	T1->Update();
// Configuration de la légende
	auto legend = new TLegend(0.4,0.7,0.6,0.9);
		//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
		legend->AddEntry(H1_pp,"Data","lep");
		legend->AddEntry(func1,"Exponential law fit","l");
		legend->AddEntry(func2,"Boltzmann law fit","l");
		legend->AddEntry(func3,"Levy-Tsallis fit","l");
		//legend->AddEntry(func4, "Power law fit from [5]", "l");
		//legend->AddEntry(func5, "Power law fit from [7]", "l");
		//legend->AddEntry(func6,"Blast-wave model","l");
		legend->Draw("SAMES");
		
	T1->Write();
	cout << " Scale after = " << H1_pp->Integral("width") << endl;
	double_t c= func3->Integral(0,20,1E-12);
	cout << " integral value = " << c << endl;
	OutputHisto->Close();
	cout << " fin " << endl ;
	*/
