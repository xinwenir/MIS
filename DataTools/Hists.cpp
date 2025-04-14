// Hists.cpp
// J. Matousek

#include "Hists.h"

Plot::Plot(TVirtualPad* targPad, int nHists, const char* plotName,
			 const char* plotTitle, const char* xTitle, const char* yTitle)
{
	pad = targPad;
	n = nHists + 1;
	nLeg   = 0;
	nLines = 0;
	nTexts = 0;
	hists  = new TH1*[n];
	graphs = new TGraphErrors*[n];
	func   = new TF1*[n];
	legLabels = new char*[n];
	flags     = new int[n];
	legend    = NULL;
	xLog      = false;
	yLog      = false;
	normalize = false;
	textSize  = 0.07;
	yMinUser  = NAN;
	yMaxUser  = NAN;
	xMinUser  = NAN;
	xMaxUser  = NAN;
	legX0     = -1;
	legY0     = -1;
	lMargin   = -1;
	rMargin   = -1;
	bMargin   = -1;
	tMargin   = -1;
	xLabelOffset = -1;
	yLabelOffset = -1;
	xTitleOffset = -1;
	yTitleOffset = -1;

	char* buffname = new char[strlen(plotName) + 1];
	sprintf(buffname, "%s", plotName);
	name = buffname;

	if (xTitle == NULL)
		xLabel = NULL;
	else
	{
		char* buffxl = new char[strlen(xTitle) + 1];
		sprintf(buffxl, "%s", xTitle);
		xLabel = buffxl;
	}
	if (yTitle == NULL)
		yLabel = NULL;
	else
	{
		char* buffyl = new char[strlen(yTitle) + 1];
		sprintf(buffyl, "%s", yTitle);
		yLabel = buffyl;
	}
	if (plotTitle == NULL)
		title = NULL;
	else
	{
		char* bufft = new char[strlen(plotTitle) + 1];
		sprintf(bufft, "%s", plotTitle);
		title = bufft;
	}

	for (int i = 0; i < n; i++)
	{
		hists[i] = NULL;
		graphs[i] = NULL;
		func[i] = NULL;
		legLabels[i] = NULL;
	}
}

Plot::~Plot()
{
	for (int i = 0; i < n; i++)
	{
		if (hists[i] != NULL)
			delete hists[i];
		if (legLabels[i] != NULL)
			delete legLabels[i];
		if (graphs[i] != NULL)
			delete graphs[i];
		if (func[i] != NULL)
			delete func[i];
	}
	delete[] flags;
	delete[] hists;
	delete[] graphs;
	delete[] func;
	delete[] legLabels;
	delete[] name;
	if (xLabel != NULL)
		delete[] xLabel;
	if (yLabel != NULL)
		delete[] yLabel;
	if (title != NULL)
		delete[] title;
}

void Plot::Add(const TH1* h, int style, int index, const char* leg)
{
	int i;
	for (i = 1; i < n; i++)
	{
		if (hists[i] == NULL && graphs[i] == NULL && func[i] == NULL)
			break;
	}
	if (i == n)
	{
		printf("Plot '%s' capacity full, cannot add more objects...\n", name);
		return;
	}
	char buff[64];
	sprintf(buff, "%s_h%i", name, i);
	hists[i] = (TH1*)h->Clone(buff);
	flags[i] = style;

	if (strcmp(leg,"") != 0)
	{
		legLabels[i] = new char[strlen(leg) + 1];
		sprintf(legLabels[i],"%s",leg);
		nLeg++;
	}

	if (index == -1)
		index = i-1;

	if (style & Line || style & Hist)
	{
		hists[i]->SetLineColor(Hists::GetMyLineColor(index));
		hists[i]->SetLineStyle(index % 2 + 1);
		hists[i]->SetLineWidth(2 + (index % 2));
	}
	if (style & Fill)
	{
		hists[i]->SetFillColor(Hists::GetMyFillColor(index));
		hists[i]->SetFillStyle(3001 + i);
		if ( ! (style & Line) )
			hists[i]->SetLineColor(Hists::GetMyFillColor(index));
	}
	if (style & Erry)
	{
		hists[i]->SetLineColor(Hists::GetMyLineColor(index));
		hists[i]->SetLineWidth(2);
	}
	if (style & Point)
	{
		hists[i]->SetMarkerColor(Hists::GetMyLineColor(index));
		hists[i]->SetMarkerSize(1.2);
		hists[i]->SetMarkerStyle(19+i);
	}
}

void Plot::Add(const TH1* h, const char* leg, int style, int index)
{
	Add(h, style, index, leg);
}

void Plot::Add(const TGraph* h, int style, int index, const char* leg)
{
	int i;
	for (i = 1; i < n; i++)
	{
		if (hists[i] == NULL && graphs[i] == NULL && func[i] == NULL)
			break;
	}
	if (i == n)
	{
		printf("Plot '%s' capacity full, cannot add more objects...\n", name);
		return;
	}
	char buff[64];
	sprintf(buff, "%s_g%i", name, i);
	graphs[i] = new TGraphErrors(h->GetN());
	for (int j = 0; j < h->GetN(); j++)
	{
		graphs[i]->SetPoint(j,h->GetX()[j],h->GetY()[j]);
		graphs[i]->SetPointError(j,h->GetEX()[j],h->GetEY()[j]);
	}
	// make a histogram, which takes care of the graphics:
	pad->cd();
	graphs[i]->Draw();
	hists[i] = graphs[i]->GetHistogram();
	flags[i] = style;

	if (strcmp(leg,"") != 0)
	{
		legLabels[i] = new char[strlen(leg) + 1];
		sprintf(legLabels[i],"%s",leg);
		nLeg++;
	}

	if (index == -1)
		index = i-1;
	graphs[i]->SetLineWidth(2 + (index % 2));
	graphs[i]->SetFillColor(0);
	graphs[i]->SetMarkerStyle(0);
	graphs[i]->SetLineColor(0);

	if (style & Line)
	{
		graphs[i]->SetLineColor(Hists::GetMyLineColor(index));
		graphs[i]->SetLineStyle(index % 2 + 1);
	}
	if (style & Point)
	{
		graphs[i]->SetMarkerColor(Hists::GetMyLineColor(index));
		graphs[i]->SetMarkerSize(1.2);
		graphs[i]->SetMarkerStyle(19+i);
	}
	if (style & Erry)
	{
		graphs[i]->SetLineColor(Hists::GetMyLineColor(index));
		graphs[i]->SetLineWidth(2);
	}
}

void Plot::Add(const TGraph* h, const char* leg, int style, int index)
{
	Add(h, style, index, leg);
}

void Plot::Add(const TF1* f, int index, const char* leg)
{
	int i;
	for (i = 1; i < n; i++)
	{
		if (hists[i] == NULL && graphs[i] == NULL && func[i] == NULL)
			break;
	}
	if (i == n)
	{
		printf("Plot '%s' capacity full, cannot add more objects...\n", name);
		return;
	}
	char buff[64];
	sprintf(buff, "%s_g%i", name, i);
	func[i] = new TF1(*f);
	// make a histogram, which takes care of the graphics:
	pad->cd();
	func[i]->Draw();
	hists[i] = func[i]->GetHistogram();

	if (index == -1)
		index = i-1;
	func[i]->SetLineWidth(2 + (index != 0));
	func[i]->SetFillColor(0);
	func[i]->SetLineColor(Hists::GetMyLineColor(index));
	func[i]->SetLineStyle(index % 2 + 1);

	if (strcmp(leg,"") != 0)
	{
		legLabels[i] = new char[strlen(leg) + 1];
		sprintf(legLabels[i],"%s",leg);
		nLeg++;
	}
}

void Plot::Add(const TF1* f, const char* leg, int index)
{
	Add(f, index, leg);
}

void Plot::AddLine(float x0, float y0, float x1, float y1)
{
	lineCoords[nLines][0] = x0;
	lineCoords[nLines][1] = y0;
	lineCoords[nLines][2] = x1;
	lineCoords[nLines][3] = y1;
	nLines++;
}

void Plot::AddText(float x0, float y0, const char* text, float size)
{
	texts[nTexts] = new TLatex(x0,y0,text);
	texts[nTexts]->SetTextFont(132);
	if (size == -1)
		texts[nTexts]->SetTextSize(textSize);
	else
		texts[nTexts]->SetTextSize(size);

	nTexts++;
}

void Plot::Draw(bool verb)
{
	if (verb)
		printf("Drawing %s into %s...\n", name,pad->GetName());
	if (hists[1] == NULL)
	{
		printf("Nothing to draw in '%s'.\n", name);
		return;
	}
	if (normalize && graphs[1] == NULL)
		hists[1]->Scale(1.0 / hists[1]->Integral());

	float xMin = hists[1]->GetXaxis()->GetXmin();
	float xMax = hists[1]->GetXaxis()->GetXmax();
	float yMin = hists[1]->GetMinimum();
	float yMax = hists[1]->GetMaximum();
	if (yLog)
		yMin = Hists::FindSmallestAbove(hists[1],0.0) * 0.9;

	float scale = 1.0;
	if (flags[1] & Erry && graphs[1] == NULL)
	{
		if (hists[1]->GetSumw2N() != 0)
			scale = 1.0 / hists[1]->GetSumw2N();
		if (yMin - sqrt(yMin) * scale > 0)
			yMin -= sqrt(yMin) * scale;
		yMax += sqrt(yMax) * scale;
	}
	for (int i = 2; i < n; i++)
	{
		if (hists[i] != NULL)
		{
			if (normalize && graphs[i] == NULL)
				hists[i]->Scale(1.0 / hists[i]->Integral());
			if (hists[i]->GetXaxis()->GetXmin() < xMin)
				xMin = hists[i]->GetXaxis()->GetXmin();
			if (hists[i]->GetXaxis()->GetXmax() > xMax)
				xMax = hists[i]->GetXaxis()->GetXmax();

			float val;
			scale = 1.0;
			if (graphs[i] == NULL && hists[i]->GetSumw2N() != 0)
				scale = 1.0 / hists[i]->GetSumw2N();
			if (yLog)
				val = Hists::FindSmallestAbove(hists[i],0.0) * 0.9;
			else
				val = hists[i]->GetMinimum();

			if (graphs[i] == NULL && flags[i] & Erry
					&& val-sqrt(val)*scale < yMin && val-sqrt(val)*scale > 0)
				yMin = val - sqrt(val) * scale;
			else if (val < yMin)
				yMin = val;

			val = hists[i]->GetMaximum();
			if (graphs[i] == NULL && flags[i] & Erry
					&& val + sqrt(val)*scale > yMax)
				yMax = val + sqrt(val) * scale;
			else if (val > yMax)
				yMax = val;
		}
	}
	if (!isnan(yMinUser))
		yMin = yMinUser;
	if (!isnan(yMaxUser))
		yMax = yMaxUser;
	if (!isnan(xMinUser))
		xMin = xMinUser;
	if (!isnan(xMaxUser))
		xMax = xMaxUser;

	hists[0] = new TH1F(name, name, hists[1]->GetXaxis()->GetNbins(),xMin,xMax);
	if (graphs[1] != NULL || flags[1] & Erry)	// put y axis min lower
	{
		if (yLog)
			hists[0]->GetYaxis()->SetRangeUser(pow(yMin,0.9),pow(yMax,1.1));
		else
			hists[0]->GetYaxis()->SetRangeUser(yMin - 0.1 * (yMax - yMin),
											   yMax + 0.1 * (yMax - yMin));
	}
	else					// histograms - put y axis minimum at hist. minimum
		hists[0]->GetYaxis()->SetRangeUser(yMin,yMax + 0.1 * (yMax - yMin));
	Hists::SetMyTextOptions(hists[0],textSize);

	// ********** set labels and margins *********************
	if (xLabel != NULL)
	{
		hists[0]->GetXaxis()->SetTitle(xLabel);
		hists[0]->GetXaxis()->SetTitleOffset(1.05);
		pad->SetBottomMargin(0.16 * textSize / 0.07);
	}
	if (yLabel != NULL)
	{
		hists[0]->GetYaxis()->SetTitle(yLabel);
		hists[0]->GetYaxis()->SetTitleOffset(1.24);
		pad->SetLeftMargin(0.2 * textSize / 0.07);
	}
	else
		pad->SetLeftMargin(0.16);

	pad->SetRightMargin(0.05);
	// ******* automatic margins set, override by user settings if any ******
	if (lMargin != -1)
		pad->SetLeftMargin(lMargin);
	if (rMargin != -1)
		pad->SetRightMargin(rMargin);
	if (bMargin != -1)
		pad->SetBottomMargin(bMargin);
	if (tMargin != -1)
		pad->SetTopMargin(tMargin);
	// same for offsets
	if (xLabelOffset != -1)
		hists[0]->GetXaxis()->SetLabelOffset(xLabelOffset);
	if (yLabelOffset != -1)
		hists[0]->GetYaxis()->SetLabelOffset(yLabelOffset);
	if (xTitleOffset != -1)
		hists[0]->GetXaxis()->SetTitleOffset(xTitleOffset);
	if (yTitleOffset != -1)
		hists[0]->GetYaxis()->SetTitleOffset(yTitleOffset);

	// ******************* pad title ********************************
	if (title != NULL)
		hists[0]->SetTitle(title);
	else
		hists[0]->SetTitle(hists[1]->GetTitle());
	// ******************** log-scale *******************************
	pad->cd();
	if (xLog)
		pad->SetLogx();
	if (yLog)
		pad->SetLogy();
	hists[0]->Draw("P0");
	lMargin = pad->GetLeftMargin();
	rMargin = pad->GetRightMargin();
	bMargin = pad->GetBottomMargin();
	tMargin = pad->GetTopMargin();
	// ********************* legend **********************************
	unsigned int legLen = 0;
	for(int i = 1; i < n; i++)
	{
		if (legLabels[i] != NULL && strlen(legLabels[i]) > legLen)
			legLen = strlen(legLabels[i]);
	}
	if (legX0 == -1)
		legX0 = 0.89 - rMargin - 0.018*legLen;
	if (legY0 == -1)
		legY0 = 1.0 - tMargin - textSize*nLeg*0.8;

	legend = new TLegend(legX0, legY0, 0.99 - rMargin, 0.99 - tMargin);
	legend->SetBorderSize(0);
	legend->SetTextSize(textSize * 0.75);
	legend->SetTextFont(132);

	char buff[128];
	for (int i = 1; i < n; i++)		// ***** D R A W I N G ****************
	{
		const char* erry = "";
		const char* points = "";
		const char* lines = "";

		if (flags[i] & Point)
			points = "P";

		if (flags[i] & Line)
			lines = "L";

		if (flags[i] & Erry && graphs[i] != NULL)	// graph with y errorbars
			erry = "E1";
		else if (flags[i] & Erry)					// histo with y errorbars
			erry = "E1X0";
		else if (graphs[i] != NULL)					// graph without y errorbars
			erry = "X";
		else										// histo without errorbars
			erry = "HIST";	// (hide them for 'histos with errors calculated')


		sprintf(buff,"%s%s%sSAME",erry,points,lines);
//		printf("Draw options: '%s'\n",buff);

		if (graphs[i] != NULL)
		{
			graphs[i]->Draw(buff);
			if (legLabels[i] != NULL)
				legend->AddEntry(graphs[i],legLabels[i]);
		}
		else if (func[i] != NULL)
		{
			func[i]->Draw("L SAME");
			if (legLabels[i] != NULL)
				legend->AddEntry(func[i],legLabels[i]);
		}
		else if (hists[i] != NULL)
		{
			hists[i]->Draw(buff);

			Hists::CheckOverflow(hists[i]);
			if (legLabels[i] != NULL)
				legend->AddEntry(hists[i],legLabels[i]);
		}
	}
	// draw lines if any
	if (nLines > 0)
	{
		line = new TLine();
		line->SetLineColor(Hists::GetMyLineColor(7));
		line->SetLineWidth(2);
		for (int i = 0; i < nLines; i++)
			line->DrawLine(	lineCoords[i][0], lineCoords[i][1],
							lineCoords[i][2], lineCoords[i][3]);
	}
	for (int i = 0; i < nTexts; i++)
		texts[i]->Draw();

	if (legend->GetNRows() > 0)
	{
		float scaleY = (1.0 - bMargin - tMargin) / (legY0 - bMargin) + 0.05;
		if (yLog)
		{
			if (graphs[1] != NULL)
				hists[0]->GetYaxis()->SetRangeUser(yMin / pow(yMax/yMin, 0.1),
									yMin * pow(yMax/yMin, scaleY));
			else
				hists[0]->GetYaxis()->SetRangeUser(yMin,
									yMin * pow(yMax/yMin, scaleY));
		}
		else
		{
			if (graphs[1] != NULL)
				hists[0]->GetYaxis()->SetRangeUser(yMin - 0.1 * (yMax - yMin),
										yMin + scaleY * (yMax - yMin));
			else
				hists[0]->GetYaxis()->SetRangeUser(yMin,
										yMin + scaleY * (yMax - yMin));
		}
		legend->Draw();
	}
	hists[0]->Draw("AXIS SAME");	// re-draw axes
}

void Plot::SetMargins(float l, float r, float b, float t)
{
	lMargin = l;
	rMargin = r;
	bMargin = b;
	tMargin = t;
}

void Plot::MeanGraph(TGraphErrors* gMean, bool wErr, TH1* hPulls)
{
	int N = 0;
	N = graphs[1]->GetN();
	printf("Calculating mean graph from all graphs in the plot '%s'.\n",name);
	printf("Assuming all the graphs have the same # points.\n");

	double wy = 1.0;
	double x,y;

	for (int j = 0; j < N; j++)
	{
		double mx  = 0.0;
		double swx = 0.0;
		double my  = 0.0;
		double swy = 0.0;
		for (int i = 1; i < n; i++)
		{
			if (graphs[i] != NULL)
			{
				graphs[i]->GetPoint(j,x,y);
				wy = 1.0;
				if (graphs[i]->GetErrorY(j) > 0)
					wy /= ( graphs[i]->GetErrorY(j) * graphs[i]->GetErrorY(j) );
				mx  += x;
				swx += 1.0;
				my  += y * wy;
				swy += wy;
			}
		}
		mx /= swx;
		my /= swy;
		gMean->SetPoint(j,mx,my);

		double schi = 0.0;
		for (int i = 1; i < n; i++)
		{
			graphs[i]->GetPoint(j,x,y);
			if (hPulls != NULL)
				hPulls->Fill(y - my);

			wy = 1.0;
			if (graphs[i]->GetErrorY(j) > 0)
				wy /= ( graphs[i]->GetErrorY(j) * graphs[i]->GetErrorY(j) );
			schi += (y - my) * (y - my) * wy;
		}

		if (wErr)	// calculate error as for a weighted mean:
		{			// (assume all y are measurements of the same quantity)
			gMean->SetPointError(j,0.0, sqrt(1.0 / swy / (n-2) * schi) );
		}
		else	// calculate error as for a sum of independent quantities:
		{		// (in this case sqrt( 1 / (sum of weights) ) )
			gMean->SetPointError(j,0.0,sqrt(1.0 / swy));
		}
	}
}

void Plot::MeanHisto(TH1* hMean, TH1* hPulls)
{
	printf("Plot::MeanHisto: Not implemented yet.\n");
}

// ***************************************************************************

PanelPlot::PanelPlot(TVirtualPad* targPad, int nCol, int nRow, int nHistsPerPad,
					const char* plotName, const char** xAxisTitles,
					const char** yAxisTitles)
{
	motherPad = targPad;
	nx = nCol;
	ny = nRow;
	nHists = nHistsPerPad;
	char* buffname = new char[strlen(plotName) + 1];
	sprintf(buffname, "%s", plotName);
	name = buffname;

	xTitles = new char*[nx];
	yTitles = new char*[ny];

	char buff[32];
	double xlow, xup, ylow, yup;
	float lMargin = 0.1;			// left margin
	float bMargin = 0.1;			// bottom margin
	float relLM = lMargin / ( lMargin + (0.98 - lMargin) / nx );
	float relBM = bMargin / ( bMargin + (0.98 - bMargin) / nx );
	xup = 0.0;
	yup = 0.0;
	plots = new Plot**[nx];
	pads  = new TPad**[nx];

	for (int i = 0; i < nx; i++)
	{
		xlow = xup;							// xlow(i) = xup(i-1)
		yup = 0.0;
		plots[i] = new Plot*[ny];
		pads[i]  = new TPad*[ny];

		if (xAxisTitles != NULL && xAxisTitles[i] != NULL)
		{
			buffname = new char[strlen(xAxisTitles[i]) + 1];
			sprintf(buffname, "%s", xAxisTitles[i]);
			xTitles[i] = buffname;
		}
		else
			xTitles[i] = NULL;

		for (int j = 0; j < ny; j++)
		{
			if (i == 0 && yAxisTitles != NULL && yAxisTitles[i] != NULL)
			{
				buffname = new char[strlen(yAxisTitles[j]) + 1];
				sprintf(buffname, "%s", yAxisTitles[j]);
				yTitles[j] = buffname;
			}
			else if (i == 0)
				yTitles[j] = NULL;

			sprintf(buff,"%s_pad_%i_%i",name,i,j);
			xup  = lMargin + (0.98 - lMargin) / nx * (i + 1);
			ylow = yup;						// ylow(j) = yup(j-1)
			yup  = bMargin + (0.98 - bMargin) / ny * (j + 1);
			motherPad->cd();
			pads[i][j] = new TPad(buff,"",xlow,ylow,xup,yup,0,0,0);
			pads[i][j]->AppendPad();

			sprintf(buff,"%s_plot_%i_%i",name,i,j);
			plots[i][j] = new Plot(pads[i][j],nHists,buff,"",
									xTitles[i],yTitles[j]);

			if (j == 0 && i == 0)	// bottom left
			{
				plots[i][j]->SetTextSize(0.086);
				plots[i][j]->SetMargins( relLM, 0.01*(i==nx-1), relBM, 0);
				plots[i][j]->SetXLabelOffset(0.025);
				plots[i][j]->SetXTitleOffset(1.3);
			}
			else if (j == 0)		// just bottom
			{
				plots[i][j]->SetTextSize(0.104);
				plots[i][j]->SetMargins( 0, 0.01*(i==nx-1), relBM, 0);
			}
			else
			{
				plots[i][j]->SetTextSize(0.125);
				plots[i][j]->SetMargins( (i==0)*relLM,	0.01*(i==nx-1),
										 0, 			0.01*(j==ny-1));
				plots[i][j]->SetYTitleOffset(0.85);
			}
		}
	}
}

PanelPlot::~PanelPlot()
{
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			if (i == 0 && yTitles[j] != NULL)	delete yTitles[j];
			delete plots[i][j]; delete pads[i][j];
		}
		delete[] pads[i]; delete[] plots[i];
		if (xTitles[i] != NULL)	delete xTitles[i];
	}
	delete[] pads; delete[] plots; delete name;
	delete[] xTitles; delete[] yTitles;
}

Plot* PanelPlot::GetPlot(int ix, int iy)
{
	if (ix >= 0 && ix < nx && iy >= 0 && iy < ny)
		return plots[ix][iy];
	printf("Plot indices %i,%i are out of range.\n",ix,iy);
	return NULL;
}

void PanelPlot::Draw(bool verb)
{
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
			plots[i][j]->Draw(verb);
	}
}

// ***************************************************************************

ComPlot::ComPlot(TVirtualPad* targPad,int nHist,const char* plotName,
			 const char* plotTitle, const char* xTitle, const char* yTitle)
{
	n    	  = nHist;
	hists     = new TH1*[n];
	for(i = 0; i < n; i++)
		hists[i] = NULL;
	motherPad = targPad;
	i		  = 0;
	iRef	  = 0;

	char* buffname = new char[strlen(plotName) + 1];
	sprintf(buffname, "%s", plotName);
	name = buffname;
	buffname = new char[strlen(plotName) + 10];

	motherPad->cd();
	sprintf(buffname,"%s_pad_D",name);
	dPad = new TPad(buffname,"",0,0.5,1.0,1.0,0,0,0);
	dPad->AppendPad();
	sprintf(buffname,"%s_D",name);
	dPlot = new Plot(dPad,n,buffname,plotTitle,xTitle,yTitle);
	dPlot->SetMargins(-1,-1,0,-1);
	dPlot->SetNormalized();
//	dPlot->SetLegY0(0.6);

	motherPad->cd();
	sprintf(buffname,"%s_pad_R",name);
	rPad = new TPad(buffname,"",0,0.0,1.0,0.5,0,0,0);
	rPad->AppendPad();
	sprintf(buffname,"%s_R",name);
	rPlot = new Plot(rPad,n,buffname,plotTitle,xTitle,NULL);
	rPlot->SetMargins(-1,-1,-1,0);
//	rPlot->SetLegY0(0.8);

	SetTextSize(0.12);
	delete buffname;
}

ComPlot::~ComPlot()
{
	for(i = 0; i < n; i++)
		{ if (hists[i] != NULL) delete hists[i]; }
	delete[] hists; delete[] name;
	delete dPlot; delete rPlot; delete dPad; delete rPad;
}

void ComPlot::Add(const TH1* h,int style,int index,const char* leg)
{
	dPlot->Add(h,style,index,leg);
	char buff[64];
	sprintf(buff, "%s_h%i", name, i);
	hists[i] = (TH1*)h->Clone(buff);
	i++;
}

void ComPlot::Add(const TH1* h, const char* leg, int style, int index)
{
	Add(h, style, index, leg);
}

void ComPlot::Draw(bool verb)
{
	dPlot->Draw(verb);
	char buff[32];
	for (int j = 0; j < n; j++)
	{
		if (j != iRef && hists[j] != NULL && hists[iRef] != NULL)
		{
			hists[j]->Divide(hists[iRef]);
			sprintf(buff,"%i. dist. / %i. dist.",j+1,iRef+1);
			rPlot->Add(hists[j],buff,Plot::Point | Plot::Erry);
		}
		else if (j != iRef)
			printf("ComPlot: hists[%i] or hists[%i] was not initialized.\n",j,iRef);
	}
	rPlot->Draw(verb);
}

// ***************************************************************************
float Hists::FindSmallestAbove(TH1* hist, float min)
{
	float m = hist->GetMaximum();
	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		float val = hist->GetBinContent(i);
		if (val > min && val < m)
			m = val;
	}
	return m;
}

void Hists::SetLogBins(TAxis* axis)
{
	int nbins = axis->GetNbins();
	Axis_t from = axis->GetXmin();
	Axis_t to = axis->GetXmax();
	Axis_t width = (to - from) / nbins;
	Axis_t* newBins = new Axis_t[nbins + 1];

	for (int i = 0; i <= nbins; i++)
		newBins[i] = pow(10, from + i * width);

	axis->Set(nbins, newBins);
	delete[] newBins;
}

void Hists::CheckOverflow(TH1* hist)
{
	float u = hist->GetBinContent(0);
	float o = hist->GetBinContent(hist->GetNbinsX() + 1);

	if (u > 0)
		printf("Histogram '%s' underflow bin: %.0f events\n",hist->GetTitle(),u);
	if (o > 0)
		printf("Histogram '%s' overflow bin: %.0f events\n",hist->GetTitle(),o);
}

void Hists::DrawBarChart(TH1* hist, int n, int index,float font)
{
	float offset = 0.1 / n;
	float width = (1.0 - offset * (n + 1)) / n;
	SetMyTextOptions(hist,font);
	hist->SetFillColor(GetMyFillColor(index));
	hist->SetBarWidth(width);
	hist->SetBarOffset(offset + (width + offset) * index);
	if (index == 0)
		hist->DrawCopy("B");
	else
		hist->DrawCopy("B SAME");
}

void Hists::DrawHBarChart(TH1* hist, int n, int index, float font)
{
	float offset = 0.1 / n;
	float width = (1.0 - offset * (n + 1)) / n;
	SetMyTextOptions(hist,font);
	hist->SetFillColor(GetMyFillColor(index));
	hist->SetBarWidth(width);
	hist->SetBarOffset(offset + (width + offset) * index);
	if (index == 0)
		hist->DrawCopy("HBAR");
	else
		hist->DrawCopy("HBAR SAME");
}

int Hists::GetMyFillColor(int i)	// i = 0..7
{
	int colors[8] = { 	kBlue-8,	kRed-9,		kGreen-6,	kMagenta-9,
						kOrange-3,	kYellow-3,	kAzure-9,	kRed-4};
	if (i >= 0 && i < 8)
		return colors[i];
	else
		return i;
}

int Hists::GetMyLineColor(int i)
{
	int colors[8] = { 	kBlue+3,	kRed+2,		kSpring-5,	kAzure-3,
						kOrange-2,	kGreen+2,	kAzure-9,	kRed-4};
	if (i >= 0 && i < 8)
		return colors[i];
	else
		return i;
}

void Hists::SetMyTextOptions(TH1* hist, float size)
{
	hist->SetNdivisions(505,"xyz");
	hist->SetLabelOffset(0.01,"xyz");
	hist->SetLabelFont(132,"xyz");
	hist->SetTitleFont(132,"xyz");
	hist->GetXaxis()->SetLabelSize(size);
	hist->GetXaxis()->SetTitleSize(size);
	hist->GetYaxis()->SetLabelSize(size);
	hist->GetYaxis()->SetTitleSize(size);
	hist->GetZaxis()->SetLabelSize(size);
	hist->GetZaxis()->SetTitleSize(size);
}

float Hists::Integral(TH1* h,float x0,float x1,bool binWidth)
{
	if (binWidth)
		return h->Integral(h->FindBin(x0),h->FindBin(x1),"width");
	else
		return h->Integral(h->FindBin(x0),h->FindBin(x1));
}

double Hists::Chi2(double *x, double* p)
{
	if (x[0] < 0)
		return 0.0;
	double f = pow(x[0], p[0] / 2.0 - 1.0) * exp(-x[0] / 2.0);
	f /= ( pow(2.0, p[0] / 2.0) * TMath::Gamma(p[0] / 2.0));
	return f;
}