//#ifndef DRAW_STYLE
//#define DRAW_STYLE
//#define DRAW_STYLE

using namespace std;
#include "TStyle.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"
#include "TLine.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "iostream"
#include "TF1.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TROOT.h"


//#define titlesize 0.06
//#define labelsize 0.05
//#define titleoffsetx 0.7
//#define titleoffsety 1.2
//#define font 42


void  drawStyle()
{
    // root -l
    // .L drawStyle.C++

    // gSystem->Load("/misc/galatyuk/tools/drawStyle_C.so");
    // drawStyle();

    cout<<" -- SET DRAW OPTIONS LOADED --"<<endl;

    // defining canvas options ------------------------
    gStyle->SetCanvasColor(10);
    gStyle->SetPadColor(10);

    // defining size of pad ---------------------------
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadBorderMode(0);

    // switch options of the histogram ----------------
//    gStyle->SetOptStat(0);         // switch off stat box
//    gStyle->SetErrorX(0);          // switch off error in x
    gStyle->SetTitleStyle(0);      // make title box transparent
    gStyle->SetTitleBorderSize(0); // make border of box disapear
}



void setCanvas(TCanvas *canvas)
{
    gStyle->SetLineStyleString(22,"80 18 12 18 12 12");
    canvas->SetFillColor(10);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetFrameLineWidth(2);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderSize(0);
    canvas->SetLeftMargin(0.1);
    canvas->SetRightMargin(0.2);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.1);
    canvas->ToggleEventStatus();
//    canvas->Range(-194.483,-10.3682,1041.38,-2.08469);
}
void setOPT_pad(Bool_t kLog=false, Bool_t kgridX=false, Bool_t kgridY=false)
{
    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    gPad->SetBottomMargin(0.18);
    gPad->SetFillColor(10);
    gPad->SetGridx(kgridX);
    gPad->SetGridy(kgridY);
    gPad->SetLogz(kLog);
    gPad->Update();                            // hello ...wake up...something changed
}

void setOPT_padY(Bool_t kLog=false, Bool_t kgridX=true, Bool_t kgridY=true)
{
    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    gPad->SetBottomMargin(0.18);
    gPad->SetFillColor(10);
    gPad->SetGridx(kgridX);
    gPad->SetGridy(kgridY);
    gPad->SetLogy(kLog);
    gPad->Update();                            // hello ...wake up...something changed
}

void setOPT_pad_logY(Bool_t kLog=false)
{
    gPad->SetTopMargin(0.1);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.2);
    gPad->SetFillColor(10);
    gPad->SetLogy(kLog);
    gPad->Update();                            // hello ...wake up...something changed
}

// this function is setting the style of the plotting historgrams
// as an argument one should give name of the histogram, title of the X and Y axis
// Ndevision for the X axis, style, size and the color of the marker
// e.g for the official PRL plot for the RAW spectra of the NOV02 it was used in the following way:
// setOPT_hists(signal,20,1.3);
// setOPT_hists(background,24,1.3,2);
// setOPT_hists(unlikesign,22,1.3,4);
void setTH1(TH1 *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=505, Int_t marker_style =20, Float_t marker_size =0.8, Int_t color=1)
{
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
    hist->GetYaxis()->SetNdivisions(Ndevision);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetTitleOffset(1.1);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    //hist->GetXaxis()->SetRangeUser(0., 650.);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(1);
    hist->Sumw2();
//    hist->SetMinimum(5e-5);
//    hist->SetMaximum(1e2);
}
void setTProfile(TProfile *prof, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=505, Int_t marker_style =20, Float_t marker_size =0.8, Int_t color=1)
{
    prof->GetXaxis()->SetTitle(xAxisTitle);
    prof->GetYaxis()->SetTitle(yAxisTitle);
    prof->GetXaxis()->SetTitleSize(0.06);
    prof->GetYaxis()->SetTitleSize(0.06);
    prof->GetXaxis()->SetTitleFont(42);
    prof->GetYaxis()->SetTitleFont(42);
    prof->GetXaxis()->SetNdivisions(Ndevision);
    prof->GetYaxis()->SetNdivisions(Ndevision);
    prof->GetYaxis()->SetTitleOffset(1.4);
    prof->GetXaxis()->SetTitleOffset(1.1);
//    prof->GetXaxis()->CenterTitle();
//    prof->GetYaxis()->CenterTitle();
    prof->GetXaxis()->SetLabelFont(42);
    prof->GetYaxis()->SetLabelFont(42);
    prof->GetXaxis()->SetLabelSize(0.05);
    prof->GetYaxis()->SetLabelSize(0.05);
    //prof->GetXaxis()->SetRangeUser(0., 650.);
    prof->SetMarkerStyle(marker_style);
    prof->SetMarkerSize(marker_size);
    prof->SetMarkerColor(color);
    prof->SetLineColor(color);
    prof->SetLineWidth(2);
//    prof->SetMinimum(5e-5);
//    prof->SetMaximum(1e2);
}

void setTGaxis(TGaxis *ax, TString AxisTitle=" ", Int_t Ndevision=510)
{
    ax->SetTitle(AxisTitle);
    ax->SetTitleSize(0.05);
    ax->SetTitleFont(42);
    ax->SetNdivisions(Ndevision);
    ax->SetLabelFont(42);
    ax->SetLabelSize(0.05);
    ax->SetTitleOffset(1.6);

//    gStyle->SetPadLeftMargin(0.17);
//    gStyle->SetPadRightMargin(0.13);

}


void setTH2(TH2 *hist2, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510)
{
    //gStyle->SetPadRightMargin(0.15);
    hist2->GetXaxis()->SetTitle(xAxisTitle);
    hist2->GetYaxis()->SetTitle(yAxisTitle);
    hist2->GetXaxis()->SetTitleSize(0.08);
    hist2->GetYaxis()->SetTitleSize(0.08);
//    hist2->GetXaxis()->SetTitleSize(0.06);
//    hist2->GetYaxis()->SetTitleSize(0.06);
    hist2->GetZaxis()->SetTitleSize(0.06);
    hist2->GetXaxis()->SetTitleFont(42);
    hist2->GetYaxis()->SetTitleFont(42);
    hist2->GetZaxis()->SetTitleFont(42);
    hist2->GetXaxis()->SetNdivisions(Ndevision);
    hist2->GetYaxis()->SetNdivisions(Ndevision);
    hist2->GetXaxis()->SetLabelFont(42);
    hist2->GetYaxis()->SetLabelFont(42);
    hist2->GetZaxis()->SetLabelFont(42);
    hist2->GetXaxis()->SetLabelSize(0.06);
    hist2->GetYaxis()->SetLabelSize(0.06);
//    hist2->GetXaxis()->SetLabelSize(0.04);
//    hist2->GetYaxis()->SetLabelSize(0.04);
    hist2->GetZaxis()->SetLabelSize(0.04);
    hist2->GetXaxis()->SetLabelOffset(0.015);
    hist2->GetYaxis()->SetLabelOffset(0.015);
    hist2->GetZaxis()->SetLabelOffset(0.015);
    hist2->GetXaxis()->SetTitleOffset(1.);
    hist2->GetZaxis()->SetTitleOffset(1.);
    hist2->GetYaxis()->SetTitleOffset(1.);

//    gStyle->SetPadLeftMargin(0.17);
//    gStyle->SetPadRightMargin(0.13);

}

void setTF1(TF1 *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1, Int_t lineWidth=2)
{
    gStyle->SetPadRightMargin(0.05);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(lineWidth);
}



void setGraph(TGraph *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1, Int_t lineWidth=2)
{
//////////////////
    //gStyle->SetPadRightMargin(0.05);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(lineWidth);
}

void setGraphE(TGraphErrors *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1, Int_t lineWidth=2)
{
//////////////////
    //gStyle->SetPadRightMargin(0.05);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(lineWidth);
}



void setText(TPaveText *p, Int_t color=1)
{
    p->SetFillColor(10);
    p->SetTextColor(color);
    p->SetBorderSize(-1);
    p->SetTextFont(42);
}



void setLabel(TPaveLabel *p)
{
    p->SetFillColor(10);
    p->SetBorderSize(-1);
    p->SetTextFont(42);
}

void setLegend(TLegend *p)
{
    p->SetFillColor(10);
    p->SetBorderSize(-1);
    p->SetTextFont(42);
}



TLegend* plotLegend(TString pos="right_top",TString Title="No Title",Float_t scaleX=0.9,Float_t scaleY=0.9,Float_t offsetX=0.0,Float_t offsetY=0.0,TString Comment="",Int_t commentcolor=1)
{
    // pos places the Legend and can be
    // right_top, mid_top, left_top, right_bottom,mid_bottom,left_bottom
    // Title gives you the Title of the legend
    // scaleX,scaleY defines the size (=1 means 1/3 pad with,1/3 pad height (according to pos))
    // offsetX , offsetY (NDC coordinates of pad -> 0.1) shift by 10% of pad with) shifts the legend by the value in x and y
    // comment is optional text below title before legend entris will start
    // commentcolor defines the text color of the comment

    Float_t left  =gPad->GetLeftMargin()*1.15;
    Float_t right =1-gPad->GetRightMargin();
    Float_t top   =1-gPad->GetTopMargin();
    Float_t bottom=gPad->GetBottomMargin()*1.15;
    Float_t mid   =gPad->GetLeftMargin() + (1-(gPad->GetRightMargin()+gPad->GetLeftMargin()))/2;
    Float_t width =(right-left)/2;
    Float_t heith =(top-bottom)/2;
    TLegend* legend;
    TLine* dummy=new TLine();
    dummy->SetLineColor(10);

    if(pos.CompareTo("left_top")==0)    legend=new TLegend(left+offsetX,top+offsetY-(scaleY*heith),left+offsetX+(scaleX*width),top+offsetY,Title.Data());                           // left top
    if(pos.CompareTo("mid_top")==0)     legend=new TLegend(mid+offsetX-((scaleX*width)/2),top+offsetY-(scaleY*heith),mid+offsetX+((scaleX*width)/2),top+offsetY,Title.Data());      // mid up
    if(pos.CompareTo("right_top")==0)   legend=new TLegend(right+offsetX-(scaleX*width),top+offsetY-(scaleY*heith),right+offsetX,top+offsetY,Title.Data());                         // right top
    if(pos.CompareTo("left_bottom")==0) legend=new TLegend(left+offsetX,bottom+offsetY,left+offsetX+(scaleX*width),bottom+offsetY+(scaleY*heith),Title.Data());                     // left bottom
    if(pos.CompareTo("mid_bottom")==0)  legend=new TLegend(mid+offsetX-((scaleX*width)/2),bottom+offsetY,mid+offsetX+((scaleX*width)/2),bottom+offsetY+(scaleY*heith),Title.Data());// mid down
    if(pos.CompareTo("right_bottom")==0)legend=new TLegend(right+offsetX-(scaleX*width),bottom+offsetY,right+offsetX,bottom+offsetY+(scaleY*heith),Title.Data());                   // right bottom

    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetMargin(0.15);



    if(Comment.CompareTo("")!=0)
    {   // comment below header
        TLegendEntry* entry=legend->AddEntry(dummy,Comment.Data());
        entry->SetTextColor(commentcolor);
        entry->SetTextFont(62);
        entry->SetOption("l");
    }
    return legend;
}
void setLegendEntry(TLegend* legend,TObject* object,char* label,Int_t col,TString opt,Float_t size=0.05)
{
    // add entry for object with label and color and option to legend
    // legend is thge pointer to the legend where you want to add the entry
    // object ist the histogram or graph which will be added
    // label is the text which will appear in the entry
    // col is th color of the text
    // opt is the option ("pL"(marker and line), "p"(marker),"l"(line))
    TLegendEntry* entry=legend->AddEntry((TObject*)object,label);
    entry->SetTextColor(col);
    entry->SetTextSize(size);
    entry->SetOption(opt.Data());
}
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.05,Int_t color=1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    if(x<0||y<0)
    {   // defaults
        x=gPad->GetLeftMargin()*1.15;
        y=(1-gPad->GetTopMargin())*1.02;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextSize(size);
    text->SetNDC();
    text->SetTextColor(color);
    text->Draw();
    return text;
}
/*
void saveToPdf(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
{
    // saves a canvas to pdf file by storing it as ps and calling
    // ps2pdf to convert it to pdf. The ps file will be removed.
    // canvas is the pointer of the canvas which you want to save
    // path is the path to the directory where you want to store the pdf
    // filename is the name of the pdf which wou create
    // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
    // separate=kFALSE as the canvas in one pdf.
    TString File    =filename;
    TString Path    =path;
    Path            =Path + "/";
    TString command ="";
    TString totalin ="";
    TString totalout="";

    File.ReplaceAll(".ps","");
    File.ReplaceAll(".pdf","");

    if(!separate)
    {
        totalin = Path + File + ".ps";
        totalout= Path + File + ".pdf";
        canvas->SaveAs(totalin.Data());
        command= "ps2pdf " + totalin + " " + totalout;
        gSystem->Exec(command.Data());
        command= "rm " + totalin;
        gSystem->Exec(command.Data());
    }else{
        TString subpad="";
        Int_t npads= canvas->GetListOfPrimitives()->GetSize();
        for(Int_t i=1;i<=npads;i++)
        {
            subpad=i;
            totalin = Path + File + "_" + subpad + ".ps";
            totalout= Path + File + "_" + subpad + ".pdf";
            canvas->cd(i);
            gPad->SaveAs(totalin.Data());
            command= "ps2pdf " + totalin + " " + totalout;
            gSystem->Exec(command.Data());
            command= "rm " + totalin;
            gSystem->Exec(command.Data());
        }
    }
}
void saveToPng(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
{
    // saves a canvas to png file by storing it as gif and calling
    // convert to convert it to png. The gif file will be removed.
    // canvas is the pointer of the canvas which you want to save
    // path is the path to the directory where you want to store the png
    // filename is the name of the png which wou create
    // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
    // separate=kFALSE as the canvas in one png.
    TString File    =filename;
    TString Path    =path;
    Path            =Path + "/";
    TString command ="";
    TString totalin ="";
    TString totalout="";


    File.ReplaceAll(".ps","");
    if(!separate)
    {
        totalin = Path + File + ".ps";
        totalout= Path + File + ".png";
        canvas->SaveAs(totalin.Data());
        command= "convert " + totalin + " " + totalout;
        gSystem->Exec(command.Data());
        command= "rm " + totalin;
        gSystem->Exec(command.Data());
    }else{
        TString subpad="";
        Int_t npads= canvas->GetListOfPrimitives()->GetSize();
        for(Int_t i=1;i<=npads;i++)
        {
            subpad=i;
            totalin = Path + File + "_" + subpad + ".ps";
            totalout= Path + File + "_" + subpad + ".png";
            canvas->cd(i);
            gPad->SaveAs(totalin.Data());
            command= "convert " + totalin + " " + totalout;
            gSystem->Exec(command.Data());
            command= "rm " + totalin;
            gSystem->Exec(command.Data());
        }
    }
}

void saveToC(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
{
    // saves a canvas to png file by storing it as gif and calling
    // convert to convert it to png. The gif file will be removed.
    // canvas is the pointer of the canvas which you want to save
    // path is the path to the directory where you want to store the png
    // filename is the name of the png which wou create
    // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
    // separate=kFALSE as the canvas in one png.
    TString File    =filename;
    TString Path    =path;
    Path            =Path + "/";
    TString command ="";
    TString totalin ="";
    TString totalout="";

    File.ReplaceAll(".C","");
    if(!separate)
    {
        totalin = Path + File + ".C";
        canvas->SaveAs(totalin.Data());
    }else{
        TString subpad="";
        Int_t npads= canvas->GetListOfPrimitives()->GetSize();
        for(Int_t i=1;i<=npads;i++)
        {
            subpad=i;
            totalin = Path + File + "_" + subpad + ".C";
            canvas->cd(i);
            gPad->SaveAs(totalin.Data());
        }
    }
}

void saveGifToPdf(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
{
    // saves a canvas to pdf file by storing it as gif and calling
    // convert to convert it to pdf. The gif file will be removed.
    // canvas is the pointer of the canvas which you want to save
    // path is the path to the directory where you want to store the pdf
    // filename is the name of the pdf which wou create
    // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
    // separate=kFALSE as the canvas in one pdf.
    TString File    =filename;
    TString Path    =path;
    Path            =Path + "/";
    TString command ="";
    TString totalin ="";
    TString totalout="";

    File.ReplaceAll(".pdf","");
    if(!separate)
    {
        totalin = Path + File + ".gif";
        totalout= Path + File + ".pdf";
        canvas->SaveAs(totalin.Data());
        command= "convert " + totalin + " " + totalout;
        gSystem->Exec(command.Data());
        command= "rm " + totalin;
        gSystem->Exec(command.Data());
    }else{
        TString subpad="";
        Int_t npads= canvas->GetListOfPrimitives()->GetSize();
        for(Int_t i=1;i<=npads;i++)
        {
            subpad=i;
            totalin = Path + File + "_" + subpad + ".gif";
            totalout= Path + File + "_" + subpad + ".pdf";
            canvas->cd(i);
            gPad->SaveAs(totalin.Data());
            command= "convert " + totalin + " " + totalout;
            gSystem->Exec(command.Data());
            command= "rm " + totalin;
            gSystem->Exec(command.Data());
        }
    }
    }
    */
Int_t getColor(Int_t i)
{
    // pre defined set of 10 colors
    // i is larger 9 it will start from the first color again
    Int_t colors[10]   ={2,4,6,8,38,46,14,1,30,43};
    return colors[i%10];
}
Int_t getMarker(Int_t i)
{
    // pre defined set of 10 marker styles
    // i is larger 9 it will start from the first marker style again
    Int_t style[10]   ={20,21,22,23,24,25,26,27,28,29};
    return style[i%10];
}
void setGraph(TGraph* graph,Int_t i,Int_t markerstyle=20,Int_t markercolor=1,Float_t markersize=1.,Int_t linecolor=1)
{
    // sets the graphics of graph (use full in a loop)
    // graph is the pointer of the graph
    //
    // i switches to modes:
    // (if i>0) i the number of the color of the marker and the markert style (getColor(),getMarker())
    // (if i<0) marker color and style will be  markercolor, markerstyle
    //
    // markersize, and linecolor defines the size of the marker and the line color
    // can be use in two ways:
    // values >=0 the style and color will be set via getColor(i),getMarker(i)
    // values < 0 the style and color will be set -markercolor -markerstyle

    if(i<0)
    {
        graph->SetMarkerColor(markercolor);
        graph->SetMarkerSize(markersize);
        graph->SetMarkerStyle(markerstyle);
        graph->SetLineColor(linecolor);
    }else{
        if     (markercolor>=0)graph->SetMarkerColor(getColor(i));
        else if(markercolor<0) graph->SetMarkerColor(getColor(-markercolor));
        graph->SetMarkerSize(markersize);
        if     (markerstyle>=0)graph->SetMarkerStyle(getMarker(i));
        else if(markerstyle<0) graph->SetMarkerStyle(getMarker(-markerstyle));
        graph->SetLineColor(linecolor);
    }
}
void setHist(TH1* hist,Int_t i,Int_t markerstyle=20,Int_t markercolor=1,Float_t markersize=1.,Int_t linecolor=1)
{
    // sets the graphics of histogram (use full in a loop)
    // hist is the pointer of the hiistogram
    //
    // i switches to modes:
    // (if i>0) i the number of the color of the marker and the markert style (getColor(),getMarker())
    // (if i<0) marker color and style will be  markercolor, markerstyle
    //
    // markersize, and linecolor defines the size of the marker and the line color
    // can be use in two ways:
    // values >=0 the style and color will be set via getColor(i),getMarker(i)
    // values < 0 the style and color will be set -markercolor -markerstyle
    if(i==-99)
    {
        hist->SetMarkerColor(markercolor);
        hist->SetMarkerSize(markersize);
        hist->SetMarkerStyle(markerstyle);
        hist->SetLineColor(linecolor);
    }else{
        if     (markercolor>=0)hist->SetMarkerColor(getColor(i));
        else if(markercolor<0) hist->SetMarkerColor(getMarker(-markercolor));
        hist->SetMarkerSize(markersize);
        if     (markerstyle>=0)hist->SetMarkerStyle(getMarker(i));
        else if(markerstyle<0) hist->SetMarkerStyle(getMarker(-markerstyle));
        hist->SetLineColor(linecolor);
      }
}

TPaveText *setOPT_text(TString cut_desc ,  Double_t x1, Double_t y1, Double_t x2, Double_t y2,Int_t color=1, Int_t textSize=0.2)
{
  TPaveText *description=new TPaveText(x1,y1,x2,y2,"NDC");
  //TPaveText *description=new TPaveText(x1,y1,x2,y2,"blNDC");
 
  description->SetLineWidth(0);
  description->AddText(cut_desc);
//   description->SetTextSize(0.046);
  description->SetTextSize(textSize);
  description->SetBorderSize(0.0);
  description->SetTextFont(42);
  description->SetTextColor(color);
  description->SetFillColor(0);
  description->SetFillStyle(0);
  description->Draw();
  return description;
}

TPaveText *setOPT_text2(TString cut_desc ,  Double_t x1, Double_t y1, Double_t x2, Double_t y2,Int_t color=1, Int_t textSize=0.2)
{
  TPaveText *description=new TPaveText(x1,y1,x2,y2,"NDC");
  //TPaveText *description=new TPaveText(x1,y1,x2,y2,"blNDC");
 
  description->SetLineWidth(0);
  description->AddText(cut_desc);
//   description->SetTextSize(0.046);
  description->SetTextSize(textSize);
  description->SetBorderSize(0.0);
  description->SetTextFont(42);
  description->SetTextColor(color);
  description->SetFillColor(0);
  description->SetFillStyle(0);
  description->SetTextAlign(13);
  description->Draw();
  return description;
}
void setOPT_padfull(Double_t left=0, Double_t right=0, Double_t top=0, Double_t bot=0)
{
    gPad->SetTopMargin(top);
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetBottomMargin(bot);
    gPad->SetFillColor(10);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
//    gPad->SetLogy(kLog);
    gPad->Update();                            // hello ...wake up...something changed
}
