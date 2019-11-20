#include "Style.h"


Style::Style() {
  styler = new TStyle("Styler", "Style");
  read_info("style/default");
  setStyle();
}

Style::Style(string infile) {
  styler = new TStyle("Styler", "Style");
  read_info(infile);
  setStyle();
}

Style::Style(const Style& old) {
  styler = (TStyle*)old.styler->Clone();
  values = old.values;
  values = old.values;
  binlimit = old.binlimit;
  padratio = old.padratio;
  heightratio = old.heightratio;
  rebinlimit = old.rebinlimit;
  dividebins = old.dividebins;

}

Style& Style::operator=(const Style& rhs) {
  if(this == &rhs) return *this;
  
  if(!styler) delete styler;
  styler = (TStyle*)rhs.styler->Clone();
  values = rhs.values;
  binlimit = rhs.binlimit;
  padratio = rhs.padratio;
  heightratio = rhs.heightratio;
  rebinlimit = rhs.rebinlimit;
  dividebins = rhs.dividebins;

  return *this;
}

Style::~Style() {
  delete styler;
}


/////// Long function for all the different styles things that can be done
///// most are unnecessary, but kept just in case.  All values are read in as
/// doubles, so if the value is a bool, need either 0 or 1, not true false
void Style::read_info(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }

  vector<string> stemp;
  string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);

    }
    if(stemp.size() == 2) {
      char* p;
      strtod(stemp[1].c_str(), &p);
      if(! *p) values[stemp[0]]=stod(stemp[1]);
    }
  }
  info_file.close();

}


TStyle* Style::getStyle() {
  return styler;
}


///// actual function that sets the values in the style object
void Style::setStyle() {
  for(map<string,double>::const_iterator it = values.begin(); it != values.end(); ++it) {

    if(it->first == "OptStat") styler->SetOptStat(it->second);
    else if(it->first == "OptTitle") styler->SetOptTitle(it->second);
    else if(it->first == "OptDate") styler->SetOptDate(it->second);
    else if(it->first == "OptFile") styler->SetOptFile(it->second);
    else if(it->first == "OptFit") styler->SetOptFit(it->second);
    else if(it->first == "OptLogx") styler->SetOptLogx(it->second);
    else if(it->first == "OptLogy") styler->SetOptLogy(it->second);
    else if(it->first == "LabelColor") styler->SetLabelColor(it->second, "xy");
    else if(it->first == "LabelOffset") styler->SetLabelOffset(it->second, "xy");
    else if(it->first == "LabelFont") styler->SetLabelFont(it->second, "xy");
    else if(it->first == "LabelSize") styler->SetLabelSize(it->second, "xy");
    else if(it->first == "LegendBorderSize") styler->SetLegendBorderSize(it->second);
    else if(it->first == "LegendFillColor") styler->SetLegendFillColor(it->second);
    else if(it->first == "LegendFont") styler->SetLegendFont(it->second);
    else if(it->first == "LegendTextSize") styler->SetLegendTextSize(it->second);
    else if(it->first == "LegoInnerR") styler->SetLegoInnerR(it->second);
    //    else if(it->first == "LineAttributes") styler->SetLineAttributes(it->second);
    else if(it->first == "LineColor") styler->SetLineColor(it->second);
    //    else if(it->first == "LineColorAlpha") styler->SetLineColorAlpha(it->second);
    else if(it->first == "LineStyle") styler->SetLineStyle(it->second);
    //    else if(it->first == "LineStyleString") styler->SetLineStyleString(it->second);
    else if(it->first == "LineWidth") styler->SetLineWidth(it->second);
    //    else if(it->first == "MarkerAttributes") styler->SetMarkerAttributes(it->second);
    else if(it->first == "MarkerColor") styler->SetMarkerColor(it->second);
    //    else if(it->first == "MarkerColorAlpha") styler->SetMarkerColorAlpha(it->second);
    else if(it->first == "MarkerSize") styler->SetMarkerSize(it->second);
    else if(it->first == "MarkerStyle") styler->SetMarkerStyle(it->second);
    //    else if(it->first == "Name") styler->SetName(it->second);
    //    else if(it->first == "NameTitle") styler->SetNameTitle(it->second);
    else if(it->first == "Ndivisions") styler->SetNdivisions(it->second);
    else if(it->first == "NumberContours") styler->SetNumberContours(it->second);
    //    else if(it->first == "TextAttributes") styler->SetTextAttributes(it->second);
    else if(it->first == "TextColor") styler->SetTextColor(it->second);
    //    else if(it->first == "TextColorAlpha") styler->SetTextColorAlpha(it->second);
    else if(it->first == "TextFont") styler->SetTextFont(it->second);
    else if(it->first == "TextSize") styler->SetTextSize(it->second);
    else if(it->first == "TickLength") styler->SetTickLength(it->second);
    //    else if(it->first == "Title") styler->SetTitle(it->second);
    else if(it->first == "TitleAlign") styler->SetTitleAlign(it->second);
    else if(it->first == "TitleBorderSize") styler->SetTitleBorderSize(it->second);
    else if(it->first == "TitleColor") styler->SetTitleColor(it->second);
    else if(it->first == "TitleFillColor") styler->SetTitleFillColor(it->second);
    else if(it->first == "TitleOffset") styler->SetTitleOffset(it->second);
    else if(it->first == "TitleFont") styler->SetTitleFont(it->second);
    else if(it->first == "TitleFontSize") styler->SetTitleFontSize(it->second);
    else if(it->first == "TitleH") styler->SetTitleH(it->second);
    else if(it->first == "TitleSize") styler->SetTitleSize(it->second);
    else if(it->first == "TitleStyle") styler->SetTitleStyle(it->second);
    else if(it->first == "TitleTextColor") styler->SetTitleTextColor(it->second);
    else if(it->first == "TitleW") styler->SetTitleW(it->second);
    else if(it->first == "TitleX") styler->SetTitleX(it->second);
    else if(it->first == "TitleXOffset") styler->SetTitleXOffset(it->second);
    else if(it->first == "TitleXSize") styler->SetTitleXSize(it->second);
    else if(it->first == "TitleY") styler->SetTitleY(it->second);
    else if(it->first == "TitleYOffset") styler->SetTitleYOffset(it->second);
    else if(it->first == "TitleYSize") styler->SetTitleYSize(it->second);
    else if(it->first == "CanvasColor") styler->SetCanvasColor(it->second);
    else if(it->first == "CanvasBorderMode") styler->SetCanvasBorderMode(it->second);
    else if(it->first == "CanvasBorderSize") styler->SetCanvasBorderSize(it->second);
    else if(it->first == "CanvasDefH") styler->SetCanvasDefH(it->second);
    else if(it->first == "CanvasDefW") styler->SetCanvasDefW(it->second);
    else if(it->first == "CanvasDefX") styler->SetCanvasDefX(it->second);
    else if(it->first == "CanvasDefY") styler->SetCanvasDefY(it->second);
    else if(it->first == "CanvasPreferGL") styler->SetCanvasPreferGL(it->second);
    else if(it->first == "DateX") styler->SetDateX(it->second);
    else if(it->first == "DateY") styler->SetDateY(it->second);
    else if(it->first == "DrawBorder") styler->SetDrawBorder(it->second);
    else if(it->first == "EndErrorSize") styler->SetEndErrorSize(it->second);
    else if(it->first == "ErrorX") styler->SetErrorX(it->second);
    //    else if(it->first == "FillAttributes") styler->SetFillAttributes(it->second);
    else if(it->first == "FillColor") styler->SetFillColor(it->second);
    //    else if(it->first == "FillColorAlpha") styler->SetFillColorAlpha(it->second);
    else if(it->first == "FillStyle") styler->SetFillStyle(it->second);
    //    else if(it->first == "FitFormat") styler->SetFitFormat(it->second);
    else if(it->first == "FrameBorderSize") styler->SetFrameBorderSize(it->second);
    else if(it->first == "FrameFillColor") styler->SetFrameFillColor(it->second);
    else if(it->first == "FrameFillStyle") styler->SetFrameFillStyle(it->second);
    else if(it->first == "FrameLineColor") styler->SetFrameLineColor(it->second);
    else if(it->first == "FrameLineStyle") styler->SetFrameLineStyle(it->second);
    else if(it->first == "FrameLineWidth") styler->SetFrameLineWidth(it->second);
    else if(it->first == "FuncColor") styler->SetFuncColor(it->second);
    else if(it->first == "FuncStyle") styler->SetFuncStyle(it->second);
    else if(it->first == "FuncWidth") styler->SetFuncWidth(it->second);
    else if(it->first == "GridColor") styler->SetGridColor(it->second);
    else if(it->first == "GridStyle") styler->SetGridStyle(it->second);
    else if(it->first == "GridWidth") styler->SetGridWidth(it->second);
    else if(it->first == "HatchesLineWidth") styler->SetHatchesLineWidth(it->second);
    else if(it->first == "HatchesSpacing") styler->SetHatchesSpacing(it->second);
    else if(it->first == "HistFillColor") styler->SetHistFillColor(it->second);
    else if(it->first == "HistFillStyle") styler->SetHistFillStyle(it->second);
    else if(it->first == "HistLineColor") styler->SetHistLineColor(it->second);
    else if(it->first == "HistLineStyle") styler->SetHistLineStyle(it->second);
    else if(it->first == "HistLineWidth") styler->SetHistLineWidth(it->second);
    else if(it->first == "HistMinimumZero") styler->SetHistMinimumZero(it->second);
    else if(it->first == "HistTopMargin") styler->SetHistTopMargin(it->second);
    else if(it->first == "IsReading") styler->SetIsReading(it->second);
    else if(it->first == "PadBorderMode") styler->SetPadBorderMode(it->second);
    else if(it->first == "PadBorderSize") styler->SetPadBorderSize(it->second);
    else if(it->first == "PadBottomMargin") styler->SetPadBottomMargin(it->second);
    else if(it->first == "PadColor") styler->SetPadColor(it->second);
    else if(it->first == "PadBottomMargin") styler->SetPadBottomMargin(it->second);
    else if(it->first == "PadTopMargin") styler->SetPadTopMargin(it->second);
    else if(it->first == "PadLeftMargin") styler->SetPadLeftMargin(it->second);
    else if(it->first == "PadRightMargin") styler->SetPadRightMargin(it->second);
    else if(it->first == "PadGridX") styler->SetPadGridX(it->second);
    else if(it->first == "PadGridY") styler->SetPadGridY(it->second);
    else if(it->first == "PadTickX") styler->SetPadTickX(it->second);
    else if(it->first == "PadTickY") styler->SetPadTickY(it->second);
    else if(it->first == "FrameBorderMode") styler->SetFrameBorderMode(it->second);
    else if(it->first == "AxisColor") styler->SetAxisColor(it->second);
    else if(it->first == "BarOffset") styler->SetBarOffset(it->second);
    //    else if(it->first == "Barwidth") styler->SetBarwidth(it->second);
    else if(it->first == "Bit") styler->SetBit(it->second);
    else if(it->first == "PadRatio") padratio = it->second;
    else if(it->first == "TopWSRatio") heightratio = it->second;
    else if(it->first == "RebinLimit") rebinlimit = it->second;
    else if(it->first == "DivideBins") dividebins = ((int)it->second != 0);
    else if(it->first == "BinLimit") binlimit = it->second;
  }

  gStyle = styler;
}

