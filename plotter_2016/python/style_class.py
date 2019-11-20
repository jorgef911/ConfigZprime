
import matplotlib
import sys

class style_container():
    ##------------------------------------------------------------------
    ## Public functions
    ##------------------------------------------------------------------
    def __init__(self, style = 'Plain', kind = 'Standard', useRoot = False, cmsPositon = "upper right", legendPosition = "upper right", content = 'Histogram', lumi = 42000, cms = 13):
        self._style = style

        if not (kind == 'Standard' or kind == 'Lines' or kind == 'Graphs' or kind == 'Linegraphs'):
            print(('\n\tThis kind of plot (' + kind + ') is not supported'))
            print('\tThe allowed values are:')
            print('\t  - \'Standard\'   (Stacked backgrounds, signal lines and data points)')
            print('\t  - \'Lines\'      (Only lines, like for a gen level plot)')
            print('\t  - \'Graphs\'     (Only graphs, like for an efficiency plot)\n')
            print('\t  - \'Linegraphs\' (Only graphs, like for an limit plot)\n')
            sys.exit(42)
        self._kind = kind

        self._useRoot = useRoot

        if not (content == 'Histogram' or content == 'Efficiencies'):
            print(('\n\tThis content for a plot (' + content + ') is not supported'))
            print('\tThe allowed values are:')
            print('\t  - \'Histogram\'    (Normal histograms)')
            print('\t  - \'Efficiencies\' (Efficiency histograms)\n')
            sys.exit(42)
        self._content = content

        #rc('text', usetex=True)
        self._additional_text  = 'Preliminary'
        self._y_label_offset   = -0.11
        self._histaxis_y_label_offset   = 1.07
        self._error_bands_ecol = ['darkmagenta','darkcyan']
        self._error_bands_fcol = ['m','cyan']
        self._error_bands_alph = 0.7
        self._error_bands_labl = ['Sys. uncert. 1','Sys. uncert. 2']
        self._error_bands_center = 'ref'
        self._error_stacking = 'No'
        self._line_width    =0.7
        self._spine_line_width = 0.5
        self._xerr = False
        self._logx = False
        self._logy = True
        self._logz = False
        self._ymin = -1
        self._ymax = -1
        self._xmin = -1
        self._xmax = -1
        self._zmin = -1
        self._zmax = -1
        self._histaxis_logy = False
        self._histaxis_ymin = -1
        self._histaxis_ymax = -1
        self._lumi_val = lumi
        self._cms_val = cms
        self._add_lumi_text = True
        self._no_legend=False
        self._legend_ncol= 1
        self._grid = False
        self._batch_mode=True
        self._legend_font_size = 0
        self._tick_font_size = 0
        self._axis_title_font_size = 0
        self._forceBinWidth=False
        self._poisson_error=True
        self._poisson_error_max_event=False
        self._do_blinding=False
        self._blind_min=False
        self._blind_max=False
        self._do_data_errors = True
        self._show_minor_tick_labels = False
        self._do_minor_ticks         = False
        self._do_overflowbin = False

        self._cmsTextPosition = position(cmsPositon, isText = True, useRoot=self._useRoot)
        self._LegendPosition = position(legendPosition, useRoot=self._useRoot)

    def __del__(self):
        pass

    def AddAxisTitle(self, hist):
        try:
            self._xaxis_title      = hist.xaxis.GetTitle()
            self._yaxis_title      = hist.yaxis.GetTitle()
        except:
            print("(Setting axis title) Unexpected error:", sys.exc_info()[0], sys.exc_info()[1])
            self._xaxis_title      = 'bla'
            self._yaxis_title      = '#epsilon'

    def AddAxisTitle_histaxis(self, hist):
        try:
            self._histaxis_xaxis_title      = hist.xaxis.GetTitle()
            self._histaxis_yaxis_title      = hist.yaxis.GetTitle()
        except:
            print("(Setting axis title) Unexpected error:", sys.exc_info()[0], sys.exc_info()[1])
            self._histaxis_xaxis_title      = 'bla'
            self._histaxis_yaxis_title      = '#epsilon'

    def InitStyle(self, addplots = ['', '', ''], addheights = [0, 0, 0], histaxis = []):
        self.addplots = addplots
        self.addheights = addheights
        self.histaxis = histaxis

        if self._useRoot:
            self._Set_Root_style()
        else:
            matplotlib.rcParams.update({'font.size': 10})
            matplotlib.rcParams.update({'lines.linewidth' : 1})
        if self._style == 'CMS':
            self._Set_CMS_style()
        elif self._style == 'Plain':
            self._Set_Plain_style()
        elif self._style == 'Cool':
            self._Set_Cool_style()
        else:
            print('\n\tNo style chosen, this can\'t be right!\n\tDo something!\n')
            sys.exit(42)

    def Get_useRoot(self):
        return self._useRoot

    def Get_xaxis_title(self):
        try:
            #use the correct math font family (at least for one case)
            if self._axis_title_font["family"]=='sans-serif':
                return self._xaxis_title.replace("mathrm","mathsf")
            else:
                return self._xaxis_title
        except(AttributeError):
            print('\n\t It seems that no x-axis title is set, we will therefore use a dummy value\n')
            return 'dummy'

    def Get_yaxis_title(self):
        try:
            return self._yaxis_title
        except(AttributeError):
            print('\n\t It seems that no y-axis title is set, we will therefore use a dummy value\n')
            return 'dummy'

    def Get_histaxis_xaxis_title(self):
        return self._histaxis_xaxis_title

    def Get_histaxis_yaxis_title(self):
        return self._histaxis_yaxis_title

    def Get_additional_text(self):
        return self._additional_text

    def Get_y_label_offset(self):
        return self._y_label_offset

    def Get_histaxis_y_label_offset(self):
        return self._histaxis_y_label_offset

    def Get_error_bands_ecol(self):
        return self._error_bands_ecol

    def Get_error_bands_fcol(self):
        return self._error_bands_fcol

    def Get_error_bands_alph(self):
        return self._error_bands_alph

    def Get_error_bands_labl(self):
        return self._error_bands_labl

    def Get_error_bands_center(self):
        return self._error_bands_center

    def Get_error_stacking(self):
        return self._error_stacking

    def Get_error_line_width(self):
        return self._line_width

    def Get_spine_line_width(self):
        return self._spine_line_width

    def Get_logx(self):
        return self._logx

    def Get_logy(self):
        return self._logy

    def Get_logz(self):
        return self._logz

    def Get_ymin(self):
        return self._ymin

    def Get_ymax(self):
        return self._ymax

    def Get_xmin(self):
        return self._xmin

    def Get_xmax(self):
        return self._xmax

    def Get_xerr(self):
        return self._xerr

    def Get_poisson_error(self):
        return self._poisson_error

    def Get_poisson_error_max_event(self):
        return self._poisson_error_max_event

    def Get_do_blinding(self):
        return self._do_blinding

    def Get_blinding_region(self):
        return self._blind_min, self._blind_max

    def Get_do_data_errors(self):
        return self._do_data_errors

    def Get_histaxis_logy(self):
        return self._histaxis_logy

    def Get_histaxis_ymin(self):
        return self._histaxis_ymin

    def Get_histaxis_ymax(self):
        return self._histaxis_ymax

    def Get_add_cms_text(self):
        return self._add_cms_text

    def Get_add_lumi_text(self):
        return self._add_lumi_text

    def Get_label_text_color(self):
        return self._label_text_color

    def Get_histaxis_label_text_color(self):
        return self._histaxis_label_text_color

    def Get_annotation_text_color(self):
        return self._annotation_text_color

    def Get_bg_color(self):
        return self._bg_color

    def Get_batch_mode(self):
        return self._batch_mode

    def Get_ref_line_color(self):
        return self._ref_line_color

    def Get_spine_color(self):
        return self._spine_color

    def Get_tick_color(self):
        return self._tick_color

    def Gets_tick_font_size(self):
        return self._tick_font_size

    def Get_marker_style(self):
        return self._marker_style

    def Get_line_style(self):
        return self._line_style

    def Get_marker_size(self):
        return self._marker_size

    def Get_marker_color(self):
        return self._marker_color

    def Get_marker_error_cap_width(self):
        return self._marker_error_cap_width

    def Get_cms_text_alignment(self):
        return self._cms_text_alignment

    def Get_show_minor_tick_labels(self):
        return self._show_minor_tick_labels

    def Get_legend_font_size(self):
        return self._legend_font_size

    def Get_axis_title_font_size(self, value):
        return self._axis_title_font_size

    def Get_cmsTextPosition(self):
        return self._cmsTextPosition

    def Get_LegendPosition(self):
        return self._LegendPosition

    def Get_no_legend(self):
        return self._no_legend

    def Get_n_legend_columns(self):
        return self._legend_ncol

    def Get_do_minor_ticks(self):
        return self._do_minor_ticks

    def Get_kind(self):
        return self._kind

    def Get_grid(self):
        return self._grid

    def Get_grid_style(self):
        return self._grid_style

    def Get_grid_width(self):
        return self._grid_width

    def Get_grid_color(self):
        return self._grid_color

    def Get_content(self):
        return self._content

    def Get_lumi_val(self):
        return self._lumi_val

    def Get_cms_val(self):
        return self._cms_val

    def Get_zmin(self):
        return self._zmin

    def Get_zmax(self):
        return self._zmax

    def Get_do_overflowbin(self):
        return self._do_overflowbin

    def Set_line_style(self, value):
        self._line_style = value

    def Set_error_bands_ecol(self, value):
        self._error_bands_ecol = value

    def Set_error_bands_fcol(self, value):
        self._error_bands_fcol = value

    def Get_axis_title_font(self):
        return self._axis_title_font

    def Get_axis_text_main_to_sub_ratio(self):
        return self._axis_text_main_to_sub_ratio

    def Set_y_label_offset(self, value):
        self._y_label_offset = value

    def Set_add_lumi_text(self, value):
        self._add_lumi_text = value

    def Set_lumi_val(self, lumi):
        self._lumi_val = lumi

    def Set_batch_mode(self,value):
        self._batch_mode=value

    def Set_legend_font_size(self, value):
        self._legend_font_size = value

    def Set_axis_title_font_size(self, value):
        self._axis_title_font_size = value

    def Set_tick_font_size(self, value):
        self._tick_font_size = value

    #to have the same syntax as root
    def SetBatchMode(self,value):
        self.Set_batch_mode(value)

    def Set_error_bands_labl(self, label):
        self._error_bands_labl = label

    def Set_error_bands_center(self, center):
        self._error_bands_center = center

    def Set_error_stacking(self, stacking):
        self._error_stacking = stacking

    def Set_error_line_width(self,width):
        self._line_width=width

    def Set_xerr(self):
        self._xerr=True

    def Set_no_legend(self):
        self._no_legend=True

    def Set_additional_text(self, text):
        self._additional_text = text

    def Set_n_legend_columns(self,n):
        self._legend_ncol=n

    def Set_poisson_error_off(self):
        self._poisson_error = False

    def Set_poisson_error_max_event(self):
        self._poisson_error_max_event = True

    def Set_do_blinding_and_region(self,xmin=0,xmax=10):
        self._do_blinding = True
        self._blind_min = xmin
        self._blind_max = xmax

    def Set_do_blinding(self):
        self._do_blinding = True
        self._blind_min = xmin
        self._blind_max = xmax

    def Set_blinding_region(self,xmin=-10000,xmax=10000):
        self._blind_min = xmin
        self._blind_max = xmax

    def Set_do_data_errors(self,do_it=True):
        self._do_data_errors = do_it

    def Set_minor_ticks(self):
        self._do_minor_ticks = True

    def Set_do_overflowbin(self):
        self._do_overflowbin = True

    ## Function to set properties of the plotting axis
    #
    # This function sets axis properties like the y-range or
    # if any axis should be logarithmic.
    # @param[in] logx Boolean if the x-axis should be logarithmic (Default = False)
    # @param[in] logy Boolean if the y-axis should be logarithmic (Default = True)
    # @param[in] ymin Minimum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] ymax Maximum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] logy Boolean if the second / additional y-axis should be logarithmic (Default = False)
    # @param[in] ymin Minimum plotting range for the second / additional y-axis (Default = -1 automatic values)
    # @param[in] ymax Maximum plotting range for the second / additional y-axis (Default = -1 automatic values)
    # @param[in] xmin Minimum plotting range for the x-axis (Default = -1 range from hist)
    # @param[in] xmax Maximum plotting range for the x-axis (Default = -1 range from hist)
    def Set_axis(self, logx = False, logy = True, logz = False, ymin = -1, ymax = -1, xmin = -1, xmax = -1, zmin = -1, zmax = -1, histaxis_logy = False, histaxis_ymin = -1, histaxis_ymax = -1, grid = False):
        self._logx = logx
        self._logy = logy
        self._logz = logz
        self._ymin = ymin
        self._ymax = ymax
        self._zmin = zmin
        self._zmax = zmax
        self._histaxis_logy = histaxis_logy
        self._histaxis_ymin = histaxis_ymin
        self._histaxis_ymax = histaxis_ymax
        self._xmin = xmin
        self._xmax = xmax
        self._grid = grid
    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    def _Set_CMS_style(self):
        self._add_cms_text           = True
        self._label_text_color       = 'black'
        self._annotation_text_color  = 'black'
        self._bg_color               = 'w'
        self._ref_line_color         = 'blue'
        self._spine_color            = 'black'
        self._grid_style             = '-.'
        self._grid_width             = 0.2
        self._grid_color             = 'black'
        self._tick_color             = 'black'
        self._marker_style           = 'o'
        self._line_style             = '-'
        self._marker_size            = 3
        self._marker_color           = 'black'
        self._marker_error_cap_width = 0
        self._cms_text_alignment     = 'row'
        if self._legend_font_size == 0:
            self._legend_font_size       = 10
        if self._axis_title_font_size == 0:
            self._axis_title_font_size       = 14
        self._axis_title_font =  {'family' : 'sans-serif',
                                  'color'  : 'black',
                                  'weight' : 'medium',
                                  'size'   : self._axis_title_font_size,
                                 }
        if self._tick_font_size == 0:
            self._tick_font_size             = 10
        self._axis_text_main_to_sub_ratio = 0.8
        if len(self.histaxis) > 0:
            self._histaxis_label_text_color= self.histaxis[0].linecolor
        else:
            self._histaxis_label_text_color='red'

    def _Set_Plain_style(self):
        self._add_cms_text                    = False
        self._label_text_color                = 'black'
        self._annotation_text_color           = 'black'
        self._bg_color                        = 'w'
        self._ref_line_color                  = 'blue'
        self._spine_color                     = 'black'
        self._grid_style                      = '-'
        self._grid_width                      = 0.2
        self._grid_color                      = 'black'
        self._tick_color                      = 'black'
        self._marker_style                    = 'o'
        self._line_style                      = '-'
        self._marker_size                     = 4
        self._marker_color                    = 'black'
        self._marker_error_cap_width          = 1
        self._cms_text_alignment              = 'row'
        if self._legend_font_size == 0:
            self._legend_font_size                = 10
        if self._axis_title_font_size == 0:
            self._axis_title_font_size       = 9
        self._axis_title_font =  {'family' : 'sans-serif',
                                  'color'  : 'black',
                                  'weight' : 'medium',
                                  'size'   : self._axis_title_font_size,
                                 }
        if self._tick_font_size == 0:
            self._tick_font_size             = 8
        self._axis_text_main_to_sub_ratio = 0.8
        if len(self.histaxis) > 0:
            self._histaxis_label_text_color= self.histaxis[0].linecolor
        else:
            self._histaxis_label_text_color='red'

    def _Set_Cool_style(self):
        self._add_cms_text           = True
        self._label_text_color       = 'white'
        self._annotation_text_color  = 'white'
        self._bg_color               = '#07000d'
        self._ref_line_color         = 'y'
        self._spine_color            = '#5998ff'
        self._grid_style             = '-'
        self._grid_width             = 0.2
        self._grid_color             = '#5998ff'
        self._tick_color             = 'w'
        self._marker_style           = 'o'
        self._line_style             = '-'
        self._marker_size            = 3
        self._marker_color           = 'lightgray'
        self._marker_error_cap_width = 0
        self._cms_text_alignment     = 'column'
        if self._legend_font_size == 0:
            self._legend_font_size       = 9
        if self._axis_title_font_size == 0:
            self._axis_title_font_size       = 12
        self._axis_title_font =  {'family' : 'serif',
                                  'color'  : 'white',
                                  'weight' : 'medium',
                                  'size'   : self._axis_title_font_size,
                                 }
        if self._tick_font_size == 0:
            self._tick_font_size             = 10
        self._axis_text_main_to_sub_ratio = 0.8
        if len(self.histaxis) > 0:
            self._histaxis_label_text_color= self.histaxis[0].linecolor
        else:
            self._histaxis_label_text_color='red'

    def _Set_Root_style(self):
        self.cmsTextFont          = 63   # Fonts
        self.lumiTextFont         = 43
        self.extraTextFont        = 53
        self.additionalTextFont   = 43
        self.cmsTextSize          = 22  #Text sizes
        self.lumiTextSize         = 20
        self.extraTextSize        = 20
        self.additionalTextSize   = 20
        self.legendTextSize       = 14
        self.lumiTextOffset       = 0.2
        self.extraTextOffset      = 2.5  # only used in outOfFrame version
        self.axisLabelTextSize    = 18
        self.axisTitleTextSize    = 20
        self.axisOffsetY          = 1.5
        self.axisOffsetX          = 1.5

class position():
    def __init__(self,positiontext="upper right", refference="", isText=False ,useRoot=False):

        self._positiontext=positiontext
        if not isinstance(positiontext,str):
            self._definedCoorinates=True
            self._valign="left"
            self._align="left"
        else:
            self._definedCoorinates=False
            self._valign=self._positiontext.split(" ")[0]
            self._align=self._positiontext.split(" ")[1]
        self.addY=0
        self.addX=0
        self._isText=isText

        self._correctcms={"left":0.,
                    "middle":0.,
                    "right":-0.15,
                    "outside":0.0,
                    "upper":-0.04,
                    "center":0.,
                    "lower":0.,
        }
        if useRoot:
            self._correctcms={"left":0.,
                    "middle":0.,
                    "right": -0.025,
                    "outside":0.0,
                    "upper": -0.05,
                    "center":0.,
                    "lower":0.08,
        }
        if self._definedCoorinates:
            self._x=self._positiontext[0]
            self._y=self._positiontext[1]

    def __eq__(self,other):
        return (self._positiontext==other._positiontext)

    def addYspace(self,y):
        if self._valign=="upper" and y<0.:
            self.addY+=y
        elif self._valign!="upper" and self.getY()+y>0.1:
            self.addY+=y

    def addXspace(self,x):
        self.addX+=x

    def setPosition(self,positiontext):
        self._positiontext=positiontext
        self.valign=self._positiontext.split(" ")[0]
        self.align=self._positiontext.split(" ")[1]

    def getText(self):
        return self._positiontext

    def getX(self):
        if self._definedCoorinates:
            return self._x
        alignDict={
                    "left":0.12,
                    "middle":0.5,
                    "right":0.95,
        }
        if self._isText:
            return self.addX+alignDict[self._align]+self._correctcms[self._align]
        return self.addX+alignDict[self._align]

    def get_positiontext(self):
        return self._positiontext.replace('middle','center')

    def getY(self):
        if self._definedCoorinates:
            return self._y
        alignDict={
                    "outside":0.96,
                    "upper":0.95,
                    "center":0.5,
                    "lower":0.12,
        }
        if self._isText:
            return self.addY+alignDict[self._valign]+self._correctcms[self._valign]
        return self.addY+alignDict[self._valign]
