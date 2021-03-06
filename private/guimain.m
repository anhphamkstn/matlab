function fig = guimain()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load guimain

h0 = figure('Units','points', ...
	'CloseRequestFcn','BNBGUICB(''quit main''); closereq;', ...
	'Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','C:\MATLABR11\toolbox\bnb\guimain.m', ...
	'HandleVisibility','callback', ...
	'MenuBar','none', ...
	'Name','main BNB GUI', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576.0000000000001 432.0000000000002], ...
	'PaperUnits','points', ...
	'Position',mat1, ...
	'Resize','off', ...
	'Tag','main BNB GUI', ...
	'ToolBar','none');
h1 = uimenu('Parent',h0, ...
	'Label','file', ...
	'Tag','file');
h2 = uimenu('Parent',h1, ...
	'Callback','bnbguicb(''save'')', ...
	'Label','save', ...
	'Tag','save');
h2 = uimenu('Parent',h1, ...
	'Callback','bnbguicb(''load'')', ...
	'Label','load', ...
	'Tag','load');
h1 = uimenu('Parent',h0, ...
	'Label','edit', ...
	'Tag','edit');
h2 = uimenu('Parent',h1, ...
	'Callback','bnbguicb(''X -> x0'')', ...
	'Label','X -> x0', ...
	'Tag','X -> x0');
h2 = uimenu('Parent',h1, ...
	'Callback','bnbguicb(''Z X t c fail -> base'')', ...
	'Label','Z X t c fail -> base', ...
	'Tag','Z X t c fail -> base');
h1 = uimenu('Parent',h0, ...
	'Label','help', ...
	'Tag','help');
h2 = uimenu('Parent',h1, ...
	'Callback','bnbguicb(''help'')', ...
	'Label','help', ...
	'Tag','help');
h2 = uimenu('Parent',h1, ...
	'Callback','bnbguicb(''copyright'')', ...
	'Label','copyright', ...
	'Tag','copyright');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'ListboxTop',0, ...
	'Position',[6.206896551724139 142.7586206896552 279.3103448275863 93.10344827586209], ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'ListboxTop',0, ...
	'Position',[148.9655172413793 68.27586206896554 136.5517241379311 68.27586206896554], ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'ListboxTop',0, ...
	'Position',[6.206896551724139 6.206896551724139 136.5517241379311 130.3448275862069], ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'FontName','Courier-LD', ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[12.41379310344828 148.9655172413793 254.4827586206897 80.68965517241381], ...
	'Style','text', ...
	'Tag','results', ...
	'UserData',mat2);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback',mat3, ...
	'FontName','Courier-LD', ...
	'ListboxTop',0, ...
	'Position',[155.1724137931035 37.24137931034483 62.0689655172414 24.82758620689656], ...
	'String','function', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback','bnbguicb(''settings'');', ...
	'FontName','Courier-LD', ...
	'ListboxTop',0, ...
	'Position',[155.1724137931035 6.206896551724139 62.0689655172414 24.82758620689656], ...
	'String','settings', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'BusyAction','cancel', ...
	'Callback','bnbguicb(''optimize'');', ...
	'FontName','Courier-LD', ...
	'ListboxTop',0, ...
	'Position',[223.448275862069 37.24137931034483 62.0689655172414 24.82758620689656], ...
	'String','optimize', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback',mat4, ...
	'FontName','Courier-LD', ...
	'ListboxTop',0, ...
	'Position',[223.448275862069 6.206896551724139 62.0689655172414 24.82758620689656], ...
	'String','quit', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','bnbguicb(''xlist'');', ...
	'FontName','Courier-LD', ...
	'Position',mat5, ...
	'Style','listbox', ...
	'Tag','xlist', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','bnbguicb(''x0'');', ...
	'FontName','Courier-LD', ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[74.48275862068967 74.48275862068967 62.0689655172414 14.27586206896552], ...
	'Style','edit', ...
	'Tag','x0');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','bnbguicb(''xlb'');', ...
	'FontName','Courier-LD', ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[74.48275862068967 55.86206896551725 62.0689655172414 14.27586206896552], ...
	'Style','edit', ...
	'Tag','xlb');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat6, ...
	'FontName','Courier-LD', ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[74.48275862068967 93.10344827586209 62.0689655172414 14.27586206896552], ...
	'Style','edit', ...
	'Tag','xub');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[55.86206896551725 93.10344827586209 18.62068965517242 12.41379310344828], ...
	'String','xub', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[55.86206896551725 55.86206896551725 18.62068965517242 12.41379310344828], ...
	'String','xlb', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[55.86206896551725 74.48275862068967 12.41379310344828 12.41379310344828], ...
	'String','x0', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'Callback','bnbguicb(''continuous'');', ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[55.86206896551725 37.24137931034483 62.0689655172414 12.41379310344828], ...
	'String','continuous', ...
	'Style','radiobutton', ...
	'Tag','continuous', ...
	'UserData','[ ]');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'Callback','bnbguicb(''integer'');', ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[55.86206896551725 24.82758620689656 55.86206896551725 12.41379310344828], ...
	'String','integer', ...
	'Style','radiobutton', ...
	'Tag','integer');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'Callback','bnbguicb(''fixed'');', ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[55.86206896551725 12.41379310344828 55.86206896551725 12.41379310344828], ...
	'String','fixed', ...
	'Style','radiobutton', ...
	'Tag','fixed');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','bnbguicb(''parlist'');', ...
	'FontName','Courier-LD', ...
	'Position',[155.1724137931035 74.48275862068967 37.24137931034483 43.44827586206898], ...
	'String',' ', ...
	'Style','listbox', ...
	'Tag','parlist', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[198.6206896551724 93.10344827586209 18.62068965517242 12.41379310344828], ...
	'String','val', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','bnbguicb(''par'');', ...
	'FontName','Courier-LD', ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[217.2413793103449 93.10344827586209 62.0689655172414 14.27586206896552], ...
	'Style','edit', ...
	'Tag','par');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback','bnbguicb(''resultsslider'');', ...
	'ListboxTop',0, ...
	'Position',[266.896551724138 148.9655172413793 12.41379310344828 80.68965517241381], ...
	'Style','slider', ...
	'Tag','resultsslider', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',mat7, ...
	'String','variable', ...
	'Style','text', ...
	'Tag','StaticText3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',mat8, ...
	'String','parameter', ...
	'Style','text', ...
	'Tag','StaticText3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'Callback',mat9, ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',mat10, ...
	'Style','edit', ...
	'Tag','xtag');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725], ...
	'Callback',mat11, ...
	'FontName','Courier-LD', ...
	'ForegroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',mat12, ...
	'Style','edit', ...
	'Tag','partag');
if nargout > 0, fig = h0; end
