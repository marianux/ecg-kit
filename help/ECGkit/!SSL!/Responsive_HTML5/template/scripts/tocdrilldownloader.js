var gTocRootFileName = "toc.js";
var gTocDataFolder = "whxdata";
var gParentDataFile = "parentdata.js";

gRootRelPath = ".";
gCommonRootRelPath = "..";

var gTocPageHtml = "{%LINK_NAME%}";
var gTocBookClosedHtml = "{%LINK_NAME%}";
var gTocBookOpenHtml = "{%LINK_NAME%}";
var gTocUrlHtml = "{%LINK_NAME%}";
var gTocPageIconSrc = "";
var gTocPageIconHoverSrc = "";
var gTocPageIconSelSrc = "";
var gTocPageIconClass = "";
var gTocPageIconStyle = "";
var gTocBookClosedIconSrc = "";
var gTocBookClosedIconHoverSrc = "";
var gTocBookClosedIconSelSrc = "";
var gTocBookClosedIconClass = "";
var gTocBookClosedIconStyle = "";
var gTocBookOpenIconSrc = "";
var gTocBookOpenIconHoverSrc = "";
var gTocBookOpenIconSelSrc = "";
var gTocBookOpenIconClass = "";
var gTocBookOpenIconStyle = "";
var gTocUrlIconSrc = "";
var gTocUrlIconHoverSrc = "";
var gTocUrlIconSelSrc = "";
var gTocUrlIconClass = "";
var gTocUrlIconStyle = "";
var gDDList = null;


LOADINGSTRING = "Loading...";

initTocPage();
function initTocPage()
{
	addRhLoadCompleteEvent(loadToc);
}
function loadToc()
{
	if(gbPreviewMode)
		displayToc(gRootRelPath, gCommonRootRelPath);
	else
		initAndLoadParentData(null, SCR_PARENT_TOC);
}
function displayToc(rootRelPath, commonRootRelPath)
{
	gDDList = new DDList(rootRelPath, gTocDataFolder, gTocRootFileName, commonRootRelPath);
	gDDList.rootId = "toc";
	gDDList.closedBookClass = "wTOCNodeCloseBook";
	gDDList.openBookClass = "wTOCNodeOpenBook";
	gDDList.pageClass = "wTOCNodePage";
	gDDList.urlClass = "wTOCNodeLink";
	gDDList.closedBookClassHover = "wTOCNodeCloseBookHover";
	gDDList.openBookClassHover = "wTOCNodeOpenBookHover";
	gDDList.pageClassHover = "wTOCNodePageHover";
	gDDList.urlClassHover = "wTOCNodeLinkHover";
	gDDList.closedBookClassClick = "";
	gDDList.openBookClassClick = "";
	gDDList.pageClassClick = "";
	gDDList.urlClassClick = "";
	gDDList.closedBookClassSelected = "wTOCNodeCloseBookSelected";
	gDDList.openBookClassSelected = "wTOCNodeOpenBookSelected";
	gDDList.pageClassSelected = "wTOCNodePageSelected";
	gDDList.urlClassSelected = "wTOCNodeLinkSelected";
	gDDList.pageHtml = gTocPageHtml;
	gDDList.bookClosedHtml = gTocBookClosedHtml;
	gDDList.bookOpenHtml = gTocBookOpenHtml;
	gDDList.urlHtml = gTocUrlHtml;
	gDDList.upBtnClass = "wTOCNodeUp";
	gDDList.upBtnClassHover = "wTOCNodeUpHover";
	gDDList.upBtnClassClick = "";

	gDDList.iconClass[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconClass;
	gDDList.iconClass[ITEMTYPEBOOKOPEN] = gTocBookOpenIconClass;
	gDDList.iconClass[ITEMTYPEPAGE] = gTocPageIconClass;
	gDDList.iconClass[ITEMTYPEURL] = gTocUrlIconClass;
	
	gDDList.iconStyle[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconStyle;
	gDDList.iconStyle[ITEMTYPEBOOKOPEN] = gTocBookOpenIconStyle;
	gDDList.iconStyle[ITEMTYPEPAGE] = gTocPageIconStyle;
	gDDList.iconStyle[ITEMTYPEURL] = gTocUrlIconStyle;
	
	gDDList.iconSrc[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconSrc;
	gDDList.iconSrc[ITEMTYPEBOOKOPEN] = gTocBookOpenIconSrc;
	gDDList.iconSrc[ITEMTYPEPAGE] = gTocPageIconSrc;
	gDDList.iconSrc[ITEMTYPEURL] = gTocUrlIconSrc;
	
	gDDList.iconHoverSrc[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconHoverSrc;
	gDDList.iconHoverSrc[ITEMTYPEBOOKOPEN] = gTocBookOpenIconHoverSrc;
	gDDList.iconHoverSrc[ITEMTYPEPAGE] = gTocPageIconHoverSrc;
	gDDList.iconHoverSrc[ITEMTYPEURL] = gTocUrlIconHoverSrc;
	
	gDDList.iconSelSrc[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconSelSrc;
	gDDList.iconSelSrc[ITEMTYPEBOOKOPEN] = gTocBookOpenIconSelSrc;
	gDDList.iconSelSrc[ITEMTYPEPAGE] = gTocPageIconSelSrc;
	gDDList.iconSelSrc[ITEMTYPEURL] = gTocUrlIconSelSrc;

	gDDList.upButtonId = "upBtn";
	gDDList.errorMsg = "Unknown error";
	gDDList.setLoadingDisplayInfo("loadingicon", "<img src='" + gRootRelPath + "/template/resources/LoadingData.gif' />", "loadingtext", LOADINGSTRING);
	gDDList.init();
	gDDList.load();
}