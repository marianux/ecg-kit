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

var gTree = null;
var gParentDataPathStack = new MhStack();

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
	gTree = new Tree(rootRelPath, gTocDataFolder, gTocRootFileName, commonRootRelPath);
	gTree.rootId = "toc";
	gTree.closedBookClass = "wTOCTreeCloseBook";
	gTree.openBookClass = "wTOCTreeOpenBook";
	gTree.pageClass = "wTOCTreePage";
	gTree.urlClass = "wTOCTreeLink";
	gTree.closedBookClassHover = "wTOCTreeCloseBookHover";
	gTree.openBookClassHover = "wTOCTreeOpenBookHover";
	gTree.pageClassHover = "wTOCTreePageHover";
	gTree.urlClassHover = "wTOCTreeLinkHover";
	gTree.closedBookClassClick = "";
	gTree.openBookClassClick = "";
	gTree.pageClassClick = "";
	gTree.urlClassClick = "";
	gTree.closedBookClassSelected = "wTOCTreeCloseBookSelected";
	gTree.openBookClassSelected = "wTOCTreeOpenBookSelected";
	gTree.pageClassSelected = "wTOCTreePageSelected";
	gTree.urlClassSelected = "wTOCTreeLinkSelected";
	gTree.pageHtml = gTocPageHtml;
	gTree.bookClosedHtml = gTocBookClosedHtml;
	gTree.bookOpenHtml = gTocBookOpenHtml;
	gTree.urlHtml = gTocUrlHtml;

	gTree.bookChildsClass = "wTOCTreeOpenBookChildBlock";
	
	gTree.iconClass[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconClass;
	gTree.iconClass[ITEMTYPEBOOKOPEN] = gTocBookOpenIconClass;
	gTree.iconClass[ITEMTYPEPAGE] = gTocPageIconClass;
	gTree.iconClass[ITEMTYPEURL] = gTocUrlIconClass;
	
	gTree.iconStyle[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconStyle;
	gTree.iconStyle[ITEMTYPEBOOKOPEN] = gTocBookOpenIconStyle;
	gTree.iconStyle[ITEMTYPEPAGE] = gTocPageIconStyle;
	gTree.iconStyle[ITEMTYPEURL] = gTocUrlIconStyle;
	
	gTree.iconSrc[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconSrc;
	gTree.iconSrc[ITEMTYPEBOOKOPEN] = gTocBookOpenIconSrc;
	gTree.iconSrc[ITEMTYPEPAGE] = gTocPageIconSrc;
	gTree.iconSrc[ITEMTYPEURL] = gTocUrlIconSrc;
	
	gTree.iconHoverSrc[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconHoverSrc;
	gTree.iconHoverSrc[ITEMTYPEBOOKOPEN] = gTocBookOpenIconHoverSrc;
	gTree.iconHoverSrc[ITEMTYPEPAGE] = gTocPageIconHoverSrc;
	gTree.iconHoverSrc[ITEMTYPEURL] = gTocUrlIconHoverSrc;
	
	gTree.iconSelSrc[ITEMTYPEBOOKCLOSED] = gTocBookClosedIconSelSrc;
	gTree.iconSelSrc[ITEMTYPEBOOKOPEN] = gTocBookOpenIconSelSrc;
	gTree.iconSelSrc[ITEMTYPEPAGE] = gTocPageIconSelSrc;
	gTree.iconSelSrc[ITEMTYPEURL] = gTocUrlIconSelSrc;

	gTree.errorMsg = "Unknown error";
	gTree.setLoadingDisplayInfo("loadingicon", "<img src='" + gRootRelPath + "/template/resources/LoadingData.gif' alt='Loading' />", "loadingtext", LOADINGSTRING);
	gTree.init();
	gTree.load();
}