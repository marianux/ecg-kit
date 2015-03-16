var gGloRootFileName = "glo.js";
var gGloDataFolder = "whxdata";

var gGloTermHtml = "{%LINK_NAME%}";
var gGloDefHtml = "{%LINK_NAME%}";
var gGloCategoryHtml = "{%LINK_NAME%}";


LOADINGSTRING = "Loading...";

initGloPage();
function initGloPage()
{
	addRhLoadCompleteEvent(loadGlo);
}
var gGloList = null;
function loadGlo()
{
	if(gbPreviewMode)
	{
		var gloRootPathsArr = new Array;
		gloRootPathsArr[0] = gRootRelPath;
		displayGlo(gloRootPathsArr);
	}
	else
		initAndLoadParentData(null, SCR_PARENT_GLO);
}

function displayGlo(gloRootPathsArr)
{
	gGloList = new GloList(gloRootPathsArr, gGloDataFolder, gGloRootFileName);
	gGloList.rootId = "glo";
	gGloList.termClass = "wGloTerm";
	gGloList.termStyle = "";
	gGloList.termClassHover = "wGloTermHover";
	gGloList.termClassClick = "";
	gGloList.defClass = "wGloDefinition";
	gGloList.defStyle = "";
	gGloList.defClassHover = "wGloDefinitionHover";
	gGloList.defClassClick = "";
	gGloList.categoryClass = "wGloAlphabet";
	gGloList.categoryStyle = "";
	gGloList.termHtml = gGloTermHtml;
	gGloList.defHtml = gGloDefHtml;
	gGloList.categoryHtml = gGloCategoryHtml;


	gGloList.filterBoxId = "gloFilterBox";
	gGloList.errorMsg = "Unknown error";
	gGloList.setLoadingDisplayInfo("loadingicon", "<img src='" + gRootRelPath + "/template/resources/LoadingData.gif' alt='Loading' />", "loadingtext", LOADINGSTRING);
	gGloList.init();
	gGloList.load();
}