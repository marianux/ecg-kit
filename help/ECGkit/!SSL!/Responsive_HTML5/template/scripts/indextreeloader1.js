var gIdxRootFileName = "idx.js";
var gIdxDataFolder = "whxdata";

var gIdxKWHtml = "{%LINK_NAME%}";
var gIdxLinkHtml = "{%LINK_NAME%}";
var gIdxCategoryHtml = "{%LINK_NAME%}";

LOADINGSTRING = "Loading...";

initIdxPage();
function initIdxPage()
{
	addRhLoadCompleteEvent(loadIdx);
}
var gIdxTree = null;
function loadIdx()
{
	if(gbPreviewMode)
	{
		var idxRootPathsArr = new Array;
		idxRootPathsArr[0] = gRootRelPath;
		displayIdx(idxRootPathsArr);
	}
	else
		initAndLoadParentData(null, SCR_PARENT_IDX);
}

function displayIdx(idxRootPathsArr)
{
	gIdxTree = new IdxTree(idxRootPathsArr, gIdxDataFolder, gIdxRootFileName);
	gIdxTree.rootId = "idx";
	gIdxTree.kWClass = "wIdxKeyword";
	gIdxTree.kWStyle = "";
	gIdxTree.kWClassHover = "wIdxKeywordHover";
	gIdxTree.kWClassClick = "";
	gIdxTree.linkClass = "wIdxLink";
	gIdxTree.linkStyle = "";
	gIdxTree.linkClassHover = "wIdxLinkHover";
	gIdxTree.linkClassClick = "";
	gIdxTree.categoryClass = "wIdxAlphabet";
	gIdxTree.categoryStyle = "";
	gIdxTree.kWHtml = gIdxKWHtml;
	gIdxTree.linkHtml = gIdxLinkHtml;
	gIdxTree.categoryHtml = gIdxCategoryHtml;

	gIdxTree.bookChildsClass = "wIdxChildBlock";
	gIdxTree.filterBoxId = "idxFilterBox";
	gIdxTree.errorMsg = "Unknown error";
	gIdxTree.setLoadingDisplayInfo("loadingicon", "<img src='" + gRootRelPath + "/template/resources/LoadingData.gif' alt='Loading' />", "loadingtext", LOADINGSTRING);
	gIdxTree.init();
	gIdxTree.load();
}