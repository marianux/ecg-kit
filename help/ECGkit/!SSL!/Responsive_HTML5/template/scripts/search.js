var gSearchMsgId = "searchMsg";
var gResultsFoundString = "%1 result(s) found for %2";
var gSearchResultHtml = "{%LINK_NAME%}\n							\n							<span  class=\"wSearchURL\">{%SEARCH_URL%}</span> \n							 <br />\n							<span  class=\"wSearchContext\">{%SEARCH_SUMMARY%}</span>";
var gSearchResClassName = "wSearchResultItem";
var gSearchResTitleClassName = "wSearchResultTitle";
var gSearchResTitleClassHover = "wSearchResultTitleHover";
var gSearchResStyle = "";
var gSearchResTitleStyle = "";
var gSearchPrevBtnId = "searchprev";
var gSearchNextBtnId = "searchnext";
var gsResultDivID="searchResList";
var gPageListBarID ="pageList";
var gPageLinkClass = "wSearchPageNumberLink";
var gPageClass = "wSearchPageNumberSelected";
var gSearchDropdownID = "searchResCount";
var gSearchPageFilePath = "";
var gSearchResultsCount = "15";
var gSearchHighlightControlID = "highlightsearch";
var gbHighLight = 1;
var gTextHighlightColor = "#000000";
var gbgHighlightColor = "#b2b4bf";

gRootRelPath = ".";

var gPageRange = 0;
function initSearchCountDropDown()
{
	readSetting(RHSEARCHCOUNT, callbackSearchCountCookieRead);
}
function initHighlightSearchControl()
{
	readSetting(RHHIGHLIGHT, callbackHighlightCookieRead);
}
function callbackSearchCountCookieRead(maxValCookie)
{
	var val;
	if(maxValCookie != "")
		val = maxValCookie;
	else 
		val = gSearchResultsCount;
	var dropdown = getElement(gSearchDropdownID);
	if(dropdown)
		dropdown.value = val;
	g_nMaxResult = val;
}
function callbackHighlightCookieRead(highlightFlag)
{
	if(highlightFlag == TRUESTR)
		gbHighLight = 1;

	else if(highlightFlag == FALSESTR)
		gbHighLight = 0;

	var highlightElem = document.getElementById(gSearchHighlightControlID);
	if(highlightElem)
	{
		if(gbHighLight)
		{
			highlightElem.checked = true;
			saveSetting(RHHIGHLIGHT, TRUESTR, true);
		}
		else
		{
			highlightElem.checked = false;
			saveSetting(RHHIGHLIGHT, FALSESTR, true);
		}
		saveSetting(RHHIGHLIGHTTEXTCOLOR, gTextHighlightColor, true);
		saveSetting(RHHIGHLIGHTBGCOLOR, gbgHighlightColor, true);
	}
}
function onToggleHighlightSearch()
{
	if(gbHighLight)
	{
		gbHighLight = 0;
		saveSetting(RHHIGHLIGHT, FALSESTR, true);
	}
	else
	{
		gbHighLight = 1;
		saveSetting(RHHIGHLIGHT, TRUESTR, true);
	}
}
function onMaxPageCountChange(maxVal)
{
	g_nMaxResult = maxVal;
	
	if(rh.model.get(rh.consts('KEY_SEARCHED_TERM')))
		onClickPage(null, 1);
	saveSetting(RHSEARCHCOUNT, maxVal, true);
}
function onClickPrevNext( btn, a_nPageNumber )
{
	onClickPage(a_nPageNumber);	
}
function updateNavigationPagesBar(nCurPage, nNumPages)
{
	var pageListBarDiv = document.getElementById(gPageListBarID);
	if(pageListBarDiv == null || pageListBarDiv == 'undefined'){
		return;
	}
	if(nNumPages <= 1) {
		pageListBarDiv.innerHTML = '';
		return;
	}	
		
	var resDiv = document.getElementById(gsResultDivID);
	if(gPageRange == 0)
		gPageRange = Math.floor(resDiv.offsetWidth/SEARCHPAGEWIDTHRATIO);
	var startPage = nCurPage - Math.floor(gPageRange/2);
	var endPage = 0;
	if(startPage < 1)
		startPage = 1;
	endPage = startPage + gPageRange -1;
	if(endPage > nNumPages)
	{
		endPage = nNumPages;
		startPage = endPage - gPageRange + 1;
		if(startPage < 1)
			startPage = 1;
	}
	var sHTML = "";
	sHTML += "<ul style='margin: 0px; padding: 0px;'>";
	for(var i=startPage; i<=endPage; i++)
	{
		if(i == nCurPage)
			sHTML += "<li class='" + gPageClass + "' style='display:inline;'>" + i.toString() + "</li>";
		else
			sHTML += "<li class='" + gPageLinkClass + " " + HLISTCLASS + " " + HANDCURSORCLASS + "' onclick=\"onClickPrevNext(this,'" + i.toString() + "')\" >" + i.toString() + "</li>";
	}
	sHTML += "</ul>";
	pageListBarDiv.innerHTML = sHTML;
}
function updatePrevNextButtons(nCurPage, nNumPages)
{
	var prevBtn = document.getElementById(gSearchPrevBtnId);
	var nextBtn = document.getElementById(gSearchNextBtnId);
	var isPrevBtn = false;
	var isNextBtn = false;
	if(prevBtn != null && prevBtn != 'undefined')
		isPrevBtn = true;
	if(nextBtn != null && nextBtn != 'undefined')
		isNextBtn = true;
	if (nNumPages > 1)
	{
		if(nCurPage > 1)
		{
			if(isPrevBtn)
			{
				prevBtn.style.display = "inline";
				prevBtn.onclick = function(){onClickPrevNext(prevBtn, (parseInt(nCurPage)-1).toString());};
			}
		}
		else if(isPrevBtn)
				prevBtn.style.display = "none";
		if(nCurPage < nNumPages)
		{
			if(isNextBtn)
			{
				nextBtn.style.display = "inline";
				nextBtn.onclick = function(){onClickPrevNext(nextBtn, (parseInt(nCurPage)+1).toString());};
			}
		}
		else if(isNextBtn)
				nextBtn.style.display = "none";
	}
	else
	{
		if(isPrevBtn)
			prevBtn.style.display = "none";
		if(isNextBtn)
			nextBtn.style.display = "none";
	}
}
function initSearchPage()
{
	initSearchCountDropDown();
	initHighlightSearchControl();
	updatePrevNextButtons(0,0);
}

function writeResult( a_strUrl, a_strTitle, a_nIndex, a_sSummary )
{
	var strTitleStyle = "";
	if(gSearchResTitleStyle != "")
		strTitleStyle = "style=\"" + gSearchResTitleStyle + "\" ";
	var strHoverEvents = "";
	if(isTouchDevice())
	{
		strHoverEvents += " ontouchstart=\"onSearchItemHover(this,'" + gSearchResTitleClassHover + "')\" ";
		strHoverEvents += " ontouchend=\"onSearchItemHoverOut(this,'" + gSearchResTitleClassName + "')\"";
		strHoverEvents += " ontouchmove=\"onSearchItemHoverOut(this,'" + gSearchResTitleClassName + "')\"";
	}
	else
	{
		strHoverEvents += " onmouseover=\"onSearchItemHover(this,'" + gSearchResTitleClassHover + "')\" ";
		strHoverEvents += " onmouseout=\"onSearchItemHoverOut(this,'" + gSearchResTitleClassName + "')\"";
	}
	var anchorStartTag = "<a class='"+ NOLINKANCHORCLASS + "' href=\"" + a_strUrl + "\" >";
	var divStartTag = "<div class='" + gSearchResTitleClassName + "' " + strTitleStyle + strHoverEvents + ">";
	var title = anchorStartTag + divStartTag + _textToHtml_nonbsp(a_strTitle) + "</div></a>";
	
	var html = gSearchResultHtml.replace(LINK_NAME_MACRO, title);
	if(a_sSummary.length > 0)
	{
		var summary = _textToHtml_nonbsp(lTrim(a_sSummary));
		html = html.replace(SEARCH_SUMMARY_MACRO, summary);
	}
	else
		html = html.replace(SEARCH_SUMMARY_MACRO, "");
	if(a_strUrl.length > 2 && a_strUrl.charAt(0) == '.' && a_strUrl.charAt(1) == '/')
		a_strUrl = a_strUrl.substring(2, a_strUrl.length);
	var pos = a_strUrl.indexOf("?");
	if(pos != -1)
		a_strUrl = a_strUrl.substring(0, pos);
	html = html.replace(SEARCH_URL_MACRO, _htmlToText(a_strUrl));
	var strStyle ="";
	if(gSearchResStyle != "")
		strStyle = "style=\"" + gSearchResStyle + "\" ";
	var	strOutput = "<div class='" + gSearchResClassName + "' " + strStyle + " >" + html + "</div>";
	return strOutput;
}
function setResultsStringHTML(results_no, searchStr)
{
	var msg = gResultsFoundString;
	msg = msg.replace("%1", results_no);
	msg = msg.replace("%2", "\'" + searchStr + "\'");
	displayMsg(msg);
}
function displayMsg(msg)
{
	var spanNode = document.getElementById(gSearchMsgId);
	if(spanNode != null && spanNode != 'undefined')
		spanNode.innerHTML = msg;
}
function onSearchItemHover(node, className)
{
	if(className != "")
		node.className = className;
}
function onSearchItemHoverOut(node, className)
{
	if(className != "")
		node.className = className;
}