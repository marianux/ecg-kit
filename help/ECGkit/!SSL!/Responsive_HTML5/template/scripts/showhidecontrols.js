//This list of tags can be extended to add more tags for Show/Hide support
var gShowHideTagsList = ["div", "a"];

addRhLoadCompleteEvent(doShowHide);

function doShowHide(evt, bIgnoreUrlParam)
{
	var showHideMode;
	var cshModeFlag = "";
	if(bIgnoreUrlParam == true)
		cshModeFlag = "";
	else
		cshModeFlag = getUrlParameter(RHCSHMODE);
	if(cshModeFlag == TRUESTR)
	{
		showHideMode = CSHMODE;
		saveSetting(RHCSHMODE, CSHMODE);
		callbackDoShowHideCSHModeCookie(showHideMode);
	}
	else if(cshModeFlag == FALSESTR)
	{
		showHideMode = NONCSHMODE;
		saveSetting(RHCSHMODE, NONCSHMODE);
		callbackDoShowHideCSHModeCookie(showHideMode);
	}
	else
		readSetting(RHCSHMODE, callbackDoShowHideCSHModeCookie);
}

function callbackDoShowHideCSHModeCookie(showHideMode)
{		
	if(showHideMode == "")
		showHideMode = NONCSHMODE;
	
	publishCSHMode(showHideMode);

	var len = gShowHideTagsList.length;
	for(var i=0; i<len; i++)
	{
		var tagName = gShowHideTagsList[i];
		var elemsList = document.getElementsByTagName(tagName);
		showHideElems(elemsList, showHideMode);
	}
}

function publishCSHMode(showHideMode) {
	if (showHideMode === CSHMODE)
		rh.model.publish(rh.consts('KEY_CSH_MODE'), true);
	else
		rh.model.publish(rh.consts('KEY_CSH_MODE'), false);
}

function showHideElems(elemsList, showHideMode)
{
	var len = elemsList.length;
	for(var i=0; i<len; i++)
	{
		var elem = elemsList[i];
		var showinAttrib = elem.getAttribute(DATASHOWIN);
		if(showinAttrib == SHOWINCSHMODE)
		{
			if(showHideMode == CSHMODE)
				elem.style.display = "";
			else
				elem.style.display = "none";
		}
		else if(showinAttrib == SHOWINNONCSHMODE)
		{
			if(showHideMode == CSHMODE)
				elem.style.display = "none";
			else
				elem.style.display = "";
		}
	}

}
function onShowHideClick()
{
	readSetting(RHCSHMODE, callbackShowHideClickCSHModeCookie);
}

function callbackShowHideClickCSHModeCookie(showHideMode)
{
	if(showHideMode == CSHMODE)
		showHideMode = NONCSHMODE;
	else if(showHideMode == "" || showHideMode == NONCSHMODE)
		showHideMode = CSHMODE;

	saveSetting(RHCSHMODE, showHideMode);
	doShowHide(null, true);
}