gDataFolder = "whxdata";
gCSHDataFile = "csh.js";
gWndDataFile = "window.js";

var gBrowserWnd = null;
var gTopicURL = null;
var gDefTopicURL = null;

var gCSHWnd = new Object();
gCSHWnd.sName=RHMSWINNAME;
gCSHWnd.nBOptions=WINRESIZABLE|WINSCROLLBARS;
gCSHWnd.sBLeft="49%";
gCSHWnd.sBTop="0";
gCSHWnd.sBWidth="49%";
gCSHWnd.sBHeight="49%";

var gMapType = null;
var gMapData = "";
var gCshRootPathsArr = null;
function loadTopic(defaultTopicURL)
{
	gDefTopicURL = defaultTopicURL;
	var mapnum = getUrlParameter(RHMAPNO);
	if(mapnum!="")
	{
		showCSHTopic(ITEMTYPEMAPNO, mapnum);
		return true;
	}
	else
	{
		var mapid = getUrlParameter(RHMAPID);
		if(mapid!="")
		{
			showCSHTopic(ITEMTYPEMAPID, mapid);
			return true;
		}
	}
	gTopicURL = gDefTopicURL;
	redirectToTopic();
}

function showCSHTopic(maptype, mapdata)
{
	gMapType = maptype;
	gMapData = mapdata;
	initAndCollectAllChildPaths(gRootRelPath, gCommonRootRelPath, SCR_CHILD_CSH);
}
function loadCSH(cshRootPathsArr)
{
	gCshRootPathsArr = cshRootPathsArr;
	var arrIndex = 0;
	loadCSHFile(arrIndex);
}
function loadCSHFile(arrIndex)
{
	xmlJsReader.loadFile(gCshRootPathsArr[arrIndex] + "/" + gDataFolder + "/" + gCSHDataFile, callBackCSHLoaded, arrIndex);
}
function callBackCSHLoaded(xmlDoc, arrIndex)
{
	var cshInfoXmlNode = xmlDoc.getElementsByTagName(CSHINFONODE)[0];
	var itemNodes = cshInfoXmlNode.getElementsByTagName(ITEMNODE);
	var len = itemNodes.length;
	for(var i=0; i<len; i++)
	{
		var itemNode = itemNodes[i];
		var mapdata;
		if(gMapType == ITEMTYPEMAPNO)
			mapdata = itemNode.getAttribute(MAPNUM);
		else if(gMapType = ITEMTYPEMAPID)
			mapdata = itemNode.getAttribute(MAPID);
		if(mapdata.toLowerCase() == gMapData.toLowerCase())
		{
			gTopicURL = gCshRootPathsArr[arrIndex] + "/" + itemNode.getAttribute(TOPICURL);
			break;
		}
	}
	if(gTopicURL == null || gTopicURL == "")
	{
		if(arrIndex < gCshRootPathsArr.length-1)
		{
			loadCSHFile(arrIndex + 1);
			return;
		}
		else
			gTopicURL = gDefTopicURL;
	}	
	var cshModeParam = "";
	var appendChar = "";
	var cshModeFlag = getUrlParameter(RHCSHMODE);
	if(cshModeFlag == TRUESTR)
	{
		cshModeParam = RHCSHMODE + "=" + cshModeFlag;
		if(gTopicURL.indexOf("?") == -1)
			appendChar = "?";
		else
			appendChar = "&";
	}
	gTopicURL = gTopicURL + appendChar + cshModeParam;
	redirectToTopic();
}
function redirectToTopic()
{
	var bNewWin = getUrlParameter(RHNEWWINDOW);
	var strWndName = getUrlParameter(RHWINDOW);
	if(bNewWin == TRUESTR || strWndName != "")
	{
		if(strWndName != "")
			loadWindow(strWndName);
		else
			showTopicWindow(gCSHWnd);
	}
	else
		window.location = gTopicURL;

}

function loadWindow(strWndName)
{
	xmlJsReader.loadFile(gRootRelPath + "/" + gDataFolder + "/" + gWndDataFile, callBackWndLoaded, strWndName);
}
function callBackWndLoaded(xmlDoc, strWndName)
{
	var winListXmlNode = xmlDoc.getElementsByTagName(WINDOWLISTNODE)[0];
	var windowNodes = winListXmlNode.getElementsByTagName(WINDOWNODE);
	var len = windowNodes.length;
	for(var i=0; i<len; i++)
	{
		var windowNode = windowNodes[i];
		var name = windowNode.getAttribute(NAME);
		if(name.toLowerCase() == strWndName.toLowerCase())
		{
			var attribVal = windowNode.getAttribute(XCOORD);
			gCSHWnd.sBLeft = attribVal;
			
			attribVal = windowNode.getAttribute(YCOORD);
			gCSHWnd.sBTop = attribVal;
			
			attribVal = windowNode.getAttribute(WIDTH);
			gCSHWnd.sBWidth = attribVal;
			
			attribVal = windowNode.getAttribute(HEIGHT);
			gCSHWnd.sBHeight = attribVal;
			
			attribVal = windowNode.getAttribute(OPTIONS);
			gCSHWnd.nBOptions = parseInt(attribVal);
		}
	}
	showTopicWindow(gCSHWnd);
}
function showTopicWindow(oWnd)
{
	var	strOpt = getBrowserOptionString(oWnd);
	var	sNewName=convertWindowName(oWnd.sName);

	if(gbNav4 || gbSafari)
	{
		if (gbNav6)
		{
			if (navigator.appVersion.indexOf("rv:11.0") > -1)
			{
				// IE 11
				gBrowserWnd = window.open(gTopicURL, sNewName, strOpt);
			}
			else
			{
			gBrowserWnd=window.open("about:blank",sNewName,strOpt);
			setTimeout("postWindowNSOpen();",100);
		}
		}
		else
		{
			window.open("about:blank",sNewName,strOpt);
			var oNewWnd=window.open(gTopicURL,sNewName);
			window.close();
			oNewWnd.focus();
			top.blur();
		}
	}
	else
	{
		if(gbIE5)
		{
			var curWnd = null;	
			curWnd = window.open("about:blank",sNewName,strOpt);
			gBrowserWnd=window.open(gTopicURL,sNewName);
		}
		else
		{
			// IE4 had hard time to handle bookmark.
			gBrowserWnd=window.open("about:blank",sNewName,strOpt);
		}
		setTimeout("postWindowOpen();",100);
	}
}
function postWindowNSOpen()
{
	if(gBrowserWnd)
	{
		if (gTopicURL)
			gBrowserWnd.document.location.href=gTopicURL;
		window.close();
		gBrowserWnd.focus();
		top.blur();
	}
}

function postWindowOpen()
{
	if(gBrowserWnd)
	{
		if (gTopicURL&&!gbIE5&&gbIE4)
			gBrowserWnd.document.location.href=gTopicURL;
		gBrowserWnd.focus();
	}
}
function getBrowserOptionString(oWnd)
{
	var strOpts="";
	if(oWnd.bUseDefault)
		return strOpts;
	if(oWnd.nBOptions&WINLOCATION)
		strOpts+="location=yes";
	else
		strOpts+="location=no";
	if(oWnd.nBOptions&WINTOOLBAR)
		strOpts+=",toolbar=yes";		
	else
		strOpts+=",toolbar=no";		
	if(oWnd.nBOptions&WINMENUBAR)
		strOpts+=",menubar=yes";		
	else
		strOpts+=",menubar=no";
	if(oWnd.nBOptions&WINSTATUS)
		strOpts+=",status=yes";		
	else
		strOpts+=",status=no";		
	if(oWnd.nBOptions&WINSCROLLBARS)
		strOpts+=",scrollbars=yes";
	else
		strOpts+=",scrollbars=no";	
	if(oWnd.nBOptions&WINRESIZABLE)
		strOpts+=",resizable=yes";
	else
		strOpts+=",resizable=no";
	if(oWnd.sBTop)
	{
		var nTop=getSValue(oWnd.sBTop,screen.height);
		strOpts+=",top="+nTop;
		strOpts+=",screenY="+nTop;
	}
	if(oWnd.sBLeft)
	{
		var nLeft=getSValue(oWnd.sBLeft,screen.width);
		strOpts+=",left="+nLeft;
		strOpts+=",screenX="+nLeft;
	}
	if(oWnd.sBWidth)
	{
		var nWidth=getSValue(oWnd.sBWidth,screen.width);
		strOpts+=",width="+nWidth;
		strOpts+=",outerWidth="+nWidth;
	}
	if(oWnd.sBHeight)
	{
		var nHeight=getSValue(oWnd.sBHeight,screen.height);
		strOpts+=",height="+nHeight;
		strOpts+=",outerHeight="+nHeight;
	}
	return strOpts;
}
function getSValue(sValue,nLength)
{
	var nValue=0;
	var nPos=sValue.indexOf("%");
	if(nPos!=-1)
	{
		if(nPos>0)
		{
			var nPart=parseInt(sValue.substring(0,nPos));
			nValue=nLength*nPart/100;
		}
	}
	else
		nValue=parseInt(sValue);
	return nValue;
}
function convertWindowName(strName)
{
	var strNewName = strName;
	var strResultName = "";
	var re=new RegExp("_","g");
	strNewName = strName.replace(re,"__");
	for (var i=0;i<strNewName.length;i++)
		if (!(strNewName[i] == '_' ||
			(strNewName[i] <= '9' && strNewName[i] >= '0') ||
			(strNewName[i] <= 'z' && strNewName[i] >= 'a') ||
			(strNewName[i] <= 'Z' && strNewName[i] >= 'A')))
		{
			strResultName += "_" + strNewName.charCodeAt(i);
		}
		else
			strResultName += strNewName[i];
	return strResultName;
}
