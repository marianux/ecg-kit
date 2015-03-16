var gCurTabId = null;
var gTabListIdArr = ["tab"];
var gTabContentIdArr = ["tocTabPane" , "idxTabPane" , "gloTabPane"];

addRhLoadCompleteEvent(loadTabContainer);

function loadTabContainer()
{
	var noTabContentPanes = gTabContentIdArr.length;
	for(var i=0; i<noTabContentPanes; i++)
	{
		var tabContentPaneId = gTabContentIdArr[i];
		var tabContentPaneElem = document.getElementById(tabContentPaneId);
		tabContentPaneElem.style.display = "none";
	}
		
	var noTabContainers = gTabListIdArr.length;
	for(var ii=0; ii<noTabContainers; ii++)
	{
		var tabId = gTabListIdArr[ii];
		var tabListElem = document.getElementById(tabId);
		selectTab(null, tabId);
	}
}

function selectTab(tabButtonElem, tabId)
{
	readSetting(TABBUTTONID + "-" + tabId, callbackSelectTab, tabButtonElem, tabId);
}
function callbackSelectTab(curTabButtonId, tabButtonElem, tabId)
{
	var tabElem = document.getElementById(tabId);
	
	if(curTabButtonId == "" || curTabButtonId == "null")
	{
		if(tabElem != null)
			curTabButtonId = tabElem.getAttribute(DATATABBUTTONID);
	}
	
	var newTabButtonId = "";
	var newTabButtonElem = null;
	if(tabButtonElem != null)
	{
		newTabButtonId = tabButtonElem.getAttribute('id');
		newTabButtonElem = tabButtonElem;
		var curTabButtonElem = document.getElementById(curTabButtonId);
		var curTabPaneId = "";
		if(curTabButtonElem != null)
		{
			var className = curTabButtonElem.getAttribute(DATACLASS);
			if(className == null)
				className = "";
			curTabButtonElem.className = className;
			curTabPaneId = curTabButtonElem.getAttribute(DATAPANEID);
		}
		if(curTabPaneId != null)
		{
			var curTabPaneElem = document.getElementById(curTabPaneId);
			if(curTabPaneElem != null)
				curTabPaneElem.style.display = "none";
		}
	}
	else
	{
		newTabButtonId = curTabButtonId;
		newTabButtonElem = document.getElementById(curTabButtonId);
		if(newTabButtonElem == null)
		{
			var tabbuttonsArr = tabElem.getElementsByTagName('li');
			if(tabbuttonsArr.length > 0)
				newTabButtonElem = tabbuttonsArr[0];
		}
	}

	gCurTabId = newTabButtonId;

	var newTabPaneId = "";
	if(newTabButtonElem != null)
	{
		var className = newTabButtonElem.getAttribute(DATACLASSSEL);
		if(className == null)
			className = "";
		newTabButtonElem.className = className;
		newTabPaneId = newTabButtonElem.getAttribute(DATAPANEID);
	}
	
	if(newTabPaneId != null)
	{
		var newTabPaneElem = document.getElementById(newTabPaneId);
		if(newTabPaneElem != null)
			newTabPaneElem.style.display = "block";
	}
	tabElem.setAttribute(DATATABBUTTONID, newTabButtonId);
	saveSetting(TABBUTTONID + "-" + tabId, newTabButtonId, true);
}

function onTabHover(tabButtonElem, tabId)
{
	if(tabButtonElem == null)
		return;
	var tabButtonId = tabButtonElem.getAttribute('id');
	if(tabButtonId == gCurTabId)
		return;
	var classHover = tabButtonElem.getAttribute(DATACLASSHOVER);
	if(classHover == null)
		classHover = "";
	tabButtonElem.className = classHover;
}

function onTabHoverOut(tabButtonElem, tabId)
{
	if(tabButtonElem == null)
		return;
	var tabButtonId = tabButtonElem.getAttribute('id');
	if(tabButtonId == gCurTabId)
		return;
	var classHover = tabButtonElem.getAttribute(DATACLASS);
	if(classHover == null)
		classHover = "";
	tabButtonElem.className = classHover;
}