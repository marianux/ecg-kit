var gParentDataFile = "parentdata.js";
gFinalScrFolder = "";
gbLoadingParentData = false;
gbParentDataLoaded = false;
gTocChildPrefixStr = "";
gTocChildOrder = "";
gFlowTypeArrParentData = new Array;

function parentDataCallbackObj()
{
	this.path = null;
	this.flowType = SCR_NONE;
}

function initAndLoadParentData(path, flowType)
{
	if(gbLoadingParentData == true)
	{
		gFlowTypeArrParentData[gFlowTypeArrParentData.length] = flowType;
		return;
	}
	else if(gbParentDataLoaded == true)
	{
		doReturnParentDataCallAction(flowType);
		return;
	}
	gbLoadingParentData = true;
	gFlowTypeArrParentData[0] = flowType;
	gFinalRootRelPath = "";
	gFinalCommonRootRelPath = "";
	loadParentData(path, flowType);

}
function loadParentData(path, flowType)
{
	var parentDataCBObj = new parentDataCallbackObj();
	parentDataCBObj.path = path;
	parentDataCBObj.flowType = flowType;
	
	if(path == null)
		parentDataFile = gCommonRootRelPath + "/" + gParentDataFile;
	else
		parentDataFile = path + "/" + gParentDataFile;
	
	xmlJsReader.loadFile(parentDataFile, callbackParentDataLoaded, parentDataCBObj);
}
function callbackParentDataLoaded(xmlDoc, parentDataCBObj)
{
	var path = parentDataCBObj.path;
	var flowType = parentDataCBObj.flowType;
	var mpXmlTags = null;	
	var mpXmlTag = null;
	var parentName = null;
	var parentUrl = null;
	
	if (xmlDoc != undefined && xmlDoc != null)
	{
		mpXmlTags = xmlDoc.getElementsByTagName(MASTERPROJECT);	
		mpXmlTag = mpXmlTags[0];
		if (mpXmlTag != null)
		{
			parentName = mpXmlTag.getAttribute(NAME);
			parentUrl = mpXmlTag.getAttribute(URL);
		}
	}
	
	if(parentUrl == null)
	{
		if(path == null)
		{
			gFinalRootRelPath = gRootRelPath;
			gFinalCommonRootRelPath = gCommonRootRelPath;
		}
		returnParentDataCall();
	}
	else
	{
		if(path == null)
			path = gCommonRootRelPath;
		var parentPath = path + "/" + parentUrl;
		loadScreenData(parentPath, flowType);
		parentDataCBObj.path = parentPath;
	}
}

function returnParentDataCall()
{		
	gbLoadingParentData = false;
	gbParentDataLoaded = true;
	gLastScreenObj = null;
	initSettings(gFinalCommonRootRelPath);
	fireRhLoadCompleteEvent();
	for(var i=0; i<gFlowTypeArrParentData.length; i++)
		doReturnParentDataCallAction(gFlowTypeArrParentData[i]);
	
}

function doReturnParentDataCallAction(flowType)
{
	if(flowType == SCR_PARENT_TOC)
		displayToc(gFinalRootRelPath, gFinalCommonRootRelPath);
	else if(flowType == SCR_PARENT_IDX)
		initAndCollectAllChildPaths(gFinalRootRelPath, gFinalCommonRootRelPath, SCR_CHILD_IDX);
	else if(flowType == SCR_PARENT_GLO)
		initAndCollectAllChildPaths(gFinalRootRelPath, gFinalCommonRootRelPath, SCR_CHILD_GLO);
	else if(flowType == SCR_PARENT_FTS)
		initAndCollectAllChildPaths(gFinalRootRelPath, gFinalCommonRootRelPath, SCR_CHILD_FTS);
}

function loadParentDataForSyncing(commonRootRelPath, flowType)
{
	if(gbPreviewMode)
	{
		returnParentDataCallForSyncing(flowType, null);
		return;
	}
	var parentDataCBObj = new parentDataCallbackObj();
	parentDataCBObj.path = commonRootRelPath;
	parentDataCBObj.flowType = flowType;
	
	var parentDataFile = commonRootRelPath + "/" + gParentDataFile;
	xmlJsReader.loadFile(parentDataFile, callbackParentDataLoadedForSyncing, parentDataCBObj);
}

function callbackParentDataLoadedForSyncing(xmlDoc, parentDataCBObj)
{
	var mpXmlTags = null;	
	var mpXmlTag = null;
	var parentName = null;
	var parentUrl = null;
	
	if (xmlDoc != undefined && xmlDoc != null)
	{
		mpXmlTags = xmlDoc.getElementsByTagName(MASTERPROJECT);
		mpXmlTag = mpXmlTags[0];
		
		if (mpXmlTag != undefined && mpXmlTag != null)
		{
			parentName = mpXmlTag.getAttribute(NAME);
			parentUrl = mpXmlTag.getAttribute(URL);
		}
	}
	
	if(parentUrl != null)
	{
		var childName = mpXmlTag.getAttribute(CHILDNAME);
		loadScreenData(parentDataCBObj.path + "/" + parentUrl, parentDataCBObj.flowType, childName);
	}
	else
		returnParentDataCallForSyncing(parentDataCBObj.flowType, null);
}

function returnParentDataCallForSyncing(flowType, childName, finalRootRelPath, finalCommonRootRelPath)
{
	if(childName != null && childName != "")
		loadProjDataForSyncing(flowType, finalRootRelPath, finalCommonRootRelPath, childName);
	else
	{
		if(flowType == SCR_PARENT_TOCSYNC)
			syncToc(gTocChildPrefixStr, gTocChildOrder);
		else if(flowType == SCR_PARENT_BC)
			returnProjDataCallForSyncing(flowType, null, null);
	}
}