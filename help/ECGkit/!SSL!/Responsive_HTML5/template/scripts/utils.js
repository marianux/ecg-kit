var gbBlockIOSScaling = 1;
var gbPreviewMode = 0;
var gRhEvtFuncsList = new Array;
var gHost = null;
var gHostPath = "/";
var gbRHLoadComplete = false;

addRhLoadCompleteEvent(initInputTextBoxes);

function blockIOSScaling()
{
	var metaTagsList = document.getElementsByTagName('meta');
	var i;
	if (navigator.userAgent.indexOf('iPad') != -1 || navigator.userAgent.indexOf('iPhone') != -1)
	{
		var contentString = "user-scalable=no, initial-scale=1.0, maximum-scale=1.0, minimum-scale=1.0";
		for (i=0; i<metaTagsList.length; i++)
		{
			if (metaTagsList[i].name == "viewport")
		        metaTagsList[i].content = contentString;
		}
		if(i == metaTagsList.length)
		{
			var metaTag = document.createElement('meta');
			metaTag.setAttribute("name", "viewport");
			metaTag.setAttribute("content", contentString);
			var headTags = document.getElementsByTagName('head');
			headTags[0].appendChild(metaTag);
		}
    }
}
function getUrlParamString(url)
{
	var paramstr = "";
	if(url == null || url == 'undefined')
		url = document.location.href;
	if(url.indexOf("?") != -1)
		paramstr = getUrlWithoutBookmark(url.substring(url.indexOf("?")));
	return paramstr;
}
function getUrlParameter(paramName, url)
{
  if(url == null || url == 'undefined')
	url = document.location.href;
  paramName = paramName.replace(/[\[]/,"\\\[").replace(/[\]]/,"\\\]");
  var regexS = "[\\?&#]"+paramName+"=([^&#]*)";
  var regex = new RegExp(regexS);
  var results = regex.exec(url);
  if( results == null )
    return "";
  else
    return results[1] && decodeURIComponent(results[1]);
}
function getUrlBookmark(url)
{
  var bookmark = "";
  if(url == null || url == 'undefined')
	url = document.location.href;
  if(url.indexOf("#") != -1)
	bookmark = url.substring(url.indexOf("#"));
  return bookmark;
}
function getUrlWithoutBookmark(url)
{
  if(url == null || url == 'undefined')
	url = document.location.href;
  
  var urlwithoutbookmark = url;
  
  if(url.indexOf("#") != -1)
	urlwithoutbookmark = url.substring(0, url.indexOf("#"));
  return urlwithoutbookmark;
}
function getUrlWithoutParameterAndBookMark(url)
{
  if(url == null || url == 'undefined')
	url = document.location.href;
  
  var urlwithoutparameter = url;

  if(url.indexOf("?") != -1)
	urlwithoutparameter = url.substring(0, url.indexOf("?"));
  return getUrlWithoutBookmark(urlwithoutparameter);
}
function GetSearchTextFromURL(bRemoveSlash)
{
	return getUrlParameter(RHSEARCHSTR);
}

function GetSynonymsFromURL()
{
	var strQuery = getUrlParameter(RHSYNSTR);
	return strQuery.split(" ");
}

function GetHighlightTextFromURL()
{
	return getUrlParameter(RHHIGHLIGHTTERM);
}

function lTrim(str)
{
	return str.replace(/^\s+/,'');
}

function MhStack()
{
    var container = new Array;

    this.getLength = function () {
        return (container.length - frontOffset);
    }
    this.isEmpty = function () {
        return (container.length == 0);
    }
    this.push = function (elem) {
        container.push(elem);
    }
    this.pop = function () {
        if (this.isEmpty()) return null;
        var elem = container[container.length-1];
		container.splice(container.length-1, 1);
        return elem;
    }
    this.peek = function () {
		if (this.isEmpty()) return null;
        var elem = container[container.length-1];
        return elem;
    }
}

function MhQueue() {
    var container = new Array;
    var frontOffset = 0;
    this.getLength = function () {
        return (container.length - frontOffset);
    }
    this.isEmpty = function () {
        return (container.length == 0);
    }
    this.enqueue = function (elem) {
        container.push(elem);
    }
    this.dequeue = function () {
        if (this.isEmpty()) return undefined;
        var elem = container[frontOffset];
        frontOffset++;
        if (frontOffset * 2 >= container.length) {
            container = container.slice(frontOffset);
            frontOffset = 0;
        }
        return elem;
    }
    this.peek = function () {
        if (container.length > 0)
            return container[frontOffset];
        return undefined;
    }
}

function getChildElementsByTagName(parentNode, tagName)
{
	if(parentNode == null || tagName == null || tagName == "" || parentNode.nodeType != JS_TAGTOKEN)
		return null;
	var childNodes = parentNode.childNodes;
	var childNodesArr = new Array;
	var len = childNodes.length;
	for(var i=0; i<len; i++)
	{
		var childNode = childNodes[i];
		if(childNode.nodeType == JS_TAGTOKEN && childNode.nodeName.toLowerCase() == tagName.toLowerCase())
			childNodesArr[childNodesArr.length] = childNode;		
	}
	return childNodesArr;
}

function trimString(str)
{
	str = str.replace(/^\s+/, '');
	for (var i = str.length - 1; i >= 0; i--) 
	{
		if (/\S/.test(str.charAt(i))) 
		{
			str = str.substring(0, i + 1);
			break;
		}
	}
	return str;
}

function objCookie(value, bPersistent)
{
	this.value = value;
	this.bPersistent = bPersistent;
}

function onTextBoxFocus()
{
	var placeholderText = this.getAttribute(DATAPH);
	if (trimString(this.value)==placeholderText)
		this.value = "";
}
function onTextBoxBlur()
{
	if (trimString(this.value).length==0)
	{
		var placeholderText = this.getAttribute(DATAPH);
		if(placeholderText != null)
			this.value = placeholderText;
	}
}
function initInputTextBoxes()
{
	var searchText = GetSearchTextFromURL(true);
	var inputs = document.getElementsByTagName('input');
	var len = inputs.length;
	var i=0;
	for(i=0; i<len; i++)
	{
		var searchAttr = inputs[i].getAttribute('data-search');
		if(searchAttr != null && searchAttr != 'undefined' && searchAttr == 'true' && searchText != "")
			inputs[i].value = searchText;		


		placeholderText = inputs[i].getAttribute(DATAPH);
		if(placeholderText != null)
		{
			if(gbIE5 && !gbIE10)
			{
				if(inputs[i].value == 'undefined' || inputs[i].value == "")
					inputs[i].value = placeholderText;
			}
			else
				inputs[i].setAttribute("placeholder", placeholderText);
		}


		if(gbIE5 && !gbIE10)
		{
			if(inputs[i].addEventListener)
			{
				var input = inputs[i];
				input.addEventListener('focus', onTextBoxFocus, false);
				input.addEventListener('blur', onTextBoxBlur, false);
			}
		}
	}
}

function isTouchDevice() {
    return "ontouchstart" in window;
}

////BreadCrumb functions Start
var gBreadCrumbInfo = new Array;
var gBCId = 0;
function BreadCrumbInfo(relHomePage, styleInfo, separator, strHome, strHomePath) {
    this.relHomePage = relHomePage;
    this.styleInfo = styleInfo;
    this.separator = separator;
    this.strHome = strHome;
    this.strHomePath = strHomePath;
    this.bcLinks = [];
}

function AddMasterBreadcrumbs(relHomePage, styleInfo, separator, strHome, strHomePath) {
    document.write("<span id=\"brseq" + gBCId + "\" ></span>");
    gBreadCrumbInfo[gBCId] = new BreadCrumbInfo(relHomePage, styleInfo, separator, strHome, strHomePath);
    gBCId++;
	addRhLoadCompleteEvent(UpdateBreadCrumbsMarker);
}

function UpdateBreadCrumbsMarker() {  
	if(gBreadCrumbInfo.length > 0)
	{
		if(gbPreviewMode)
			writeBreadCrumbs();
		else
			loadParentDataForSyncing(gCommonRootRelPath, SCR_PARENT_BC);
	}
}

function writeBreadCrumbs() {
    for(var i=0;i<gBCId;i++) {  
		var bHomeFound = false;
        var strTrail = "";
        if(gBreadCrumbInfo[i].bcLinks.length == 0)
        {   
	        if(gBreadCrumbInfo[i].styleInfo == "breadcrumbs")
		        strTrail = "<a class=\""+ gBreadCrumbInfo[i].styleInfo + "\"" + " href=\"" + gBreadCrumbInfo[i].strHomePath + "\">" + gBreadCrumbInfo[i].strHome + "</a> " + ((gBreadCrumbInfo[i].strHome == "")? "":gBreadCrumbInfo[i].separator) + " ";
	        else
	            strTrail = "<a style=\""+ gBreadCrumbInfo[i].styleInfo + "\"" + " href=\"" + gBreadCrumbInfo[i].strHomePath + "\">" + gBreadCrumbInfo[i].strHome + "</a> " + ((gBreadCrumbInfo[i].strHome == "")? "":gBreadCrumbInfo[i].separator) + " ";
        }
        else{
            var len = gBreadCrumbInfo[i].bcLinks.length;
			var bcName = "";
            for(var j=len-1;j>=0;j--)
            { 
				if(gBreadCrumbInfo[i].bcLinks[j].firstEntry == true)
				{
					if(bHomeFound)
						continue;
					else
						bHomeFound = true;
				}					

				bcName = gBreadCrumbInfo[i].bcLinks[j].name;
				
                if(gBreadCrumbInfo[i].bcLinks[j].strLink == "")
                {
                    strTrail += bcName + " " + gBreadCrumbInfo[i].separator + " ";
                }
                else{
                    if(gBreadCrumbInfo[i].styleInfo == "breadcrumbs")
 			            strTrail += "<a class=\""+ gBreadCrumbInfo[i].styleInfo + "\"" + " href=\"" + gBreadCrumbInfo[i].bcLinks[j].strLink + "\">" + bcName + "</a> " + gBreadCrumbInfo[i].separator + " ";
 			        else
 			            strTrail += "<a style=\""+ gBreadCrumbInfo[i].styleInfo + "\"" + " href=\"" + gBreadCrumbInfo[i].bcLinks[j].strLink + "\">" + bcName + "</a> " + gBreadCrumbInfo[i].separator + " ";
                }
            }
        }
        var brselem = document.getElementById("brseq"+i);
        brselem.innerHTML = strTrail;
    }
}

////BreadCrumb functions End

function addEvent(obj, type, func)
{
	if(obj == null || obj == 'undefined')
		return;
	if(obj.addEventListener)
		obj.addEventListener(type, func, false);
}
function removeEvent(obj, type, func)
{
	if(obj == null || obj == 'undefined')
		return;
	if(obj.removeEventListener)
		obj.removeEventListener(type, func, false);
}
function fireEvent(obj, type)
{
	if(obj == null || obj == 'undefined')
		return;
    if(document.createEventObject && obj!=window)
	{
		var evt = document.createEventObject();
		return obj.fireEvent('on'+type,evt)
    }
    else if(document.createEvent && obj.dispatchEvent)
	{
		var evt = document.createEvent("HTMLEvents");
		evt.initEvent(type, true, true ); // event type,bubbling,cancelable
		return !obj.dispatchEvent(evt);
    }
}
function preventEvent(e)
{
	if(e != null)
	{
		if (e.preventDefault)
			e.preventDefault(); 
		else
			e.returnValue = false;
	}
}
function addRhLoadCompleteEvent(func)
{
	if(gbRHLoadComplete)
		func();
	else
		gRhEvtFuncsList.push(func);
}
function fireRhLoadCompleteEvent()
{
	gbRHLoadComplete = true;
	var len = gRhEvtFuncsList.length;
	for(var i=0; i<len; i++)
		gRhEvtFuncsList[i]();
	gRhEvtFuncsList.splice(0, gRhEvtFuncsList.length);
}
function patchInputForSubmit(input, func)
{
	var formElem = document.createElement('form');
	formElem.setAttribute('method', 'POST');
	if(formElem.addEventListener)
		formElem.addEventListener('submit', func, false);
	else if(formElem.attachEvent)
		formElem.attachEvent('onsubmit', func, false);
	var parent = input.parentNode;
	parent.insertBefore(formElem, input);
	if(input.onkeydown)
		input.onkeydown = null;
	parent.removeChild(input);
	formElem.appendChild(input);
}